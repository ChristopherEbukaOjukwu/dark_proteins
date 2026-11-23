#!/usr/bin/env python3
"""Step 5: every PubMed paper on our genes, by species, with year and R03 link.

Run from the dark_proteins folder, after steps 3 and 4:
    python3 scripts/05_gene_papers.py --email you@colorado.edu [--api-key KEY] [--refresh-r03]

Paper lists: JensenLab text-mining mention files (download.jensenlab.org, Oct 2026), one per species.
Each line = one protein, then every PubMed paper (PMID) that mentions it.

Inputs
    data/processed/r03_awards/gene_summary.tsv                 our 388 genes
    data/processed/r03_awards/foa_genes_with_awards.tsv        which genes each award funded, award start
    data/raw/jensenlab/mentions/<species>_textmining_mentions.tsv   human mouse rat zebrafish fly worm yeast pig
    data/raw/jensenlab/<taxid>.protein.aliases.v12.0.txt.gz    STRING: protein ID -> gene ID (per species)
    data/raw/pubtator3/human_all_hcop_sixteen_column.txt.gz    HCOP: human gene <-> gene in other species

How papers reach our genes
    human:  Jensen protein (ENSP) -> HGNC ID via STRING "HGNC:" lines (CILK1, ERN2 matched by hand)
    others: Jensen protein -> species gene (STRING "Ensembl_gene") -> human gene (HCOP)
            Papers stay under their own species. Each match is labelled, never dropped:
              support     how many orthology databases agree (1-12)
              match_type  decided only from matches with support >= MIN_SUPPORT:
                one-to-one     1 animal gene <-> 1 human gene
                animal copies  1 human gene <-> several animal genes (e.g. zebrafish duplicates)
                shared animal  1 animal gene <-> several human genes (e.g. 1 fly gene for a family)
                many-to-many   both
                weak           this match itself has support < MIN_SUPPORT

R03 papers
    RePORTER's paper list for all 98 awards (papers that list the R03 grant number). A paper counts
    as R03 for a gene only if (1) it cites an award that funded that gene, (2) the gene is found in it
    (by Jensen, any species, or confirmed by hand in the review file), and (3) it was published in or
    after the award's start year. Earlier papers are kept in r03_before_start.
    Manual review (optional): data/raw/r03_awards/r03_untied_papers_check.csv
        about_funded_gene = family -> status family-level (manual review), not tied
        about_funded_gene = yes  -> tied; gene_named = the funded gene(s) it is about (blank = all of
        the award's funded genes). A named gene the award did not fund is not counted (printed).

Saved along the way (reruns only fetch what is new)
    data/raw/pubmed/pmid_years.tsv                       PMID -> year from PubMed
    data/raw/r03_awards/reporter_r03_publications.json   RePORTER papers per award (rebuild: --refresh-r03)

Outputs (data/processed/attention/)
    gene_papers.tsv               one row per gene x species x paper
    gene_year_species_papers.tsv  papers per gene, year, species, match type; and how many are R03
    r03_papers.tsv                every paper that cites an R03, one row per paper x award, with status:
                                  tied | tied (manual review) | before award start |
                                  cites, no funded gene found | award has no confirmed gene
"""
import argparse
import csv
import gzip
import json
import re
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

GENES = Path("data/processed/r03_awards/gene_summary.tsv")
AWARDS = Path("data/processed/r03_awards/foa_genes_with_awards.tsv")
MENTIONS = Path("data/raw/jensenlab/mentions")
ALIASES = Path("data/raw/jensenlab")
HCOP = Path("data/raw/pubtator3/human_all_hcop_sixteen_column.txt.gz")
YEARS = Path("data/raw/pubmed/pmid_years.tsv")
R03_PUBS = Path("data/raw/r03_awards/reporter_r03_publications.json")
AWARD_LIST = Path("data/raw/r03_awards/R03_Award_List.csv")
R03_PROJECTS = Path("data/raw/r03_awards/reporter_r03_projects.json")
REVIEW = Path("data/raw/r03_awards/r03_untied_papers_check.csv")
TITLES = Path("data/raw/r03_awards/r03_paper_titles.tsv")
OUT_DIR = Path("data/processed/attention")
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
REPORTER_PUBS = "https://api.reporter.nih.gov/v2/publications/search"
SPECIES = {"9606": "human", "10090": "mouse", "10116": "rat", "7955": "zebrafish",
           "7227": "fly", "6239": "worm", "4932": "yeast", "9823": "pig"}
MIN_SUPPORT = 3
# No "HGNC:" line in STRING for these; their protein was found by gene name (checked by hand, step 4).
MANUAL_ENSP = {"CILK1": "ENSP00000349458", "ERN2": "ENSP00000256797"}


def read_tsv(path):
    if not path.is_file():
        raise SystemExit(f"Missing {path}")
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, cols):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def split_ids(v):
    return [x.strip() for x in (v or "").split(",") if x.strip() not in ("", "-")]


# ---------------------------------------------------------------- 1. which proteins are our genes
def string_aliases(taxid):
    """protein -> HGNC IDs (human) or protein -> species gene IDs (others); also all protein IDs."""
    path = ALIASES / f"{taxid}.protein.aliases.v12.0.txt.gz"
    if not path.is_file():
        raise SystemExit(f"Missing {path}")
    link, proteins = defaultdict(set), set()
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p, alias, src = line.rstrip("\n").split("\t")
            prot = p.split(".", 1)[1]
            proteins.add(prot)
            if taxid == "9606" and alias.startswith("HGNC:"):
                link[prot].add(alias)
            elif taxid != "9606" and src == "Ensembl_gene":
                link[prot].add(alias)
    return link, proteins


def hcop_matches(ours):
    """(taxid, species gene ID) -> list of match labels for our genes."""
    pairs, h_partners, a_partners = {}, defaultdict(set), defaultdict(set)
    with gzip.open(HCOP, "rt", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            t = r["ortholog_species"]
            if t not in SPECIES:
                continue
            n = len([s for s in r["support"].split(",") if s.strip()])
            hkey = r["hgnc_id"] if r["hgnc_id"] not in ("", "-") else r["human_ensembl_gene"]
            for ag in split_ids(r["ortholog_species_ensembl_gene"]):
                if n >= MIN_SUPPORT:
                    h_partners[(t, ag)].add(hkey)      # strong human partners of this animal gene
                    a_partners[(t, hkey)].add(ag)      # strong animal partners of this human gene
                if r["hgnc_id"] in ours:
                    k = (t, r["hgnc_id"], ag)
                    pairs[k] = max(pairs.get(k, 0), n)
    out = defaultdict(list)
    for (t, h, ag), n in pairs.items():
        nh, na = len(h_partners[(t, ag)]), len(a_partners[(t, h)])
        if n < MIN_SUPPORT:
            mt = "weak"
        else:
            mt = ("one-to-one" if nh == 1 and na == 1 else "animal copies" if nh == 1
                  else "shared animal" if na == 1 else "many-to-many")
        out[(t, ag)].append({"hgnc_id": h, "match_type": mt, "support": n,
                             "n_human_partners": nh, "n_animal_partners": na})
    return out


# ---------------------------------------------------------------- 2. papers from Jensen
def read_mentions(taxid, wanted, proteins):
    """Pass 1: papers for the wanted proteins. Pass 2: proteins per paper (for n_genes_in_paper)."""
    path = MENTIONS / f"{SPECIES[taxid]}_textmining_mentions.tsv"
    if not path.is_file():
        raise SystemExit(f"Missing {path}")
    papers = defaultdict(set)                       # protein -> PMIDs
    with path.open(encoding="utf-8") as f:
        for line in f:
            eid, _, rest = line.rstrip("\n").partition("\t")
            if eid in wanted:
                papers[eid].update(rest.split())
    our_pmids = set().union(*papers.values()) if papers else set()
    n_prot = defaultdict(int)
    with path.open(encoding="utf-8") as f:
        for line in f:
            eid, _, rest = line.rstrip("\n").partition("\t")
            if eid in proteins:
                for p in set(rest.split()) & our_pmids:
                    n_prot[p] += 1
    return papers, n_prot


# ---------------------------------------------------------------- 3. years
def fetch_years(pmids, email, api_key):
    known = {}
    if YEARS.is_file():
        known = {r["pmid"]: r["year"] for r in read_tsv(YEARS)}
    todo = sorted(set(pmids) - set(known), key=int)
    print(f"Years: {len(known):,} cached, {len(todo):,} to fetch from PubMed")
    YEARS.parent.mkdir(parents=True, exist_ok=True)
    new_file = not YEARS.is_file()
    pause = 0.11 if api_key else 0.35
    with YEARS.open("a", encoding="utf-8", newline="") as out:
        if new_file:
            out.write("pmid\tyear\n")
        for b in range(0, len(todo), 200):
            batch = todo[b:b + 200]
            params = {"db": "pubmed", "id": ",".join(batch), "retmode": "json", "tool": "dark_proteins", "email": email}
            if api_key:
                params["api_key"] = api_key
            data = None
            for attempt in range(4):
                try:
                    req = urllib.request.Request(ESUMMARY, data=urllib.parse.urlencode(params).encode())
                    data = json.load(urllib.request.urlopen(req, timeout=60))["result"]
                    break
                except Exception as e:
                    print(f"   retry {attempt + 1} for batch at {b}: {e}")
                    time.sleep(3 * (attempt + 1))
            if data is None:
                raise SystemExit("PubMed did not answer after 4 tries; rerun later (finished batches are kept).")
            for p in batch:
                rec = data.get(p, {})
                m = re.match(r"(\d{4})", rec.get("pubdate", "")) or re.match(r"(\d{4})", rec.get("sortpubdate", ""))
                known[p] = m.group(1) if m else ""
                out.write(f"{p}\t{known[p]}\n")
            out.flush()
            if (b // 200) % 50 == 0:
                print(f"   {min(b + 200, len(todo)):,} of {len(todo):,}")
            time.sleep(pause)
    return known


def fetch_titles(pmids, email, api_key):
    """PMID -> (journal, title) for the R03 papers, cached."""
    known = {}
    if TITLES.is_file():
        known = {r["pmid"]: (r["journal"], r["title"]) for r in read_tsv(TITLES)}
    todo = sorted(set(pmids) - set(known), key=int)
    for b in range(0, len(todo), 200):
        params = {"db": "pubmed", "id": ",".join(todo[b:b + 200]), "retmode": "json",
                  "tool": "dark_proteins", "email": email}
        if api_key:
            params["api_key"] = api_key
        req = urllib.request.Request(ESUMMARY, data=urllib.parse.urlencode(params).encode())
        data = json.load(urllib.request.urlopen(req, timeout=60))["result"]
        for p in todo[b:b + 200]:
            known[p] = (data.get(p, {}).get("source", ""), data.get(p, {}).get("title", ""))
        time.sleep(0.35)
    write_tsv(TITLES, [{"pmid": p, "journal": j, "title": t} for p, (j, t) in sorted(known.items())],
              ["pmid", "journal", "title"])
    return known


# ---------------------------------------------------------------- 4. R03 papers
def r03_papers(cores, refresh):
    cached = json.loads(R03_PUBS.read_text(encoding="utf-8")) if R03_PUBS.is_file() else {}
    if refresh or set(cached.get("cores", [])) != set(cores):
        print(f"Downloading RePORTER paper lists for {len(cores)} awards ...")
        results, offset = [], 0
        while True:
            body = {"criteria": {"core_project_nums": sorted(cores)}, "offset": offset, "limit": 500}
            req = urllib.request.Request(REPORTER_PUBS, data=json.dumps(body).encode(),
                                         headers={"Content-Type": "application/json"})
            data = json.load(urllib.request.urlopen(req, timeout=120))
            results += data.get("results", [])
            offset += 500
            if offset >= data.get("meta", {}).get("total", len(results)):
                break
        R03_PUBS.parent.mkdir(parents=True, exist_ok=True)
        R03_PUBS.write_text(json.dumps({"cores": sorted(cores), "results": results}, indent=1), encoding="utf-8")
    by_pmid = defaultdict(set)
    for r in json.loads(R03_PUBS.read_text(encoding="utf-8"))["results"]:
        by_pmid[str(r.get("pmid"))].add(r.get("coreproject", ""))
    return by_pmid


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True)
    ap.add_argument("--api-key", default="")
    ap.add_argument("--refresh-r03", action="store_true", help="download the R03 paper list again")
    a = ap.parse_args()

    ours = {r["hgnc_id"]: r["gene"] for r in read_tsv(GENES)}
    by_animal_gene = hcop_matches(ours)

    rows_raw = []          # (taxid, protein, hgnc, match dict, pmid, n_genes)
    per = defaultdict(lambda: defaultdict(set))
    for taxid, name in SPECIES.items():
        link, proteins = string_aliases(taxid)
        prot_matches = defaultdict(list)
        if taxid == "9606":
            for prot, hs in link.items():
                for h in hs & set(ours):
                    prot_matches[prot].append({"hgnc_id": h, "match_type": "same gene", "support": "",
                                               "n_human_partners": 1, "n_animal_partners": 1})
            for h, g in ours.items():
                if g in MANUAL_ENSP and not any(m["hgnc_id"] == h for ms in prot_matches.values() for m in ms):
                    prot_matches[MANUAL_ENSP[g]].append({"hgnc_id": h, "match_type": "same gene", "support": "",
                                                         "n_human_partners": 1, "n_animal_partners": 1})
        else:
            for prot, genes in link.items():
                for ag in genes:
                    for m in by_animal_gene.get((taxid, ag), []):
                        prot_matches[prot].append({**m, "animal_gene": ag})
        for ms in prot_matches.values():
            for m in ms:
                per[name][m["match_type"]].add(m["hgnc_id"])
        print(f"Reading {name} papers ...")
        papers, n_prot = read_mentions(taxid, set(prot_matches), proteins)
        for prot, pmids in papers.items():
            for m in prot_matches[prot]:
                for p in pmids:
                    rows_raw.append((name, prot, m, p, n_prot.get(p, 0)))

    print("\n===== Matches: our genes with at least one match of each type =====")
    types = ["same gene", "one-to-one", "animal copies", "shared animal", "many-to-many", "weak"]
    print(f"   {'species':<11}" + "".join(f"{t:>14}" for t in types))
    for name in SPECIES.values():
        print(f"   {name:<11}" + "".join(f"{len(per[name][t]):>14}" for t in types))
    print(f"   (strong = at least {MIN_SUPPORT} databases agree; a gene can appear under more than one type)")

    core_of = lambda pn: "R03" + re.search(r"[A-Z]{2}\d{6}", pn).group(0)
    with AWARD_LIST.open(encoding="utf-8-sig", errors="replace", newline="") as f:
        all_cores = {core_of(r["Project number"]) for r in csv.DictReader(f)}
    award_genes = defaultdict(set)                    # award -> confirmed listed genes it funded (step 3)
    for r in read_tsv(AWARDS):
        if r["project_number"]:
            award_genes[core_of(r["project_number"])].add(r["hgnc_id"])
    award_start = {}                                  # award -> start year (RePORTER)
    for r in json.loads(R03_PROJECTS.read_text(encoding="utf-8"))["results"]:
        c, d = r.get("core_project_num") or "", (r.get("project_start_date") or "")[:4]
        if c in all_cores and d:
            award_start[c] = min(award_start.get(c, "9999"), d)
    r03 = r03_papers(all_cores, a.refresh_r03)
    years = fetch_years({x[3] for x in rows_raw} | set(r03), a.email, a.api_key)
    titles = fetch_titles(set(r03), a.email, a.api_key)

    rows, seen = [], set()
    for name, prot, m, p, ng in rows_raw:
        h = m["hgnc_id"]
        key = (h, name, p)
        if key in seen:                       # same paper via two proteins of one gene in one species
            continue
        seen.add(key)
        y = years.get(p, "")
        cores = sorted(r03.get(p, ()))
        own = [c for c in cores if h in award_genes.get(c, ())]
        early = [c for c in own if y and y < award_start.get(c, "0000")]
        rows.append({"gene": ours[h], "hgnc_id": h, "species": name, "jensen_protein": prot,
                     "species_gene": m.get("animal_gene", ""), "match_type": m["match_type"],
                     "support": m["support"], "n_human_partners": m["n_human_partners"],
                     "n_animal_partners": m["n_animal_partners"], "pmid": p, "year": y,
                     "n_genes_in_paper": ng, "r03_cited": ";".join(cores),
                     "r03_for_this_gene": ";".join(c for c in own if c not in early),
                     "r03_before_start": ";".join(early)})

    # manual review: papers confirmed by hand to be about a funded gene
    sym2h = {g: h for h, g in ours.items()}
    review, manual_tied = {}, set()
    if REVIEW.is_file():
        with REVIEW.open(encoding="utf-8-sig", newline="") as f:
            for r in csv.DictReader(f):
                review[(r["pmid"].strip(), r["award"].strip())] = r
        by_key = defaultdict(list)
        for r in rows:
            by_key[(r["hgnc_id"], r["pmid"])].append(r)
        for (p, c), r in review.items():
            if r.get("about_funded_gene", "").strip().lower() != "yes":
                continue
            named = [x.strip() for x in r.get("gene_named", "").split(";") if x.strip()]
            hs = [sym2h.get(x, "") for x in named] if named else sorted(award_genes.get(c, ()))
            for x, h in zip(named or [ours[h] for h in hs], hs):
                if not h or h not in award_genes.get(c, ()):
                    print(f"   review: PMID {p} names {x}, which {c} did not fund (as a listed gene); not counted")
                    continue
                y = years.get(p, "")
                if y and y < award_start.get(c, "0000"):
                    continue
                manual_tied.add((p, c))
                if by_key[(h, p)]:
                    for row in by_key[(h, p)]:
                        row["r03_for_this_gene"] = ";".join(sorted(set(filter(None, row["r03_for_this_gene"].split(";"))) | {c}))
                else:
                    rows.append({"gene": ours[h], "hgnc_id": h, "species": "manual review", "jensen_protein": "",
                                 "species_gene": "", "match_type": "manual review", "support": "",
                                 "n_human_partners": "", "n_animal_partners": "", "pmid": p, "year": y,
                                 "n_genes_in_paper": "", "r03_cited": ";".join(sorted(r03.get(p, ()))),
                                 "r03_for_this_gene": c, "r03_before_start": ""})

    summ = defaultdict(lambda: [set(), set()])
    for r in rows:
        k = (r["hgnc_id"], r["gene"], r["year"], r["species"], r["match_type"])
        summ[k][0].add(r["pmid"])
        if r["r03_for_this_gene"]:
            summ[k][1].add(r["pmid"])
    srows = [{"gene": k[1], "hgnc_id": k[0], "year": k[2], "species": k[3], "match_type": k[4],
              "n_papers": len(v[0]), "n_r03_papers": len(v[1])} for k, v in sorted(summ.items())]

    print("\n===== Papers =====")
    print(f"Gene x species x paper rows: {len(rows):,} | distinct papers: {len({r['pmid'] for r in rows}):,}")
    print(f"Papers with no year from PubMed: {sum(1 for p in {r['pmid'] for r in rows} if not years.get(p))}")
    print("\nGene-paper pairs per year, our genes added together (check that recent years are full):")
    print("   year    human   other species, one-to-one   other species, any match")
    for y in range(2015, 2027):
        ys = str(y)
        hu = sum(1 for r in rows if r["year"] == ys and r["species"] == "human")
        oo = sum(1 for r in rows if r["year"] == ys and r["species"] != "human" and r["match_type"] == "one-to-one")
        an = sum(1 for r in rows if r["year"] == ys and r["species"] != "human")
        print(f"   {y}  {hu:>7}  {oo:>26}  {an:>24}")

    # one row per R03-citing paper x award, with a status
    tied_by = defaultdict(set)                         # (pmid, award) -> genes found
    for r in rows:
        for c in filter(None, r["r03_for_this_gene"].split(";")):
            tied_by[(r["pmid"], c)].add(r["gene"])
    r03_rows = []
    for p in sorted(r03, key=int):
        for c in sorted(r03[p]):
            y = years.get(p, "")
            if not award_genes.get(c):
                status = "award has no confirmed gene"
            elif (p, c) in tied_by:
                status = "tied (manual review)" if (p, c) in manual_tied else "tied"
            elif y and y < award_start.get(c, "0000"):
                status = "before award start"
            elif review.get((p, c), {}).get("about_funded_gene", "").strip().lower() == "family":
                status = "family-level (manual review)"
            else:
                status = "cites, no funded gene found"
            rv = review.get((p, c), {})
            r03_rows.append({"award": c, "award_genes": "; ".join(sorted(ours[h] for h in award_genes.get(c, ()))),
                             "award_start": award_start.get(c, ""), "pmid": p, "year": y, "status": status,
                             "funded_genes_found": "; ".join(sorted(tied_by.get((p, c), ()))),
                             "journal": titles.get(p, ("", ""))[0], "title": titles.get(p, ("", ""))[1],
                             "review_about_funded_gene": rv.get("about_funded_gene", ""),
                             "review_notes": rv.get("notes", "")})

    print("\n===== R03 papers (papers that list an R03 grant number) =====")
    print(f"Distinct papers: {len(r03)} | paper x award rows: {len(r03_rows)} | awards with papers: "
          f"{len({r['award'] for r in r03_rows})} of {len(all_cores)}")
    for st, n in Counter(r["status"] for r in r03_rows).most_common():
        print(f"   {st}: {n}")
    print("Tied papers per award (first 10 by count):")
    per_award = Counter(r["award"] for r in r03_rows if r["status"].startswith("tied"))
    for c, n in per_award.most_common(10):
        print(f"   {c}: {n}")

    cols = ["gene", "hgnc_id", "species", "jensen_protein", "species_gene", "match_type", "support",
            "n_human_partners", "n_animal_partners", "pmid", "year", "n_genes_in_paper",
            "r03_cited", "r03_for_this_gene", "r03_before_start"]
    write_tsv(OUT_DIR / "gene_papers.tsv", rows, cols)
    write_tsv(OUT_DIR / "gene_year_species_papers.tsv", srows, list(srows[0]) if srows else ["gene"])
    write_tsv(OUT_DIR / "r03_papers.tsv", r03_rows, list(r03_rows[0]) if r03_rows else ["award"])
    print(f"\nSaved: {OUT_DIR / 'gene_papers.tsv'} ({len(rows):,} rows)")
    print(f"Saved: {OUT_DIR / 'gene_year_species_papers.tsv'} ({len(srows):,} rows)")
    print(f"Saved: {OUT_DIR / 'r03_papers.tsv'} ({len(r03_rows):,} rows)")


if __name__ == "__main__":
    main()
