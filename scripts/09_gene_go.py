#!/usr/bin/env python3
"""Step 9: Gene Ontology annotations for our genes, by the year of the supporting paper, with the R03 tie.

Run from the dark_proteins folder, after steps 3 and 5:
    python3 scripts/09_gene_go.py --email you@colorado.edu

Inputs
    data/raw/go/goa_human.gaf.gz                          GO annotations for human (GAF 2.2)
    data/processed/r03_awards/gene_summary.tsv            our 388 genes
    data/raw/hgnc/hgnc_complete_set.txt                   HGNC ID -> UniProt ID(s)
    data/processed/r03_awards/foa_genes_with_awards.tsv   which genes each R03 funded
    data/raw/r03_awards/reporter_r03_projects.json        R03 start dates
    data/raw/r03_awards/reporter_r03_publications.json    papers citing each R03
    data/raw/pubmed/pmid_years.tsv                        PMID -> year (new PMIDs are fetched and added)

Labels kept on every annotation (nothing filtered):
    evidence_group  experimental (EXP IDA IPI IMP IGI IEP) | high-throughput (HTP HDA HMP HGI HEP) |
                    phylogenetic (IBA IBD IKR IRD) | computational (ISS ISO ISA ISM IGC RCA) |
                    author/curator (TAS NAS IC ND) | electronic (IEA)
    aspect          F function | P process | C location
    negated         1 if the qualifier starts with NOT (the gene does NOT do this)
    informative     0 for "protein binding" (GO:0005515) and "identical protein binding" (GO:0042802)
    year            publication year of the supporting paper (blank if the reference is not a PMID);
                    annotation_date is GO's last-update date, not when the annotation was first made
    r03             the paper cites an R03 that funded this gene, in or after its start year

Outputs (data/processed/go/)
    gene_go_annotations.tsv  one row per gene x annotation
    gene_year_go.tsv         per gene x year (paper year):
        n_experimental, n_high_throughput, n_other (annotations, not negated)
        n_new_terms_experimental   GO terms with experimental support for this gene for the first time
        n_new_terms_informative    the same, without the two protein-binding terms
        n_new_F / n_new_P / n_new_C  the same, by aspect
        n_papers_experimental, n_r03_annotations, n_r03_papers
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

GAF = Path("data/raw/go/goa_human.gaf.gz")
GENES = Path("data/processed/r03_awards/gene_summary.tsv")
HGNC = Path("data/raw/hgnc/hgnc_complete_set.txt")
AWARDS = Path("data/processed/r03_awards/foa_genes_with_awards.tsv")
R03_PROJECTS = Path("data/raw/r03_awards/reporter_r03_projects.json")
R03_PUBS = Path("data/raw/r03_awards/reporter_r03_publications.json")
YEARS = Path("data/raw/pubmed/pmid_years.tsv")
OUT_DIR = Path("data/processed/go")
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
GROUPS = {"experimental": "EXP IDA IPI IMP IGI IEP", "high-throughput": "HTP HDA HMP HGI HEP",
          "phylogenetic": "IBA IBD IKR IRD", "computational": "ISS ISO ISA ISM IGC RCA",
          "author/curator": "TAS NAS IC ND", "electronic": "IEA"}
GROUP_OF = {code: g for g, codes in GROUPS.items() for code in codes.split()}
# Terms that only say "binds some protein": mostly from interaction screens, already covered by BioGRID.
UNINFORMATIVE = {"GO:0005515": "protein binding", "GO:0042802": "identical protein binding"}


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


def fetch_years(pmids, email):
    known = {r["pmid"]: r["year"] for r in read_tsv(YEARS)} if YEARS.is_file() else {}
    todo = sorted(set(pmids) - set(known), key=int)
    print(f"Paper years: {len(known):,} cached, {len(todo):,} to fetch from PubMed")
    new_file = not YEARS.is_file()
    YEARS.parent.mkdir(parents=True, exist_ok=True)
    with YEARS.open("a", encoding="utf-8") as out:
        if new_file:
            out.write("pmid\tyear\n")
        for b in range(0, len(todo), 200):
            batch = todo[b:b + 200]
            q = {"db": "pubmed", "id": ",".join(batch), "retmode": "json", "tool": "dark_proteins", "email": email}
            data = None
            for attempt in range(4):
                try:
                    req = urllib.request.Request(ESUMMARY, data=urllib.parse.urlencode(q).encode())
                    data = json.load(urllib.request.urlopen(req, timeout=60))["result"]
                    break
                except Exception as e:
                    print(f"   retry {attempt + 1}: {e}")
                    time.sleep(3 * (attempt + 1))
            if data is None:
                raise SystemExit("PubMed did not answer; rerun later (finished batches are kept).")
            for p in batch:
                rec = data.get(p, {})
                m = re.match(r"(\d{4})", rec.get("pubdate", "")) or re.match(r"(\d{4})", rec.get("sortpubdate", ""))
                known[p] = m.group(1) if m else ""
                out.write(f"{p}\t{known[p]}\n")
            time.sleep(0.35)
    return known


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True)
    a = ap.parse_args()

    ours = {r["hgnc_id"]: r["gene"] for r in read_tsv(GENES)}
    by_uniprot, by_symbol = {}, {}
    for r in read_tsv(HGNC):
        if r["hgnc_id"] in ours:
            for u in (r.get("uniprot_ids") or "").replace('"', "").split("|"):
                if u:
                    by_uniprot[u] = r["hgnc_id"]
            by_symbol[r["symbol"]] = r["hgnc_id"]

    header, raw = [], []
    with gzip.open(GAF, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("!"):
                if line.startswith("!date-generated") and not header:
                    header.append(line.strip())
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 15:
                continue
            h = by_uniprot.get(p[1]) or by_symbol.get(p[2])
            if h:
                raw.append((h, p))
    print(f"GO file: {header[0] if header else '(no date line)'} | annotation lines for our genes: {len(raw):,}")

    pmids = {ref.split(":", 1)[1] for _, p in raw for ref in p[5].split("|") if ref.startswith("PMID:")}
    years = fetch_years(pmids, a.email)

    core = lambda pn: "R03" + re.search(r"[A-Z]{2}\d{6}", pn).group(0)
    award_genes = defaultdict(set)
    for r in read_tsv(AWARDS):
        if r["project_number"]:
            award_genes[core(r["project_number"])].add(r["hgnc_id"])
    award_start = {}
    for r in json.loads(R03_PROJECTS.read_text(encoding="utf-8"))["results"]:
        c, d = r.get("core_project_num") or "", (r.get("project_start_date") or "")[:4]
        if d:
            award_start[c] = min(award_start.get(c, "9999"), d)
    r03 = defaultdict(set)
    for r in json.loads(R03_PUBS.read_text(encoding="utf-8"))["results"]:
        r03[str(r.get("pmid"))].add(r.get("coreproject", ""))

    rows = []
    for h, p in raw:
        refs = p[5].split("|")
        pmid = next((x.split(":", 1)[1] for x in refs if x.startswith("PMID:")), "")
        year = years.get(pmid, "") if pmid else ""
        own = [c for c in sorted(r03.get(pmid, ())) if h in award_genes.get(c, ())
               and year and year >= award_start.get(c, "9999")] if pmid else []
        rows.append({"gene": ours[h], "hgnc_id": h, "go_id": p[4], "aspect": p[8], "qualifier": p[3],
                     "negated": int(p[3].upper().startswith("NOT")), "informative": int(p[4] not in UNINFORMATIVE),
                     "evidence": p[6],
                     "evidence_group": GROUP_OF.get(p[6], "other"), "reference": p[5], "pmid": pmid, "year": year,
                     "annotation_date": p[13], "assigned_by": p[14], "r03_for_this_gene": ";".join(own)})

    # first year a term gets experimental support for a gene
    first = {}
    for r in rows:
        if r["evidence_group"] == "experimental" and not r["negated"] and r["year"]:
            k = (r["hgnc_id"], r["go_id"])
            first[k] = min(first.get(k, "9999"), r["year"])
    agg = defaultdict(lambda: defaultdict(set))
    cnt = defaultdict(Counter)
    for r in rows:
        if r["negated"] or not r["year"]:
            continue
        k = (r["hgnc_id"], r["gene"], r["year"])
        grp = r["evidence_group"]
        cnt[k]["n_experimental" if grp == "experimental" else "n_high_throughput" if grp == "high-throughput" else "n_other"] += 1
        if grp == "experimental":
            agg[k]["papers"].add(r["pmid"])
            if first.get((r["hgnc_id"], r["go_id"])) == r["year"]:
                agg[k]["new"].add(r["go_id"])
                agg[k]["new_" + r["aspect"]].add(r["go_id"])
                if r["informative"]:
                    agg[k]["new_inf"].add(r["go_id"])
        if r["r03_for_this_gene"]:
            cnt[k]["n_r03_annotations"] += 1
            agg[k]["r03p"].add(r["pmid"])
    yrows = []
    for k in sorted(set(cnt) | set(agg), key=lambda k: (k[1], k[2])):
        c, s = cnt[k], agg[k]
        yrows.append({"gene": k[1], "hgnc_id": k[0], "year": k[2], "n_experimental": c["n_experimental"],
                      "n_high_throughput": c["n_high_throughput"], "n_other": c["n_other"],
                      "n_new_terms_experimental": len(s["new"]), "n_new_terms_informative": len(s["new_inf"]),
                      "n_new_F": len(s["new_F"]),
                      "n_new_P": len(s["new_P"]), "n_new_C": len(s["new_C"]),
                      "n_papers_experimental": len(s["papers"]), "n_r03_annotations": c["n_r03_annotations"],
                      "n_r03_papers": len(s["r03p"])})

    # report
    print(f"\nAnnotations for our genes: {len(rows):,} | genes with any: {len({r['hgnc_id'] for r in rows})} of {len(ours)}")
    print("By evidence group:", dict(Counter(r["evidence_group"] for r in rows).most_common()))
    pb = [r for r in rows if not r["informative"] and r["evidence_group"] == "experimental"]
    print(f"Experimental annotations that are only 'protein binding': {len(pb):,} "
          f"(from {len({r['pmid'] for r in pb}):,} papers)")
    exp_genes = {r["hgnc_id"] for r in rows if r["evidence_group"] == "experimental" and not r["negated"]}
    print(f"Genes with at least one experimental annotation: {len(exp_genes)} of {len(ours)}")
    print("\nNew experimentally supported terms per paper year, all our genes (and papers behind them):")
    for y in range(2010, 2027):
        v = [x for x in yrows if x["year"] == str(y)]
        print(f"   {y}  new terms {sum(x['n_new_terms_experimental'] for x in v):>4} "
              f"(informative {sum(x['n_new_terms_informative'] for x in v):>4})   "
              f"papers {sum(x['n_papers_experimental'] for x in v):>4}   "
              f"high-throughput annotations {sum(x['n_high_throughput'] for x in v):>5}")
    tied = [r for r in rows if r["r03_for_this_gene"]]
    print(f"\nR03-tied annotations: {len(tied)} | papers: {len({r['pmid'] for r in tied})} | "
          f"genes: {sorted({r['gene'] for r in tied})}")
    print("\nCheck three genes (paper year: experimental annotations / new terms):")
    for g in ("GPR68", "PRKCQ", "ADGRB3"):
        v = [x for x in yrows if x["gene"] == g and x["year"] >= "2012" and x["n_experimental"]]
        print(f"   {g}: " + ("  ".join(f"{x['year']}: {x['n_experimental']}/{x['n_new_terms_experimental']}" for x in v) or "none"))

    cols = ["gene", "hgnc_id", "go_id", "aspect", "qualifier", "negated", "informative", "evidence", "evidence_group", "reference",
            "pmid", "year", "annotation_date", "assigned_by", "r03_for_this_gene"]
    write_tsv(OUT_DIR / "gene_go_annotations.tsv", rows, cols)
    write_tsv(OUT_DIR / "gene_year_go.tsv", yrows, list(yrows[0]) if yrows else ["gene"])
    print(f"\nSaved: {OUT_DIR / 'gene_go_annotations.tsv'} ({len(rows):,} rows)")
    print(f"Saved: {OUT_DIR / 'gene_year_go.tsv'} ({len(yrows):,} rows)")


if __name__ == "__main__":
    main()
