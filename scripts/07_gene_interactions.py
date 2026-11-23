#!/usr/bin/env python3
"""Step 7: BioGRID interactions for our genes, by species and year, with the R03 tie.

Run from the dark_proteins folder, after steps 3 and 5:
    python3 scripts/07_gene_interactions.py

Inputs
    data/raw/biogrid/BIOGRID-ALL-LATEST.tab3.zip               BioGRID (release in the file name inside)
    data/processed/r03_awards/gene_summary.tsv                 our 388 genes
    data/processed/r03_awards/foa_genes_with_awards.tsv        which genes each R03 funded
    data/raw/r03_awards/reporter_r03_projects.json             R03 start dates
    data/raw/r03_awards/reporter_r03_publications.json         papers citing each R03 (from step 5)
    data/raw/hgnc/hgnc_complete_set.txt                        HGNC ID -> human NCBI Gene ID
    data/raw/pubtator3/human_all_hcop_sixteen_column.txt.gz    HCOP: human gene <-> gene in other species
    data/raw/pubmed/pmid_years.tsv                             PMID -> year (fallback for the year)

Which records
    Any BioGRID record where one partner is one of our genes: human = the gene's own NCBI Gene ID;
    other species = its HCOP ortholog (NCBI Gene ID), labelled exactly as for the papers:
      support (databases agreeing) and match_type (one-to-one, animal copies, shared animal,
      many-to-many; weak = fewer than MIN_SUPPORT databases).
Labels kept on every record (nothing is filtered): physical/genetic, Low/High throughput,
experimental system, partner, PMID, year (from BioGRID's author field, else PubMed), R03 tie.
R03 tie (same rule as the papers): the record's paper cites an R03 that funded this gene,
published in or after the award's start year.

Outputs (data/processed/interactions/)
    gene_interactions.tsv       one row per our gene x BioGRID record
    gene_year_interactions.tsv  per gene x year x species x match type:
        n_records, n_low, n_high, n_physical, n_genetic, n_papers,
        n_new_partners      partners seen for the first time that year (for this gene/species/match type)
        n_new_partners_low  same, counting only low-throughput evidence
        n_r03_records, n_r03_papers
"""
import csv
import gzip
import io
import json
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

BIOGRID = Path("data/raw/biogrid/BIOGRID-ALL-LATEST.tab3.zip")
GENES = Path("data/processed/r03_awards/gene_summary.tsv")
AWARDS = Path("data/processed/r03_awards/foa_genes_with_awards.tsv")
R03_PROJECTS = Path("data/raw/r03_awards/reporter_r03_projects.json")
R03_PUBS = Path("data/raw/r03_awards/reporter_r03_publications.json")
HGNC = Path("data/raw/hgnc/hgnc_complete_set.txt")
HCOP = Path("data/raw/pubtator3/human_all_hcop_sixteen_column.txt.gz")
YEARS = Path("data/raw/pubmed/pmid_years.tsv")
OUT_DIR = Path("data/processed/interactions")
MIN_SUPPORT = 3


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


def matches(ours):
    """NCBI Gene ID -> list of {hgnc_id, match_type, support} for our genes (human + orthologs)."""
    entrez = {r["hgnc_id"]: r.get("entrez_id", "") for r in read_tsv(HGNC)}
    by_id = defaultdict(list)
    for h in ours:
        if entrez.get(h):
            by_id[entrez[h]].append({"hgnc_id": h, "match_type": "same gene", "support": ""})
    pairs, h_partners, a_partners = {}, defaultdict(set), defaultdict(set)
    with gzip.open(HCOP, "rt", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            t, oe = r["ortholog_species"], r["ortholog_species_entrez_gene"]
            if oe in ("", "-"):
                continue
            n = len([s for s in r["support"].split(",") if s.strip()])
            hkey = r["hgnc_id"] if r["hgnc_id"] not in ("", "-") else r["human_entrez_gene"]
            if n >= MIN_SUPPORT:
                h_partners[(t, oe)].add(hkey)
                a_partners[(t, hkey)].add(oe)
            if r["hgnc_id"] in ours:
                k = (t, r["hgnc_id"], oe)
                pairs[k] = max(pairs.get(k, 0), n)
    for (t, h, oe), n in pairs.items():
        nh, na = len(h_partners[(t, oe)]), len(a_partners[(t, h)])
        mt = ("weak" if n < MIN_SUPPORT else "one-to-one" if nh == 1 and na == 1
              else "animal copies" if nh == 1 else "shared animal" if na == 1 else "many-to-many")
        by_id[oe].append({"hgnc_id": h, "match_type": mt, "support": n})
    return by_id


def main():
    ours = {r["hgnc_id"]: r["gene"] for r in read_tsv(GENES)}
    by_id = matches(ours)

    # R03 tie inputs
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
    years = {r["pmid"]: r["year"] for r in read_tsv(YEARS)} if YEARS.is_file() else {}

    zf = zipfile.ZipFile(BIOGRID)
    name = zf.namelist()[0]
    print(f"BioGRID file: {name}")
    rows, n_lines = [], 0
    with io.TextIOWrapper(zf.open(name), encoding="utf-8", errors="replace") as f:
        header = f.readline().lstrip("#").rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        for line in f:
            n_lines += 1
            p = line.rstrip("\n").split("\t")
            a, b = p[col["Entrez Gene Interactor A"]], p[col["Entrez Gene Interactor B"]]
            if a not in by_id and b not in by_id:
                continue
            src = p[col["Publication Source"]]
            pmid = src.split(":", 1)[1] if src.upper().startswith("PUBMED:") else ""
            m = re.search(r"\((\d{4})\)", p[col["Author"]])
            year = m.group(1) if m else years.get(pmid, "")
            thr = p[col["Throughput"]]
            cores = sorted(r03.get(pmid, ())) if pmid else []
            base = {"biogrid_id": p[col["BioGRID Interaction ID"]],
                    "experimental_system": p[col["Experimental System"]],
                    "type": p[col["Experimental System Type"]],
                    "throughput": "low+high" if "Low" in thr and "High" in thr else "low" if "Low" in thr
                                  else "high" if "High" in thr else "",
                    "pmid": pmid, "source": src if not pmid else "", "year": year}
            sides = [(a, b, "A"), (b, a, "B")] if a != b else [(a, b, "A")]
            for me, partner, side in sides:
                for mt in by_id.get(me, []):
                    h = mt["hgnc_id"]
                    own = [c for c in cores if h in award_genes.get(c, ())
                           and year and year >= award_start.get(c, "9999")]
                    rows.append({"gene": ours[h], "hgnc_id": h,
                                 "species": p[col[f"Organism Name Interactor {side}"]],
                                 "gene_id_in_record": me, "match_type": mt["match_type"], "support": mt["support"],
                                 "partner_gene_id": partner,
                                 "partner_symbol": p[col["Official Symbol Interactor " + ("B" if side == "A" else "A")]],
                                 "partner_species": p[col["Organism Name Interactor " + ("B" if side == "A" else "A")]],
                                 "self_interaction": int(a == b), **base,
                                 "r03_cited": ";".join(cores), "r03_for_this_gene": ";".join(own)})
    print(f"BioGRID records read: {n_lines:,} | gene x record rows for our genes: {len(rows):,}")

    # first year each partner is seen (any evidence / low-throughput evidence)
    first, first_low = {}, {}
    for r in rows:
        if not r["year"]:
            continue
        k = (r["hgnc_id"], r["species"], r["match_type"], r["partner_gene_id"])
        first[k] = min(first.get(k, "9999"), r["year"])
        if "low" in r["throughput"]:
            first_low[k] = min(first_low.get(k, "9999"), r["year"])
    agg = defaultdict(lambda: Counter())
    sets = defaultdict(lambda: defaultdict(set))
    for r in rows:
        k = (r["hgnc_id"], r["gene"], r["year"], r["species"], r["match_type"])
        c, s = agg[k], sets[k]
        c["n_records"] += 1
        c["n_low"] += "low" in r["throughput"]
        c["n_high"] += "high" in r["throughput"]
        c["n_physical"] += r["type"] == "physical"
        c["n_genetic"] += r["type"] == "genetic"
        if r["pmid"]:
            s["papers"].add(r["pmid"])
        pk = (r["hgnc_id"], r["species"], r["match_type"], r["partner_gene_id"])
        if r["year"] and first.get(pk) == r["year"]:
            s["new"].add(r["partner_gene_id"])
        if r["year"] and first_low.get(pk) == r["year"]:
            s["new_low"].add(r["partner_gene_id"])
        if r["r03_for_this_gene"]:
            c["n_r03_records"] += 1
            s["r03_papers"].add(r["pmid"])
    yrows = []
    for k in sorted(agg, key=lambda k: (k[1], k[2], k[3], k[4])):
        c, s = agg[k], sets[k]
        yrows.append({"gene": k[1], "hgnc_id": k[0], "year": k[2], "species": k[3], "match_type": k[4],
                      "n_records": c["n_records"], "n_low": c["n_low"], "n_high": c["n_high"],
                      "n_physical": c["n_physical"], "n_genetic": c["n_genetic"], "n_papers": len(s["papers"]),
                      "n_new_partners": len(s["new"]), "n_new_partners_low": len(s["new_low"]),
                      "n_r03_records": c["n_r03_records"], "n_r03_papers": len(s["r03_papers"])})

    # report
    print("\nRecords by species and match type (top 12):")
    sp = Counter((r["species"], r["match_type"]) for r in rows)
    for (s, m), n in sp.most_common(12):
        print(f"   {s:<32} {m:<14} {n:>8,}")
    hum = [r for r in rows if r["match_type"] == "same gene"]
    print(f"\nHuman: our genes with any interaction: {len({r['hgnc_id'] for r in hum})} of {len(ours)} | "
          f"records {len(hum):,} | low {sum('low' in r['throughput'] for r in hum):,} | "
          f"high {sum('high' in r['throughput'] for r in hum):,}")
    print("Human records per year (low / high throughput) and new partners (any / low):")
    for y in range(2012, 2027):
        ys = str(y)
        rr = [r for r in hum if r["year"] == ys]
        nn = sum(v["n_new_partners"] for v in yrows if v["year"] == ys and v["match_type"] == "same gene")
        nl = sum(v["n_new_partners_low"] for v in yrows if v["year"] == ys and v["match_type"] == "same gene")
        print(f"   {y}  records {len(rr):>6} ({sum('low' in r['throughput'] for r in rr):>5} / "
              f"{sum('high' in r['throughput'] for r in rr):>6})  new partners {nn:>5} / {nl:>4}")
    tied = [r for r in rows if r["r03_for_this_gene"]]
    print(f"\nR03-tied records: {len(tied)} | papers: {len({r['pmid'] for r in tied})} | "
          f"genes: {sorted({r['gene'] for r in tied})}")
    print("\nCheck three genes (human, by year: records / new partners):")
    for g in ("GPR68", "PRKCQ", "ADGRB3"):
        v = [x for x in yrows if x["gene"] == g and x["match_type"] == "same gene" and x["year"] >= "2015"]
        print(f"   {g}: " + "  ".join(f"{x['year']}: {x['n_records']}/{x['n_new_partners']}" for x in v))

    cols = ["gene", "hgnc_id", "species", "gene_id_in_record", "match_type", "support", "partner_gene_id",
            "partner_symbol", "partner_species", "self_interaction", "biogrid_id", "experimental_system", "type",
            "throughput", "pmid", "source", "year", "r03_cited", "r03_for_this_gene"]
    write_tsv(OUT_DIR / "gene_interactions.tsv", rows, cols)
    write_tsv(OUT_DIR / "gene_year_interactions.tsv", yrows, list(yrows[0]) if yrows else ["gene"])
    print(f"\nSaved: {OUT_DIR / 'gene_interactions.tsv'} ({len(rows):,} rows)")
    print(f"Saved: {OUT_DIR / 'gene_year_interactions.tsv'} ({len(yrows):,} rows)")


if __name__ == "__main__":
    main()
