#!/usr/bin/env python3
"""Step 4: yearly paper scores (Jensen and PubTator) for each of our genes.

Run from the dark_proteins folder, after step 3:
    python3 scripts/04_gene_year_scores.py

Inputs
    data/processed/r03_awards/gene_summary.tsv               our 388 genes (step 3)
    data/raw/hgnc/hgnc_complete_set.txt                      HGNC ID -> NCBI Gene ID (entrez_id)
    data/raw/jensenlab/protein_counts_2026-10.tsv            Jensen PubMed score: ENSP, year, score
                                                             (download.jensenlab.org, file dated 2026-10-26)
    data/raw/jensenlab/pubtator_counts.tsv                   PubTator score: NCBI Gene ID, year, score
    data/raw/jensenlab/9606.protein.aliases.v12.0.txt.gz     STRING: ENSP -> HGNC ID

How genes are matched
    PubTator: HGNC entrez_id = NCBI Gene ID. Human gene IDs only, so human papers only.
    Jensen:   STRING alias lines "HGNC:..." link each ENSP to a gene.
              CILK1 and ERN2 have no HGNC line in STRING; they are matched by name (MANUAL_ENSP).
              Jensen matches gene names in text, so it can include papers on other species.

Scores are copied exactly as the files give them (fractional: a paper's point is split across
the genes it mentions). Nothing is smoothed or changed.

Output
    data/processed/attention/gene_year_scores.tsv
        one row per gene x year, 1990-2026
        jensen_score / pubtator_score: that year's score; 0 = in the file but no score that year;
                                       blank = gene not in that file at all
        note: which values not to trust that year
              PubTator: totals jump in 2020-21 and collapse from 2022 -> use through 2021 only
              Jensen:   2026 is a partial year (file built Oct 2026)
"""
import csv
import gzip
from collections import defaultdict
from pathlib import Path

GENES = Path("data/processed/r03_awards/gene_summary.tsv")
HGNC = Path("data/raw/hgnc/hgnc_complete_set.txt")
JENSEN = Path("data/raw/jensenlab/protein_counts_2026-10.tsv")
PUBTATOR = Path("data/raw/jensenlab/pubtator_counts.tsv")
STRING = Path("data/raw/jensenlab/9606.protein.aliases.v12.0.txt.gz")
OUT = Path("data/processed/attention/gene_year_scores.tsv")
YEARS = range(1990, 2027)
PUBTATOR_LAST_GOOD = 2021
JENSEN_PARTIAL = {2026}


def note_for(year):
    parts = []
    if year > PUBTATOR_LAST_GOOD:
        parts.append("pubtator unreliable")
    if year in JENSEN_PARTIAL:
        parts.append("jensen partial year")
    return "; ".join(parts)
# No "HGNC:" line in STRING for these; their STRING protein was found by gene name (checked by hand).
MANUAL_ENSP = {"CILK1": "ENSP00000349458", "ERN2": "ENSP00000256797"}


def read_tsv(path):
    if not path.is_file():
        raise SystemExit(f"Missing {path}")
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def scores_by_year(path, wanted):
    """{id: {year: score}} for the ids in `wanted` only."""
    out = defaultdict(lambda: defaultdict(float))
    with path.open(encoding="utf-8") as f:
        for line in f:
            i, y, s = line.rstrip("\n").split("\t")
            if i in wanted:
                out[i][int(y)] += float(s)
    return out


def main():
    genes = read_tsv(GENES)
    entrez = {r["hgnc_id"]: r.get("entrez_id", "") for r in read_tsv(HGNC)}

    # Jensen: ENSP -> HGNC from STRING, plus the two manual matches
    hgnc2ensp = defaultdict(set)
    with gzip.open(STRING, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p, alias, _ = line.rstrip("\n").split("\t")
            if alias.startswith("HGNC:"):
                hgnc2ensp[alias].add(p.replace("9606.", ""))
    for g in genes:
        if g["gene"] in MANUAL_ENSP and not hgnc2ensp.get(g["hgnc_id"]):
            hgnc2ensp[g["hgnc_id"]].add(MANUAL_ENSP[g["gene"]])

    want_ensp = {p for g in genes for p in hgnc2ensp.get(g["hgnc_id"], ())}
    want_entrez = {entrez[g["hgnc_id"]] for g in genes if entrez.get(g["hgnc_id"])}
    print("Reading Jensen file ...")
    jensen = scores_by_year(JENSEN, want_ensp)
    print("Reading PubTator file (12 million rows) ...")
    pubtator = scores_by_year(PUBTATOR, want_entrez)

    rows, no_jensen, no_pubtator, two_proteins = [], [], [], []
    for g in genes:
        hid = g["hgnc_id"]
        ensps = [p for p in sorted(hgnc2ensp.get(hid, ())) if p in jensen]
        if len(ensps) > 1:
            two_proteins.append(f"{g['gene']}: {', '.join(ensps)} (scores added)")
        eid = entrez.get(hid, "")
        in_pt = eid in pubtator
        if not ensps:
            no_jensen.append(g["gene"])
        if not in_pt:
            no_pubtator.append(g["gene"])
        for y in YEARS:
            js = sum(jensen[p].get(y, 0.0) for p in ensps)
            ps = pubtator[eid].get(y, 0.0) if in_pt else 0.0
            rows.append({"gene": g["gene"], "hgnc_id": hid, "year": y,
                         "jensen_score": f"{js:.6f}" if ensps else "",
                         "pubtator_score": f"{ps:.6f}" if in_pt else "",
                         "note": note_for(y)})

    print(f"\nGenes: {len(genes)}")
    print(f"With a Jensen score: {len(genes) - len(no_jensen)} | without: {no_jensen or 'none'}")
    print(f"With a PubTator score: {len(genes) - len(no_pubtator)} | without: {no_pubtator or 'none'}")
    print(f"Genes with more than one Jensen protein: {two_proteins or 'none'}")

    print("\nAll our genes added together, per year (to see where the files run out):")
    print("   year   jensen  pubtator")
    tot = defaultdict(lambda: [0.0, 0.0])
    for r in rows:
        tot[r["year"]][0] += float(r["jensen_score"] or 0)
        tot[r["year"]][1] += float(r["pubtator_score"] or 0)
    for y in range(2015, 2027):
        print(f"   {y}  {tot[y][0]:7.1f}  {tot[y][1]:8.1f}" + (f"   <- {note_for(y)}" if note_for(y) else ""))

    print("\nCheck three genes, 2016-2026 (jensen / pubtator):")
    for name in ("ADGRB3", "GPR68", "PRKCQ"):
        rr = [r for r in rows if r["gene"] == name and r["year"] >= 2016]
        if rr:
            print(f"   {name}: " + "  ".join(f"{r['year']}: {float(r['jensen_score'] or 0):.1f}/"
                                            f"{float(r['pubtator_score'] or 0):.1f}" for r in rr))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"\nSaved: {OUT} ({len(rows):,} rows = {len(genes)} genes x {len(YEARS)} years)")


if __name__ == "__main__":
    main()
