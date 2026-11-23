#!/usr/bin/env python3
"""Step 11: gene-disease claims (GenCC) for our genes, by year, with the R03 tie.

Run from the dark_proteins folder, after steps 3 and 5:
    python3 scripts/11_gene_disease.py --email you@colorado.edu

Inputs
    data/raw/gencc/gencc-submissions.tsv                  GenCC submissions export
    data/processed/r03_awards/gene_summary.tsv            our 388 genes
    data/processed/r03_awards/foa_genes_with_awards.tsv   which genes each R03 funded
    data/raw/r03_awards/reporter_r03_projects.json        R03 start dates
    data/raw/r03_awards/reporter_r03_publications.json    papers citing each R03
    data/raw/pubmed/pmid_years.tsv                        PMID -> year (new PMIDs are fetched and added)

Each GenCC submission = one submitter's grade for one gene-disease link.
    grade         Definitive | Strong | Moderate | Supportive | Limited | Disputed | Refuted | No Known ...
    moderate_plus 1 for Definitive, Strong or Moderate
    submitted_year  when the submitter made the claim (official recognition; often years after discovery)
    first_paper_year  year of the earliest paper the submission cites (closer to discovery; blank if none)
    r03           a cited paper cites an R03 that funded this gene, in or after its start year

Outputs (data/processed/disease/)
    gene_disease_submissions.tsv  one row per submission for our genes
    gene_disease_links.tsv        one row per gene x disease: best grade, first submitted year,
                                  first Moderate+ year, earliest paper year, submitters
    gene_year_disease.tsv         per gene x year:
        n_new_links_submitted       links submitted for the first time that year
        n_new_moderate_plus         links reaching Moderate or stronger for the first time that year
        n_new_links_by_paper_year   links whose earliest supporting paper is from that year
        n_r03_submissions
"""
import argparse
import csv
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

GENCC = Path("data/raw/gencc/gencc-submissions.tsv")
GENES = Path("data/processed/r03_awards/gene_summary.tsv")
AWARDS = Path("data/processed/r03_awards/foa_genes_with_awards.tsv")
R03_PROJECTS = Path("data/raw/r03_awards/reporter_r03_projects.json")
R03_PUBS = Path("data/raw/r03_awards/reporter_r03_publications.json")
YEARS = Path("data/raw/pubmed/pmid_years.tsv")
OUT_DIR = Path("data/processed/disease")
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
RANK = ["Definitive", "Strong", "Moderate", "Supportive", "Limited", "Disputed", "Refuted", "No Known"]
MODERATE_PLUS = {"Definitive", "Strong", "Moderate"}


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


def grade(title):
    for g in RANK:
        if (title or "").lower().startswith(g.lower()):
            return g
    return title or ""


def fetch_years(pmids, email):
    known = {r["pmid"]: r["year"] for r in read_tsv(YEARS)} if YEARS.is_file() else {}
    todo = sorted(set(pmids) - set(known), key=int)
    print(f"Paper years: {len(known):,} cached, {len(todo):,} to fetch from PubMed")
    with YEARS.open("a", encoding="utf-8") as out:
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
    csv.field_size_limit(sys.maxsize)

    ours = {r["hgnc_id"]: r["gene"] for r in read_tsv(GENES)}
    with GENCC.open(encoding="utf-8", newline="") as f:
        allrows = list(csv.DictReader(f, delimiter="\t"))
    subs = [r for r in allrows if r["gene_curie"] in ours]
    print(f"GenCC submissions: {len(allrows):,} | for our genes: {len(subs):,} | "
          f"genes with any: {len({r['gene_curie'] for r in subs})} of {len(ours)}")

    pm = lambda s: re.findall(r"\d{5,9}", s or "")
    years = fetch_years({p for r in subs for p in pm(r["submitted_as_pmids"])}, a.email)

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
    for r in subs:
        h = r["gene_curie"]
        pmids = pm(r["submitted_as_pmids"])
        pyears = sorted(years.get(p) for p in pmids if years.get(p))
        own = sorted({c for p in pmids for c in r03.get(p, ()) if h in award_genes.get(c, ())
                      and years.get(p, "") >= award_start.get(c, "9999")})
        g = grade(r["classification_title"])
        rows.append({"gene": ours[h], "hgnc_id": h, "disease_curie": r["disease_curie"], "disease_title": r["disease_title"],
                     "grade": g, "moderate_plus": int(g in MODERATE_PLUS), "moi": r["moi_title"],
                     "submitter": r["submitter_title"], "submitted_year": (r["submitted_as_date"] or "")[:4],
                     "pmids": ";".join(pmids), "first_paper_year": pyears[0] if pyears else "",
                     "r03_for_this_gene": ";".join(own)})

    links = defaultdict(list)
    for r in rows:
        links[(r["hgnc_id"], r["disease_curie"])].append(r)
    lrows = []
    for (h, d), rs in links.items():
        sub_y = [r["submitted_year"] for r in rs if r["submitted_year"]]
        mod_y = [r["submitted_year"] for r in rs if r["moderate_plus"] and r["submitted_year"]]
        pap_y = [r["first_paper_year"] for r in rs if r["first_paper_year"]]
        best = min((r["grade"] for r in rs), key=lambda g: RANK.index(g) if g in RANK else 99)
        lrows.append({"gene": ours[h], "hgnc_id": h, "disease_curie": d, "disease_title": rs[0]["disease_title"],
                      "best_grade": best, "n_submissions": len(rs),
                      "first_submitted_year": min(sub_y) if sub_y else "",
                      "first_moderate_plus_year": min(mod_y) if mod_y else "",
                      "earliest_paper_year": min(pap_y) if pap_y else "",
                      "submitters": "; ".join(sorted({r["submitter"] for r in rs})),
                      "r03_for_this_gene": ";".join(sorted({c for r in rs for c in r["r03_for_this_gene"].split(";") if c}))})
    lrows.sort(key=lambda r: (r["gene"], r["disease_title"]))

    agg = defaultdict(Counter)
    for l in lrows:
        if l["first_submitted_year"]:
            agg[(l["hgnc_id"], l["gene"], l["first_submitted_year"])]["n_new_links_submitted"] += 1
        if l["first_moderate_plus_year"]:
            agg[(l["hgnc_id"], l["gene"], l["first_moderate_plus_year"])]["n_new_moderate_plus"] += 1
        if l["earliest_paper_year"]:
            agg[(l["hgnc_id"], l["gene"], l["earliest_paper_year"])]["n_new_links_by_paper_year"] += 1
    for r in rows:
        if r["r03_for_this_gene"] and r["submitted_year"]:
            agg[(r["hgnc_id"], r["gene"], r["submitted_year"])]["n_r03_submissions"] += 1
    yrows = [{"gene": g, "hgnc_id": h, "year": y, "n_new_links_submitted": c["n_new_links_submitted"],
              "n_new_moderate_plus": c["n_new_moderate_plus"], "n_new_links_by_paper_year": c["n_new_links_by_paper_year"],
              "n_r03_submissions": c["n_r03_submissions"]}
             for (h, g, y), c in sorted(agg.items(), key=lambda kv: (kv[0][1], kv[0][2]))]

    # report
    print("Grades (submissions):", dict(Counter(r["grade"] for r in rows).most_common()))
    print(f"Gene-disease links: {len(lrows):,} | genes with any link: {len({l['hgnc_id'] for l in lrows})} | "
          f"genes with a Moderate+ link: {len({l['hgnc_id'] for l in lrows if l['first_moderate_plus_year']})}")
    print("Top submitters:", dict(Counter(r["submitter"] for r in rows).most_common(6)))
    print("\nPer year, all our genes: new links submitted / new Moderate+ / new links by earliest paper year")
    for y in range(2005, 2027):
        v = [x for x in yrows if x["year"] == str(y)]
        print(f"   {y}  {sum(x['n_new_links_submitted'] for x in v):>4}  {sum(x['n_new_moderate_plus'] for x in v):>4}  "
              f"{sum(x['n_new_links_by_paper_year'] for x in v):>4}")
    first_mod = {}
    for l in lrows:
        if l["first_moderate_plus_year"]:
            first_mod[l["hgnc_id"]] = min(first_mod.get(l["hgnc_id"], "9999"), l["first_moderate_plus_year"])
    print("Year genes first got a Moderate+ link:", dict(sorted(Counter(first_mod.values()).items())))
    tied = [r for r in rows if r["r03_for_this_gene"]]
    print(f"\nR03-tied submissions: {len(tied)} | genes: {sorted({r['gene'] for r in tied})}")
    print("\nCheck three genes (disease | best grade | first submitted | earliest paper):")
    for g in ("GPR68", "PRKCQ", "ADGRB3"):
        v = [l for l in lrows if l["gene"] == g]
        print(f"   {g}: " + (" || ".join(f"{l['disease_title'][:40]} | {l['best_grade']} | {l['first_submitted_year']} | "
                                     f"{l['earliest_paper_year']}" for l in v) or "none"))

    cols = ["gene", "hgnc_id", "disease_curie", "disease_title", "grade", "moderate_plus", "moi", "submitter",
            "submitted_year", "pmids", "first_paper_year", "r03_for_this_gene"]
    write_tsv(OUT_DIR / "gene_disease_submissions.tsv", rows, cols)
    write_tsv(OUT_DIR / "gene_disease_links.tsv", lrows, list(lrows[0]) if lrows else ["gene"])
    write_tsv(OUT_DIR / "gene_year_disease.tsv", yrows, list(yrows[0]) if yrows else ["gene"])
    print(f"\nSaved: {OUT_DIR / 'gene_disease_submissions.tsv'} ({len(rows):,} rows)")
    print(f"Saved: {OUT_DIR / 'gene_disease_links.tsv'} ({len(lrows):,} rows)")
    print(f"Saved: {OUT_DIR / 'gene_year_disease.tsv'} ({len(yrows):,} rows)")


if __name__ == "__main__":
    main()
