#!/usr/bin/env python3
"""Step 8: ChEMBL compound activities for our genes, by year, with IDG's potency cutoffs and the R03 tie.

Run from the dark_proteins folder, after steps 3 and 5:
    python3 scripts/08_gene_chemistry.py [--refresh]

Source: ChEMBL web services (https://www.ebi.ac.uk/chembl/api/data), release printed at the start.
Human targets only: ChEMBL targets are found by each gene's UniProt ID (from HGNC).

Inputs
    data/processed/r03_awards/gene_summary.tsv            our 388 genes and their family
    data/raw/hgnc/hgnc_complete_set.txt                   HGNC ID -> UniProt ID(s)
    data/processed/r03_awards/foa_genes_with_awards.tsv   which genes each R03 funded
    data/raw/r03_awards/reporter_r03_projects.json        R03 start dates
    data/raw/r03_awards/reporter_r03_publications.json    papers citing each R03

Labels kept on every activity (nothing filtered):
    target_type     SINGLE PROTEIN (main counts) | PROTEIN COMPLEX | PROTEIN FAMILY | ... (kept, labelled)
    assay_type      B binding | F functional | A ADME | T toxicity | P physchem | U unclassified
    valid           1 unless ChEMBL flags the value (data_validity_comment) or marks a potential duplicate
    potent          pChEMBL at or above IDG's Tchem cutoff for the family:
                    kinase >= 7.5 (30 nM) | GPCR >= 7.0 (100 nM) | ion channel >= 5.0 (10 uM)
    document        ChEMBL document, its year, PMID, type (publication or deposited dataset)
    r03             the document's paper cites an R03 that funded this gene, in or after its start year

Saved along the way (reruns are fast; --refresh downloads again)
    data/raw/chembl/targets.json          targets per UniProt ID
    data/raw/chembl/activities.jsonl      one line per target: all its activities
    data/raw/chembl/documents.json        document details (PMID, type, year)

Outputs (data/processed/chemistry/)
    gene_activities.tsv       one row per gene x activity
    gene_year_chemistry.tsv   per gene x year (single-protein targets, valid activities):
        n_activities, n_compounds, n_new_compounds, n_potent_compounds, n_new_potent_compounds,
        n_documents, n_r03_activities, n_r03_documents, tchem_by_this_year (1 once a potent compound exists)
"""
import argparse
import csv
import json
import re
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

API = "https://www.ebi.ac.uk/chembl/api/data"
GENES = Path("data/processed/r03_awards/gene_summary.tsv")
HGNC = Path("data/raw/hgnc/hgnc_complete_set.txt")
AWARDS = Path("data/processed/r03_awards/foa_genes_with_awards.tsv")
R03_PROJECTS = Path("data/raw/r03_awards/reporter_r03_projects.json")
R03_PUBS = Path("data/raw/r03_awards/reporter_r03_publications.json")
RAW = Path("data/raw/chembl")
OUT_DIR = Path("data/processed/chemistry")
CUTOFF = {"Kinase": 7.5, "GPCR": 7.0, "Ion channel": 5.0}       # IDG Tchem potency cutoffs (pChEMBL)
ONLY = ("activity_id,molecule_chembl_id,parent_molecule_chembl_id,standard_type,standard_relation,"
        "standard_value,standard_units,pchembl_value,document_chembl_id,document_year,assay_type,"
        "data_validity_comment,potential_duplicate,src_id,target_chembl_id")


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


def get(path):
    for attempt in range(5):
        try:
            req = urllib.request.Request(API + path, headers={"Accept": "application/json"})
            data = json.load(urllib.request.urlopen(req, timeout=120))
            time.sleep(0.2)
            return data
        except Exception as e:
            print(f"   retry {attempt + 1}: {path[:80]} ({e})")
            time.sleep(5 * (attempt + 1))
    raise SystemExit("ChEMBL did not answer after 5 tries; rerun later (finished parts are kept).")


def all_pages(path, key):
    out, nxt = [], path
    while nxt:
        d = get(nxt)
        out += d.get(key, [])
        nxt = d.get("page_meta", {}).get("next")
        if nxt:
            nxt = nxt.replace("/chembl/api/data", "")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)

    genes = read_tsv(GENES)
    ours = {r["hgnc_id"]: r["gene"] for r in genes}
    fam = {r["hgnc_id"]: r.get("family", "") for r in genes}
    uni = {r["hgnc_id"]: [u for u in (r.get("uniprot_ids") or "").replace('"', "").split("|") if u]
           for r in read_tsv(HGNC) if r["hgnc_id"] in ours}
    release = get("/status.json").get("chembl_db_version")
    print(f"ChEMBL release: {release} | genes with a UniProt ID: {sum(1 for h in ours if uni.get(h))} of {len(ours)}")

    # targets per UniProt ID
    tpath = RAW / "targets.json"
    targets = json.loads(tpath.read_text()) if tpath.is_file() and not a.refresh else {}
    todo = sorted({u for us in uni.values() for u in us} - set(targets))
    for i, u in enumerate(todo, 1):
        targets[u] = [{"id": t["target_chembl_id"], "type": t["target_type"], "organism": t["organism"],
                       "name": t["pref_name"]}
                      for t in all_pages(f"/target.json?target_components__accession={u}&limit=100", "targets")]
        if i % 50 == 0:
            print(f"   targets: {i} of {len(todo)} UniProt IDs")
            tpath.write_text(json.dumps(targets))
    tpath.write_text(json.dumps(targets))
    gene_targets = defaultdict(dict)                       # hgnc -> target id -> info (human only)
    for h, us in uni.items():
        for u in us:
            for t in targets.get(u, []):
                if t["organism"] == "Homo sapiens":
                    gene_targets[h][t["id"]] = t
    all_t = {tid for d in gene_targets.values() for tid in d}
    print(f"Human ChEMBL targets linked to our genes: {len(all_t)} | genes with any target: {len(gene_targets)}")
    print("   by target type:", dict(Counter(t["type"] for d in gene_targets.values() for t in d.values())))

    # activities per target
    apath = RAW / "activities.jsonl"
    acts = {}
    if apath.is_file() and not a.refresh:
        for line in apath.open(encoding="utf-8"):
            d = json.loads(line)
            acts[d["target"]] = d["activities"]
    elif apath.is_file():
        apath.unlink()
    todo = sorted(all_t - set(acts))
    print(f"Activities: {len(acts)} targets cached, {len(todo)} to download")
    with apath.open("a", encoding="utf-8") as f:
        for i, tid in enumerate(todo, 1):
            acts[tid] = all_pages(f"/activity.json?target_chembl_id={tid}&only={ONLY}&limit=1000", "activities")
            f.write(json.dumps({"target": tid, "activities": acts[tid]}) + "\n")
            f.flush()
            if i % 25 == 0:
                print(f"   {i} of {len(todo)} targets")

    # documents (PMID, type, year)
    dpath = RAW / "documents.json"
    docs = json.loads(dpath.read_text()) if dpath.is_file() and not a.refresh else {}
    need = sorted({x["document_chembl_id"] for tid in all_t for x in acts.get(tid, []) if x.get("document_chembl_id")} - set(docs))
    print(f"Documents: {len(docs)} cached, {len(need)} to download")
    for b in range(0, len(need), 50):
        ids = ",".join(need[b:b + 50])
        for d in all_pages(f"/document.json?document_chembl_id__in={ids}&only=document_chembl_id,pubmed_id,doc_type,year,journal,src_id&limit=50", "documents"):
            docs[d["document_chembl_id"]] = {"pmid": str(d.get("pubmed_id") or ""), "type": d.get("doc_type") or "",
                                             "year": d.get("year"), "journal": d.get("journal") or ""}
    dpath.write_text(json.dumps(docs))

    # R03 tie
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
    for h, tmap in gene_targets.items():
        cut = CUTOFF.get(fam.get(h, ""))
        for tid, t in tmap.items():
            for x in acts.get(tid, []):
                doc = docs.get(x.get("document_chembl_id") or "", {})
                year = str(x.get("document_year") or doc.get("year") or "")
                p = x.get("pchembl_value")
                pv = float(p) if p not in (None, "") else None
                valid = int(not x.get("data_validity_comment") and not x.get("potential_duplicate"))
                pmid = doc.get("pmid", "")
                own = [c for c in sorted(r03.get(pmid, ())) if h in award_genes.get(c, ())
                       and year and year >= award_start.get(c, "9999")] if pmid else []
                rows.append({"gene": ours[h], "hgnc_id": h, "family": fam.get(h, ""), "target_chembl_id": tid,
                             "target_type": t["type"], "activity_id": x.get("activity_id"),
                             "compound": x.get("parent_molecule_chembl_id") or x.get("molecule_chembl_id"),
                             "standard_type": x.get("standard_type"), "standard_relation": x.get("standard_relation"),
                             "standard_value": x.get("standard_value"), "standard_units": x.get("standard_units"),
                             "pchembl_value": p or "", "assay_type": x.get("assay_type"), "valid": valid,
                             "potent": int(bool(valid and pv is not None and cut is not None and pv >= cut)),
                             "document_chembl_id": x.get("document_chembl_id"), "document_type": doc.get("type", ""),
                             "journal": doc.get("journal", ""), "pmid": pmid, "year": year,
                             "r03_for_this_gene": ";".join(own)})

    # per gene x year (single protein, valid)
    main_rows = [r for r in rows if r["target_type"] == "SINGLE PROTEIN" and r["valid"] and r["year"]]
    first, first_potent, tchem = {}, {}, {}
    for r in main_rows:
        k = (r["hgnc_id"], r["compound"])
        first[k] = min(first.get(k, "9999"), r["year"])
        if r["potent"]:
            first_potent[k] = min(first_potent.get(k, "9999"), r["year"])
            tchem[r["hgnc_id"]] = min(tchem.get(r["hgnc_id"], "9999"), r["year"])
    agg = defaultdict(lambda: defaultdict(set))
    for r in main_rows:
        s = agg[(r["hgnc_id"], r["gene"], r["year"])]
        s["acts"].add(r["activity_id"])
        s["cmp"].add(r["compound"])
        s["docs"].add(r["document_chembl_id"])
        k = (r["hgnc_id"], r["compound"])
        if first[k] == r["year"]:
            s["new"].add(r["compound"])
        if r["potent"]:
            s["potent"].add(r["compound"])
            if first_potent.get(k) == r["year"]:
                s["new_potent"].add(r["compound"])
        if r["r03_for_this_gene"]:
            s["r03a"].add(r["activity_id"])
            s["r03d"].add(r["document_chembl_id"])
    yrows = []
    for (h, g, y), s in sorted(agg.items(), key=lambda kv: (kv[0][1], kv[0][2])):
        yrows.append({"gene": g, "hgnc_id": h, "year": y, "n_activities": len(s["acts"]), "n_compounds": len(s["cmp"]),
                      "n_new_compounds": len(s["new"]), "n_potent_compounds": len(s["potent"]),
                      "n_new_potent_compounds": len(s["new_potent"]), "n_documents": len(s["docs"]),
                      "n_r03_activities": len(s["r03a"]), "n_r03_documents": len(s["r03d"]),
                      "tchem_by_this_year": int(h in tchem and y >= tchem[h])})

    # report
    print(f"\nGene x activity rows: {len(rows):,} | single-protein, valid, dated: {len(main_rows):,}")
    print(f"Genes with any single-protein activity: {len({r['hgnc_id'] for r in main_rows})} of {len(ours)}")
    print("By family: genes with any activity / genes reaching IDG's potency cutoff (Tchem):")
    for f_ in ("Kinase", "GPCR", "Ion channel"):
        gs = {h for h in ours if fam.get(h) == f_}
        print(f"   {f_:<12} {len(gs & {r['hgnc_id'] for r in main_rows}):>4} / {len(gs & set(tchem)):>4} of {len(gs)}")
    print("Year genes first reached the cutoff:", dict(sorted(Counter(tchem.values()).items())))
    print("\nActivities per year, all our genes (single protein, valid):")
    c = Counter(r["year"] for r in main_rows)
    print("   " + "  ".join(f"{y}:{c.get(str(y), 0)}" for y in range(2010, 2027)))
    tied = [r for r in rows if r["r03_for_this_gene"]]
    print(f"\nR03-tied activities: {len(tied)} | documents: {len({r['document_chembl_id'] for r in tied})} | "
          f"genes: {sorted({r['gene'] for r in tied})}")
    print("\nCheck three genes (year: activities / new compounds / potent):")
    for g in ("GPR68", "PRKCQ", "ADGRB3"):
        v = [x for x in yrows if x["gene"] == g and x["year"] >= "2012"]
        print(f"   {g}: " + ("  ".join(f"{x['year']}: {x['n_activities']}/{x['n_new_compounds']}/{x['n_potent_compounds']}" for x in v) or "none"))

    cols = ["gene", "hgnc_id", "family", "target_chembl_id", "target_type", "activity_id", "compound", "standard_type",
            "standard_relation", "standard_value", "standard_units", "pchembl_value", "assay_type", "valid", "potent",
            "document_chembl_id", "document_type", "journal", "pmid", "year", "r03_for_this_gene"]
    write_tsv(OUT_DIR / "gene_activities.tsv", rows, cols)
    write_tsv(OUT_DIR / "gene_year_chemistry.tsv", yrows, list(yrows[0]) if yrows else ["gene"])
    (RAW / "release.txt").write_text(f"{release}\n")
    print(f"\nSaved: {OUT_DIR / 'gene_activities.tsv'} ({len(rows):,} rows)")
    print(f"Saved: {OUT_DIR / 'gene_year_chemistry.tsv'} ({len(yrows):,} rows)")


if __name__ == "__main__":
    main()
