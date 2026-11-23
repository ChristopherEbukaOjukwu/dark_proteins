#!/usr/bin/env python3
"""Step 10: experimentally solved PDB structures containing our proteins, by year, with the R03 tie.

Run from the dark_proteins folder, after steps 3 and 5:
    python3 scripts/10_gene_structures.py [--refresh]

Source: RCSB PDB search service (entries containing a UniProt protein; experimental structures only,
no computed models) and the RCSB data service (entry details).

Inputs
    data/processed/r03_awards/gene_summary.tsv            our 388 genes
    data/raw/hgnc/hgnc_complete_set.txt                   HGNC ID -> UniProt ID(s)
    data/processed/r03_awards/foa_genes_with_awards.tsv   which genes each R03 funded
    data/raw/r03_awards/reporter_r03_projects.json        R03 start dates
    data/raw/r03_awards/reporter_r03_publications.json    papers citing each R03

Labels kept on every entry:
    year            year the structure was released (made public); deposit year also kept
    method          X-ray | cryo-EM | NMR | other
    resolution      in angstroms (lower = sharper)
    ligand          1 if a bound small molecule is present that is not a common additive
                    (ions, water, buffers, detergents, sugars from expression, lipids used for packing)
    pmid            the entry's main paper (blank if unpublished)
    r03             main paper cites an R03 that funded this gene, in or after its start year
    An entry counts even if it contains only part of the protein (e.g. one domain).

Saved along the way (reruns are fast; --refresh downloads again)
    data/raw/pdb/search.json     UniProt ID -> PDB entry IDs
    data/raw/pdb/entries.json    PDB entry -> details

Outputs (data/processed/structure/)
    gene_pdb_entries.tsv        one row per gene x PDB entry
    gene_year_structure.tsv     per gene x release year: n_entries, n_xray, n_cryoem, n_nmr, n_with_ligand,
                                first_structure_year (the gene's first structure ever), n_r03_entries
"""
import argparse
import csv
import json
import re
import time
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"
GENES = Path("data/processed/r03_awards/gene_summary.tsv")
HGNC = Path("data/raw/hgnc/hgnc_complete_set.txt")
AWARDS = Path("data/processed/r03_awards/foa_genes_with_awards.tsv")
R03_PROJECTS = Path("data/raw/r03_awards/reporter_r03_projects.json")
R03_PUBS = Path("data/raw/r03_awards/reporter_r03_publications.json")
RAW = Path("data/raw/pdb")
OUT_DIR = Path("data/processed/structure")
# Common crystallisation/cryo-EM additives, ions and packing lipids: not counted as a ligand.
ADDITIVES = set("""HOH DOD NA CL K MG CA ZN MN FE CO NI CU CD IOD BR SO4 PO4 PO3 NO3 ACT ACY EDO GOL PEG PGE PG4 1PE
P6G MPD DMS TRS EPE MES BME DTT FMT IPA EOH MOH CIT TLA MLI NAG NDG BMA MAN FUC GAL GLC BGC SUC LMT LMN DDM
OLC OLA PLM MYR STE CLR Y01 PCW POV PEE LDA BOG C8E HTG UNX UNL CO3 BCT AZI SCN OXL IMD""".split())


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


def post(url, body):
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, data=json.dumps(body).encode(), headers={"Content-Type": "application/json"})
            resp = urllib.request.urlopen(req, timeout=120)
            raw = resp.read()
            time.sleep(0.2)
            return json.loads(raw) if raw else {}          # RCSB search returns an empty body for no hits
        except Exception as e:
            print(f"   retry {attempt + 1}: {e}")
            time.sleep(5 * (attempt + 1))
    raise SystemExit("RCSB did not answer after 5 tries; rerun later (finished parts are kept).")


def search_entries(uniprot):
    ids, start = [], 0
    while True:
        q = {"query": {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match", "value": uniprot}},
             "return_type": "entry", "request_options": {"paginate": {"start": start, "rows": 1000}}}
        r = post(SEARCH, q)
        ids += [x["identifier"] for x in r.get("result_set", [])]
        start += 1000
        if start >= r.get("total_count", 0):
            return ids


def entry_details(ids):
    gql = ('{ entries(entry_ids: %s) { rcsb_id struct { title } rcsb_accession_info { initial_release_date deposit_date } '
           'exptl { method } rcsb_entry_info { resolution_combined nonpolymer_bound_components } '
           'rcsb_primary_citation { pdbx_database_id_PubMed year } } }') % json.dumps(ids)
    return post(GRAPHQL, {"query": gql}).get("data", {}).get("entries") or []


def method_label(methods):
    m = " ".join(methods).upper()
    return ("cryo-EM" if "ELECTRON MICROSCOPY" in m else "X-ray" if "X-RAY" in m
            else "NMR" if "NMR" in m else "other")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)

    ours = {r["hgnc_id"]: r["gene"] for r in read_tsv(GENES)}
    uni = {r["hgnc_id"]: [u for u in (r.get("uniprot_ids") or "").replace('"', "").split("|") if u]
           for r in read_tsv(HGNC) if r["hgnc_id"] in ours}

    spath, epath = RAW / "search.json", RAW / "entries.json"
    found = json.loads(spath.read_text()) if spath.is_file() and not a.refresh else {}
    todo = sorted({u for us in uni.values() for u in us} - set(found))
    print(f"UniProt IDs: {sum(len(v) for v in uni.values())} | already searched: {len(found)} | to search: {len(todo)}")
    for i, u in enumerate(todo, 1):
        found[u] = search_entries(u)
        if i % 50 == 0:
            print(f"   {i} of {len(todo)}")
            spath.write_text(json.dumps(found))
    spath.write_text(json.dumps(found))

    entries = json.loads(epath.read_text()) if epath.is_file() and not a.refresh else {}
    need = sorted({e for ids in found.values() for e in ids} - set(entries))
    print(f"PDB entries: {len(entries)} cached, {len(need)} to download")
    for b in range(0, len(need), 100):
        for e in entry_details(need[b:b + 100]):
            entries[e["rcsb_id"]] = e
    epath.write_text(json.dumps(entries))

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
    for h, us in uni.items():
        ids = sorted({e for u in us for e in found.get(u, [])})
        for eid in ids:
            e = entries.get(eid)
            if not e:
                continue
            acc = e.get("rcsb_accession_info") or {}
            year = (acc.get("initial_release_date") or "")[:4]
            info = e.get("rcsb_entry_info") or {}
            comps = info.get("nonpolymer_bound_components") or []
            cit = e.get("rcsb_primary_citation") or {}
            pmid = str(cit.get("pdbx_database_id_PubMed") or "")
            res = info.get("resolution_combined") or []
            own = [c for c in sorted(r03.get(pmid, ())) if h in award_genes.get(c, ())
                   and year and year >= award_start.get(c, "9999")] if pmid else []
            rows.append({"gene": ours[h], "hgnc_id": h, "pdb_id": eid, "title": ((e.get("struct") or {}).get("title") or "")[:150],
                         "year": year, "deposit_year": (acc.get("deposit_date") or "")[:4],
                         "method": method_label([x.get("method", "") for x in e.get("exptl") or []]),
                         "resolution": res[0] if res else "", "bound_components": ";".join(comps),
                         "ligand": int(any(c not in ADDITIVES for c in comps)), "pmid": pmid,
                         "r03_for_this_gene": ";".join(own)})

    first = {}
    for r in rows:
        if r["year"]:
            first[r["hgnc_id"]] = min(first.get(r["hgnc_id"], "9999"), r["year"])
    agg = defaultdict(Counter)
    for r in rows:
        if not r["year"]:
            continue
        c = agg[(r["hgnc_id"], r["gene"], r["year"])]
        c["n_entries"] += 1
        c[{"X-ray": "n_xray", "cryo-EM": "n_cryoem", "NMR": "n_nmr"}.get(r["method"], "n_other")] += 1
        c["n_with_ligand"] += r["ligand"]
        c["n_r03_entries"] += bool(r["r03_for_this_gene"])
    yrows = [{"gene": g, "hgnc_id": h, "year": y, "n_entries": c["n_entries"], "n_xray": c["n_xray"],
              "n_cryoem": c["n_cryoem"], "n_nmr": c["n_nmr"], "n_other": c["n_other"],
              "n_with_ligand": c["n_with_ligand"], "first_structure_year": first.get(h, ""),
              "n_r03_entries": c["n_r03_entries"]}
             for (h, g, y), c in sorted(agg.items(), key=lambda kv: (kv[0][1], kv[0][2]))]

    # report
    print(f"\nGene x entry rows: {len(rows):,} | distinct entries: {len({r['pdb_id'] for r in rows}):,} | "
          f"genes with any structure: {len(first)} of {len(ours)}")
    print("By method:", dict(Counter(r["method"] for r in rows).most_common()))
    print("First-ever structure per gene, by year:", dict(sorted(Counter(first.values()).items())))
    print("\nNew entries per release year, all our genes (cryo-EM share, with ligand):")
    for y in range(2012, 2027):
        v = [r for r in rows if r["year"] == str(y)]
        print(f"   {y}  entries {len(v):>5}  cryo-EM {sum(r['method'] == 'cryo-EM' for r in v):>5}  "
              f"with ligand {sum(r['ligand'] for r in v):>5}  genes {len({r['hgnc_id'] for r in v}):>4}")
    tied = [r for r in rows if r["r03_for_this_gene"]]
    print(f"\nR03-tied entries: {len(tied)} | papers: {len({r['pmid'] for r in tied})} | genes: {sorted({r['gene'] for r in tied})}")
    print("\nCheck three genes (release year: entries):")
    for g in ("GPR68", "PRKCQ", "ADGRB3"):
        v = [x for x in yrows if x["gene"] == g]
        print(f"   {g}: " + ("  ".join(f"{x['year']}: {x['n_entries']}" for x in v) or "none"))

    cols = ["gene", "hgnc_id", "pdb_id", "title", "year", "deposit_year", "method", "resolution", "bound_components",
            "ligand", "pmid", "r03_for_this_gene"]
    write_tsv(OUT_DIR / "gene_pdb_entries.tsv", rows, cols)
    write_tsv(OUT_DIR / "gene_year_structure.tsv", yrows, list(yrows[0]) if yrows else ["gene"])
    print(f"\nSaved: {OUT_DIR / 'gene_pdb_entries.tsv'} ({len(rows):,} rows)")
    print(f"Saved: {OUT_DIR / 'gene_year_structure.tsv'} ({len(yrows):,} rows)")


if __name__ == "__main__":
    main()
