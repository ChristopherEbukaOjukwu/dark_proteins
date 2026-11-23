#!/usr/bin/env python3
"""Step 2: map every printed symbol to its official gene ID, then compare the lists year by year.

Run from the dark_proteins folder, after step 1 has been run on all five pages:
    python3 scripts/02_map_and_compare_foa_lists.py

Part A - mapping (uses data/raw/hgnc/hgnc_complete_set.txt from HGNC).
For each symbol printed on an announcement, in this order:
    1. it is a current official symbol          -> that gene
    2. it is an old (previous) symbol of exactly one gene -> that gene
    3. it is an alias of exactly one gene       -> that gene
    4. otherwise                                -> NOT MAPPED (listed for you to check by hand)
Each gene is then identified by its HGNC ID, so a renamed gene (e.g. FAM26D and
CALHM4) counts as one gene.

Part B - comparison, by HGNC ID. For each step (2018->2019, ...), each gene is:
    stayed    on the previous list and on this one
    left      on the previous list, not on this one
    added     on this list, never on any earlier list
    re-added  on this list, not the previous one, but on some earlier list
Genes printed under a different symbol than the year before are reported too.

Outputs (data/processed/foa_genes/):
    foa_symbols_mapped.tsv   every printed symbol, per announcement, with its gene ID and how it was matched
    foa_list_changes.tsv     one row per gene per step, with its label
    foa_gene_history.tsv     one row per gene: on/off each year, symbols used, family, pattern
"""
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

FOLDER = Path("data/processed/foa_genes")
HGNC = Path("data/raw/hgnc/hgnc_complete_set.txt")


# ---------------------------------------------------------------- Part A: mapping
def load_hgnc():
    if not HGNC.is_file():
        raise SystemExit(f"Missing {HGNC}. Download it first (see instructions).")
    rows = list(csv.DictReader(HGNC.open(encoding="utf-8"), delimiter="\t"))
    current, previous, alias, info = {}, defaultdict(set), defaultdict(set), {}
    for r in rows:
        if r.get("status", "Approved") != "Approved":
            continue
        gid = r["hgnc_id"]
        info[gid] = {"symbol": r["symbol"], "locus_type": r.get("locus_type", ""), "name": r.get("name", "")}
        current[r["symbol"].upper()] = gid
        for s in (r.get("prev_symbol") or "").replace('"', "").split("|"):
            if s.strip():
                previous[s.strip().upper()].add(gid)
        for s in (r.get("alias_symbol") or "").replace('"', "").split("|"):
            if s.strip():
                alias[s.strip().upper()].add(gid)
    print(f"HGNC genes loaded: {len(info):,}")
    return current, previous, alias, info


def map_symbol(sym, current, previous, alias):
    key = sym.upper()
    if key in current:
        return current[key], "current symbol"
    if len(previous.get(key, ())) == 1:
        return next(iter(previous[key])), "old symbol"
    if len(alias.get(key, ())) == 1:
        return next(iter(alias[key])), "alias"
    if previous.get(key) or alias.get(key):
        return "", f"ambiguous ({len(previous.get(key, set()) | alias.get(key, set()))} genes)"
    return "", "not found"


def year_of(foa):
    return 2000 + int(re.search(r"RM-(\d\d)-", foa).group(1))


def load_lists(current, previous, alias, info):
    files = sorted(FOLDER.glob("foa_genes_RFA-RM-*.tsv"))
    if not files:
        raise SystemExit(f"No step-1 files found in {FOLDER}")
    mapped_rows, lists = [], {}
    for f in files:
        rows = list(csv.DictReader(f.open(encoding="utf-8"), delimiter="\t"))
        foa = rows[0]["foa"]
        year = year_of(foa)
        genes = {}
        for r in rows:
            gid, how = map_symbol(r["symbol"], current, previous, alias)
            key = gid or f"UNMAPPED:{r['symbol']}"
            genes[key] = {"printed": r["symbol"], "family": r["family"]}
            mapped_rows.append({"foa": foa, "year": year, "family": r["family"], "printed_symbol": r["symbol"],
                                "hgnc_id": gid, "current_symbol": info[gid]["symbol"] if gid else "",
                                "matched_by": how, "locus_type": info[gid]["locus_type"] if gid else ""})
        lists[year] = {"foa": foa, "genes": genes}
    return dict(sorted(lists.items())), mapped_rows


def label(key, info):
    if key.startswith("UNMAPPED:"):
        return key
    return info[key]["symbol"]


# ---------------------------------------------------------------- main
def main():
    current, previous, alias, info = load_hgnc()
    lists, mapped = load_lists(current, previous, alias, info)
    years = list(lists)

    # ---- Part A report
    print("\n===== Mapping printed symbols to genes =====")
    for y in years:
        rows = [m for m in mapped if m["year"] == y]
        print(f"{y}: {len(rows)} symbols -> {dict(Counter(m['matched_by'] for m in rows))}")
    renamed = sorted({(m["printed_symbol"], m["current_symbol"]) for m in mapped
                      if m["matched_by"] in ("old symbol", "alias")})
    print(f"\nPrinted under an old name or alias ({len(renamed)}):")
    for printed, cur in renamed:
        print(f"   {printed} -> {cur}")
    bad = sorted({(m["printed_symbol"], m["matched_by"]) for m in mapped if not m["hgnc_id"]})
    print(f"\nNOT MAPPED - please check by hand ({len(bad)}): {bad or 'none'}")
    non_coding = sorted({(m["current_symbol"], m["locus_type"]) for m in mapped
                         if m["hgnc_id"] and m["locus_type"] != "gene with protein product"})
    print(f"Mapped, but not a protein-coding gene ({len(non_coding)}): {non_coding or 'none'}")
    for y in years:
        per_gene = Counter(m["hgnc_id"] for m in mapped if m["year"] == y and m["hgnc_id"])
        twice = [(label(g, info), n) for g, n in per_gene.items() if n > 1]
        if twice:
            print(f"   {y}: two printed symbols map to the same gene: {twice}")

    # ---- Part B: year-by-year comparison by gene ID
    change_rows = []
    for prev, cur in zip(years, years[1:]):
        before, now = lists[prev]["genes"], lists[cur]["genes"]
        earlier = set().union(*(lists[y]["genes"] for y in years if y < prev))
        groups = {"stayed": [], "left": [], "added": [], "re-added": []}
        moved, renamed_here = [], []
        for g in sorted(set(before) | set(now), key=lambda k: label(k, info)):
            if g in before and g in now:
                lab = "stayed"
                if before[g]["family"] != now[g]["family"]:
                    moved.append(f"{label(g, info)} ({before[g]['family']} -> {now[g]['family']})")
                if before[g]["printed"] != now[g]["printed"]:
                    renamed_here.append(f"{before[g]['printed']} -> {now[g]['printed']}")
            elif g in before:
                lab = "left"
            elif g in earlier:
                lab = "re-added"
            else:
                lab = "added"
            fam = (now.get(g) or before.get(g))["family"]
            groups[lab].append((label(g, info), fam))
            change_rows.append({"from_year": prev, "to_year": cur, "hgnc_id": g if not g.startswith("UNMAPPED") else "",
                                "gene": label(g, info), "family": fam, "change": lab,
                                "printed_before": before.get(g, {}).get("printed", ""),
                                "printed_now": now.get(g, {}).get("printed", "")})
        print(f"\n===== {prev} -> {cur} =====")
        print(f"stayed {len(groups['stayed'])} | left {len(groups['left'])} | "
              f"added {len(groups['added'])} | re-added {len(groups['re-added'])}")
        for lab in ("left", "added", "re-added"):
            items = groups[lab]
            if items:
                for fam in sorted({f for _, f in items}):
                    names = [g for g, f in items if f == fam]
                    print(f"  {lab} {fam} ({len(names)}): {', '.join(names)}")
        if renamed_here:
            print(f"  same gene, printed differently: {', '.join(renamed_here)}")
        if moved:
            print(f"  same gene, family changed: {', '.join(moved)}")

    # ---- one row per gene
    keys = sorted(set().union(*(lists[y]["genes"] for y in years)), key=lambda k: label(k, info))
    history = []
    for g in keys:
        on = [g in lists[y]["genes"] for y in years]
        printed = sorted({lists[y]["genes"][g]["printed"] for y in years if g in lists[y]["genes"]})
        fams = sorted({lists[y]["genes"][g]["family"] for y in years if g in lists[y]["genes"]})
        history.append({"hgnc_id": g if not g.startswith("UNMAPPED") else "", "gene": label(g, info),
                        "locus_type": info[g]["locus_type"] if g in info else "",
                        "family": "/".join(fams), "printed_as": "/".join(printed),
                        **{f"on_{y}": int(x) for y, x in zip(years, on)},
                        "pattern": "".join("1" if x else "0" for x in on), "n_lists": sum(on)})
    print("\n===== Whole period =====")
    print(f"Distinct genes on any list: {len(history)}")
    print(f"On every list: {sum(h['n_lists'] == len(years) for h in history)}")
    print("Most common on/off patterns (" + ", ".join(map(str, years)) + "):")
    for pat, n in Counter(h["pattern"] for h in history).most_common(10):
        print(f"   {pat}: {n}")

    for name, rows in (("foa_symbols_mapped.tsv", mapped), ("foa_list_changes.tsv", change_rows),
                       ("foa_gene_history.tsv", history)):
        with (FOLDER / name).open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        print(f"Saved: {FOLDER / name}")


if __name__ == "__main__":
    main()
