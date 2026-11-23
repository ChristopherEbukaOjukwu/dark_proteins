#!/usr/bin/env python3
"""Step 3: for each year's list, every gene on it and the R03 awards tied to it.

Run from the dark_proteins folder, after step 2:
    python3 scripts/03_map_funded_genes.py

Year = the FOA list year. Round 1 used the 2018 list, round 2 the 2019 list, ... round 5 the 2022 list.
An award is tied to a gene in a year only if the gene was on that year's list.
Only genes confirmed by name are used (no inferred genes).

Inputs
    data/processed/foa_genes/foa_symbols_mapped.tsv   from step 2: genes on each year's list
    data/raw/r03_awards/R03_Award_List.csv           IDG award list (98 awards)
    data/raw/r03_awards/r03_unclear_genes_manual.csv the 11 awards checked by hand
    data/raw/r03_awards/reporter_r03_projects.json   NIH RePORTER records for the 98 awards
    data/raw/hgnc/hgnc_complete_set.txt              official gene names
    data/raw/foa_pages/RFA-RM-*.html                 the 5 announcements (posted and due dates)

Output
    data/processed/r03_awards/foa_genes_with_awards.tsv
        one row per year x gene x award; a gene with no award that year has one row with the award columns empty
        list_posted / list_due: when that year's list was posted and applications were due
        award_years: project end minus project start, in years (includes any extension)
    data/processed/r03_awards/gene_summary.tsv
        one row per gene (made from the file above): gene, hgnc_id, family,
        listed_years, funded_years, awards_per_year (e.g. 2018:2;2020:1), total_awards
"""
import csv
import json
import re
from collections import defaultdict
from datetime import date, datetime
from pathlib import Path

MAPPED = Path("data/processed/foa_genes/foa_symbols_mapped.tsv")
AWARDS = Path("data/raw/r03_awards/R03_Award_List.csv")
MANUAL = Path("data/raw/r03_awards/r03_unclear_genes_manual.csv")
REPORTER = Path("data/raw/r03_awards/reporter_r03_projects.json")
HGNC = Path("data/raw/hgnc/hgnc_complete_set.txt")
FOA_PAGES = Path("data/raw/foa_pages")
OUT = Path("data/processed/r03_awards/foa_genes_with_awards.tsv")
SUMMARY = Path("data/processed/r03_awards/gene_summary.tsv")
ROUND_YEAR = {"1": 2018, "2": 2019, "3": 2020, "4": 2021, "5": 2022}

# Printed name on the IDG award list -> gene symbol(s), checked against each award's title. [] = not a gene.
TOKEN_FIX = {
    "PKCQ": ["PRKCQ"], "PKR1": ["PROKR1"], "PKR2": ["PROKR2"], "PKACBETA": ["PRKACB"],
    "BPBWR1": ["NPBWR1"], "KIR4.2": ["KCNJ15"], "KIR5.1": ["KCNJ16"],
    "MRCKALPHA": ["CDC42BPA"], "MRCKBETA": ["CDC42BPB"], "MRCKGAMMA": ["CDC42BPG"],
    "CDC42BP": [], "TMEM16": [], "KCNG2-4": ["KCNG2", "KCNG3", "KCNG4"], "KV1.6": ["KCNA6"],
}
# Whole award replaced, because the title names the gene(s), or no gene is confirmed.
AWARD_FIX = {
    "MH122261": ["ADGRB1", "ADGRB2", "ADGRB3"],   # title: BAIs/ADGRBs; "ADGRB4" is not a gene
    "TR004471": ["KCNN1"],                         # blank on list; title names KCNN1
    "TR004452": ["MAP3K10"],                       # listed as MAP3K; title says MLK2 = MAP3K10
    "TR004204": [],                                # blank; gene never named (CLCC1 only inferred)
    "TR003700": [],                                # "PKA-C"; isoforms never named (only inferred)
    "TR002929": ["NPY2R", "NPY5R"],                # IDG list says NPY2R, NPY4R; abstract says "Y2 and Y5 receptors"
}


def read(path, delimiter=","):
    if not path.is_file():
        raise SystemExit(f"Missing {path}")
    with path.open(encoding="utf-8-sig", errors="replace", newline="") as f:
        return [{(k or "").strip(): (v or "").strip() for k, v in r.items()}
                for r in csv.DictReader(f, delimiter=delimiter)]


def core_of(project):
    return re.search(r"[A-Z]{2}\d{6}", project).group(0)


def day(value):
    return (value or "")[:10]


DATE = r"([A-Z][a-z]+ \d{1,2}, \d{4})"


def foa_dates():
    """year -> (posted date, application due date), read from the saved announcement pages."""
    out = {}
    for f in sorted(FOA_PAGES.glob("RFA-RM-*.html")):
        year = 2000 + int(re.search(r"RM-(\d\d)-", f.name).group(1))
        text = " ".join(re.sub(r"<[^>]+>", " ", f.read_text(encoding="utf-8", errors="replace")).split())
        posted = re.search(r"Posted Date\s*:?\s*" + DATE, text)
        app = re.search(r"Application Due Date", text)
        due = re.search(DATE, text[app.end():app.end() + 600]) if app else None
        iso = lambda m: datetime.strptime(m.group(1), "%B %d, %Y").date().isoformat() if m else ""
        out[year] = (iso(posted), iso(due))
    return out


def years_between(start, end):
    if not (start and end):
        return ""
    d = (date.fromisoformat(end) - date.fromisoformat(start)).days + 1
    return round(d / 365.25, 1)


def load_hgnc():
    current, previous, alias = {}, defaultdict(set), defaultdict(set)
    for r in read(HGNC, "\t"):
        if r.get("status", "Approved") != "Approved":
            continue
        current[r["symbol"].upper()] = r["hgnc_id"]
        for col, d in (("prev_symbol", previous), ("alias_symbol", alias)):
            for s in (r.get(col) or "").replace('"', "").split("|"):
                if s.strip():
                    d[s.strip().upper()].add(r["hgnc_id"])
    return current, previous, alias


def map_symbol(sym, current, previous, alias):
    key = sym.upper()
    if key in current:
        return current[key]
    if len(previous.get(key, ())) == 1:
        return next(iter(previous[key]))
    if len(alias.get(key, ())) == 1:
        return next(iter(alias[key]))
    return ""


def award_genes(core, listed, manual):
    """Gene symbols named for one award. Inferred (?), family: and unclear entries are dropped."""
    if core in manual:
        parts = [p.strip() for p in manual[core]["genes_found"].split(";")]
        return [p for p in parts if p and not p.endswith("?") and ":" not in p and p.lower() != "unclear"]
    if core in AWARD_FIX:
        return AWARD_FIX[core]
    out = []
    for tok in (t.strip() for t in listed.split(",")):
        if tok:
            out += TOKEN_FIX.get(tok.upper(), [tok])
    return out


def main():
    current, previous, alias = load_hgnc()
    list_dates = foa_dates()
    manual = {core_of(m["project_number"]): m for m in read(MANUAL)}

    # NIH RePORTER records, by full project number
    recs = json.loads(REPORTER.read_text(encoding="utf-8"))["results"]
    by_num = defaultdict(list)
    for r in recs:
        by_num[r.get("project_num", "")].append(r)

    # genes on each year's list (from step 2)
    on_list = defaultdict(dict)                      # year -> hgnc_id -> row
    for m in read(MAPPED, "\t"):
        on_list[int(m["year"])][m["hgnc_id"]] = m

    # awards -> (year, hgnc_id)
    tied = defaultdict(list)                         # (year, hgnc_id) -> [award info]
    not_on_list, not_mapped, no_gene, extra_recs = [], [], [], []
    for a in read(AWARDS):
        pnum, rnd = a["Project number"], a["Round of R03"]
        year, core = ROUND_YEAR[rnd], core_of(pnum)
        rr = by_num.get(pnum, [])
        if not rr:
            raise SystemExit(f"{pnum}: not found in {REPORTER}")
        nih = rr[0]
        others = sorted({r.get("project_num") for r in recs if core in (r.get("core_project_num") or "")} - {pnum})
        if len(rr) > 1 or others:
            extra_recs.append(f"{pnum}: {len(rr)} record(s) with this number; other numbers: {others or 'none'}")
        info = {"round": rnd, "project_number": pnum, "pi": a["PI Name (Last name, First name)"],
                "institution": a["Institution Name"], "title": a["R03 Title"],
                "fiscal_year": nih.get("fiscal_year", ""), "award_notice_date": day(nih.get("award_notice_date")),
                "project_start": day(nih.get("project_start_date")), "project_end": day(nih.get("project_end_date")),
                "budget_start": day(nih.get("budget_start")), "budget_end": day(nih.get("budget_end")),
                "award_years": years_between(day(nih.get("project_start_date")), day(nih.get("project_end_date"))),
                "award_amount": nih.get("award_amount", ""),
                "reporter_link": f"https://reporter.nih.gov/project-details/{nih.get('appl_id', '')}"}
        kept = set()
        for sym in award_genes(core, a["Specific Protein (if listed)"], manual):
            gid = map_symbol(sym, current, previous, alias)
            if not gid:
                not_mapped.append(f"{pnum}: {sym}")
            elif gid not in on_list[year]:
                not_on_list.append(f"{year} {pnum}: {sym}")
            elif gid not in kept:
                kept.add(gid)
                tied[(year, gid)].append(info)
        if not kept:
            no_gene.append(f"{year} {pnum}: {a['Specific Protein (if listed)'] or '(blank)'}")

    # one row per year x gene x award
    award_cols = ["round", "project_number", "pi", "institution", "title", "fiscal_year", "award_notice_date",
                  "project_start", "project_end", "award_years", "budget_start", "budget_end", "award_amount",
                  "reporter_link"]
    rows = []
    for year in sorted(on_list):
        for gid, m in sorted(on_list[year].items(), key=lambda kv: (kv[1]["family"], kv[1]["current_symbol"])):
            base = {"year": year, "list_posted": list_dates.get(year, ("", ""))[0],
                    "list_due": list_dates.get(year, ("", ""))[1], "family": m["family"], "gene": m["current_symbol"], "hgnc_id": gid,
                    "printed_on_list": m["printed_symbol"], "n_awards": len(tied.get((year, gid), []))}
            for info in tied.get((year, gid), [{}]):
                rows.append({**base, **{c: info.get(c, "") for c in award_cols}})

    # report
    print("Year | list posted | applications due")
    for year in sorted(on_list):
        posted, due = list_dates.get(year, ("", ""))
        print(f"{year} | {posted or 'NOT FOUND'} | {due or 'NOT FOUND'}")
    print("\nYear | genes on list | genes with an award | awards tied")
    for year in sorted(on_list):
        g = [gid for gid in on_list[year] if (year, gid) in tied]
        n_aw = len({i["project_number"] for gid in g for i in tied[(year, gid)]})
        print(f"{year} | {len(on_list[year])} | {len(g)} | {n_aw}")
    print(f"\nAwards with no confirmed gene on their year's list ({len(no_gene)}):")
    for x in no_gene:
        print("   " + x)
    print(f"\nNamed genes not on that year's list, left out ({len(not_on_list)}):")
    for x in not_on_list:
        print("   " + x)
    print(f"\nNames that did not map to a gene ({len(not_mapped)}): {not_mapped or 'none'}")
    print("\nYear | award length in years: number of awards")
    for year in sorted(on_list):
        lens = defaultdict(set)
        for gid in on_list[year]:
            for i in tied.get((year, gid), []):
                lens[i["award_years"]].add(i["project_number"])
        print(f"{year} | " + ", ".join(f"{k}: {len(v)}" for k, v in sorted(lens.items())))
    print(f"\nGrants with more than one RePORTER record ({len(extra_recs)}):")
    for x in extra_recs:
        print("   " + x)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"\nSaved: {OUT} ({len(rows)} rows)")

    # one row per gene, built from the rows above
    per_gene = {}
    for r in rows:
        g = per_gene.setdefault(r["hgnc_id"], {"gene": r["gene"], "hgnc_id": r["hgnc_id"], "family": set(),
                                               "listed": set(), "awards": defaultdict(set)})
        g["family"].add(r["family"])
        g["listed"].add(r["year"])
        if r["project_number"]:
            g["awards"][r["year"]].add(r["project_number"])
    summary = []
    for g in sorted(per_gene.values(), key=lambda g: (sorted(g["family"]), g["gene"])):
        fy = sorted(g["awards"])
        summary.append({"gene": g["gene"], "hgnc_id": g["hgnc_id"], "family": "/".join(sorted(g["family"])),
                        "listed_years": ";".join(map(str, sorted(g["listed"]))),
                        "funded_years": ";".join(map(str, fy)),
                        "awards_per_year": ";".join(f"{y}:{len(g['awards'][y])}" for y in fy),
                        "total_awards": sum(len(v) for v in g["awards"].values())})
    with SUMMARY.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(summary)
    funded = sum(1 for s in summary if s["total_awards"])
    print(f"Saved: {SUMMARY} ({len(summary)} genes: {funded} funded, {len(summary) - funded} not funded)")


if __name__ == "__main__":
    main()
