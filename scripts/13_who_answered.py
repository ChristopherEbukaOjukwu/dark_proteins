#!/usr/bin/env python3
"""Analysis 2: who answered the call? The R03 awards that funded genes on the 2018 list (any round).

Run from the dark_proteins folder, after scripts 03 and 06:
    python3 scripts/13_who_answered.py
It asks NIH RePORTER for each PI's grants once (about 1 query per PI, ~2 min); answers are saved and
reused, so reruns are offline.

Awards: every IDG R03 with at least one confirmed gene that was on the 2018 list (rounds 2018-2022).
"Before" = grants whose first fiscal year is earlier than the R03's first fiscal year.
PIs are matched by RePORTER profile ID, not by name.

Per award (one row each):
    round, award, PIs, institution, all genes on the award, which of them were on the 2018 list
    pi_prior_gene_grants     NIH grants naming the award's 2018-list gene(s) held by the award's PI(s)
                             before the R03 (gene-specific matches; IDG R03s and family names excluded)
    pi_prior_gene_first_fy   earliest such grant
    gene_prior_grants_any_pi NIH grants naming the gene(s) before the R03, by anyone
    pi_status                the most senior of the award's PIs, from all their NIH grants as PI before the R03:
                               R01-equivalent   R01, R37, RF1, R35, DP1, DP2
                               other research   any other research grant (R21, R03, R15, R00, R56, U, P, VA ...)
                               training only    only F, K or T awards
                               none as PI       no earlier NIH grant as PI
    pi_first_nih_fy          first fiscal year any of the award's PIs held an NIH grant as PI
    pi_awards_all_rounds     how many of the 98 IDG R03s the award's PI(s) held (repeat PIs)
    state                    the institution's state (RePORTER)
    inst_nih_rank            the institution's rank by total NIH funding in the R03's first fiscal year
    inst_band                top 25 / 26-100 / 101+ / not ranked (needs the ExPORTER files; blank without them)
    idg_centre_institution   the institution led an IDG phase-2 centre (IDG_CENTRES below)
Reference for inst_band: all new (type 1) NIH R03 awards in the same fiscal years, IDG R03s left out.
Limits: NIH grants only (no NSF, foundation, industry or non-US funding). Co-investigators are not PIs.

Inputs
    data/processed/r03_awards/foa_genes_with_awards.tsv, gene_summary.tsv   (script 03)
    data/raw/r03_awards/reporter_r03_projects.json                           (the 98 R03s in RePORTER)
    data/processed/funding/gene_grants.tsv                                   (script 06)
    data/raw/exporter/*.zip or *.csv   NIH ExPORTER project files, one per fiscal year (FY2019-FY2023),
                                       downloaded by hand from https://reporter.nih.gov/exporter (Projects).
                                       Optional: without them the rank columns stay blank.
Saved along the way
    data/raw/reporter_pis/pi_projects.jsonl    one line per PI: their RePORTER records as PI
Outputs
    data/processed/analysis/r03_awards_2018genes_who.tsv
    data/processed/analysis/r03_institution_bands_reference.tsv   IDG vs all new NIH R03s, % per band
                                                                   (only with the ExPORTER files)
"""
import csv
import json
import re
import time
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

AWARDS = Path("data/processed/r03_awards/foa_genes_with_awards.tsv")
GENES = Path("data/processed/r03_awards/gene_summary.tsv")
R03_PROJECTS = Path("data/raw/r03_awards/reporter_r03_projects.json")
GRANTS = Path("data/processed/funding/gene_grants.tsv")
CACHE = Path("data/raw/reporter_pis/pi_projects.jsonl")
OUT = Path("data/processed/analysis/r03_awards_2018genes_who.tsv")
REF_OUT = Path("data/processed/analysis/r03_institution_bands_reference.tsv")
API = "https://api.reporter.nih.gov/v2/projects/search"
FIELDS = ["ApplId", "ProjectNum", "CoreProjectNum", "ActivityCode", "FiscalYear", "ProjectStartDate",
          "PrincipalInvestigators", "Organization"]
EXPORTER = Path("data/raw/exporter")
R01_EQUIV = {"R01", "R37", "RF1", "R35", "DP1", "DP2"}
# Lead institutions of the IDG phase-2 centres (Oprea et al. 2018, Nat Rev Drug Discov, PMC6309548):
# GPCR and kinase DRGCs (UNC Chapel Hill), ion channel DRGC (UCSF), KMC (UNM), RDOC (Miami).
# Partner sites of these centres are not included. Patterns match NIH organisation names.
IDG_CENTRES = {r"NORTH CAROLINA.*CHAPEL HILL": "UNC Chapel Hill (GPCR and kinase DRGCs)",
               r"CALIFORNIA,? SAN FRANCISCO": "UCSF (ion channel DRGC)",
               r"NEW MEXICO": "UNM (KMC)",
               r"UNIVERSITY OF MIAMI": "Miami (RDOC)"}
BANDS = [("top 25", 1, 25), ("26-100", 26, 100), ("101+", 101, 10 ** 9)]
STATUS_ORDER = ["R01-equivalent", "other research", "training only", "none as PI"]


def read_tsv(path):
    if not path.is_file():
        raise SystemExit(f"Missing {path}")
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def core_of(project_number):
    m = re.search(r"[A-Z]{2}\d{6}", project_number or "")
    return "R03" + m.group(0) if m else ""


def fy(v):
    try:
        return int(str(v)[:4])
    except ValueError:
        return None


def kind(code):
    code = (code or "").upper()
    if code in R01_EQUIV:
        return "R01-equivalent"
    if code[:1] in ("F", "K", "T"):
        return "training only"
    return "other research"


# ---------------------------------------------------------------- RePORTER, by PI profile ID
def load_cache():
    out = {}
    if CACHE.is_file():
        for line in CACHE.open(encoding="utf-8"):
            d = json.loads(line)
            out[str(d["pi_profile_id"])] = d["records"]
    return out


def fetch_pi(pid, state):
    records, offset = [], 0
    while True:
        body = {"criteria": {"pi_profile_ids": [int(pid)], "fiscal_years": list(range(1985, 2027))},
                "include_fields": FIELDS, "offset": offset, "limit": 500}
        data = None
        for attempt in range(4):
            wait = 1.1 - (time.monotonic() - state["last"])
            if wait > 0:
                time.sleep(wait)
            try:
                req = urllib.request.Request(API, data=json.dumps(body).encode(),
                                             headers={"Content-Type": "application/json"})
                data = json.load(urllib.request.urlopen(req, timeout=120))
                state["last"] = time.monotonic()
                break
            except Exception as e:
                state["last"] = time.monotonic()
                print(f"   retry {attempt + 1} for PI {pid}: {e}")
                time.sleep(5 * (attempt + 1))
        if data is None:
            raise SystemExit("RePORTER did not answer after 4 tries; rerun later (finished PIs are kept).")
        for r in data.get("results", []):
            records.append({"core_project_num": r.get("core_project_num"), "activity_code": r.get("activity_code"),
                            "fiscal_year": r.get("fiscal_year"), "project_start": r.get("project_start_date"),
                            "org": (r.get("organization") or {}).get("org_name", "")})
        total = data.get("meta", {}).get("total", 0)
        offset += 500
        if offset >= total or offset >= 15000:
            break
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with CACHE.open("a", encoding="utf-8") as f:
        f.write(json.dumps({"pi_profile_id": str(pid), "records": records}) + "\n")
    return records


# ---------------------------------------------------------------- NIH ExPORTER: institution ranks
def exporter_rows():
    """Yield dict rows from every ExPORTER project file (zip or csv) in EXPORTER."""
    import io
    import zipfile
    for path in sorted(EXPORTER.glob("*")):
        if path.suffix.lower() == ".zip":
            with zipfile.ZipFile(path) as z:
                for name in z.namelist():
                    if name.lower().endswith(".csv"):
                        with z.open(name) as f:
                            yield from csv.DictReader(io.TextIOWrapper(f, encoding="latin-1", newline=""))
        elif path.suffix.lower() == ".csv":
            with path.open(encoding="latin-1", newline="") as f:
                yield from csv.DictReader(f)


def institution_ranks(years, idg_cores):
    """{fy: {org: rank}} by total NIH cost, and the new-R03 reference: {fy: Counter(org)}."""
    total, ref = defaultdict(lambda: defaultdict(float)), defaultdict(Counter)
    n = 0
    for r in exporter_rows():
        r = {k.strip().upper(): (v or "").strip() for k, v in r.items() if k}
        y = fy(r.get("FY"))
        if y not in years:
            continue
        n += 1
        org = r.get("ORG_NAME", "")
        try:
            total[y][org] += float(r.get("TOTAL_COST") or 0)   # sub-project rows carry no TOTAL_COST
        except ValueError:
            pass
        if r.get("ACTIVITY") == "R03" and r.get("APPLICATION_TYPE") == "1" \
                and r.get("CORE_PROJECT_NUM", "") not in idg_cores:
            ref[y][org] += 1
    ranks = {}
    for y, d in total.items():
        order = sorted(d, key=lambda o: -d[o])
        ranks[y] = {o: i + 1 for i, o in enumerate(order)}
    return ranks, ref, n


def band_of(rank):
    if not rank:
        return "not ranked"
    return next(b for b, lo, hi in BANDS if lo <= rank <= hi)


# ---------------------------------------------------------------- main
def main():
    listed_2018 = {g["hgnc_id"] for g in read_tsv(GENES) if "2018" in g["listed_years"].split(";")}

    # awards: genes (confirmed, on that year's list), round, PI/institution as listed by IDG
    aw = defaultdict(lambda: {"genes": {}, "years": set(), "pi_list": set(), "inst_list": set(), "pnums": set()})
    for r in read_tsv(AWARDS):
        if not r["project_number"]:
            continue
        c = core_of(r["project_number"])
        a = aw[c]
        a["genes"][r["hgnc_id"]] = r["gene"]
        a["years"].add(r["year"])
        a["pnums"].add(r["project_number"])
        a["pi_list"].add(r["pi"])
        a["inst_list"].add(r["institution"])

    # RePORTER records of the R03s: PI profile IDs, organisation, first fiscal year
    rep = defaultdict(lambda: {"pis": {}, "org": set(), "fy": set(), "state": set()})
    for r in json.loads(R03_PROJECTS.read_text(encoding="utf-8"))["results"]:
        c = r.get("core_project_num") or core_of(r.get("project_num", ""))
        for p in r.get("principal_investigators") or []:
            if p.get("profile_id"):
                rep[c]["pis"][str(p["profile_id"])] = p.get("full_name", "")
        if r.get("organization"):
            rep[c]["org"].add(r["organization"].get("org_name", ""))
            rep[c]["state"].add(r["organization"].get("org_state", "") or "")
        if fy(r.get("fiscal_year")):
            rep[c]["fy"].add(fy(r.get("fiscal_year")))

    # repeat PIs: across all 98 awards
    pi_awards = Counter(pid for c in rep for pid in rep[c]["pis"])

    chosen = sorted(c for c, a in aw.items() if set(a["genes"]) & listed_2018)
    print(f"IDG R03 awards: {len(aw)} with a confirmed gene | {len(chosen)} funded a gene on the 2018 list")

    # grants naming each gene: first fiscal year and PIs per core grant
    gg = defaultdict(lambda: defaultdict(lambda: {"fy": 9999, "pis": set(), "code": ""}))
    for r in read_tsv(GRANTS):
        if r["grant_group"] == "IDG R03" or r.get("family_name_only") == "1" or not fy(r["fiscal_year"]):
            continue
        d = gg[r["hgnc_id"]][r["core_project_num"]]
        d["fy"] = min(d["fy"], fy(r["fiscal_year"]))
        d["pis"] |= {p.strip() for p in r["pi_profile_ids"].split(";") if p.strip()}
        d["code"] = r.get("activity_code", "")

    # PI histories from RePORTER (cached)
    cache, state = load_cache(), {"last": 0.0}
    need = sorted({pid for c in chosen for pid in rep[c]["pis"]} - set(cache))
    if need:
        print(f"Asking RePORTER for {len(need)} PIs (about {len(need) * 1.2 / 60:.0f} min) ...")
    for k, pid in enumerate(need, 1):
        cache[pid] = fetch_pi(pid, state)
        if k % 20 == 0:
            print(f"   {k} of {len(need)}")

    # institution ranks (optional ExPORTER files)
    ranks, ref, n_exp = {}, {}, 0
    if EXPORTER.is_dir() and any(EXPORTER.glob("*")):
        years = {min(rep[c]["fy"]) for c in chosen if rep.get(c) and rep[c]["fy"]}
        print(f"Reading ExPORTER files for FY{min(years)}-FY{max(years)} ...")
        ranks, ref, n_exp = institution_ranks(years, set(rep))
        missing_fy = sorted(years - set(ranks))
        print(f"   {n_exp:,} project records" + (f" | no file for FY {missing_fy}" if missing_fy else ""))
    else:
        print(f"No ExPORTER files in {EXPORTER}: institution rank columns left blank")

    rows, missing = [], []
    for c in chosen:
        a, rp = aw[c], rep.get(c)
        r03_fy = min(rp["fy"]) if rp and rp["fy"] else int(min(a["years"])) + 1
        pids = set(rp["pis"]) if rp else set()
        if not pids:
            missing.append(c)
        g18 = {h: n for h, n in a["genes"].items() if h in listed_2018}

        # prior NIH grants naming the award's 2018-list genes
        pi_prior, any_prior = {}, set()
        for h in g18:
            for core, d in gg[h].items():
                if d["fy"] < r03_fy:
                    any_prior.add(core)
                    if d["pis"] & pids:
                        pi_prior[core] = min(pi_prior.get(core, 9999), d["fy"])

        # PI status from all their NIH grants as PI before the R03
        best, first_fy, prior_cores = "none as PI", None, set()
        for pid in pids:
            firsts = defaultdict(lambda: (9999, ""))
            for rec in cache.get(pid, []):
                y = fy(rec["fiscal_year"])
                core = rec["core_project_num"] or ""
                if not y or core_of(core) == c or core == c:
                    continue
                if y < firsts[core][0]:
                    firsts[core] = (y, rec["activity_code"])
            for core, (y, code) in firsts.items():
                first_fy = y if first_fy is None else min(first_fy, y)
                if y < r03_fy:
                    prior_cores.add(core)
                    k = kind(code)
                    if STATUS_ORDER.index(k) < STATUS_ORDER.index(best):
                        best = k

        rows.append({
            "round": min(a["years"]), "award": c, "project_numbers": "; ".join(sorted(a["pnums"])),
            "r03_first_fy": r03_fy,
            "pis": "; ".join(sorted(rp["pis"].values())) if rp and rp["pis"] else "; ".join(sorted(a["pi_list"])),
            "n_pis": len(pids) or len(a["pi_list"]),
            "institution": "; ".join(sorted(rp["org"])) if rp and rp["org"] else "; ".join(sorted(a["inst_list"])),
            "genes_all": "; ".join(sorted(a["genes"].values())), "n_genes": len(a["genes"]),
            "genes_2018_list": "; ".join(sorted(g18.values())),
            "pi_prior_gene_grants": len(pi_prior),
            "pi_prior_gene_first_fy": min(pi_prior.values()) if pi_prior else "",
            "pi_prior_gene_grant_ids": "; ".join(sorted(pi_prior)),
            "gene_prior_grants_any_pi": len(any_prior),
            "pi_status": best if pids else "not in RePORTER",
            "pi_prior_nih_grants": len(prior_cores),
            "pi_first_nih_fy": first_fy or "",
            "years_since_first_nih": (r03_fy - first_fy) if first_fy and first_fy < r03_fy else "",
            "pi_awards_all_rounds": max((pi_awards[p] for p in pids), default=1),
            "state": "; ".join(sorted(x for x in rp["state"] if x)) if rp else "",
            "inst_nih_rank": "", "inst_band": "",
            "idg_centre_institution": "; ".join(lab for pat, lab in IDG_CENTRES.items()
                                                for o in (rp["org"] if rp else a["inst_list"])
                                                if re.search(pat, o, re.I)),
        })
        if ranks:
            orgs = sorted(rp["org"]) if rp else []
            rk = [ranks.get(r03_fy, {}).get(o) for o in orgs]
            rk = [x for x in rk if x]
            rows[-1]["inst_nih_rank"] = min(rk) if rk else ""
            rows[-1]["inst_band"] = band_of(min(rk) if rk else None)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    # ---- plain summary
    n = len(rows)
    pct = lambda k: f"{k} ({100 * k / n:.0f}%)"
    print(f"\nAwards: {n} | by round: {dict(sorted(Counter(r['round'] for r in rows).items()))}")
    print(f"Awards with more than one PI: {sum(r['n_pis'] > 1 for r in rows)}")
    insts = Counter(i for r in rows for i in r["institution"].split("; ") if i)
    print(f"Institutions: {len(insts)} | with 2+ awards: "
          + (", ".join(f"{k} ({v})" for k, v in insts.most_common() if v > 1) or "none"))
    rep_pis = [r for r in rows if r["pi_awards_all_rounds"] > 1]
    print(f"Awards whose PI held 2+ IDG R03s (any round): {pct(len(rep_pis))}")
    print(f"\nPI already had an NIH grant naming the gene before the R03: {pct(sum(r['pi_prior_gene_grants'] > 0 for r in rows))}")
    print(f"The gene had an NIH grant by anyone before the R03:          {pct(sum(r['gene_prior_grants_any_pi'] > 0 for r in rows))}")
    print("\nPI status before the R03 (most senior PI on the award):")
    for s in STATUS_ORDER + ["not in RePORTER"]:
        k = sum(r["pi_status"] == s for r in rows)
        if k:
            print(f"   {s:16} {pct(k)}")
    yrs = sorted(r["years_since_first_nih"] for r in rows if r["years_since_first_nih"] != "")
    if yrs:
        print(f"Years from the PI's first NIH grant to the R03: median {yrs[len(yrs) // 2]} (range {yrs[0]}-{yrs[-1]})")
    if missing:
        print(f"Not in RePORTER (no PI ID): {', '.join(missing)}")
    idg = [r for r in rows if r["idg_centre_institution"]]
    print(f"\nAwards at an IDG centre's lead institution: {pct(len(idg))}"
          + (": " + ", ".join(f"{r['award']} ({r['idg_centre_institution']})" for r in idg) if idg else ""))
    if ranks:
        print("Institution's NIH funding rank in the R03's fiscal year (IDG awards vs all new NIH R03s, same years):")
        years = Counter(r["r03_first_fy"] for r in rows)
        ref_rows = []
        for b, _, _ in BANDS + [("not ranked", 0, 0)]:
            k = sum(r["inst_band"] == b for r in rows)
            # reference weighted to the IDG awards' fiscal years
            refshare = 0.0
            for y, w in years.items():
                tot = sum(ref[y].values())
                if tot:
                    refshare += w / n * sum(v for o, v in ref[y].items()
                                            if band_of(ranks.get(y, {}).get(o)) == b) / tot
            ref_rows.append({"band": b, "idg_awards": k, "idg_pct": round(100 * k / n, 1),
                             "all_new_r03_pct": round(100 * refshare, 1)})
            if k or refshare:
                print(f"   {b:10}  IDG {pct(k):10}   all new R03s {100 * refshare:4.0f}%")
        print("   (rank by total NIH cost from ExPORTER; reference weighted to the IDG awards' fiscal years)")
        with REF_OUT.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(ref_rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(ref_rows)
        print(f"Saved: {REF_OUT}")

    print(f"\nAwards where the PI already had a grant naming the gene:")
    for r in rows:
        if r["pi_prior_gene_grants"]:
            print(f"   {r['round']} {r['award']:13} {r['genes_2018_list']:22} {r['pis'][:40]:40} since FY{r['pi_prior_gene_first_fy']}")
    print(f"Saved: {OUT} ({n} rows)")


if __name__ == "__main__":
    main()
