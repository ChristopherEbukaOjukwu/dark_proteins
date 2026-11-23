#!/usr/bin/env python3
"""Analysis 5: what changed after the 2018 announcement? The 328 genes on the 2018 list,
later funded (an IDG R03 in any round) vs never funded, each gene against its own earlier pace.

Run from the dark_proteins folder, after scripts 04-12:
    python3 scripts/14_what_changed.py

Per gene and measure:
    before   yearly average over 2013-2017 (the five years before the announcement)
    after    yearly average from 2019 to the measure's last reliable year (2018 = transition year, left out)
    change   after - before   (0 = kept its earlier pace)
Measures (yearly counts of NEW items unless said otherwise):
    5a  IDG's criteria   jensen (score per year, to 2025), pubtator (score per year, to 2021 only),
                         r01 (new R01-equivalent grants naming the gene: R01, R37, RF1; to FY2026)
    5b  ours             focused papers (<=10 genes, human + one-to-one, to 2025), grants (new NIH grants,
                         any kind, IDG R03s left out, to FY2026), partners (new small-scale human partners,
                         to 2023), go (new informative experimental GO terms, to 2024), compounds (new
                         compounds tested, single protein, valid, to 2024), structures (new PDB entries, to 2026)
Milestones (firsts): among genes WITHOUT it by the end of 2017, % that first reached it from 2018 to the
    measure's last reliable year: first R01, first structure, first potent compound, first small-scale
    partner, first informative GO term, first Moderate+ disease link (GenCC, by earliest paper year).
From the program: for funded genes, how many of the new items after the announcement are tied to the
    gene's own R03 (the source paper cites an R03 that funded this gene; scripts 05-10).

Tests, funded vs never funded (and, as a check, round-1-funded vs never funded):
    change      Mann-Whitney; effect = P(funded change higher); also % of genes whose pace rose
    milestones  Fisher's exact; difference in % with a 95% range (Newcombe)
    all p       Benjamini-Hochberg together (p_adj)

Expected vs actual (ANCOVA), per measure, on log10(value + c) so the before-after band is a straight line:
    band       after ~ before + family, all 328 genes  -> each gene's "vs_expected" = how far above (+) or
               below (-) the band it sits (log10 units; +0.30 = about 2x its expected pace)
    effect     after ~ before + family + funded          -> the funded coefficient, 95% range, p (normal approx.);
               shown as "% above expected". Also round-1 funded vs never funded. BH within this set.
    took off   the funded genes furthest above the band (Jensen, focused papers, grants), with what ties them
               to their R03 (R03-tied papers, a later same-PI grant, a first structure from 2018)
Timing, each funded gene against itself: year 0 = the year its funding started (first fiscal year of its
    first R03); +1, +2 ... = the years after. Baseline = the gene's own average over the 5 years before funding
    (at least 3 of them inside the measure's reliable years). For year 0 to +5: the gene's value that year
    minus its baseline, in the measure's own units; averaged over funded genes, 95% range from resampling
    genes (2,000 times).
    Never-funded genes, same analysis: they have no award, so each borrows the start years of the funded genes
    in its family, in the same proportions (only start years that reach year k). For each never-funded gene:
    its own change from its own 5-year baseline before each borrowed start year, weighted by how many funded
    genes started that year -> one number per gene per year, as for funded genes. Line = average over genes,
    weighted to the funded genes' family mix; 95% range from resampling genes (within family), as for funded.
Checks
    round 1 by family   round-1 funded vs never funded within GPCRs, kinases and ion channels (first structure,
                        new structures, new grants): is the round-1 signal one family?
    R03 -> later grant  per award (awards that funded a 2018-list gene): did the award's PI later hold an NIH
                        grant naming that gene, first funded AFTER the R03's first fiscal year? R01 group
                        (R01, R37, RF1) and any grant. PIs matched by RePORTER profile ID. The FOA's stated aim
                        was preliminary data for a later R01. No comparison group: description only.

Outputs (data/processed/analysis/)
    gene_change.tsv            one row per gene: before / after / change per measure, milestones, R03-tied counts
    gene_change_tests.tsv      one row per comparison x measure
    r03_followon_grants.tsv    one row per award: later grants by the same PI naming the award's gene(s)
    gene_change_ancova.tsv     one row per comparison x measure: funded effect after allowing for before + family
    gene_event_time.tsv        one row per measure x year from funding: mean change from own baseline, 95% range, n
"""
import csv
import json
import math
import random
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path

P = Path("data/processed")
OUT = P / "analysis"
BEFORE = range(2013, 2018)
AFTER_START = 2019
LAST = {"jensen": 2025, "pubtator": 2021, "r01": 2026, "papers": 2025, "grants": 2026,
        "partners": 2023, "go": 2024, "compounds": 2024, "structures": 2026}
PART_A = ["jensen", "pubtator", "r01"]
PART_B = ["papers", "grants", "partners", "go", "compounds", "structures"]
MILESTONES = ["first_r01", "first_structure", "first_potent", "first_partner", "first_go", "first_disease"]
# offset c in log10(value + c): lets 0 be used; about the smallest non-zero yearly value
OFFSET = {"jensen": 0.1, "pubtator": 0.1, "r01": 0.05, "papers": 0.2, "grants": 0.05,
          "partners": 0.05, "go": 0.05, "compounds": 0.2, "structures": 0.05}
EVENT_K = range(0, 6)          # years from the start of funding
PRE_K = range(-5, 0)           # the 5 years before: each gene's own baseline
EVENT_MEASURES = ["jensen", "pubtator", "r01", "papers", "grants", "partners", "go", "compounds", "structures"]
MILESTONE_LAST = {"first_r01": 2026, "first_structure": 2026, "first_potent": 2024, "first_partner": 2023,
                  "first_go": 2024, "first_disease": 2025}


def read_tsv(path):
    if not path.is_file():
        raise SystemExit(f"Missing {path}")
    with path.open(encoding="utf-8", newline="") as f:
        yield from csv.DictReader(f, delimiter="\t")


def yr(v):
    try:
        return int(str(v)[:4])
    except ValueError:
        return None


# ---------------------------------------------------------------- tests (same as script 12)
def norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2))


def mann_whitney(x, y):
    n1, n2 = len(x), len(y)
    allv = sorted((v, i) for i, v in enumerate(x + y))
    ranks, i, ties = [0.0] * (n1 + n2), 0, 0.0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        for k in range(i, j + 1):
            ranks[allv[k][1]] = (i + j) / 2 + 1
        t = j - i + 1
        ties += t ** 3 - t
        i = j + 1
    u1 = sum(ranks[:n1]) - n1 * (n1 + 1) / 2
    n = n1 + n2
    sd = math.sqrt(n1 * n2 / 12 * ((n + 1) - ties / (n * (n - 1)))) if n > 1 else 0
    if sd == 0:
        return 1.0, u1 / (n1 * n2)
    z = (abs(u1 - n1 * n2 / 2) - 0.5) / sd
    return min(1.0, 2 * norm_sf(max(z, 0))), u1 / (n1 * n2)


def fisher_exact(a, b, c, d):
    r1, c1, n = a + b, a + c, a + b + c + d
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)
    lp = lambda k: (math.lgamma(r1 + 1) - math.lgamma(k + 1) - math.lgamma(r1 - k + 1)
                    + math.lgamma(n - r1 + 1) - math.lgamma(c1 - k + 1) - math.lgamma(n - r1 - c1 + k + 1)
                    - math.lgamma(n + 1) + math.lgamma(c1 + 1) + math.lgamma(n - c1 + 1))
    obs = lp(a)
    return min(1.0, sum(math.exp(lp(k)) for k in range(lo, hi + 1) if lp(k) <= obs + 1e-7))


def wilson(x, n, z=1.959964):
    if n == 0:
        return 0.0, 1.0
    p = x / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


def newcombe(x1, n1, x2, n2):
    p1, p2 = x1 / n1, x2 / n2
    l1, u1 = wilson(x1, n1)
    l2, u2 = wilson(x2, n2)
    d = p1 - p2
    return d - math.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2), d + math.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)


def bh(ps):
    m = len(ps)
    order = sorted(range(m), key=lambda i: ps[i])
    adj, prev = [0.0] * m, 1.0
    for rank, i in reversed(list(enumerate(order, 1))):
        prev = min(prev, ps[i] * m / rank)
        adj[i] = prev
    return adj


# ---------------------------------------------------------------- least squares (no extra packages)
def solve(a, b):
    """Solve a x = b (small square system) by Gauss-Jordan elimination."""
    n = len(a)
    m = [row[:] + [bb] for row, bb in zip(a, b)]
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(m[r][i]))
        m[i], m[piv] = m[piv], m[i]
        d = m[i][i]
        m[i] = [v / d for v in m[i]]
        for r in range(n):
            if r != i and m[r][i]:
                f = m[r][i]
                m[r] = [vr - f * vi for vr, vi in zip(m[r], m[i])]
    return [m[i][n] for i in range(n)]


def ols(X, y):
    """Coefficients, their standard errors and residuals."""
    k = len(X[0])
    xtx = [[sum(r[i] * r[j] for r in X) for j in range(k)] for i in range(k)]
    xty = [sum(r[i] * v for r, v in zip(X, y)) for i in range(k)]
    beta = solve(xtx, xty)
    res = [v - sum(b * x for b, x in zip(beta, r)) for r, v in zip(X, y)]
    s2 = sum(e * e for e in res) / (len(y) - k)
    inv = [solve(xtx, [1.0 if i == j else 0.0 for i in range(k)]) for j in range(k)]
    se = [math.sqrt(s2 * inv[i][i]) for i in range(k)]
    return beta, se, res


# ---------------------------------------------------------------- collect yearly series
def main():
    genes = {r["hgnc_id"]: r for r in read_tsv(P / "r03_awards/gene_summary.tsv")
             if "2018" in r["listed_years"].split(";")}
    print(f"Genes on the 2018 list: {len(genes)}")
    series = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))   # measure -> gene -> year -> value
    tied = defaultdict(lambda: defaultdict(int))                              # measure -> gene -> tied new items after
    has_data = defaultdict(set)                                               # measure -> genes with a record at all
    first = defaultdict(dict)                                                 # milestone -> gene -> first year

    def first_year(m, h, y):
        if y and (h not in first[m] or y < first[m][h]):
            first[m][h] = y

    # Jensen and PubTator scores
    for r in read_tsv(P / "attention/gene_year_scores.tsv"):
        h, y = r["hgnc_id"], int(r["year"])
        if h not in genes:
            continue
        for m, col in (("jensen", "jensen_score"), ("pubtator", "pubtator_score")):
            if r[col] != "":
                has_data[m].add(h)
                series[m][h][y] += float(r[col])

    # focused papers (<=10 genes, human + one-to-one)
    seen = set()
    for r in read_tsv(P / "attention/gene_papers.tsv"):
        h, y = r["hgnc_id"], yr(r["year"])
        if h not in genes or not y or not (r["species"] == "human" or r["match_type"] == "one-to-one"):
            continue
        has_data["papers"].add(h)
        if not r["n_genes_in_paper"] or int(r["n_genes_in_paper"]) > 10 or (h, r["pmid"]) in seen:
            continue
        seen.add((h, r["pmid"]))
        series["papers"][h][y] += 1
        if r["r03_for_this_gene"] and AFTER_START <= y <= LAST["papers"]:
            tied["papers"][h] += 1

    # grants: new cores per fiscal year (gene-specific, IDG R03s out); R01-equivalents separately
    cores = defaultdict(dict)
    for r in read_tsv(P / "funding/gene_grants.tsv"):
        h, y = r["hgnc_id"], yr(r["fiscal_year"])
        if h not in genes or not y or r["grant_group"] == "IDG R03" or r.get("family_name_only") == "1":
            continue
        c = r["core_project_num"]
        if c not in cores[h] or y < cores[h][c][0]:
            cores[h][c] = (y, r["grant_group"])
    for h in genes:
        has_data["grants"].add(h)
        has_data["r01"].add(h)
        for c, (y, group) in cores[h].items():
            series["grants"][h][y] += 1
            if group == "R01":                       # script 06's R01 group = R01, R37, RF1
                series["r01"][h][y] += 1
                first_year("first_r01", h, y)

    # small-scale human partners: new partner per year
    pfirst = defaultdict(dict)
    for r in read_tsv(P / "interactions/gene_interactions.tsv"):
        h, y = r["hgnc_id"], yr(r["year"])
        if (h not in genes or not y or r["match_type"] != "same gene" or r["self_interaction"] != "0"
                or "low" not in r["throughput"].lower()):
            continue
        p = r["partner_gene_id"]
        old = pfirst[h].get(p)
        if old is None or y < old[0]:
            pfirst[h][p] = (y, bool(r["r03_for_this_gene"]))
        elif y == old[0] and r["r03_for_this_gene"]:
            pfirst[h][p] = (y, True)
    for h in genes:
        has_data["partners"].add(h)
        for p, (y, t) in pfirst[h].items():
            series["partners"][h][y] += 1
            first_year("first_partner", h, y)
            if t and AFTER_START <= y <= LAST["partners"]:
                tied["partners"][h] += 1

    # compounds and potent
    cfirst = defaultdict(dict)
    for r in read_tsv(P / "chemistry/gene_activities.tsv"):
        h, y = r["hgnc_id"], yr(r["year"])
        if h not in genes or not y or r["target_type"] != "SINGLE PROTEIN" or r["valid"] != "1":
            continue
        c = r["compound"]
        if c not in cfirst[h] or y < cfirst[h][c][0]:
            cfirst[h][c] = (y, bool(r["r03_for_this_gene"]))
        if r["potent"] == "1":
            first_year("first_potent", h, y)
    for h in genes:
        has_data["compounds"].add(h)
        for c, (y, t) in cfirst[h].items():
            series["compounds"][h][y] += 1
            if t and AFTER_START <= y <= LAST["compounds"]:
                tied["compounds"][h] += 1

    # GO: new informative experimental terms
    gfirst = defaultdict(dict)
    for r in read_tsv(P / "go/gene_go_annotations.tsv"):
        h, y = r["hgnc_id"], yr(r["year"])
        if (h not in genes or not y or r["evidence_group"] != "experimental" or r["informative"] != "1"
                or r["negated"] != "0"):
            continue
        g = r["go_id"]
        if g not in gfirst[h] or y < gfirst[h][g][0]:
            gfirst[h][g] = (y, bool(r["r03_for_this_gene"]))
    for h in genes:
        has_data["go"].add(h)
        for g, (y, t) in gfirst[h].items():
            series["go"][h][y] += 1
            first_year("first_go", h, y)
            if t and AFTER_START <= y <= LAST["go"]:
                tied["go"][h] += 1

    # structures
    sseen = set()
    for r in read_tsv(P / "structure/gene_pdb_entries.tsv"):
        h, y = r["hgnc_id"], yr(r["year"])
        if h not in genes or not y or (h, r["pdb_id"]) in sseen:
            continue
        sseen.add((h, r["pdb_id"]))
        series["structures"][h][y] += 1
        first_year("first_structure", h, y)
        if r["r03_for_this_gene"] and AFTER_START <= y <= LAST["structures"]:
            tied["structures"][h] += 1
    for h in genes:
        has_data["structures"].add(h)

    # disease: first Moderate+ link by earliest paper year
    for r in read_tsv(P / "disease/gene_disease_links.tsv"):
        h = r["hgnc_id"]
        if h in genes and r["best_grade"] in ("Definitive", "Strong", "Moderate"):
            first_year("first_disease", h, yr(r["earliest_paper_year"]))

    # ---------------------------------------------------------------- per gene
    rows = []
    for h, g in sorted(genes.items(), key=lambda kv: (kv[1]["family"], kv[1]["gene"])):
        fy = sorted(int(y) for y in g["funded_years"].split(";") if y)
        row = {"gene": g["gene"], "hgnc_id": h, "family": g["family"], "funded_years": g["funded_years"],
               "funded_any": int(bool(fy)), "first_funded_round": fy[0] if fy else "",
               "funded_round1": int(2018 in fy)}
        for m in PART_A + PART_B:
            if h not in has_data[m]:
                row[f"{m}_before"] = row[f"{m}_after"] = row[f"{m}_change"] = ""
                continue
            d = series[m][h]
            after = range(AFTER_START, LAST[m] + 1)
            b = sum(d.get(y, 0) for y in BEFORE) / len(BEFORE)
            a = sum(d.get(y, 0) for y in after) / len(after)
            row[f"{m}_before"], row[f"{m}_after"], row[f"{m}_change"] = round(b, 3), round(a, 3), round(a - b, 3)
        for m in MILESTONES:
            f = first[m].get(h)
            had = f is not None and f <= 2017
            row[f"{m}_had_by_2017"] = int(had)
            row[m] = "" if had else int(f is not None and 2018 <= f <= MILESTONE_LAST[m])
        for m in PART_B:
            row[f"{m}_r03_tied_after"] = tied[m].get(h, 0)
        rows.append(row)

    # ---------------------------------------------------------------- tests
    never = [r for r in rows if not r["funded_any"]]
    comparisons = [("funded vs never", [r for r in rows if r["funded_any"]]),
                   ("round-1 funded vs never", [r for r in rows if r["funded_round1"]])]
    tests = []
    for label, grp in comparisons:
        for m in PART_A + PART_B:
            col = f"{m}_change"
            x = [float(r[col]) for r in grp if r[col] != ""]
            y = [float(r[col]) for r in never if r[col] != ""]
            p, ps = mann_whitney(x, y)
            up = lambda v: 100 * sum(t > 0 for t in v) / len(v)
            tests.append({"comparison": label, "part": "5a" if m in PART_A else "5b", "measure": m,
                          "test": "Mann-Whitney", "n_funded": len(x), "n_never": len(y),
                          "funded": f"median change {statistics.median(x):+.3g}; {up(x):.0f}% rose",
                          "never": f"median change {statistics.median(y):+.3g}; {up(y):.0f}% rose",
                          "effect": f"P(funded change higher) = {ps:.2f}", "range_95": "", "p": p})
        for m in MILESTONES:
            x = [r[m] for r in grp if r[m] != ""]
            y = [r[m] for r in never if r[m] != ""]
            if not x or not y:
                continue
            a, c = sum(x), sum(y)
            p = fisher_exact(a, len(x) - a, c, len(y) - c)
            lo, hi = newcombe(a, len(x), c, len(y))
            tests.append({"comparison": label, "part": "milestone", "measure": m, "test": "Fisher",
                          "n_funded": len(x), "n_never": len(y),
                          "funded": f"{a}/{len(x)} = {100 * a / len(x):.0f}%",
                          "never": f"{c}/{len(y)} = {100 * c / len(y):.0f}%",
                          "effect": f"{100 * (a / len(x) - c / len(y)):+.1f} points",
                          "range_95": f"{100 * lo:+.1f} to {100 * hi:+.1f}", "p": p})
    for t, pa in zip(tests, bh([t["p"] for t in tests])):
        t["p_adj"] = pa

    for label, _ in comparisons:
        print(f"\n{label.upper()}  (p_adj = Benjamini-Hochberg across all {len(tests)} tests)")
        for t in tests:
            if t["comparison"] == label:
                print(f"   {t['part']:9} {t['measure']:16} {t['funded']:34} | {t['never']:34} "
                      f"{t['effect']:28} {t['range_95']:17} p={t['p']:.3f} p_adj={t['p_adj']:.3f}")

    # from the program
    funded = [r for r in rows if r["funded_any"]]
    print("\nNew items after the announcement in funded genes, and how many are tied to the gene's own R03:")
    for m in PART_B:
        if m == "grants":
            continue
        span = LAST[m] - AFTER_START + 1
        total = sum(float(r[f"{m}_after"] or 0) * span for r in funded)
        t = sum(r[f"{m}_r03_tied_after"] for r in funded)
        g = sum(r[f"{m}_r03_tied_after"] > 0 for r in funded)
        print(f"   {m:11} {total:8.0f} new (2019-{LAST[m]})   R03-tied {t:5d} ({100 * t / total if total else 0:4.1f}%)"
              f"   on {g} of {len(funded)} funded genes")

    # ---------------------------------------------------------------- check 1: round 1 by family
    print("\nCheck: round-1 funded vs never funded, within each family")
    print(f"   {'family':12} {'n r1/never':>11}  {'first structure':>27}  {'structures change: P(r1 higher)':>32}  "
          f"{'grants change: P(r1 higher)':>28}")
    for fam in sorted({r["family"] for r in rows}):
        r1 = [r for r in rows if r["family"] == fam and r["funded_round1"]]
        nv = [r for r in rows if r["family"] == fam and not r["funded_any"]]
        if not r1:
            continue
        a = [r["first_structure"] for r in r1 if r["first_structure"] != ""]
        c = [r["first_structure"] for r in nv if r["first_structure"] != ""]
        ms = f"{sum(a)}/{len(a)} vs {sum(c)}/{len(c)}" if a and c else "-"
        pf = fisher_exact(sum(a), len(a) - sum(a), sum(c), len(c) - sum(c)) if a and c else 1.0
        out = []
        for m in ("structures", "grants"):
            x = [float(r[f"{m}_change"]) for r in r1 if r[f"{m}_change"] != ""]
            y = [float(r[f"{m}_change"]) for r in nv if r[f"{m}_change"] != ""]
            p, e = mann_whitney(x, y)
            out.append(f"{e:.2f} (p={p:.2f})")
        print(f"   {fam:12} {len(r1):>4}/{len(nv):<6}  {ms:>18} (p={pf:.2f})  {out[0]:>32}  {out[1]:>28}")

    # ---------------------------------------------------------------- check 2: R03 -> later grants by the same PI
    def core_of(pn):
        mm = re.search(r"[A-Z]{2}\d{6}", pn or "")
        return "R03" + mm.group(0) if mm else ""

    award_genes = defaultdict(set)
    for r in read_tsv(P / "r03_awards/foa_genes_with_awards.tsv"):
        if r["project_number"] and r["hgnc_id"] in genes:
            award_genes[core_of(r["project_number"])].add(r["hgnc_id"])
    pis, r03fy = defaultdict(set), {}
    for r in json.loads(Path("data/raw/r03_awards/reporter_r03_projects.json").read_text(encoding="utf-8"))["results"]:
        c = r.get("core_project_num") or core_of(r.get("project_num", ""))
        for pp in r.get("principal_investigators") or []:
            if pp.get("profile_id"):
                pis[c].add(str(pp["profile_id"]))
        if yr(r.get("fiscal_year")):
            r03fy[c] = min(r03fy.get(c, 9999), yr(r.get("fiscal_year")))
    later = defaultdict(lambda: defaultdict(dict))        # award -> core grant -> info
    gstart = defaultdict(dict)
    for r in read_tsv(P / "funding/gene_grants.tsv"):
        if r["grant_group"] == "IDG R03" or r.get("family_name_only") == "1" or not yr(r["fiscal_year"]):
            continue
        k = (r["hgnc_id"], r["core_project_num"])
        y = yr(r["fiscal_year"])
        old = gstart[k]
        if not old or y < old["fy"]:
            gstart[k] = {"fy": y, "group": r["grant_group"], "gene": r["gene"],
                         "pis": {x.strip() for x in r["pi_profile_ids"].split(";") if x.strip()}}
        else:
            old["pis"] |= {x.strip() for x in r["pi_profile_ids"].split(";") if x.strip()}
    for c, hs in award_genes.items():
        if c not in r03fy:
            continue
        for (h, core), d in gstart.items():
            if h in hs and d["pis"] & pis[c] and d["fy"] > r03fy[c]:
                later[c][core] = d
    fo = []
    for c in sorted(award_genes):
        if c not in r03fy:
            continue
        g = later.get(c, {})
        r01 = {k: d for k, d in g.items() if d["group"] == "R01"}
        fo.append({"award": c, "r03_first_fy": r03fy[c],
                   "genes_2018_list": "; ".join(sorted(genes[h]["gene"] for h in award_genes[c])),
                   "later_r01": len(r01), "later_any_grant": len(g),
                   "later_grants": "; ".join(f"{k} ({d['group']}, FY{d['fy']}, {d['gene']})"
                                             for k, d in sorted(g.items(), key=lambda kv: kv[1]["fy"]))})
    n = len(fo)
    print(f"\nCheck: did the R03 PI later hold an NIH grant naming the award's gene? ({n} awards with a RePORTER record)")
    print(f"   later R01 (R01/R37/RF1): {sum(f['later_r01'] > 0 for f in fo)} of {n} awards "
          f"({100 * sum(f['later_r01'] > 0 for f in fo) / n:.0f}%)")
    print(f"   any later NIH grant:     {sum(f['later_any_grant'] > 0 for f in fo)} of {n} awards "
          f"({100 * sum(f['later_any_grant'] > 0 for f in fo) / n:.0f}%)")
    by = defaultdict(list)
    for f in fo:
        by[f["r03_first_fy"]].append(f)
    print("   by the R03's first fiscal year (later R01 / any later grant / awards):")
    for y in sorted(by):
        v = by[y]
        print(f"      FY{y}: {sum(f['later_r01'] > 0 for f in v)} / {sum(f['later_any_grant'] > 0 for f in v)} / {len(v)}")
    for f in fo:
        if f["later_any_grant"]:
            print(f"   {f['award']:13} FY{f['r03_first_fy']}  {f['genes_2018_list'][:30]:30} {f['later_grants']}")


    # ---------------------------------------------------------------- timing: each funded gene against itself
    gene_start = {}
    for c, hs in award_genes.items():
        if c in r03fy:
            for h in hs:
                gene_start[h] = min(gene_start.get(h, 9999), r03fy[c])
    ev_rows = []
    random.seed(2018)

    def own_change(m, h, t0):
        """{k: value(t0 + k) - own average over t0-5..t0-1}, or None if fewer than 3 baseline years."""
        pre = [series[m][h].get(t0 + k, 0) for k in PRE_K if 1990 <= t0 + k <= LAST[m]]
        if len(pre) < 3:
            return None, None
        base = sum(pre) / len(pre)
        return base, {k: series[m][h].get(t0 + k, 0) - base for k in EVENT_K if t0 + k <= LAST[m]}

    def boot(vals):
        b = sorted(sum(random.choice(vals) for _ in vals) / len(vals) for _ in range(2000))
        return b[49], b[1949]

    never_by_fam = defaultdict(list)
    for r in rows:
        if not r["funded_any"]:
            never_by_fam[r["family"]].append(r["hgnc_id"])
    for m in EVENT_MEASURES:
        diffs, bases = defaultdict(list), []
        starts = defaultdict(lambda: defaultdict(Counter))   # k -> family -> Counter(start year) of funded genes used
        for r in rows:
            h = r["hgnc_id"]
            if not r["funded_any"] or h not in gene_start or h not in has_data[m]:
                continue
            t0 = gene_start[h]
            base, ch = own_change(m, h, t0)
            if ch is None:
                continue
            bases.append(base)
            for k, v in ch.items():
                diffs[k].append(v)
                starts[k][r["family"]][t0] += 1
        # never-funded: one number per gene per year, using the funded genes' start years in its family
        nchange = defaultdict(lambda: defaultdict(list))     # k -> family -> per-gene values
        memo = {}
        for k in EVENT_K:
            for fam, cnt in starts.get(k, {}).items():
                for g in never_by_fam[fam]:
                    if g not in has_data[m]:
                        continue
                    num = den = 0.0
                    for t0, w in cnt.items():
                        if (g, t0) not in memo:
                            memo[(g, t0)] = own_change(m, g, t0)[1]
                        ch = memo[(g, t0)]
                        if ch is not None and k in ch:
                            num += w * ch[k]
                            den += w
                    if den:
                        nchange[k][fam].append(num / den)

        def never_line(k, resample=False):
            tot = wsum = 0.0
            for fam, cnt in starts[k].items():
                vals = nchange[k].get(fam, [])
                if not vals:
                    continue
                if resample:
                    vals = [random.choice(vals) for _ in vals]
                w = sum(cnt.values())                        # funded genes of this family at year k
                tot += w * sum(vals) / len(vals)
                wsum += w
            return tot / wsum if wsum else None

        for k in EVENT_K:
            vals = diffs.get(k, [])
            if not vals:
                continue
            lo, hi = boot(vals)
            nm = never_line(k)
            if nm is not None:
                nb = sorted(never_line(k, True) for _ in range(2000))
                nlo, nhi = nb[49], nb[1949]
            ev_rows.append({"measure": m, "k": k, "n_genes": len(vals),
                            "n_never": sum(len(v) for v in nchange[k].values()),
                            "mean_baseline": round(sum(bases) / len(bases), 4),
                            "mean_change": round(sum(vals) / len(vals), 4),
                            "lo95": round(lo, 4), "hi95": round(hi, 4),
                            "pct_genes_above_baseline": round(100 * sum(v > 0 for v in vals) / len(vals), 1),
                            "never_mean_change": round(nm, 4) if nm is not None else "",
                            "never_lo95": round(nlo, 4) if nm is not None else "",
                            "never_hi95": round(nhi, 4) if nm is not None else ""})
    print("\nTiming: each funded gene against its own average in the 5 years before funding (year 0 = funding starts)")
    print("   value that year minus the gene's own baseline (measure's own units); never funded = same calendar years,")
    print("   same families, each never-funded gene against its own past")
    for m in EVENT_MEASURES:
        e = [x for x in ev_rows if x["measure"] == m]
        if not e:
            continue
        print(f"   {m} (average baseline {e[0]['mean_baseline']:g} a year):")
        for x in e:
            nv = f"{x['never_mean_change']:+8.3f}" if x["never_mean_change"] != "" else "       -"
            nr = (f"({x['never_lo95']:+.3f} to {x['never_hi95']:+.3f})" if x["never_mean_change"] != "" else "")
            print(f"      year {x['k']:+d}  funded (n={x['n_genes']:3d}) {x['mean_change']:+8.3f} "
                  f"({x['lo95']:+.3f} to {x['hi95']:+.3f})   never funded (n={x['n_never']:3d}) {nv} {nr}")

    # ---------------------------------------------------------------- expected vs actual (ANCOVA)
    fams = sorted({r["family"] for r in rows})
    lg = lambda v, m: math.log10(float(v) + OFFSET[m])

    def design(rs, m, extra=None):
        X = [[1.0, lg(r[f"{m}_before"], m)] + [1.0 if r["family"] == f else 0.0 for f in fams[1:]]
             + ([float(r[extra])] if extra else []) for r in rs]
        return X, [lg(r[f"{m}_after"], m) for r in rs]

    anc = []
    for m in PART_A + PART_B:
        rs = [r for r in rows if r[f"{m}_before"] != ""]
        X, y = design(rs, m)
        _, _, res = ols(X, y)
        for r, e in zip(rs, res):
            r[f"{m}_vs_expected"] = round(e, 3)
        for r in rows:
            r.setdefault(f"{m}_vs_expected", "")
        for label, subset, flag in (("funded vs never", rs, "funded_any"),
                                    ("round-1 funded vs never",
                                     [r for r in rs if r["funded_round1"] or not r["funded_any"]], "funded_round1")):
            X, y = design(subset, m, flag)
            beta, se, _ = ols(X, y)
            b, sb = beta[-1], se[-1]
            p = min(1.0, 2 * norm_sf(abs(b) / sb)) if sb > 0 else 1.0
            pct = lambda v: 100 * (10 ** v - 1)
            anc.append({"comparison": label, "part": "5a" if m in PART_A else "5b", "measure": m,
                        "n_funded": sum(1 for r in subset if r[flag]), "n_never": sum(1 for r in subset if not r[flag]),
                        "slope_before": round(beta[1], 3), "funded_coef_log10": round(b, 4),
                        "pct_above_expected": round(pct(b), 1),
                        "range_95": f"{pct(b - 1.96 * sb):+.1f}% to {pct(b + 1.96 * sb):+.1f}%", "p": p})
    for t, pa in zip(anc, bh([t["p"] for t in anc])):
        t["p_adj"] = pa
    for label in ("funded vs never", "round-1 funded vs never"):
        print(f"\nExpected vs actual ({label}): after allowing for each gene's before and family "
              f"(BH across {len(anc)} tests)")
        print(f"   {'measure':11} {'slope':>6}  {'funded: % above expected':>25}  {'95% range':>22}       p   p_adj")
        for t in anc:
            if t["comparison"] == label:
                print(f"   {t['measure']:11} {t['slope_before']:6.2f}  {t['pct_above_expected']:+24.1f}%  "
                      f"{t['range_95']:>22}  {t['p']:6.3f}  {t['p_adj']:6.3f}")
    print("   slope = how closely 'after' follows 'before' (1 = fully; 0 = not at all)")

    later_gene = {g.strip() for f in fo if f["later_any_grant"] for g in f["genes_2018_list"].split(";")}
    for m, lab in (("jensen", "Jensen score"), ("papers", "focused papers"), ("grants", "NIH grants")):
        top = sorted((r for r in rows if r["funded_any"] and r[f"{m}_vs_expected"] != ""),
                     key=lambda r: -r[f"{m}_vs_expected"])[:10]
        print(f"\nFunded genes furthest above the band: {lab} (vs_expected, log10; before -> after per year)")
        for r in top:
            ties = []
            if r["papers_r03_tied_after"]:
                ties.append(f"{r['papers_r03_tied_after']} R03-tied papers")
            if r["gene"] in later_gene:
                ties.append("later same-PI grant")
            if r["first_structure"] == 1:
                ties.append("first structure")
            print(f"   {r['gene']:10} {r['family']:12} round {r['first_funded_round']}  {r[f'{m}_vs_expected']:+.2f}  "
                  f"({r[f'{m}_before']:g} -> {r[f'{m}_after']:g})  {'; '.join(ties) or '-'}")

    for t in tests:
        t["p"], t["p_adj"] = f"{t['p']:.4g}", f"{t['p_adj']:.4g}"
    OUT.mkdir(parents=True, exist_ok=True)
    for t in anc:
        t["p"], t["p_adj"] = f"{t['p']:.4g}", f"{t['p_adj']:.4g}"
    for name, data in (("gene_change.tsv", rows), ("gene_change_tests.tsv", tests), ("r03_followon_grants.tsv", fo),
                       ("gene_change_ancova.tsv", anc), ("gene_event_time.tsv", ev_rows)):
        with (OUT / name).open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(data)
        print(f"Saved: {OUT / name} ({len(data)} rows)")


if __name__ == "__main__":
    main()
