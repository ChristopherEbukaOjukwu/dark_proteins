"""
15_list_changes.py - genes removed from the IDG list vs genes that stayed

Question
    Did genes removed from IDG eligibility gain more R01s, Jensen score and PubTator score
    than genes that stayed on the list?
    Descriptive only: we know what changed before a gene left, not why NIH removed it.

List changes and windows (the 3 full calendar years before the new list was posted)
    2018 -> 2019   posted 2019-07   window 2016-2018
    2019 -> 2020   identical lists  (checked, nothing to compare)
    2020 -> 2021   posted 2021-04   window 2018-2020
    2021 -> 2022   posted 2022-04   window 2019-2021

Measures, per gene on the earlier list
    R01        new R01 projects (R01, R37, RF1 = script 06's R01 group) starting in the window
               got_r01 = at least one; first_r01 = none ever before the window, at least one in it
               (IDG R03s are not R01s, so they never count)
    Jensen     gain in window   = score added in the 3-year window (sum of yearly scores)
                                  = how much the gene's total Jensen score went up
               change           = gain in window - gain in the 3 years before it (e.g. 2016-18 minus 2013-15)
                                  = did attention speed up, not just continue
               before window    = total score from all years before the window (starting level)
    PubTator   same as Jensen
    A gene with no Jensen (or PubTator) record at all is missing (blank), not 0, and is left out of that test.
    Within a gene with records, a year with no record counts as 0 (no papers that year).

Tests (removed vs stayed, per list change)
    got R01 / first R01   Fisher exact, difference in % with Newcombe 95% range
    Jensen / PubTator     Mann-Whitney P(removed higher), medians: gain, change, starting level,
                          growth in proportion = log10((end + 1) / (start + 1)) over the window
    Allowing for start    OLS: log10(end + 1) ~ log10(start + 1) + removed + family; the removed coefficient
                          as % above/below the end level expected from start and family (normal-approx. p)
    Benjamini-Hochberg across all tests

Inputs  (data/processed)
    r03_awards/gene_summary.tsv        gene, hgnc_id, family, listed_years, funded_years
    attention/gene_year_scores.tsv     yearly Jensen and PubTator scores
    funding/gene_grants.tsv            grants per gene (script 06)

Outputs (data/processed/analysis)  [+ list_change_family_tests.tsv: same tests within each family, where >= 5 removed and >= 5 stayed]
    list_change_genes.tsv   one row per gene per list change (removed or stayed), with all measures
    list_change_tests.tsv   one row per list change x measure
"""
import csv
import math
from collections import defaultdict
from pathlib import Path

P = Path("data/processed")
OUT = P / "analysis"
CHANGES = [  # (old list, new list, window)
    (2018, 2019, range(2016, 2019)),
    (2019, 2020, range(2017, 2020)),
    (2020, 2021, range(2018, 2021)),
    (2021, 2022, range(2019, 2022)),
]


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


# ---------------------------------------------------------------- tests (same as scripts 12 and 14)
def norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2))


def mann_whitney(x, y):
    """Two-sided p and P(x > y) (ties count half)."""
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


def solve(a, b):
    """Solve a x = b (small square system) by Gauss-Jordan elimination (as in script 14)."""
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
    """Coefficients and their standard errors (as in script 14)."""
    k = len(X[0])
    xtx = [[sum(r[i] * r[j] for r in X) for j in range(k)] for i in range(k)]
    xty = [sum(r[i] * v for r, v in zip(X, y)) for i in range(k)]
    beta = solve(xtx, xty)
    res = [v - sum(b * x for b, x in zip(beta, r)) for r, v in zip(X, y)]
    s2 = sum(e * e for e in res) / (len(y) - k)
    inv = [solve(xtx, [1.0 if i == j else 0.0 for i in range(k)]) for j in range(k)]
    return beta, [math.sqrt(s2 * inv[i][i]) for i in range(k)]


def med(v):
    v = sorted(v)
    n = len(v)
    return (v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2) if n else float("nan")


# ---------------------------------------------------------------- main
def main():
    genes = {r["hgnc_id"]: r for r in read_tsv(P / "r03_awards/gene_summary.tsv")}
    listed = {h: {int(y) for y in g["listed_years"].split(";") if y} for h, g in genes.items()}
    funded = {h: sorted(int(y) for y in g["funded_years"].split(";") if y) for h, g in genes.items()}

    jensen, pubtator = defaultdict(dict), defaultdict(dict)        # gene -> year -> yearly score
    for r in read_tsv(P / "attention/gene_year_scores.tsv"):
        h, y = r["hgnc_id"], int(r["year"])
        if h in genes:
            if r["jensen_score"] != "":
                jensen[h][y] = jensen[h].get(y, 0.0) + float(r["jensen_score"])
            if r["pubtator_score"] != "":
                pubtator[h][y] = pubtator[h].get(y, 0.0) + float(r["pubtator_score"])

    r01_start = defaultdict(dict)                                   # gene -> R01 core -> first fiscal year
    for r in read_tsv(P / "funding/gene_grants.tsv"):
        h, y = r["hgnc_id"], yr(r["fiscal_year"])
        if h not in genes or not y or r["grant_group"] != "R01" or r.get("family_name_only") == "1":
            continue
        c = r["core_project_num"]
        if c not in r01_start[h] or y < r01_start[h][c]:
            r01_start[h][c] = y

    rows, tests, groups = [], [], []
    for old, new, win in CHANGES:
        on_old = [h for h in genes if old in listed[h]]
        removed = [h for h in on_old if new not in listed[h]]
        stayed = [h for h in on_old if new in listed[h]]
        added = [h for h in genes if new in listed[h] and old not in listed[h]]
        label = f"{old}->{new}"
        print(f"\n{label}  window {win.start}-{win.stop - 1}:  on {old} list {len(on_old)}, "
              f"removed {len(removed)}, stayed {len(stayed)}, added {len(added)}")
        if not removed or not stayed:
            print("   nothing to compare (no genes removed, or none stayed)")
            continue

        rec = {}
        for h in on_old:
            ys = list(r01_start[h].values())
            r01_before = sum(y < win.start for y in ys)
            r01_win = sum(y in win for y in ys)
            back = sorted(y for y in listed[h] if y > new)
            rec[h] = {
                "change": label, "status": "removed" if h in removed else "stayed",
                "gene": genes[h]["gene"], "hgnc_id": h, "family": genes[h]["family"],
                "funded_before_change": "yes" if any(y <= old for y in funded[h]) else "no",
                "funded_years": genes[h]["funded_years"],
                "back_on_list_later": ";".join(map(str, back)) if h in removed else "",
                "r01_before_window": r01_before, "r01_in_window": r01_win,
                "got_r01": int(r01_win > 0), "first_r01": int(r01_before == 0 and r01_win > 0),
            }
            earlier = range(win.start - 3, win.start)
            for m, d in (("jensen", jensen), ("pubtator", pubtator)):
                if h not in d:                               # no record at all: missing, not 0
                    for k in ("before_window", "earlier_3y", "gain_in_window", "change", "growth_log"):
                        rec[h][f"{m}_{k}"] = None
                    continue
                gain = sum(v for y, v in d[h].items() if y in win)
                prev = sum(v for y, v in d[h].items() if y in earlier)
                rec[h][f"{m}_before_window"] = round(sum(v for y, v in d[h].items() if y < win.start), 3)
                rec[h][f"{m}_earlier_3y"] = round(prev, 3)
                rec[h][f"{m}_gain_in_window"] = round(gain, 3)
                rec[h][f"{m}_change"] = round(gain - prev, 3)
                start = rec[h][f"{m}_before_window"]
                # growth in proportion over the window: log10 of (end + 1) / (start + 1); +1 keeps zeros usable
                rec[h][f"{m}_growth_log"] = round(math.log10((start + gain + 1) / (start + 1)), 5)
            rows.append(rec[h])
        R, S = [rec[h] for h in removed], [rec[h] for h in stayed]

        def prop_test(col, what):
            a, c = sum(r[col] for r in R), sum(r[col] for r in S)
            p = fisher_exact(a, len(R) - a, c, len(S) - c)
            lo, hi = newcombe(a, len(R), c, len(S))
            tests.append({"change": label, "measure": what, "n_removed": len(R), "n_stayed": len(S),
                          "removed": f"{a}/{len(R)} ({100 * a / len(R):.0f}%)",
                          "stayed": f"{c}/{len(S)} ({100 * c / len(S):.0f}%)",
                          "effect": f"{100 * (a / len(R) - c / len(S)):+.0f} pts ({100 * lo:+.0f} to {100 * hi:+.0f})",
                          "p": p})

        def rank_test(col, what):
            x = [r[col] for r in R if r[col] is not None]
            y = [r[col] for r in S if r[col] is not None]
            p, pr = mann_whitney(x, y)
            tests.append({"change": label, "measure": what, "n_removed": len(x), "n_stayed": len(y),
                          "removed": f"median {med(x):.2f}", "stayed": f"median {med(y):.2f}",
                          "effect": f"P(removed higher) = {pr:.2f}", "p": p})

        prop_test("got_r01", "got an R01 in window")
        prop_test("first_r01", "first-ever R01 in window")
        rank_test("jensen_gain_in_window", "Jensen gained in window")
        rank_test("pubtator_gain_in_window", "PubTator gained in window")
        rank_test("jensen_change", "Jensen change (window - earlier 3y)")
        rank_test("pubtator_change", "PubTator change (window - earlier 3y)")
        rank_test("jensen_growth_log", "Jensen growth in proportion")
        rank_test("pubtator_growth_log", "PubTator growth in proportion")

        def ancova(m, what):
            """log10(end + 1) ~ log10(start + 1) + removed + family: is a removed gene's end level above what
            its start level (and family) predicts? Effect shown as % above or below expected, 95% range."""
            use = [r for r in R + S if r[f"{m}_before_window"] is not None]
            fams = sorted({r["family"] for r in use})[1:]
            X, y = [], []
            for r in use:
                a = r[f"{m}_before_window"]
                X.append([1.0, math.log10(a + 1), 1.0 if r["status"] == "removed" else 0.0]
                         + [1.0 if r["family"] == f else 0.0 for f in fams])
                y.append(math.log10(a + r[f"{m}_gain_in_window"] + 1))
            beta, se = ols(X, y)
            b, e = beta[2], se[2]
            pct = lambda v: 100 * (10 ** v - 1)
            tests.append({"change": label, "measure": what,
                          "n_removed": sum(r["status"] == "removed" for r in use),
                          "n_stayed": sum(r["status"] == "stayed" for r in use),
                          "removed": "", "stayed": "",
                          "effect": f"removed {pct(b):+.1f}% vs expected ({pct(b - 1.96 * e):+.1f} to {pct(b + 1.96 * e):+.1f})",
                          "p": min(1.0, 2 * norm_sf(abs(b / e)))})

        ancova("jensen", "Jensen end level, allowing for start + family")
        ancova("pubtator", "PubTator end level, allowing for start + family")
        rank_test("jensen_before_window", "Jensen before window (start level)")
        rank_test("pubtator_before_window", "PubTator before window (start level)")
        rank_test("r01_before_window", "R01s before window (start level)")

        print(f"   removed genes: {', '.join(sorted(r['gene'] for r in R))}")
        groups.append((label, R, S))

    for t, a in zip(tests, bh([t["p"] for t in tests])):
        t["p"], t["p_adj"] = f"{t['p']:.3g}", f"{a:.3g}"

    print(f"\n{'change':9} {'measure':38} {'removed':>16} {'stayed':>16}  {'effect':32} {'p':>7} {'p_adj':>7}")
    for t in tests:
        print(f"{t['change']:9} {t['measure']:38} {t['removed']:>16} {t['stayed']:>16}  "
              f"{t['effect']:32} {t['p']:>7} {t['p_adj']:>7}")

    f1 = lambda v: "NA" if v is None else f"{v:.1f}"
    print("\nRemoved genes, one line each\n   R01s before/in window;  Jensen and PubTator: before / earlier 3y / window 3y")
    for r in sorted((r for r in rows if r["status"] == "removed"), key=lambda r: (r["change"], r["family"], r["gene"])):
        print(f"   {r['change']}  {r['gene']:10} {r['family']:12} funded={r['funded_before_change']:3} "
              f"R01 {r['r01_before_window']}/{r['r01_in_window']}  "
              f"Jensen {f1(r['jensen_before_window'])}/{f1(r['jensen_earlier_3y'])}/{f1(r['jensen_gain_in_window'])}  "
              f"PubTator {f1(r['pubtator_before_window'])}/{f1(r['pubtator_earlier_3y'])}/{f1(r['pubtator_gain_in_window'])}"
              + (f"  back on list {r['back_on_list_later']}" if r["back_on_list_later"] else ""))

    # ------------------------------------------------------------ within each family
    # Removals came partly in family blocks (e.g. many ion channels in 2019). Same comparison inside each
    # family, so a family mix cannot make the overall gap. Tested only where a family has >= MIN_N removed
    # and >= MIN_N stayed genes; smaller groups are counted but not tested.
    MIN_N = 5
    fam_tests = []
    print("\nWithin each family")
    for label, R, S in groups:
        for fam in sorted({r["family"] for r in R + S}):
            Rf = [r for r in R if r["family"] == fam]
            Sf = [r for r in S if r["family"] == fam]
            n_on = len(Rf) + len(Sf)
            print(f"   {label}  {fam:12} on list {n_on:3}, removed {len(Rf):3} ({100 * len(Rf) / n_on:.0f}%)")
            if len(Rf) < MIN_N or len(Sf) < MIN_N:
                continue
            for col, what in (("jensen_before_window", "Jensen before window (start level)"),
                              ("pubtator_before_window", "PubTator before window (start level)"),
                              ("jensen_gain_in_window", "Jensen gained in window"),
                              ("pubtator_gain_in_window", "PubTator gained in window"),
                              ("jensen_change", "Jensen change (window - earlier 3y)"),
                              ("pubtator_change", "PubTator change (window - earlier 3y)"),
                              ("r01_before_window", "R01s before window (start level)")):
                x = [r[col] for r in Rf if r[col] is not None]
                y = [r[col] for r in Sf if r[col] is not None]
                p, pr = mann_whitney(x, y)
                fam_tests.append({"change": label, "family": fam, "measure": what,
                                  "n_removed": len(x), "n_stayed": len(y),
                                  "removed": f"median {med(x):.2f}", "stayed": f"median {med(y):.2f}",
                                  "effect": f"P(removed higher) = {pr:.2f}", "p": p})
            a, c = sum(r["got_r01"] for r in Rf), sum(r["got_r01"] for r in Sf)
            lo, hi = newcombe(a, len(Rf), c, len(Sf))
            fam_tests.append({"change": label, "family": fam, "measure": "got an R01 in window",
                              "n_removed": len(Rf), "n_stayed": len(Sf),
                              "removed": f"{a}/{len(Rf)} ({100 * a / len(Rf):.0f}%)",
                              "stayed": f"{c}/{len(Sf)} ({100 * c / len(Sf):.0f}%)",
                              "effect": f"{100 * (a / len(Rf) - c / len(Sf)):+.0f} pts ({100 * lo:+.0f} to {100 * hi:+.0f})",
                              "p": fisher_exact(a, len(Rf) - a, c, len(Sf) - c)})
    for t, a in zip(fam_tests, bh([t["p"] for t in fam_tests])):
        t["p"], t["p_adj"] = f"{t['p']:.3g}", f"{a:.3g}"
    print(f"\n{'change':9} {'family':12} {'measure':38} {'removed':>16} {'stayed':>16}  {'effect':32} {'p':>7} {'p_adj':>7}")
    for t in fam_tests:
        print(f"{t['change']:9} {t['family']:12} {t['measure']:38} {t['removed']:>16} {t['stayed']:>16}  "
              f"{t['effect']:32} {t['p']:>7} {t['p_adj']:>7}")

    OUT.mkdir(parents=True, exist_ok=True)
    for name, data in (("list_change_genes.tsv", rows), ("list_change_tests.tsv", tests),
                       ("list_change_family_tests.tsv", fam_tests)):
        if not data:
            continue
        with (OUT / name).open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0]), delimiter="\t")
            w.writeheader()
            w.writerows(data)
        print(f"Wrote {OUT / name}  ({len(data)} rows)")


if __name__ == "__main__":
    main()
