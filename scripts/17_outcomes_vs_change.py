#!/usr/bin/env python3
"""17_outcomes_vs_change.py - did funded genes whose R03 found something grow more afterwards?

Run from the dark_proteins folder, after scripts 14 and 16 and the outcome labels:
    python3 scripts/17_outcomes_vs_change.py

Genes   the 98 funded genes on the 2018 list (script 14's gene_change.tsv, funded_any = 1).
Labels  data/raw/r03_awards/r03_outcome_genes.csv, one row per gene (checked by hand):
            found_something  yes = an R03 reported a finding about THIS gene
                             (on multi-gene awards, only the genes the finding was about)
            if_yes_what      role, disease link, mechanism, location, binding partners, new ligand, structure
            if_no_why        tools only, negative result, result on another gene of the award, other gene, ...
        Every yes and no counts; the "why" is kept as context (printed per reason, not tested).

Part 1  growth after 2018, yes vs no (each gene's yearly pace, from script 14):
        before = yearly average 2013-2017; after = yearly average 2019 to the measure's last reliable year.
        Items tied to the gene's own R03 (papers, partners, compounds, GO terms, structures citing it) are
        taken out of "after", so a finding's own outputs do not count as growth.
        Regression within the 98 genes, on log10(value + c) as in script 14:
            after ~ before + family + first funding round + found_something
        -> found_something as % above/below the pace expected for a gene with the same before, family and
           round; 95% range; normal-approx. p; Benjamini-Hochberg across the 9 measures.
        Joint test (as script 12): each gene's distance from expected (same regression without
        found_something); P(yes higher) per measure, averaged over the 9; labels shuffled 10,000 times.
Part 2  matching checks, gene by gene: does the kind of finding show up where it should, after the
        R03 started (year = the finding award's round + 1, to the measure's last reliable year)?
            structure         -> new PDB entries
            new ligand        -> new compounds tested (ChEMBL, single protein, valid)
            binding partners  -> new small-scale human interaction partners
            role / mechanism / location -> new informative experimental GO terms
        For each gene: new items in its window (how many tied to the R03), and its chance of any new item if
        the finding made no difference, from other funded genes (without this kind of finding) over the SAME
        calendar years: (1) same family; (2) same family AND same starting study level (third of Jensen pace
        2013-2017), falling back to same level in any family if fewer than 5 genes, and shrunk toward the
        rate of all reference genes (weight 2 genes) so a small 0 is not read as 0%. Summed over genes =
        expected count; observed vs expected with an exact (Poisson-binomial) p on (2), a guide only.

Outputs (data/processed/analysis)
    outcome_vs_change_genes.tsv    one row per gene: label, why/what, distance from expected per measure
    outcome_vs_change_tests.tsv    one row per measure: yes vs no
    outcome_matching.tsv           one row per gene x finding type: new items after the R03 started, chance
    outcome_matching_tests.tsv     one row per finding type: observed vs expected
"""
import csv
import math
import random
import statistics
from collections import Counter, defaultdict
from pathlib import Path

P = Path("data/processed")
OUT = P / "analysis"
LABELS = Path("data/raw/r03_awards/r03_outcome_genes.csv")
AFTER_START = 2019
LAST = {"jensen": 2025, "pubtator": 2021, "r01": 2026, "papers": 2025, "grants": 2026,
        "partners": 2023, "go": 2024, "compounds": 2024, "structures": 2026}
MEASURES = ["jensen", "pubtator", "r01", "papers", "grants", "partners", "go", "compounds", "structures"]
OFFSET = {"jensen": 0.1, "pubtator": 0.1, "r01": 0.05, "papers": 0.2, "grants": 0.05,       # as script 14
          "partners": 0.05, "go": 0.05, "compounds": 0.2, "structures": 0.05}
MATCH = [("structure", "structures", "new PDB entries"),
         ("new ligand", "compounds", "new compounds tested"),
         ("binding partners", "partners", "new small-scale partners"),
         ("role / mechanism / location", "go", "new experimental GO terms")]
N_PERM = 10000


def read_tsv(path):
    if not path.is_file():
        raise SystemExit(f"Missing {path}")
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_labels(path):
    if not path.is_file():
        raise SystemExit(f"Missing {path}")
    for enc in ("utf-8-sig", "cp1252"):          # Excel may save it in Windows encoding
        try:
            with path.open(encoding=enc, newline="") as f:
                return list(csv.DictReader(f))
        except UnicodeDecodeError:
            continue
    raise SystemExit(f"Cannot read {path}: save it as CSV UTF-8")


def yr(v):
    try:
        return int(str(v)[:4])
    except ValueError:
        return None


# ---------------------------------------------------------------- statistics (as scripts 12 and 14)
def norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2))


def midranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r, i = [0.0] * len(v), 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def p_higher(x, y):
    """P(a value from x is higher than one from y), ties half."""
    r = midranks(x + y)
    return (sum(r[:len(x)]) - len(x) * (len(x) + 1) / 2) / (len(x) * len(y))


def solve(a, b):
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
    k = len(X[0])
    xtx = [[sum(r[i] * r[j] for r in X) for j in range(k)] for i in range(k)]
    xty = [sum(r[i] * v for r, v in zip(X, y)) for i in range(k)]
    beta = solve(xtx, xty)
    res = [v - sum(b * x for b, x in zip(beta, r)) for r, v in zip(X, y)]
    s2 = sum(e * e for e in res) / (len(y) - k)
    inv = [solve(xtx, [1.0 if i == j else 0.0 for i in range(k)]) for j in range(k)]
    return beta, [math.sqrt(s2 * inv[i][i]) for i in range(k)], res


def bh(ps):
    m = len(ps)
    order = sorted(range(m), key=lambda i: ps[i])
    adj, prev = [0.0] * m, 1.0
    for rank, i in reversed(list(enumerate(order, 1))):
        prev = min(prev, ps[i] * m / rank)
        adj[i] = prev
    return adj


def med(v):
    return statistics.median(v) if v else float("nan")


# ---------------------------------------------------------------- main
def main():
    change = [r for r in read_tsv(P / "analysis/gene_change.tsv") if r["funded_any"] == "1"]
    labels = {r["gene"].strip(): r for r in read_labels(LABELS)}
    missing = [r["gene"] for r in change if r["gene"] not in labels]
    if missing:
        raise SystemExit(f"No outcome label for: {', '.join(missing)}")
    genes = []
    for r in change:
        lab = labels[r["gene"]]
        g = {"gene": r["gene"], "hgnc_id": r["hgnc_id"], "family": r["family"],
             "round": r["first_funded_round"], "found": lab["found_something"].strip().lower() == "yes",
             "what": lab["if_yes_what"].strip(), "why": lab["if_no_why"].strip(),
             "yes_round": min((int(x) for x, a in zip(lab["rounds"].split(";"), lab["awards"].split(";"))
                               if a.strip() in lab["yes_award"]), default=None)}
        for m in MEASURES:
            if r[f"{m}_before"] == "":
                g[f"{m}_before"] = g[f"{m}_after"] = None
                continue
            span = LAST[m] - AFTER_START + 1
            tied = float(r.get(f"{m}_r03_tied_after") or 0)
            g[f"{m}_before"] = float(r[f"{m}_before"])
            g[f"{m}_after"] = max(0.0, float(r[f"{m}_after"]) - tied / span)     # R03-tied items taken out
            g[f"{m}_tied"] = tied
        genes.append(g)
    yes = [g for g in genes if g["found"]]
    no = [g for g in genes if not g["found"]]
    print(f"Funded genes on the 2018 list: {len(genes)}  (found something {len(yes)}, nothing {len(no)})")
    print("   why not:", ", ".join(f"{k} {v}" for k, v in
                                    Counter(t.strip() for g in no for t in g["why"].split(";") if t.strip()).most_common()))
    print("   what:   ", ", ".join(f"{k} {v}" for k, v in
                                    Counter(t.strip() for g in yes for t in g["what"].split(";") if t.strip()).most_common()))

    # ------------------------------------------------------------ part 1: growth, yes vs no
    fams = sorted({g["family"] for g in genes})[1:]
    rounds = sorted({g["round"] for g in genes})[1:]
    lg = lambda v, m: math.log10(v + OFFSET[m])

    def design(rs, m, with_found):
        X = [[1.0, lg(g[f"{m}_before"], m)] + [1.0 if g["family"] == f else 0.0 for f in fams]
             + [1.0 if g["round"] == k else 0.0 for k in rounds] + ([1.0 if g["found"] else 0.0] if with_found else [])
             for g in rs]
        return X, [lg(g[f"{m}_after"], m) for g in rs]

    tests, resid = [], {}
    pct = lambda v: 100 * (10 ** v - 1)
    for m in MEASURES:
        rs = [g for g in genes if g[f"{m}_before"] is not None]
        X, y = design(rs, m, True)
        beta, se, _ = ols(X, y)
        b, s = beta[-1], se[-1]
        _, _, res = ols(*design(rs, m, False))
        for g, e in zip(rs, res):
            g[f"{m}_vs_expected"] = round(e, 4)
        resid[m] = rs
        ry = [g for g in rs if g["found"]]
        rn = [g for g in rs if not g["found"]]
        tests.append({"measure": m, "n_yes": len(ry), "n_no": len(rn),
                      "yes_before": round(med([g[f"{m}_before"] for g in ry]), 3),
                      "yes_after": round(med([g[f"{m}_after"] for g in ry]), 3),
                      "no_before": round(med([g[f"{m}_before"] for g in rn]), 3),
                      "no_after": round(med([g[f"{m}_after"] for g in rn]), 3),
                      "p_yes_higher": round(p_higher([g[f"{m}_vs_expected"] for g in ry],
                                                     [g[f"{m}_vs_expected"] for g in rn]), 3),
                      "yes_vs_expected_pct": round(pct(b), 1),
                      "lo95": round(pct(b - 1.96 * s), 1), "hi95": round(pct(b + 1.96 * s), 1),
                      "p": min(1.0, 2 * norm_sf(abs(b / s))) if s > 0 else 1.0,
                      "r03_tied_items_removed": round(sum(g.get(f"{m}_tied", 0) for g in rs))})
    for t, a in zip(tests, bh([t["p"] for t in tests])):
        t["p"], t["p_adj"] = f"{t['p']:.3g}", f"{a:.3g}"

    # joint test: mean P(yes higher) over the 9 measures, labels shuffled
    def joint(flags):
        ps = []
        for m in MEASURES:
            rs = resid[m]
            x = [g[f"{m}_vs_expected"] for g in rs if flags[g["gene"]]]
            y = [g[f"{m}_vs_expected"] for g in rs if not flags[g["gene"]]]
            ps.append(p_higher(x, y))
        return sum(ps) / len(ps), sum(p > 0.5 for p in ps)

    obs = {g["gene"]: g["found"] for g in genes}
    j_obs, k_obs = joint(obs)
    random.seed(17)
    names, vals = list(obs), list(obs.values())
    hits = 0
    for _ in range(N_PERM):
        random.shuffle(vals)
        jm, _ = joint(dict(zip(names, vals)))
        hits += abs(jm - 0.5) >= abs(j_obs - 0.5) - 1e-12
    p_joint = (hits + 1) / (N_PERM + 1)

    print("\nPart 1. Yearly pace after 2018 (R03-tied items taken out): found something vs nothing")
    print("   'vs expected' = yes genes' pace compared with no genes of the same before, family and round\n")
    print(f"   {'':11} {'':8}   {'median yearly pace, before -> after':^62}")
    print(f"   {'measure':11} {'yes/no':>7}   {'yes':^28}   {'no':^18}   "
          f"{'P(yes higher)':>13}   {'yes vs expected (95%)':<26}{'p_adj':>6}")
    for t in tests:
        print(f"   {t['measure']:11} {t['n_yes']:3}/{t['n_no']:<3}   "
              f"{t['yes_before']:12.2f} -> {t['yes_after']:<12.2f}   {t['no_before']:7.2f} -> {t['no_after']:<7.2f}   "
              f"{t['p_yes_higher']:13.2f}   {t['yes_vs_expected_pct']:+7.1f}% ({t['lo95']:+.0f} to {t['hi95']:+.0f})"
              f"{'':6}{t['p_adj']:>6}")
    print(f"\n   Joint test (9 measures): mean P(yes higher) = {j_obs:.3f}, higher in {k_obs} of 9; "
          f"label-shuffle p = {p_joint:.3f} ({N_PERM:,} shuffles)")
    tests.append({"measure": "JOINT (mean P(yes higher), 9 measures)", "n_yes": len(yes), "n_no": len(no),
                  "p_yes_higher": round(j_obs, 3), "p": f"{p_joint:.3g}", "p_adj": ""})

    print("\n   Context, the 'no' genes by reason: median distance from expected (log10; + = faster than expected)")
    reasons = sorted({t.strip() for g in no for t in g["why"].split(";") if t.strip()})
    print(f"   {'reason':38} {'n':>3}  " + "  ".join(f"{m[:8]:>8}" for m in MEASURES))
    toks = lambda g: {t.strip() for t in g["why"].split(";") if t.strip()}      # whole labels, not substrings
    for grp, rows_ in [("found something (all yes)", yes)] + [(r, [g for g in no if r in toks(g)]) for r in reasons]:
        cells = []
        for m in MEASURES:
            v = [g[f"{m}_vs_expected"] for g in rows_ if g.get(f"{m}_vs_expected") is not None]
            cells.append(f"{med(v):+8.2f}" if v else f"{'-':>8}")
        print(f"   {grp:38} {len(rows_):3}  " + "  ".join(cells))

    # ------------------------------------------------------------ part 2: matching checks
    ids = {g["hgnc_id"]: g for g in genes}
    first_round = {g["hgnc_id"]: int(g["round"]) for g in genes if g["round"]}
    items = defaultdict(lambda: defaultdict(dict))        # measure -> gene -> item -> (first year, tied)

    def add(m, h, key, y, t):
        old = items[m][h].get(key)
        if old is None or y < old[0]:
            items[m][h][key] = (y, t)
        elif y == old[0] and t:
            items[m][h][key] = (y, True)

    for r in read_tsv(P / "structure/gene_pdb_entries.tsv"):
        h, y = r["hgnc_id"], yr(r["year"])
        if h in ids and y:
            add("structures", h, r["pdb_id"], y, bool(r["r03_for_this_gene"]))
    for r in read_tsv(P / "chemistry/gene_activities.tsv"):
        h, y = r["hgnc_id"], yr(r["year"])
        if h in ids and y and r["target_type"] == "SINGLE PROTEIN" and r["valid"] == "1":
            add("compounds", h, r["compound"], y, bool(r["r03_for_this_gene"]))
    for r in read_tsv(P / "interactions/gene_interactions.tsv"):
        h, y = r["hgnc_id"], yr(r["year"])
        if (h in ids and y and r["match_type"] == "same gene" and r["self_interaction"] == "0"
                and "low" in r["throughput"].lower()):
            add("partners", h, r["partner_gene_id"], y, bool(r["r03_for_this_gene"]))
    for r in read_tsv(P / "go/gene_go_annotations.tsv"):
        h, y = r["hgnc_id"], yr(r["year"])
        if (h in ids and y and r["evidence_group"] == "experimental" and r["informative"] == "1"
                and r["negated"] == "0"):
            add("go", h, r["go_id"], y, bool(r["r03_for_this_gene"]))

    def new_since(m, h, start):
        got = [(y, t) for y, t in items[m][h].values() if start <= y <= LAST[m]]
        return len(got), sum(t for _, t in got)

    def poisson_binomial_p(ps, k):
        """Two-sided exact p for observing k successes when each gene succeeds with its own chance."""
        dist = [1.0]
        for p_ in ps:
            nxt = [0.0] * (len(dist) + 1)
            for i, d in enumerate(dist):
                nxt[i] += d * (1 - p_)
                nxt[i + 1] += d * p_
            dist = nxt
        return min(1.0, sum(d for d in dist if d <= dist[k] * (1 + 1e-9)))

    # starting study level: thirds of the Jensen pace before 2018 (2013-2017), across the 98 genes
    jb = sorted(g["jensen_before"] for g in genes if g["jensen_before"] is not None)
    cut1, cut2 = jb[len(jb) // 3], jb[2 * len(jb) // 3]
    level = lambda g: ("?" if g["jensen_before"] is None else
                       "low" if g["jensen_before"] < cut1 else "mid" if g["jensen_before"] < cut2 else "high")
    MIN_REF, PRIOR_W = 5, 2      # smallest reference group; weight (in genes) of the overall rate, so 0 of 5 is not "0%"

    match_rows, match_tests = [], []
    print("\nPart 2. Matching checks: does the kind of finding show up after the R03 started?")
    print("   Each gene's chance of any new item if the finding made no difference, from other funded genes")
    print("   (without this kind of finding) over the SAME calendar years:")
    print("     family        = same family (as before)")
    print(f"     family+level  = same family AND same third of Jensen pace before 2018 (low < {cut1:.2f} <= mid < "
          f"{cut2:.2f} <= high);")
    print(f"                     if fewer than {MIN_REF} such genes, same level in any family; shrunk toward the rate of")
    print(f"                     all reference genes over those years with a weight of {PRIOR_W} genes, so a small 0 is not 0%.")
    print("   Expected = chances summed over genes. Exact p uses the family+level chances (guide only).")
    for kind, m, what in MATCH:
        kinds = kind.split(" / ")
        with_kind = [g for g in genes if g["found"] and any(k in g["what"] for k in kinds)]
        pool = [g for g in genes if g not in with_kind]
        print(f"\n   {kind} -> {what} (from the finding award's round + 1 to {LAST[m]})")
        print(f"      {'gene':10} {'family':12} {'level':5} {'years':9}  {'new':>4}  {'tied':>4}   "
              f"{'chance: family':>14}   {'family+level':>12}")
        hit, ch_fam, ch_lev = 0, [], []
        for g in sorted(with_kind, key=lambda g: g["gene"]):
            start = (g["yes_round"] or first_round[g["hgnc_id"]]) + 1
            n, t = new_since(m, g["hgnc_id"], start)
            got = lambda o: new_since(m, o["hgnc_id"], start)[0] > 0
            same = [o for o in pool if o["family"] == g["family"]] or pool
            c_fam = sum(map(got, same)) / len(same)
            lev = [o for o in same if level(o) == level(g)]
            note = g["family"]
            if len(lev) < MIN_REF:
                lev, note = [o for o in pool if level(o) == level(g)], "any family"
            prior = sum(map(got, pool)) / len(pool)
            c_lev = (sum(map(got, lev)) + PRIOR_W * prior) / (len(lev) + PRIOR_W)
            hit += n > 0
            ch_fam.append(c_fam)
            ch_lev.append(c_lev)
            print(f"      {g['gene']:10} {g['family']:12} {level(g):5} {start}-{str(LAST[m])[2:]:<4}  {n:4d}  "
                  f"{(str(t) if t else ''):>4}   {100 * c_fam:13.0f}%   {100 * c_lev:11.0f}%  ({len(lev)} {note})")
            match_rows.append({"finding": kind, "check": what, "gene": g["gene"], "family": g["family"],
                               "jensen_level": level(g), "from_year": start, "to_year": LAST[m], "new_items": n,
                               "r03_tied": t, "any_new": int(n > 0), "chance_family": round(c_fam, 3),
                               "chance_family_level": round(c_lev, 3), "reference_genes_level": len(lev),
                               "reference_note": note})
        if with_kind:
            p = poisson_binomial_p(ch_lev, hit)
            print(f"      -> observed {hit} of {len(with_kind)} with any new item;  expected {sum(ch_fam):.1f} (family), "
                  f"{sum(ch_lev):.1f} (family+level);  exact p = {p:.2g} (guide only)")
            match_tests.append({"finding": kind, "check": what, "n_genes": len(with_kind), "observed": hit,
                                "expected_family": round(sum(ch_fam), 2),
                                "expected_family_level": round(sum(ch_lev), 2), "p_exact_family_level": f"{p:.3g}"})

    # ------------------------------------------------------------ outputs
    OUT.mkdir(parents=True, exist_ok=True)
    gene_rows = [{"gene": g["gene"], "hgnc_id": g["hgnc_id"], "family": g["family"], "first_round": g["round"],
                  "found_something": "yes" if g["found"] else "no", "if_yes_what": g["what"], "if_no_why": g["why"],
                  **{f"{m}_vs_expected": g.get(f"{m}_vs_expected", "") for m in MEASURES}} for g in genes]
    for name, data in (("outcome_vs_change_genes.tsv", gene_rows), ("outcome_vs_change_tests.tsv", tests),
                       ("outcome_matching.tsv", match_rows), ("outcome_matching_tests.tsv", match_tests)):
        if not data:
            continue
        cols = list(dict.fromkeys(k for d in data for k in d))
        with (OUT / name).open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", restval="")
            w.writeheader()
            w.writerows(data)
        print(f"Wrote {OUT / name}  ({len(data)} rows)")


if __name__ == "__main__":
    main()
