#!/usr/bin/env python3
"""Analysis 1: where did the 2018-list genes stand before the announcement (up to the end of 2017)?

Run from the dark_proteins folder, after steps 3-11:
    python3 scripts/12_state_2017.py

Genes: the 328 on the 2018 list, marked later funded (any round) or never funded.
"Before" = up to and including 2017 (the 2018 list was posted in July 2018).

Per gene (one row each):
    attention    jensen_total_2017, jensen_mean_2013_2017
    papers       papers_total_2017, papers_mean_2013_2017   (human + one-to-one animal matches, any paper)
                 papers_le10_total_2017, papers_le10_mean_2013_2017  (only papers tagging 10 or fewer genes)
    funding      grants_total_2017, r01_total_2017           (gene-specific matches; IDG R03s excluded)
    interactions partners_all_2017, partners_small_2017     (human; small = low-throughput evidence)
    chemistry    compounds_2017, potent_by_2017             (single-protein, valid; IDG Tchem cutoffs)
    function     go_terms_2017                               (experimental, informative, not negated)
    structure    structures_2017
    disease      disease_modplus_submitted_by_2017, disease_modplus_paper_by_2017
    IDG rule     idg_jensen_under_50 (Jensen total < 50), idg_no_r01 (no R01 to 2017), idg_meets_both
Blank = no data for that gene in that source (e.g. pseudogenes have no Jensen score).

Tests (later funded vs never funded), overall and within each family:
    counts/scores  Mann-Whitney U (two-sided, normal approximation with tie and continuity correction);
                   effect = chance a random funded gene is higher than a random unfunded one (0.5 = no difference)
    yes/no         Fisher's exact test (two-sided); effect = difference in % with a 95% range (Newcombe)
    all p-values   adjusted together with Benjamini-Hochberg (p_adj)
    joint test     the 8 measures of Figure 2 together (JOINT): does the funded group lean one way across them?
                   effect per measure = P(funded higher) (for yes/no measures too: ties count half)
                   statistics = mean of the 8 effects, and how many of the 8 are above 0.5
                   p = permutation: funded labels shuffled across genes 10,000 times; each gene keeps all its
                   values, so the measures stay correlated as they really are. Two-sided. Not part of the BH set.

Outputs (data/processed/analysis/)
    gene_state_2017.tsv          one row per gene
    gene_state_2017_summary.tsv  per group (later funded / never funded), overall and by family:
                                 median of each count, % of genes with at least one
    gene_state_2017_tests.tsv    one row per group x measure: values, effect, 95% range, p, p_adj
"""
import csv
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path

P = Path("data/processed")
OUT = P / "analysis"
END = 2017
WINDOW = range(2013, 2018)
JOINT = ["papers_le10_total_2017", "grants_total_2017", "partners_small_2017", "go_terms_2017",
         "compounds_2017", "potent_by_2017", "structures_2017", "disease_modplus_paper_by_2017"]
N_PERM = 10000


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


# ---------------------------------------------------------------- tests (no extra packages needed)
def norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2))


def mann_whitney(x, y):
    """Two-sided p (normal approx., tie + continuity correction) and P(random x > random y) + 0.5 * ties."""
    n1, n2 = len(x), len(y)
    allv = sorted((v, i) for i, v in enumerate(x + y))
    ranks, i = [0.0] * (n1 + n2), 0
    ties = 0.0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[allv[k][1]] = r
        t = j - i + 1
        ties += t ** 3 - t
        i = j + 1
    u1 = sum(ranks[:n1]) - n1 * (n1 + 1) / 2
    mu = n1 * n2 / 2
    n = n1 + n2
    sd = math.sqrt(n1 * n2 / 12 * ((n + 1) - ties / (n * (n - 1)))) if n > 1 else 0
    if sd == 0:
        return 1.0, u1 / (n1 * n2)
    z = (abs(u1 - mu) - 0.5) / sd
    return min(1.0, 2 * norm_sf(max(z, 0))), u1 / (n1 * n2)


def fisher_exact(a, b, c, d):
    """Two-sided Fisher p for [[a, b], [c, d]] (sum of tables no more likely than the observed one)."""
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
    """95% range for p1 - p2 (Newcombe hybrid score)."""
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


def midranks(vals):
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    r, i = [0.0] * len(vals), 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def joint_test(rows):
    """Mean P(funded higher) over JOINT, and count above 0.5, with a label-shuffling permutation p."""
    flags = [1 if r["funded_any"] else 0 for r in rows]
    prep = []                                   # per measure: (gene indices with a value, their midranks)
    for m in JOINT:
        idx = [i for i, r in enumerate(rows) if r[m] != ""]
        prep.append((idx, midranks([float(rows[i][m]) for i in idx])))

    def effects(fl):
        out = []
        for idx, rk in prep:
            n1 = sum(fl[i] for i in idx)
            n0 = len(idx) - n1
            s1 = sum(r for i, r in zip(idx, rk) if fl[i])
            out.append((s1 - n1 * (n1 + 1) / 2) / (n1 * n0))
        return out

    obs = effects(flags)
    mean_obs = sum(obs) / len(obs)
    cnt_obs = sum(e > 0.5 for e in obs)
    random.seed(2017)
    fl = flags[:]
    hit_mean = hit_cnt = 0
    for _ in range(N_PERM):
        random.shuffle(fl)
        e = effects(fl)
        mn, ct = sum(e) / len(e), sum(x > 0.5 for x in e)
        hit_mean += abs(mn - 0.5) >= abs(mean_obs - 0.5) - 1e-12
        hit_cnt += max(ct, len(e) - ct) >= max(cnt_obs, len(obs) - cnt_obs)
    p_mean = (hit_mean + 1) / (N_PERM + 1)
    p_cnt = (hit_cnt + 1) / (N_PERM + 1)
    print(f"\nJoint test, the {len(JOINT)} Figure 2 measures together ({N_PERM:,} label shuffles):")
    for m, e in zip(JOINT, obs):
        print(f"   {m:<32} P(funded higher) = {e:.3f}")
    print(f"   mean P(funded higher) = {mean_obs:.3f}  (0.5 = no lean)   p = {p_mean:.4f}")
    print(f"   measures leaning to funded: {cnt_obs} of {len(obs)}            p = {p_cnt:.4f}")
    return {"group": "all", "measure": "joint: " + ", ".join(JOINT), "test": "permutation (label shuffle)",
            "n_funded": sum(flags), "n_unfunded": len(flags) - sum(flags), "funded": "", "unfunded": "",
            "effect": f"mean P(funded higher) = {mean_obs:.3f}; {cnt_obs} of {len(obs)} above 0.5",
            "range_95": f"p(count) = {p_cnt:.4g}", "p": p_mean, "p_adj": ""}


def main():
    genes = {r["hgnc_id"]: r for r in read_tsv(P / "r03_awards/gene_summary.tsv")
             if "2018" in r["listed_years"].split(";")}
    print(f"Genes on the 2018 list: {len(genes)}")
    S = defaultdict(dict)                       # hgnc -> column -> value

    # attention: Jensen score
    js = defaultdict(dict)
    for r in read_tsv(P / "attention/gene_year_scores.tsv"):
        if r["hgnc_id"] in genes and r["jensen_score"] != "":
            js[r["hgnc_id"]][int(r["year"])] = float(r["jensen_score"])
    for h, d in js.items():
        S[h]["jensen_total_2017"] = round(sum(v for y, v in d.items() if y <= END), 2)
        S[h]["jensen_mean_2013_2017"] = round(sum(d.get(y, 0) for y in WINDOW) / len(WINDOW), 2)

    # papers: human + one-to-one animal matches
    pp = defaultdict(lambda: defaultdict(set))
    p10 = defaultdict(lambda: defaultdict(set))
    for r in read_tsv(P / "attention/gene_papers.tsv"):
        h, y = r["hgnc_id"], yr(r["year"])
        if h in genes and y and y <= END and (r["species"] == "human" or r["match_type"] == "one-to-one"):
            pp[h][y].add(r["pmid"])
            if r["n_genes_in_paper"] and int(r["n_genes_in_paper"]) <= 10:
                p10[h][y].add(r["pmid"])
    for h in genes:
        if h in js or h in pp:
            for name, d in (("papers", pp[h]), ("papers_le10", p10[h])):
                allp = set().union(*d.values()) if d else set()
                S[h][f"{name}_total_2017"] = len(allp)
                S[h][f"{name}_mean_2013_2017"] = round(sum(len(d.get(y, ())) for y in WINDOW) / len(WINDOW), 2)

    # funding
    gr, r01 = defaultdict(set), defaultdict(set)
    for r in read_tsv(P / "funding/gene_grants.tsv"):
        h, y = r["hgnc_id"], yr(r["fiscal_year"])
        if h in genes and y and y <= END and r["grant_group"] != "IDG R03" and r.get("family_name_only", "0") != "1":
            gr[h].add(r["core_project_num"])
            if r["grant_group"] == "R01":
                r01[h].add(r["core_project_num"])
    for h in genes:
        S[h]["grants_total_2017"] = len(gr[h])
        S[h]["r01_total_2017"] = len(r01[h])

    # interactions (human)
    pa, ps = defaultdict(set), defaultdict(set)
    for r in read_tsv(P / "interactions/gene_interactions.tsv"):
        h, y = r["hgnc_id"], yr(r["year"])
        if h in genes and r["match_type"] == "same gene" and y and y <= END and r["self_interaction"] == "0":
            pa[h].add(r["partner_gene_id"])
            if "low" in r["throughput"]:
                ps[h].add(r["partner_gene_id"])
    for h in genes:
        S[h]["partners_all_2017"] = len(pa[h])
        S[h]["partners_small_2017"] = len(ps[h])

    # chemistry
    cp, pot = defaultdict(set), set()
    for r in read_tsv(P / "chemistry/gene_activities.tsv"):
        h, y = r["hgnc_id"], yr(r["year"])
        if h in genes and r["target_type"] == "SINGLE PROTEIN" and r["valid"] == "1" and y and y <= END:
            cp[h].add(r["compound"])
            if r["potent"] == "1":
                pot.add(h)
    for h in genes:
        S[h]["compounds_2017"] = len(cp[h])
        S[h]["potent_by_2017"] = int(h in pot)

    # function
    go = defaultdict(set)
    for r in read_tsv(P / "go/gene_go_annotations.tsv"):
        h, y = r["hgnc_id"], yr(r["year"])
        if (h in genes and r["evidence_group"] == "experimental" and r["informative"] == "1"
                and r["negated"] == "0" and y and y <= END):
            go[h].add(r["go_id"])
    for h in genes:
        S[h]["go_terms_2017"] = len(go[h])

    # structure
    st = defaultdict(set)
    for r in read_tsv(P / "structure/gene_pdb_entries.tsv"):
        h, y = r["hgnc_id"], yr(r["year"])
        if h in genes and y and y <= END:
            st[h].add(r["pdb_id"])
    for h in genes:
        S[h]["structures_2017"] = len(st[h])

    # disease
    dsub, dpap = set(), set()
    for r in read_tsv(P / "disease/gene_disease_links.tsv"):
        h = r["hgnc_id"]
        if h not in genes:
            continue
        if yr(r["first_moderate_plus_year"]) and yr(r["first_moderate_plus_year"]) <= END:
            dsub.add(h)
        if r["best_grade"] in ("Definitive", "Strong", "Moderate") and yr(r["earliest_paper_year"]) \
                and yr(r["earliest_paper_year"]) <= END:
            dpap.add(h)
    for h in genes:
        S[h]["disease_modplus_submitted_by_2017"] = int(h in dsub)
        S[h]["disease_modplus_paper_by_2017"] = int(h in dpap)

    # IDG rule check
    for h in genes:
        j = S[h].get("jensen_total_2017")
        S[h]["idg_jensen_under_50"] = "" if j is None else int(j < 50)
        S[h]["idg_no_r01"] = int(S[h]["r01_total_2017"] == 0)
        S[h]["idg_meets_both"] = "" if j is None else int(j < 50 and S[h]["r01_total_2017"] == 0)

    measures = ["jensen_total_2017", "jensen_mean_2013_2017", "papers_total_2017", "papers_mean_2013_2017",
                "papers_le10_total_2017", "papers_le10_mean_2013_2017",
                "grants_total_2017", "r01_total_2017", "partners_all_2017", "partners_small_2017", "compounds_2017",
                "potent_by_2017", "go_terms_2017", "structures_2017", "disease_modplus_submitted_by_2017",
                "disease_modplus_paper_by_2017", "idg_jensen_under_50", "idg_no_r01", "idg_meets_both"]
    rows = []
    for h, g in sorted(genes.items(), key=lambda kv: (kv[1]["family"], kv[1]["gene"])):
        rows.append({"gene": g["gene"], "hgnc_id": h, "family": g["family"], "listed_years": g["listed_years"],
                     "funded_years": g["funded_years"], "funded_any": int(bool(g["funded_years"])),
                     "funded_round1": int("2018" in g["funded_years"].split(";")),
                     **{m: S[h].get(m, "") for m in measures}})

    # summary
    binary = {"potent_by_2017", "disease_modplus_submitted_by_2017", "disease_modplus_paper_by_2017",
              "idg_jensen_under_50", "idg_no_r01", "idg_meets_both"}
    summ = []

    def describe(label, rs):
        out = {"group": label, "n_genes": len(rs)}
        for m in measures:
            vals = [float(r[m]) for r in rs if r[m] != ""]
            if not vals:
                out[m] = ""
            elif m in binary:
                out[m] = f"{100 * sum(vals) / len(vals):.0f}%"
            else:
                out[m] = f"median {statistics.median(vals):g} | {100 * sum(v > 0 for v in vals) / len(vals):.0f}% >0"
        summ.append(out)

    for fam in ["all"] + sorted({r["family"] for r in rows}):
        sub = rows if fam == "all" else [r for r in rows if r["family"] == fam]
        describe(f"{fam} | later funded", [r for r in sub if r["funded_any"]])
        describe(f"{fam} | never funded", [r for r in sub if not r["funded_any"]])

    print("\nAll 2018 genes, later funded vs never funded (median | % with at least one; or % yes):")
    show = ["jensen_total_2017", "papers_total_2017", "papers_le10_total_2017", "grants_total_2017", "r01_total_2017", "partners_small_2017",
            "compounds_2017", "potent_by_2017", "go_terms_2017", "structures_2017",
            "disease_modplus_paper_by_2017", "idg_meets_both"]
    for s in summ[:2]:
        print(f"\n   {s['group']} (n = {s['n_genes']})")
        for m in show:
            print(f"      {m:<32} {s[m]}")
    print("\nBy family (n later funded / never funded):")
    for s in summ[2:]:
        print(f"   {s['group']:<32} n = {s['n_genes']}")

    # tests: later funded vs never funded, overall and within each family
    tests = []
    for fam in ["all"] + sorted({r["family"] for r in rows}):
        sub = rows if fam == "all" else [r for r in rows if r["family"] == fam]
        f_, u_ = [r for r in sub if r["funded_any"]], [r for r in sub if not r["funded_any"]]
        for m in measures:
            x = [float(r[m]) for r in f_ if r[m] != ""]
            y = [float(r[m]) for r in u_ if r[m] != ""]
            if not x or not y:
                continue
            if m in binary:
                a, c = int(sum(x)), int(sum(y))
                p = fisher_exact(a, len(x) - a, c, len(y) - c)
                lo, hi = newcombe(a, len(x), c, len(y))
                tests.append({"group": fam, "measure": m, "test": "Fisher", "n_funded": len(x), "n_unfunded": len(y),
                              "funded": f"{100 * a / len(x):.1f}%", "unfunded": f"{100 * c / len(y):.1f}%",
                              "effect": f"{100 * (a / len(x) - c / len(y)):+.1f} points",
                              "range_95": f"{100 * lo:+.1f} to {100 * hi:+.1f}", "p": p})
            else:
                p, ps = mann_whitney(x, y)
                tests.append({"group": fam, "measure": m, "test": "Mann-Whitney", "n_funded": len(x),
                              "n_unfunded": len(y), "funded": f"median {statistics.median(x):g}",
                              "unfunded": f"median {statistics.median(y):g}",
                              "effect": f"P(funded higher) = {ps:.2f}", "range_95": "", "p": p})
    for t, pa in zip(tests, bh([t["p"] for t in tests])):
        t["p_adj"] = pa
    print(f"\nTests ({len(tests)}; p_adj = Benjamini-Hochberg across all of them). All genes:")
    for t in tests:
        if t["group"] == "all":
            print(f"   {t['measure']:<32} {t['funded']:>14} vs {t['unfunded']:<14} {t['effect']:<26} "
                  f"{t['range_95']:<17} p={t['p']:.3f}  p_adj={t['p_adj']:.3f}")
    print("Within families, tests with p_adj < 0.10:")
    hits = [t for t in tests if t["group"] != "all" and t["p_adj"] < 0.10]
    for t in hits:
        print(f"   {t['group']:<12} {t['measure']:<32} {t['funded']} vs {t['unfunded']}  {t['effect']}  p_adj={t['p_adj']:.3f}")
    if not hits:
        print("   none")
    for t in tests:
        t["p"], t["p_adj"] = f"{t['p']:.4g}", f"{t['p_adj']:.4g}"
    jt = joint_test(rows)
    jt["p"] = f"{jt['p']:.4g}"
    tests.append(jt)

    OUT.mkdir(parents=True, exist_ok=True)
    for name, data in (("gene_state_2017.tsv", rows), ("gene_state_2017_summary.tsv", summ),
                       ("gene_state_2017_tests.tsv", tests)):
        with (OUT / name).open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(data)
        print(f"Saved: {OUT / name} ({len(data)} rows)")


if __name__ == "__main__":
    main()
