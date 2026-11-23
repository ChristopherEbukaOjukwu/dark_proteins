#!/usr/bin/env python3
"""Figure 8: R03 outcomes and what happened to the gene afterwards.

Run from the dark_proteins folder, after script 17:
    python3 scripts/17b_fig8_outcomes.py

Panel a (3 x 3, one per measure)  did genes whose R03 reported a positive gene-specific finding grow faster?
    One dot per funded gene on the 2018 list (98), two rows:
        Positive gene-specific finding reported / No positive gene-specific finding reported.
    x = the gene's yearly pace after 2018 against the pace expected for a gene with the same pace before,
        family and first funding round (script 17's regression without the finding; R03-tied items taken out).
        1x = as expected; 2x = twice the expected pace; log scale.
    Dark line = row average on the log scale (geometric mean), the quantity the regression compares.
    Under each title: finding vs no finding, adjusted (% difference, 95% range, p_adj).
    Vertical stripes in count measures (R01s, partners, GO, structures) are genes with the same counts
    before and after (often 0), which get the same expected value.
Panel b  does the kind of finding show up later? Observed genes with any new item (dark) vs expected (grey)
    if the finding made no difference (same family, same starting study level, same calendar years; script 17).
    Genes with items tied directly to their R03 are named.

Inputs (data/processed/analysis, from script 17)
    outcome_vs_change_genes.tsv, outcome_vs_change_tests.tsv, outcome_matching.tsv, outcome_matching_tests.tsv
Outputs (results/figures/, PNG to view, PDF for the paper)
    fig8_outcomes
"""
import csv
import math
import random
import statistics
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullLocator
    from matplotlib.transforms import offset_copy
except ImportError:
    raise SystemExit("matplotlib is needed: pip install matplotlib")

A = Path("data/processed/analysis")
OUT = Path("results/figures")
FOUND, NONE = "#7a5195", "#a8a8a8"          # purple, light grey (checked: distinct in normal and colour-blind vision)
OBS, EXP = "#2b2b2a", "#bdbcb6"             # observed, expected
INK, INK2, GRID, SURFACE = "#1f1f1e", "#5f5e5a", "#e6e5e0", "#fcfcfb"
YES_LAB = "Positive gene-specific finding reported"
NO_LAB = "No positive gene-specific finding reported"
MEASURES = [("jensen", "Jensen score"), ("pubtator", "PubTator score"), ("r01", "New R01s"),
            ("papers", "Focused papers"), ("grants", "New NIH grants"), ("partners", "New small-scale partners"),
            ("go", "New GO terms (experimental)"), ("compounds", "New compounds tested"),
            ("structures", "New PDB structures")]
CHECKS = {"structure": "Structure\n→ new PDB entries", "new ligand": "New ligand\n→ new compounds",
          "binding partners": "Binding partners\n→ new partners",
          "role / mechanism / location": "Role, mechanism, location\n→ new GO terms"}
XT = [(-1.0, "¹⁄₁₀×"), (math.log10(0.25), "¼×"), (math.log10(0.5), "½×"),
      (0.0, "1×"), (math.log10(2), "2×"), (math.log10(4), "4×"), (1.0, "10×")]


def read_tsv(path):
    if not path.is_file():
        raise SystemExit(f"Missing {path} (run script 17 first)")
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5, length=0)
    ax.grid(color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)


def under_title(ax, text):
    t = ax.text(0, 1, text, transform=offset_copy(ax.transAxes, fig=ax.figure, y=4, units="points"),
                fontsize=7.8, color=INK2, va="bottom")
    t.set_in_layout(False)


def fmt_p(v):
    v = float(v)
    return "< 0.001" if v < 0.001 else f"{v:.3f}" if v < 0.01 else f"{v:.2f}"


def main():
    genes = read_tsv(A / "outcome_vs_change_genes.tsv")
    tests = {t["measure"]: t for t in read_tsv(A / "outcome_vs_change_tests.tsv")}
    mtests = read_tsv(A / "outcome_matching_tests.tsv")
    mrows = read_tsv(A / "outcome_matching.tsv")
    yes = [g for g in genes if g["found_something"] == "yes"]
    no = [g for g in genes if g["found_something"] != "yes"]

    fig = plt.figure(figsize=(15, 15.5), facecolor=SURFACE)
    gs = fig.add_gridspec(4, 3, height_ratios=[1, 1, 1, 1.35], hspace=0.62, wspace=0.1, top=0.87, bottom=0.05)

    # ------------------------------------------------------------ panel a
    allv = [float(g[f"{m}_vs_expected"]) for g in genes for m, _ in MEASURES if g[f"{m}_vs_expected"] != ""]
    lo, hi = min(allv), max(allv)
    xlim = (min(lo, -0.35) - 0.08, max(hi, 0.35) + 0.08)
    random.seed(8)
    print("Panel a: average pace vs expected (geometric mean, x), finding reported vs not")
    for k, (m, name) in enumerate(MEASURES):
        ax = fig.add_subplot(gs[k // 3, k % 3])
        style(ax)
        meds = []
        for y0, grp, col in ((1, yes, FOUND), (0, no, NONE)):
            v = [float(g[f"{m}_vs_expected"]) for g in grp if g[f"{m}_vs_expected"] != ""]
            ax.scatter(v, [y0 + random.uniform(-0.26, 0.26) for _ in v], s=14, color=col, alpha=0.8,
                       edgecolor=SURFACE, linewidth=0.3, zorder=2)
            md = statistics.fmean(v)
            ax.plot([md, md], [y0 - 0.4, y0 + 0.4], color=INK, linewidth=1.8, zorder=3, solid_capstyle="round")
            meds.append(10 ** md)
        ax.axvline(0, color=INK2, linewidth=0.8, linestyle=(0, (3, 3)), zorder=1)
        ax.set_xlim(*xlim)
        ax.set_xticks([p for p, _ in XT if xlim[0] <= p <= xlim[1]])
        ax.set_xticklabels([t for p, t in XT if xlim[0] <= p <= xlim[1]])
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_ylim(-0.6, 1.6)
        ax.set_yticks([1, 0])
        if k % 3 == 0:
            ax.set_yticklabels([f"Positive gene-specific\nfinding reported ({len(yes)})",
                                f"No positive gene-specific\nfinding reported ({len(no)})"], color=INK, fontsize=8.5)
        else:
            ax.tick_params(labelleft=False)
        ax.grid(axis="y", visible=False)
        if k // 3 == 2:
            ax.set_xlabel("yearly pace after 2018 vs expected (log scale)", color=INK2, fontsize=8.5)
        letter = "a. " if k == 0 else ""
        ax.set_title(f"{letter}{name}", loc="left", fontsize=10, color=INK, fontweight="bold", pad=18)
        t = tests.get(m)
        if t:
            under_title(ax, f"Finding vs no finding, adjusted: {float(t['yes_vs_expected_pct']):+.0f}% "
                            f"({float(t['lo95']):+.0f} to {float(t['hi95']):+.0f}), p_adj {fmt_p(t['p_adj'])}")
        print(f"   {name:30} average {meds[0]:.2f}x vs {meds[1]:.2f}x")

    # ------------------------------------------------------------ panel b
    ax = fig.add_subplot(gs[3, :])
    style(ax)
    order = [c for c in CHECKS if any(t["finding"] == c for t in mtests)]
    w = 0.36
    top = 0
    print("\nPanel b: observed vs expected (family + level)")
    for i, c in enumerate(order):
        t = next(t for t in mtests if t["finding"] == c)
        obs, exp, n = int(t["observed"]), float(t["expected_family_level"]), int(t["n_genes"])
        ax.bar(i - w / 2, obs, width=w, color=OBS, zorder=2)
        ax.bar(i + w / 2, exp, width=w, color=EXP, zorder=2)
        ax.text(i - w / 2, obs + 0.3, f"{obs} of {n}", ha="center", va="bottom", fontsize=8.5, color=INK)
        ax.text(i + w / 2, exp + 0.3, f"{exp:.1f}", ha="center", va="bottom", fontsize=8.5, color=INK2)
        ax.text(i, -0.2, f"p = {fmt_p(t['p_exact_family_level'])}", transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=8.5, color=INK2)
        tied = sorted({r["gene"] for r in mrows if r["finding"] == c and r["r03_tied"] not in ("", "0")})
        if tied:
            ax.text(i, -0.3, "tied to the R03: " + ", ".join(tied), transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=8, color=INK)
        top = max(top, obs, exp)
        print(f"   {c:28} observed {obs} of {n}, expected {exp:.1f}, p {t['p_exact_family_level']}"
              + (f"; tied: {', '.join(tied)}" if tied else ""))
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([CHECKS[c] for c in order], color=INK, fontsize=9)
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_ylim(0, top * 1.25 + 1)
    ax.set_ylabel("genes with any new item", color=INK2, fontsize=9)
    ax.grid(axis="x", visible=False)
    ax.set_title("b. Does the kind of finding show up later?", loc="left", fontsize=10.5, color=INK,
                 fontweight="bold", pad=40)
    t = ax.text(0, 1, "Dark: genes with any new item after the R03 started. Grey: expected if the finding made "
                      "no difference (other funded genes; same family, study level and years).",
                transform=offset_copy(ax.transAxes, fig=fig, y=5, units="points"), fontsize=8.5, color=INK2,
                va="bottom")
    t.set_in_layout(False)
    h = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (OBS, EXP)]
    ax.legend(h, ["Observed", "Expected"], loc="lower left", bbox_to_anchor=(-0.01, 1),
              bbox_transform=offset_copy(ax.transAxes, fig=fig, y=19, units="points"),
              frameon=False, fontsize=8.5, ncol=2, labelcolor=INK, borderaxespad=0, handletextpad=0.3)

    # ------------------------------------------------------------ titles
    joint = next((t for k_, t in tests.items() if k_.startswith("JOINT")), None)
    fig.suptitle("Figure 8. R03 outcomes and what happened to the gene afterwards",
                 x=0.01, y=0.985, ha="left", fontsize=12.5, color=INK, fontweight="bold")
    sub = ("a. One dot per funded gene on the 2018 list. x = its yearly pace after 2018 against the pace expected for a "
           "gene with the same pace before, family and funding round (R03-tied items taken out); 1× = as expected.\n"
           "Dark line: group average (geometric mean). Under each title: genes with a finding vs without, adjusted "
           "for the same factors (% difference, 95% range, Benjamini-Hochberg p_adj).")
    if joint:
        sub += (f" Across all 9 measures: finding higher on average, P = {float(joint['p_yes_higher']):.2f}, "
                f"label-shuffle p = {fmt_p(joint['p'])}.")
    fig.text(0.01, 0.962, sub, fontsize=9, color=INK2, va="top")
    h = [plt.Line2D([], [], marker="o", linestyle="", markersize=6, color=c) for c in (FOUND, NONE)]
    fig.legend(h, [YES_LAB, NO_LAB], loc="upper left", bbox_to_anchor=(0.01, 0.928), frameon=False, fontsize=9,
               ncol=2, labelcolor=INK)

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig8_outcomes.{ext}", dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print(f"\nSaved: {OUT / 'fig8_outcomes'}.png / .pdf")


if __name__ == "__main__":
    main()
