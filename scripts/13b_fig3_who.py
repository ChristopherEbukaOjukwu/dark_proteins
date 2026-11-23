#!/usr/bin/env python3
"""Figure 3: who answered the call. The R03 awards that funded a 2018-list gene (any round).

Run from the dark_proteins folder, after script 13:
    python3 scripts/13b_fig3_who.py

One bar per round (and all rounds), each bar = 100% of that round's awards:
    a  PI status before the R03 (most senior PI on the award, NIH grants as PI)
    b  the PI and the gene before the R03: PI already held an NIH grant naming the gene /
       someone else did / no NIH grant named the gene
    c  years from the PI's first NIH grant to the R03, one dot per award; open dots = the R03 was
       the PI's first NIH grant
    d  the institution's rank by total NIH funding in the R03's first fiscal year (top 25 / 26-100 / 101+),
       with a last bar for all new NIH R03s in the same years (reference, weighted to the IDG years)
Inputs  data/processed/analysis/r03_awards_2018genes_who.tsv         (script 13)
        data/processed/analysis/r03_institution_bands_reference.tsv  (script 13, with the ExPORTER files)
Output  results/figures/fig3_who_answered.png / .pdf
"""
import csv
import random
import statistics
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    raise SystemExit("matplotlib is needed: pip install matplotlib")

IN = Path("data/processed/analysis/r03_awards_2018genes_who.tsv")
REF = Path("data/processed/analysis/r03_institution_bands_reference.tsv")
OUT = Path("results/figures")
INK, INK2, GRID, SURFACE = "#1f1f1e", "#5f5e5a", "#e6e5e0", "#fcfcfb"
# ordered categories -> one hue, dark = most (seniority / closeness to the gene)
STATUS = [("R01-equivalent", "#1f3a5f"), ("other research", "#4f7196"), ("training only", "#93abc4"),
          ("none as PI", "#d5dee8")]
GENE = [("PI already funded on the gene", "#0f5257"), ("someone else funded on it", "#5e9e98"),
        ("no NIH grant on the gene", "#d3e6e3")]
BANDS = [("top 25", "#4b2e63"), ("26-100", "#8a6aa6"), ("101+", "#d8cbe4"), ("not ranked", "#ecebe8")]
DOT = "#4f7196"


def read_tsv(path):
    if not path.is_file():
        raise SystemExit(f"Missing {path}; run scripts/13_who_answered.py first")
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def gene_cat(r):
    if int(r["pi_prior_gene_grants"]) > 0:
        return GENE[0][0]
    if int(r["gene_prior_grants_any_pi"]) > 0:
        return GENE[1][0]
    return GENE[2][0]


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5, length=0)
    ax.set_axisbelow(True)


def stacked_pct(ax, labels, pcts, cats, title, ypos=None):
    """Horizontal 100% bars. pcts[i] = list of % per category (same order as cats)."""
    style(ax)
    ypos = ypos or list(range(len(labels)))
    for y, row in zip(ypos, pcts):
        left = 0
        for (name, color), w in zip(cats, row):
            if w <= 0:
                continue
            ax.barh(y, w, left=left, height=0.62, color=color, edgecolor=SURFACE, linewidth=2, zorder=2)
            if w >= 9:
                dark = color in (cats[0][1], cats[1][1])
                ax.text(left + w / 2, y, f"{w:.0f}%", ha="center", va="center", fontsize=8,
                        color="white" if dark else INK)
            left += w
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, color=INK, fontsize=8.5)
    ax.set_ylim(max(ypos) + 0.6, min(ypos) - 0.6)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(["0", "25", "50", "75", "100%"])
    ax.grid(axis="x", color=GRID, linewidth=0.7)
    ax.set_title(title, loc="left", fontsize=10, color=INK, fontweight="bold", pad=34)
    h = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c in cats]
    ax.legend(h, [n for n, _ in cats], loc="lower left", bbox_to_anchor=(0, 1.0), ncol=len(cats) if len(cats) > 3 else 2,
              frameon=False, fontsize=8, labelcolor=INK, borderaxespad=0.2, handlelength=1.2, columnspacing=1.2)


def stacked(ax, groups, cats, catfun, title):
    pcts = [[100 * sum(catfun(r) == name for r in rows) / len(rows) for name, _ in cats] for _, rows in groups]
    stacked_pct(ax, [lab for lab, _ in groups], pcts, cats, title)


def main():
    rows = read_tsv(IN)
    rounds = sorted({r["round"] for r in rows})
    groups = [(f"All ({len(rows)})", rows)] + \
             [(f"{y} ({sum(r['round'] == y for r in rows)})", [r for r in rows if r["round"] == y]) for y in rounds]

    fig, ax4 = plt.subplots(2, 2, figsize=(14, 10), facecolor=SURFACE)
    axes = ax4.flat
    stacked(ax4[0, 0], groups, STATUS, lambda r: r["pi_status"], "a. PI's NIH record before the R03")
    stacked(ax4[0, 1], groups, GENE, gene_cat, "b. The PI and the gene before the R03")

    # c: years from first NIH grant, one dot per award
    ax = ax4[1, 0]
    style(ax)
    random.seed(3)
    for y, (lab, grp) in enumerate(groups):
        yrs = [int(r["years_since_first_nih"]) for r in grp if r["years_since_first_nih"] != ""]
        first = [r for r in grp if r["years_since_first_nih"] == ""]
        ax.scatter(yrs, [y + random.uniform(-0.22, 0.22) for _ in yrs], s=16, color=DOT, alpha=0.75,
                   edgecolor=SURFACE, linewidth=0.4, zorder=2)
        ax.scatter([0] * len(first), [y + random.uniform(-0.22, 0.22) for _ in first], s=16,
                   facecolor="none", edgecolor=DOT, linewidth=1, zorder=2)
        if yrs:
            m = statistics.median(yrs)
            ax.plot([m, m], [y - 0.32, y + 0.32], color=INK, linewidth=1.8, zorder=1)
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels([lab for lab, _ in groups], color=INK, fontsize=8.5)
    ax.invert_yaxis()
    ax.grid(axis="x", color=GRID, linewidth=0.7)
    ax.set_xlabel("years from the PI's first NIH grant to the R03", color=INK, fontsize=9)
    ax.set_title("c. How long the PI had held NIH grants", loc="left", fontsize=10, color=INK,
                 fontweight="bold", pad=34)
    h = [plt.Line2D([], [], marker="o", linestyle="", markersize=6, color=DOT),
         plt.Line2D([], [], marker="o", linestyle="", markersize=6, markerfacecolor="none", markeredgecolor=DOT),
         plt.Line2D([], [], color=INK, linewidth=1.8)]
    ax.legend(h, ["one award", "R03 was the first NIH grant", "median"], loc="lower left",
              bbox_to_anchor=(0, 1.0), ncol=3, frameon=False, fontsize=8, labelcolor=INK, borderaxespad=0.2,
              handlelength=1.2, columnspacing=1.0)

    # d: institution's NIH funding rank, with all new NIH R03s as reference
    ax = ax4[1, 1]
    have_rank = any(r.get("inst_band") for r in rows)
    if have_rank and REF.is_file():
        ref = {r["band"]: float(r["all_new_r03_pct"]) for r in read_tsv(REF)}
        pcts = [[100 * sum(r["inst_band"] == b for r in grp) / len(grp) for b, _ in BANDS] for _, grp in groups]
        pcts.append([ref.get(b, 0.0) for b, _ in BANDS])
        labels = [lab for lab, _ in groups] + ["All new NIH R03s\n(same years)"]
        ypos = list(range(len(groups))) + [len(groups) + 0.5]
        stacked_pct(ax, labels, pcts, BANDS, "d. The institution's NIH funding rank", ypos)
        ax.axhline(len(groups) - 0.25, color=GRID, linewidth=1)
    else:
        style(ax)
        ax.text(0.5, 0.5, "No institution ranks: add the ExPORTER files and rerun script 13",
                ha="center", va="center", color=INK2, fontsize=9, transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(f"Figure 3. Who answered: the {len(rows)} IDG R03 awards that funded a gene on the 2018 list",
                 x=0.01, y=1.0, ha="left", fontsize=12, color=INK, fontweight="bold")
    fig.text(0.01, 0.962, "Each bar = 100% of that round's awards (n in brackets). Most senior PI on the award; "
                          "NIH grants as PI only. Round = the FOA year of the award.", fontsize=9, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.95), w_pad=3, h_pad=3)
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig3_who_answered.{ext}", dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    print(f"Awards: {len(rows)}")
    print(f"   {'group':12} " + "  ".join(f"{n[:14]:>14}" for n, _ in STATUS) + "  | PI on gene  others  none")
    for lab, grp in groups:
        st = "  ".join(f"{100 * sum(r['pi_status'] == n for r in grp) / len(grp):>13.0f}%" for n, _ in STATUS)
        gc = [100 * sum(gene_cat(r) == n for r in grp) / len(grp) for n, _ in GENE]
        print(f"   {lab:12} {st}  | {gc[0]:>9.0f}%  {gc[1]:>5.0f}%  {gc[2]:>4.0f}%")
    print(f"Saved: {OUT / 'fig3_who_answered'}.png / .pdf")


if __name__ == "__main__":
    main()
