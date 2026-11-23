#!/usr/bin/env python3
"""Figures 1 and 2: the anchor. Where each gene on the 2018 list stood at the end of 2017.

The 2018 FOA (RFA-RM-18-021): "IDG-eligible proteins have low Jensen and PubTator Scores and
minimal or no R01s." No cutoffs are given, so "low" is simply what each gene had at this point.

Run from the dark_proteins folder, after script 12:
    python3 scripts/12b_fig1_anchor.py

One dot per gene (328 genes on the 2018 list), two panels in each row:
    left   x = Jensen score, summed over all years to 2017
    right  x = PubTator score, summed over all years to 2017
    both   y = R01s naming the gene (title or abstract), any year to FY2017
Fig 1a  all genes, one colour. Dashed line = median score.
Fig 1b  the same genes: later funded (green; an R03 in any round, 2018-2022) vs never funded (magenta).
        Dashed lines = each group's median score.
Dots in the same R01 row are spread up and down a little so they don't hide each other.
Genes with no score are listed under the figure.

Inputs
    data/processed/analysis/gene_state_2017.tsv     (script 12: jensen_total_2017, r01_total_2017, funded_any)
    data/processed/attention/gene_year_scores.tsv   (script 04: pubtator_score per year)
Fig 2   the other measures, same genes, later funded (green) vs never funded (magenta).
        Counts: one row of dots per group (x = the count, log scale), short dark line = group median.
        Yes/no measures: % of each group.
            focused papers   human + one-to-one animal papers tagging 10 or fewer genes
            NIH grants       any NIH grant naming the gene (gene-specific matches; IDG R03s excluded)
            partners         small-scale (low-throughput) human interaction partners, BioGRID
            GO terms         informative experimental GO terms
            compounds        compounds tested on the protein, ChEMBL (single-protein, valid)
            potent compound  any compound reaching IDG's Tchem cutoff
            structures       experimental PDB structures
            disease link     GenCC Moderate-or-stronger link, dated by its earliest paper

Outputs (results/figures/, PNG to view, PDF for the paper)
    fig1_anchor_2017            1a and 1b together
    fig1a_anchor_2017           1a alone
    fig1b_anchor_2017_funded    1b alone
    fig2_state2017_funded       Figure 2
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
    from matplotlib.ticker import MaxNLocator, NullLocator
except ImportError:
    raise SystemExit("matplotlib is needed: pip install matplotlib")

STATE = Path("data/processed/analysis/gene_state_2017.tsv")
SCORES = Path("data/processed/attention/gene_year_scores.tsv")
OUT = Path("results/figures")
END = 2017
ALL, FUNDED, NEVER = "#5b6b7f", "#008300", "#d55181"     # neutral slate, green, magenta (colour-blind checked)
INK, INK2, GRID, SURFACE = "#1f1f1e", "#5f5e5a", "#e6e5e0", "#fcfcfb"
N_LABEL = 4          # name the top genes on each axis
R01_TICKS = [0, 1, 2, 5, 10, 20, 50]
SCORES_SHOWN = (("jensen", "Jensen"), ("pubtator", "PubTator"))
# Figure 2: (column in gene_state_2017.tsv, title, yes/no?)
FIG2 = [("papers_le10_total_2017", "Focused papers (10 or fewer genes)", False),
        ("grants_total_2017", "NIH grants naming the gene", False),
        ("partners_small_2017", "Interaction partners (small-scale)", False),
        ("go_terms_2017", "Known functions (GO, experimental)", False),
        ("compounds_2017", "Compounds tested (ChEMBL)", False),
        ("structures_2017", "Experimental structures (PDB)", False),
        ("potent_by_2017", "Potent compound (IDG cutoff)", True),
        ("disease_modplus_paper_by_2017", "Solid disease link (GenCC)", True)]


def read_tsv(path):
    if not path.is_file():
        raise SystemExit(f"Missing {path}")
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def ypos(n):
    """R01 count -> vertical position (log2(n + 1): 0, 1, 2 ... stay readable, 20 doesn't dominate)."""
    return math.log2(n + 1)


def load():
    genes = read_tsv(STATE)
    ids = {g["hgnc_id"] for g in genes}
    pub, has_pub = {}, set()
    for r in read_tsv(SCORES):
        if r["hgnc_id"] in ids and int(r["year"]) <= END and r["pubtator_score"] != "":
            has_pub.add(r["hgnc_id"])
            pub[r["hgnc_id"]] = pub.get(r["hgnc_id"], 0.0) + float(r["pubtator_score"])
    rows = []
    random.seed(7)
    for g in genes:
        h = g["hgnc_id"]
        n = int(g["r01_total_2017"] or 0)
        rows.append({"gene": g["gene"], "r01": n, "funded": g["funded_any"] == "1",
                     "jensen": float(g["jensen_total_2017"]) if g["jensen_total_2017"] != "" else None,
                     "pubtator": pub[h] if h in has_pub else None,
                     "y": ypos(n) + random.uniform(-0.17, 0.17)})
    return rows


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8.5, length=0)
    ax.grid(color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)


def name_genes(ax, pts, key, top):
    """Name the top genes on this score and on R01s; labels pushed apart so they don't overlap."""
    named, seen = [], set()
    for r in sorted(pts, key=lambda r: -r[key])[:N_LABEL] + sorted(pts, key=lambda r: -r["r01"])[:N_LABEL]:
        if r["gene"] not in seen:
            seen.add(r["gene"])
            named.append({"r": r, "lx": math.log10(r[key] + 1), "ty": r["y"] + 0.4})
    named.sort(key=lambda d: d["ty"])
    for i, d in enumerate(named):
        for e in named[:i]:
            if abs(d["lx"] - e["lx"]) < 0.45 and d["ty"] - e["ty"] < 0.32:
                d["ty"] = e["ty"] + 0.32
    right = math.log10(top + 1) * 0.66
    for d in named:
        r = d["r"]
        ax.annotate(r["gene"], (r[key], r["y"]), xytext=(r[key], d["ty"]), fontsize=7.5, color=INK,
                    ha="right" if d["lx"] > right else "left", va="center",
                    arrowprops=dict(arrowstyle="-", color=INK2, linewidth=0.5, shrinkA=0, shrinkB=2))


def draw_row(axes, rows, split, letter):
    """One row: Jensen panel and PubTator panel. split=False: all genes; True: funded vs never funded."""
    top_r01 = max(r["r01"] for r in rows)
    for ax, (key, name) in zip(axes, SCORES_SHOWN):
        style(ax)
        pts = [r for r in rows if r[key] is not None]
        top, low = max(r[key] for r in pts), min(r[key] for r in pts)
        groups = [(pts, ALL, "")] if not split else \
                 [([r for r in pts if not r["funded"]], NEVER, "never funded"),
                  ([r for r in pts if r["funded"]], FUNDED, "later funded")]
        for grp, color, _ in groups:          # never funded drawn first, funded on top
            ax.scatter([r[key] for r in grp], [r["y"] for r in grp], s=16, color=color,
                       alpha=0.6 if not split else 0.7, edgecolor=SURFACE, linewidth=0.4, zorder=2)

        # medians: one line per group, labelled at the top
        ytxt = ypos(top_r01) + 0.55
        for k, (grp, color, _) in enumerate(groups):
            m = statistics.median(r[key] for r in grp)
            ax.axvline(m, color=color if split else INK2, linewidth=1.1, linestyle=(0, (4, 3)), zorder=1)
            ax.text(m, ytxt - 0.32 * k, f" median {m:.1f}", color=color if split else INK2, fontsize=8, va="top")

        name_genes(ax, pts, key, top)

        ax.set_xscale("symlog", linthresh=1)
        ax.xaxis.set_minor_locator(NullLocator())
        xt = [t for t in (0, 1, 10, 100, 1000, 10000) if t <= top * 1.5]
        ax.set_xticks(xt)
        ax.set_xticklabels([f"{t:,}" for t in xt])
        ax.set_xlim(-0.08 if low < 1 else low * 0.6, top * 2.2)
        ax.set_xlabel(f"{name} score (all years to 2017)", color=INK, fontsize=9.5)

    yt = [t for t in R01_TICKS if t <= max(top_r01, 2)]
    for ax in axes:
        ax.set_yticks([ypos(t) for t in yt])
        ax.set_ylim(-0.45, ypos(top_r01) + 0.6)
    axes[0].set_yticklabels([str(t) for t in yt])
    axes[1].tick_params(labelleft=False)
    axes[0].set_ylabel("R01s naming the gene (any year to 2017)", color=INK, fontsize=9.5)

    n = len(rows)
    if not split:
        n0 = sum(r["r01"] == 0 for r in rows)
        head = f"{letter}. All {n} genes on the 2018 list"
        sub = f"One dot per gene. Dashed line: median score. {n0} of {n} ({100 * n0 / n:.0f}%) had no R01 naming them."
    else:
        f = [r for r in rows if r["funded"]]
        u = [r for r in rows if not r["funded"]]
        z = lambda g: f"{100 * sum(r['r01'] == 0 for r in g) / len(g):.0f}%"
        head = f"{letter}. Later funded ({len(f)}) vs never funded ({len(u)})"
        sub = (f"Later funded = an IDG R03 in any round (2018-2022). Dashed lines: each group's median. "
               f"No R01: {z(f)} of funded, {z(u)} of never funded.")
    axes[0].set_title(head, loc="left", x=0, fontsize=10.5, color=INK, fontweight="bold", pad=24)
    t = axes[0].text(0, 1.035, sub, transform=axes[0].transAxes, fontsize=8.5, color=INK2)
    t.set_in_layout(False)          # long line: don't let it push the panels apart
    if split:
        h = [plt.Line2D([], [], marker="o", linestyle="", markersize=6, color=c) for c in (FUNDED, NEVER)]
        axes[1].legend(h, ["Later funded", "Never funded"], loc="lower right", bbox_to_anchor=(1.0, 1.0),
                       frameon=False, fontsize=8.5, ncol=2, labelcolor=INK, borderaxespad=0.2)


def finish(fig, rows, path, title):
    fig.suptitle(title, x=0.01, ha="left", fontsize=12, color=INK, fontweight="bold")
    notes = "; ".join(f"no {name} score: {', '.join(sorted(r['gene'] for r in rows if r[key] is None))}"
                      for key, name in SCORES_SHOWN if any(r[key] is None for r in rows))
    if notes:
        fig.text(0.01, 0.005, notes, fontsize=7.5, color=INK2)
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{path}.{ext}", dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print(f"Saved: {OUT / path}.png / .pdf")


def figure2():
    genes = read_tsv(STATE)
    groups = [("Later funded", FUNDED, [g for g in genes if g["funded_any"] == "1"]),
              ("Never funded", NEVER, [g for g in genes if g["funded_any"] != "1"])]
    random.seed(11)
    fig, axes = plt.subplots(2, 4, figsize=(16, 7.6), facecolor=SURFACE)
    print("\nFigure 2 (to the end of 2017): later funded vs never funded")
    for ax, (col, title, yesno) in zip(axes.flat, FIG2):
        style(ax)
        vals = [[float(g[col]) for g in grp if g[col] != ""] for _, _, grp in groups]
        if yesno:
            pct = [100 * sum(v) / len(v) for v in vals]
            ax.bar([0, 1], pct, width=0.6, color=[c for _, c, _ in groups], zorder=2)
            for x, p in enumerate(pct):
                ax.text(x, p + 1.5, f"{p:.0f}%", ha="center", va="bottom", fontsize=8.5, color=INK)
            ax.set_xticks([0, 1])
            ax.set_xticklabels([lab for lab, _, _ in groups], color=INK, fontsize=8.5)
            ax.set_ylim(0, max(pct) * 1.35 + 5)
            ax.set_ylabel("% of genes", color=INK2, fontsize=8.5)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
            ax.grid(axis="x", visible=False)
            sub = f"{pct[0]:.0f}% vs {pct[1]:.0f}% yes"
            print(f"   {title:38} yes: {pct[0]:.0f}% vs {pct[1]:.0f}%")
        else:
            meds, zeros = [], []
            for k, ((lab, color, _), v) in enumerate(zip(groups, vals)):
                y0 = 1 - k                      # funded on top
                ax.scatter(v, [y0 + random.uniform(-0.28, 0.28) for _ in v], s=12, color=color, alpha=0.6,
                           edgecolor=SURFACE, linewidth=0.3, zorder=2)
                m = statistics.median(v)
                # median drawn behind the dots and a bit taller than the band, so dots on it stay visible
                ax.plot([m, m], [y0 - 0.42, y0 + 0.42], color=INK, linewidth=1.8, zorder=1, solid_capstyle="round")
                meds.append(m)
                zeros.append(100 * sum(x == 0 for x in v) / len(v))
            ax.set_xscale("symlog", linthresh=1)
            ax.xaxis.set_minor_locator(NullLocator())
            top = max(max(v) for v in vals)
            xt = [t for t in (0, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000) if t <= top * 1.2]
            if len(xt) > 7:                     # wide range: decades only, so labels don't collide
                xt = [t for t in xt if t in (0, 1, 10, 100, 1000, 10000)]
            ax.set_xticks(xt)
            ax.set_xticklabels([f"{t:,}" for t in xt])
            ax.set_xlim(-0.15, top * 1.8)
            ax.set_ylim(-0.6, 1.6)
            ax.set_yticks([1, 0])
            ax.set_yticklabels(["funded", "never"], color=INK, fontsize=8.5)
            ax.grid(axis="y", visible=False)
            sub = f"median {meds[0]:g} vs {meds[1]:g}  ·  none: {zeros[0]:.0f}% vs {zeros[1]:.0f}%"
            print(f"   {title:38} median {meds[0]:g} vs {meds[1]:g} | none: {zeros[0]:.0f}% vs {zeros[1]:.0f}%")
        ax.set_title(title, loc="left", fontsize=9.5, color=INK, fontweight="bold", pad=18)
        ax.text(0, 1.03, sub, transform=ax.transAxes, fontsize=8, color=INK2)
    h = [plt.Line2D([], [], marker="o", linestyle="", markersize=6, color=c) for _, c, _ in groups]
    fig.legend(h, [f"{lab} ({len(grp)} genes)" for lab, _, grp in groups], loc="upper right",
               bbox_to_anchor=(0.995, 1.0), frameon=False, fontsize=9, ncol=2, labelcolor=INK)
    fig.suptitle("Figure 2. The 2018-list genes at the end of 2017: everything else we measure",
                 x=0.01, y=1.0, ha="left", fontsize=12, color=INK, fontweight="bold")
    fig.text(0.01, 0.94, "One dot per gene (log scale; 0 shown). Dark line: group median. "
                         "\"none\" = % of the group with zero. Later funded = an IDG R03 in any round (2018-2022).",
             fontsize=9, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.92), h_pad=3, w_pad=2.5)
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig2_state2017_funded.{ext}", dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print(f"Saved: {OUT / 'fig2_state2017_funded'}.png / .pdf")


def main():
    rows = load()
    title = "Figure 1. The genes on the 2018 list at the end of 2017: papers and R01s"

    fig, axes = plt.subplots(2, 2, figsize=(13, 11), facecolor=SURFACE)
    draw_row(axes[0], rows, False, "a")
    draw_row(axes[1], rows, True, "b")
    fig.tight_layout(rect=(0, 0.02, 1, 0.97), w_pad=3, h_pad=4)
    finish(fig, rows, "fig1_anchor_2017", title)

    for split, letter, path in ((False, "a", "fig1a_anchor_2017"), (True, "b", "fig1b_anchor_2017_funded")):
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), facecolor=SURFACE)
        draw_row(axes, rows, split, letter)
        fig.tight_layout(rect=(0, 0.035, 1, 0.94), w_pad=3)
        finish(fig, rows, path, title)

    print(f"\nGenes on the 2018 list: {len(rows)}")
    groups = (("all", rows), ("later funded", [r for r in rows if r["funded"]]),
              ("never funded", [r for r in rows if not r["funded"]]))
    print(f"   {'group':14} {'n':>4}  {'Jensen median':>13}  {'PubTator median':>15}  {'no R01':>7}")
    for lab, g in groups:
        jm = statistics.median(r["jensen"] for r in g if r["jensen"] is not None)
        pm = statistics.median(r["pubtator"] for r in g if r["pubtator"] is not None)
        z = 100 * sum(r["r01"] == 0 for r in g) / len(g)
        print(f"   {lab:14} {len(g):>4}  {jm:>13.1f}  {pm:>15.1f}  {z:>6.0f}%")

    figure2()


if __name__ == "__main__":
    main()
