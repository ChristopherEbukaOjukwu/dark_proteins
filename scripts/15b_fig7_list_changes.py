#!/usr/bin/env python3
"""Figure 7: did genes removed from the IDG list grow more than genes that stayed?

Run from the dark_proteins folder, after script 15:
    python3 scripts/15b_fig7_list_changes.py

Rows a, b, c = Jensen, PubTator, R01; columns = the three list changes.
    Each panel goes from "through the year before the 3-year window" to "through the window's last year":
        2018 -> 2019   through 2015 -> through 2018
        2020 -> 2021   through 2017 -> through 2020
        2021 -> 2022   through 2018 -> through 2021
    Jensen, PubTator  one thin line per gene = its total score at the two points (grey stayed, orange removed);
                      thick line = group median, labelled with its growth (end median / start median - 1).
                      Log scale (0 shown), so a steeper line = faster growth in proportion, not just more papers.
                      Under each title: growth in proportion, removed vs stayed (Mann-Whitney), and the end
                      level allowing for start level + family (regression), both from script 15.
    R01               % of each group with any R01 (R01, R37, RF1) naming the gene, at the two points;
                      under the title: % with a new R01 in the window, removed vs stayed (Fisher, script 15).
    Genes with no Jensen or PubTator score are left out of that row and named under the figure.
Panel d  every list change: genes removed (down) and added (up), by family. 2019 -> 2020: no change.

Inputs
    data/processed/analysis/list_change_genes.tsv   (script 15)
    data/processed/analysis/list_change_tests.tsv   (script 15)
    data/processed/r03_awards/gene_summary.tsv      (family, listed_years: for genes added)
Outputs (results/figures/, PNG to view, PDF for the paper)
    fig7_list_changes
"""
import csv
import math
import statistics
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullLocator
    from matplotlib.transforms import offset_copy, blended_transform_factory
except ImportError:
    raise SystemExit("matplotlib is needed: pip install matplotlib")

GENES = Path("data/processed/analysis/list_change_genes.tsv")
TESTS = Path("data/processed/analysis/list_change_tests.tsv")
SUMMARY = Path("data/processed/r03_awards/gene_summary.tsv")
OUT = Path("results/figures")
REMOVED, STAYED, ADDED = "#c4651a", "#5b6b7f", "#2f6fa8"   # orange, neutral slate, blue (colour-blind checked)
INK, INK2, GRID, SURFACE = "#1f1f1e", "#5f5e5a", "#e6e5e0", "#fcfcfb"
FAMILIES = ["GPCR", "Ion channel", "Kinase"]
CHANGES = [(2018, 2019), (2020, 2021), (2021, 2022)]
WINDOW = {2019: (2016, 2018), 2021: (2018, 2020), 2022: (2019, 2021)}
MEASURES = [("jensen", "Jensen score"), ("pubtator", "PubTator score")]
FEW = 15               # draw removed genes' lines stronger when a change has this many or fewer


def read_tsv(path):
    if not path.is_file():
        raise SystemExit(f"Missing {path}")
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


def sl(v):
    """Position on the symlog axis (linthresh 1), used only to keep labels apart."""
    return v if v < 1 else 1 + math.log10(v)


def inv_sl(p):
    return p if p < 1 else 10 ** (p - 1)


def spread(items, gap):
    """items: [(position, ...)]; sorted and pushed apart so neighbours are at least `gap` apart."""
    out = []
    for p, *rest in sorted(items, key=lambda d: d[0]):
        if out and p - out[-1][0] < gap:
            p = out[-1][0] + gap
        out.append([p, *rest])
    return out


def spread_frac(items, gap, lo=0.02, hi=0.98):
    """items: [(axes fraction, ...)]; pushed apart by `gap` and kept inside [lo, hi]."""
    out = spread(items, gap)
    over = out[-1][0] - hi if out else 0
    if over > 0:                                  # slide the block down, then re-check from the bottom
        for d in out:
            d[0] -= over
        for i in range(len(out) - 2, -1, -1):
            if out[i + 1][0] - out[i][0] < gap:
                out[i][0] = out[i + 1][0] - gap
    return out


def score_panel(ax, rows, m, change, first_col):
    """Thin line per gene, thick median per group; log scale."""
    style(ax)
    old, new = change
    a, b = WINDOW[new]
    rows = [r for r in rows if r["change"] == f"{old}->{new}"]
    missing = sorted(r["gene"] for r in rows if r[f"{m}_before_window"] == "")
    rows = [r for r in rows if r[f"{m}_before_window"] != ""]
    pts = {}
    for status in ("stayed", "removed"):
        grp = [r for r in rows if r["status"] == status]
        s0 = [float(r[f"{m}_before_window"]) for r in grp]
        s1 = [x + float(r[f"{m}_gain_in_window"]) for x, r in zip(s0, grp)]
        pts[status] = (grp, s0, s1)
    few = len(pts["removed"][0]) <= FEW
    for status, color, alpha, lw in (("stayed", STAYED, 0.12, 0.6),
                                     ("removed", REMOVED, 0.85 if few else 0.45, 1.0 if few else 0.7)):
        grp, s0, s1 = pts[status]
        for x0, x1 in zip(s0, s1):
            ax.plot([0, 1], [x0, x1], color=color, alpha=alpha, linewidth=lw, zorder=2 if status == "stayed" else 3)
    labels, out = [], {}
    for status, color in (("stayed", STAYED), ("removed", REMOVED)):
        grp, s0, s1 = pts[status]
        m0, m1 = statistics.median(s0), statistics.median(s1)
        g = 100 * (m1 / m0 - 1) if m0 > 0 else float("nan")
        per_gene = [100 * (y / x - 1) for x, y in zip(s0, s1) if x > 0]
        ax.plot([0, 1], [m0, m1], color=color, linewidth=3.2, zorder=5, solid_capstyle="round")
        ax.scatter([0, 1], [m0, m1], s=36, color=color, edgecolor=SURFACE, linewidth=0.8, zorder=6)
        labels.append((sl(m1), f"{m1:.1f} ({g:+.0f}%)", color, True, sl(m1)))
        out[status] = (m0, m1, g, statistics.median(per_gene) if per_gene else float("nan"))
    # axis limits from the data only (labels never stretch the axis)
    allv = [v for p in pts.values() for v in p[1] + p[2]]
    top, low = max(allv), min(allv)
    ax.set_yscale("symlog", linthresh=1)
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_ylim(-0.1 if low < 1 else low * 0.7, top * 1.6)
    yt = [t for t in (0, 1, 10, 100, 1000, 10000) if ax.get_ylim()[0] <= t <= top * 1.6]
    ax.set_yticks(yt)
    ax.set_yticklabels([f"{t:,}" for t in yt])
    # labels placed in axes-fraction space, pushed apart, with a thin leader to their point
    frac = lambda v: ax.transAxes.inverted().transform(ax.transData.transform((0, v)))[1]
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    for p, t, c in spread_frac([(frac(out[s][0]), f"{out[s][0]:.1f}", col)
                                for s, col in (("stayed", STAYED), ("removed", REMOVED))], 0.06):
        ax.text(-0.05, p, t, transform=tr, ha="right", va="center", fontsize=8, color=c, fontweight="bold")
    for p, t, c, bold, p0 in spread_frac([(frac(inv_sl(q)), t, c, bold, frac(inv_sl(q)))
                                          for q, t, c, bold, _ in labels], 0.048):
        ax.annotate(t, (1, p0), xycoords=tr, xytext=(1.07, p), textcoords=tr, ha="left", va="center",
                    fontsize=8 if bold else 7.5, color=c, fontweight="bold" if bold else "normal",
                    arrowprops=dict(arrowstyle="-", color=INK2, linewidth=0.4, alpha=0.6, shrinkA=0, shrinkB=1))
    ax.set_xlim(-0.32, 1.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"through {a - 1}", f"through {b}"], color=INK, fontsize=8.5)
    ax.grid(axis="x", visible=False)
    if not first_col:
        ax.tick_params(labelleft=False)
    print(f"   {old}->{new}  {m:9} removed {out['removed'][0]:7.1f} -> {out['removed'][1]:7.1f} "
          f"({out['removed'][2]:+.0f}%; per-gene median {out['removed'][3]:+.0f}%)   "
          f"stayed {out['stayed'][0]:7.1f} -> {out['stayed'][1]:7.1f} "
          f"({out['stayed'][2]:+.0f}%; per-gene median {out['stayed'][3]:+.0f}%)")
    return missing


def r01_panel(ax, rows, change, first_col):
    """% of each group with any R01, at the two points."""
    style(ax)
    old, new = change
    a, b = WINDOW[new]
    rows = [r for r in rows if r["change"] == f"{old}->{new}"]
    vals = {}
    for status, color in (("stayed", STAYED), ("removed", REMOVED)):
        grp = [r for r in rows if r["status"] == status]
        p0 = 100 * sum(int(r["r01_before_window"]) > 0 for r in grp) / len(grp)
        p1 = 100 * sum(int(r["r01_before_window"]) + int(r["r01_in_window"]) > 0 for r in grp) / len(grp)
        ax.plot([0, 1], [p0, p1], color=color, linewidth=3.2, zorder=5, solid_capstyle="round")
        ax.scatter([0, 1], [p0, p1], s=36, color=color, edgecolor=SURFACE, linewidth=0.8, zorder=6)
        vals[status] = (p0, p1, color)
    for side, x, ha in ((0, -0.05, "right"), (1, 1.07, "left")):
        for p, t, c in spread([(v[side], f"{v[side]:.0f}%", v[2]) for v in vals.values()], 8):
            ax.text(x, p, t, ha=ha, va="center", fontsize=8, color=c, fontweight="bold")
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlim(-0.32, 1.8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"through {a - 1}", f"through {b}"], color=INK, fontsize=8.5)
    ax.grid(axis="x", visible=False)
    if not first_col:
        ax.tick_params(labelleft=False)
    print(f"   {old}->{new}  any R01   removed {vals['removed'][0]:.0f}% -> {vals['removed'][1]:.0f}%   "
          f"stayed {vals['stayed'][0]:.0f}% -> {vals['stayed'][1]:.0f}%")


def panel_counts(ax, summary, letter):
    """Genes removed (down) and added (up) by family, at each list change."""
    style(ax)
    listed = {r["hgnc_id"]: {int(y) for y in r["listed_years"].split(";") if y} for r in summary}
    fam = {r["hgnc_id"]: r["family"] for r in summary}
    width, gap = 0.8, 1.2
    xs, xl, groups = [], [], []
    print("\nRemoved / added by family")
    x, lim = 0, 0
    for old, new in CHANGES:
        start = x
        for f in FAMILIES:
            rem = sum(1 for h in listed if fam[h] == f and old in listed[h] and new not in listed[h])
            add = sum(1 for h in listed if fam[h] == f and new in listed[h] and old not in listed[h])
            ax.bar(x, add, width=width, color=ADDED, zorder=2)
            ax.bar(x, -rem, width=width, color=REMOVED, zorder=2)
            if add:
                ax.text(x, add + 1, f"+{add}", ha="center", va="bottom", fontsize=8, color=INK)
            if rem:
                ax.text(x, -rem - 1, f"−{rem}", ha="center", va="top", fontsize=8, color=INK)
            if not add and not rem:
                ax.text(x, 1, "0", ha="center", va="bottom", fontsize=8, color=INK2)
            xs.append(x)
            xl.append({"GPCR": "GPCR", "Ion channel": "Ion\nch.", "Kinase": "Kin."}[f])
            lim = max(lim, add, rem)
            print(f"   {old}->{new}  {f:12} removed {rem:3}  added {add:3}")
            x += 1
        groups.append(((start + x - 1) / 2, f"{old} → {new}"))
        x += gap
    ax.axhline(0, color=INK2, linewidth=0.9, zorder=3)
    ax.set_xticks(xs)
    ax.set_xticklabels(xl, color=INK, fontsize=8)
    for cx, lab in groups:
        ax.text(cx, -0.2, lab, transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=9, color=INK, fontweight="bold")
    ax.set_ylim(-lim * 1.25 - 3, lim * 1.25 + 3)
    ax.set_ylabel("genes", color=INK2, fontsize=9)
    ax.grid(axis="x", visible=False)
    ax.set_title(f"{letter}. Genes removed and added, by family", loc="left", fontsize=10.5,
                 color=INK, fontweight="bold", pad=40)
    t = ax.text(0, 1, "Down: removed. Up: added. 2019 → 2020: no change.",
                transform=offset_copy(ax.transAxes, fig=ax.figure, y=5, units="points"),
                fontsize=8.5, color=INK2, va="bottom")
    t.set_in_layout(False)
    h = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (REMOVED, ADDED)]
    ax.legend(h, ["Removed", "Added"], loc="lower left", bbox_to_anchor=(-0.01, 1),
              bbox_transform=offset_copy(ax.transAxes, fig=ax.figure, y=19, units="points"),
              frameon=False, fontsize=8.5, ncol=2, labelcolor=INK, borderaxespad=0, handletextpad=0.3)


def stats_line(ax, lines):
    """Test results under the panel title (fixed distance in points)."""
    t = ax.text(0, 1, "\n".join(lines), transform=offset_copy(ax.transAxes, fig=ax.figure, y=4, units="points"),
                fontsize=7.8, color=INK2, va="bottom", linespacing=1.35)
    t.set_in_layout(False)


def fmt_p(v):
    v = float(v)
    return "< 0.001" if v < 0.001 else f"{v:.3f}" if v < 0.01 else f"{v:.2f}"


def main():
    rows = read_tsv(GENES)
    summary = read_tsv(SUMMARY)
    tests = {(t["change"], t["measure"]): t for t in read_tsv(TESTS)}
    fig = plt.figure(figsize=(15, 18), facecolor=SURFACE)
    gs = fig.add_gridspec(4, 3, height_ratios=[1.2, 1.2, 0.75, 0.85], hspace=0.45, wspace=0.14,
                          top=0.87, bottom=0.06)
    missing = {}
    print("\nMedian totals at the two points (growth = end median / start median - 1)")
    for c, change in enumerate(CHANGES):
        old, new = change
        sub = [r for r in rows if r["change"] == f"{old}->{new}"]
        n_r = sum(r["status"] == "removed" for r in sub)
        n_s = sum(r["status"] == "stayed" for r in sub)
        for k, (m, name) in enumerate(MEASURES):
            ax = fig.add_subplot(gs[k, c])
            miss = score_panel(ax, rows, m, change, c == 0)
            if miss:
                missing.setdefault(name, set()).update(miss)
            ax.set_title(f"{'ab'[k]}. {name}" if c == 0 else name, loc="left", fontsize=10, color=INK, fontweight="bold", pad=34)
            short = name.split()[0]
            tg = tests.get((f"{old}->{new}", f"{short} growth in proportion"))
            ta = tests.get((f"{old}->{new}", f"{short} end level, allowing for start + family"))
            if tg and ta:
                pr = tg["effect"].split("= ")[1]
                stats_line(ax, [f"Growth: P(removed faster) {pr}, p_adj {fmt_p(tg['p_adj'])}",
                                f"Removed vs expected*: "
                                f"{ta['effect'].replace('removed ', '').replace(' vs expected', '')}, "
                                f"p_adj {fmt_p(ta['p_adj'])}"])
            if c == 0:
                ax.set_ylabel(f"{name}, running total (log)", color=INK, fontsize=9)
            if k == 0:
                t1 = ax.text(0, 1, f"{old} → {new}", va="bottom",
                             transform=offset_copy(ax.transAxes, fig=fig, y=70, units="points"),
                             fontsize=11.5, color=INK, fontweight="bold")
                t2 = ax.text(0, 1, f"{n_r} removed, {n_s} stayed  ·  window "
                                      f"{WINDOW[new][0]}-{WINDOW[new][1]}",
                             transform=offset_copy(ax.transAxes, fig=fig, y=56, units="points"),
                             va="bottom", fontsize=8.5, color=INK2)
        ax = fig.add_subplot(gs[2, c])
        r01_panel(ax, rows, change, c == 0)
        ax.set_title("c. Any R01 naming the gene" if c == 0 else "Any R01 naming the gene", loc="left", fontsize=10, color=INK,
                     fontweight="bold", pad=22)
        tr = tests.get((f"{old}->{new}", "got an R01 in window"))
        if tr:
            stats_line(ax, [f"New R01 in window: {tr['removed'].split(' (')[1].rstrip(')')} vs "
                            f"{tr['stayed'].split(' (')[1].rstrip(')')}, p_adj {fmt_p(tr['p_adj'])}"])
        if c == 0:
            ax.set_ylabel("% of genes", color=INK, fontsize=9)

    panel_counts(fig.add_subplot(gs[3, :]), summary, "d")

    fig.suptitle("Figure 7. Did genes removed from the IDG list grow more than genes that stayed?",
                 x=0.01, y=0.99, ha="left", fontsize=12.5, color=INK, fontweight="bold")
    fig.text(0.01, 0.97, "Thin lines: one gene's running total before the 3-year window and at its end "
                         "(grey stayed, orange removed). Thick lines: group medians, growth in brackets.\n"
                         "Log scale, so a steeper line = faster growth in proportion. Growth: Mann-Whitney on each gene's growth in proportion. "
                         "*Removed vs expected: end level compared with genes of the same start level and family (regression),\n"
                         "with 95% range. R01: Fisher on a new R01 in the window. All p_adj: Benjamini-Hochberg across script 15's tests.",
             fontsize=9, color=INK2, va="top")
    h = [plt.Line2D([], [], color=c, linewidth=3) for c in (REMOVED, STAYED)]
    fig.legend(h, ["Removed from the next list", "Stayed"], loc="upper right", bbox_to_anchor=(0.99, 0.995),
               frameon=False, fontsize=9, ncol=2, labelcolor=INK)
    if missing:
        note = "; ".join(f"no {n}: {', '.join(sorted(g))}" for n, g in missing.items())
        fig.text(0.01, 0.01, f"Left out of the score panels ({note}).", fontsize=7.5, color=INK2)
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig7_list_changes.{ext}", dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print(f"\nSaved: {OUT / 'fig7_list_changes'}.png / .pdf")


if __name__ == "__main__":
    main()
