#!/usr/bin/env python3
"""Figures 4, 5 and 6: what changed after the 2018 announcement, gene by gene, and when.

Run from the dark_proteins folder, after script 14:
    python3 scripts/14b_fig45_change.py

Every panel: one dot per gene on the 2018 list.
    x = yearly average before (2013-2017), y = yearly average after (2019 to the measure's last reliable year)
    diagonal = same pace as before; above it = faster after, below = slower
    green = later funded (an IDG R03 in any round), magenta = never funded
Both axes use log(value + c), so 0 can be shown; dots are spread very slightly so identical values stay visible.
Fig 4  IDG's criteria: Jensen score, PubTator score (after = 2019-2021 only), new R01s
Fig 5  our measures: focused papers, NIH grants, small-scale partners, GO terms, compounds, structures,
       plus milestones: % of ALL genes in each group whose first-ever one came from 2018 on
       (a gene that already had it by 2017 counts as no: its first came earlier)
Fig 6  funded genes only, each against itself: year 0 = the year funding started (first R03 fiscal year);
       y = the gene's value that year minus its own average over the 5 years before; line = average over
       funded genes, band = 95% range. Field-wide growth is included (Figures 4-5 handle the comparison).
Input   data/processed/analysis/gene_change.tsv, gene_change_tests.tsv, gene_event_time.tsv   (script 14)
Output  results/figures/fig4_change_idg_criteria, fig5_change_our_measures, fig6_timing  (.png/.pdf)
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
except ImportError:
    raise SystemExit("matplotlib is needed: pip install matplotlib")

IN = Path("data/processed/analysis/gene_change.tsv")
TESTS = Path("data/processed/analysis/gene_change_tests.tsv")
EVENT = Path("data/processed/analysis/gene_event_time.tsv")
OUT = Path("results/figures")
FUNDED, NEVER = "#008300", "#d55181"
INK, INK2, GRID, SURFACE = "#1f1f1e", "#5f5e5a", "#e6e5e0", "#fcfcfb"
# measure, title, offset c for log(value + c), after-window label
FIG4 = [("jensen", "Jensen score per year", 0.1, "2019-2025"),
        ("pubtator", "PubTator score per year", 0.1, "2019-2021"),
        ("r01", "New R01s per year", 0.05, "FY2019-2026")]
FIG5 = [("papers", "New focused papers per year", 0.2, "2019-2025"),
        ("grants", "New NIH grants per year", 0.05, "FY2019-2026"),
        ("partners", "New small-scale partners per year", 0.05, "2019-2023"),
        ("go", "New GO terms per year", 0.05, "2019-2024"),
        ("compounds", "New compounds tested per year", 0.2, "2019-2024"),
        ("structures", "New PDB structures per year", 0.05, "2019-2026")]
MILESTONES = [("first_r01", "first\nR01"), ("first_structure", "first\nstructure"),
              ("first_potent", "first potent\ncompound"), ("first_partner", "first small-\nscale partner"),
              ("first_go", "first GO\nfunction"), ("first_disease", "first\ndisease link")]
TICKS = [0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]


def read_tsv(path):
    if not path.is_file():
        raise SystemExit(f"Missing {path}; run scripts/14_what_changed.py first")
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8, length=0)
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def padj_lookup():
    out = {}
    if TESTS.is_file():
        for t in read_tsv(TESTS):
            if t["comparison"] == "funded vs never":
                out[t["measure"]] = float(t["p_adj"])
    return out


def scatter_panel(ax, rows, m, title, c, window, padj):
    style(ax)
    f = lambda v: math.log10(v + c)
    pts = [r for r in rows if r[f"{m}_before"] != ""]
    top = max(max(float(r[f"{m}_before"]), float(r[f"{m}_after"])) for r in pts)
    lo, hi = f(0) - 0.08, f(top) + 0.12
    ax.plot([f(0), f(top)], [f(0), f(top)], color=INK2, linewidth=1, linestyle=(0, (4, 3)), zorder=1)
    for flag, color in (("0", NEVER), ("1", FUNDED)):          # never funded first, funded on top
        g = [r for r in pts if r["funded_any"] == flag]
        xs = [f(float(r[f"{m}_before"])) + random.uniform(-0.025, 0.025) for r in g]
        ys = [f(float(r[f"{m}_after"])) + random.uniform(-0.025, 0.025) for r in g]
        ax.scatter(xs, ys, s=13, color=color, alpha=0.6, edgecolor=SURFACE, linewidth=0.3, zorder=2)
    ticks = [t for t in TICKS if t <= top * 1.3]
    if len(ticks) > 7:
        ticks = [t for t in ticks if t in (0, 0.1, 1, 10, 100, 1000)]
    kept = []                                  # drop ticks too close to the previous one (e.g. 0 and 0.1)
    for t in ticks:
        if not kept or f(t) - f(kept[-1]) >= 0.15:
            kept.append(t)
    ticks = kept
    for setter, lab in ((ax.set_xticks, ax.set_xticklabels), (ax.set_yticks, ax.set_yticklabels)):
        setter([f(t) for t in ticks])
        lab([f"{t:g}" for t in ticks])
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel("before (2013-2017)", color=INK2, fontsize=8)
    ax.set_ylabel(f"after ({window})", color=INK2, fontsize=8)
    up = lambda flag: [r for r in pts if r["funded_any"] == flag]
    pr = lambda g: 100 * sum(float(r[f"{m}_change"]) > 0 for r in g) / len(g)
    sub = f"faster after: {pr(up('1')):.0f}% funded vs {pr(up('0')):.0f}% never"
    p = padj.get(m)
    if p is not None:
        sub += "  ·  " + ("clear difference" if p < 0.05 else "no clear difference")
    ax.set_title(title, loc="left", fontsize=9.5, color=INK, fontweight="bold", pad=16)
    ax.text(0, 1.015, sub, transform=ax.transAxes, fontsize=7.6, color=INK2)


def legend(fig, y):
    h = [plt.Line2D([], [], marker="o", linestyle="", markersize=6, color=c) for c in (FUNDED, NEVER)] + \
        [plt.Line2D([], [], color=INK2, linestyle=(0, (4, 3)))]
    fig.legend(h, ["Later funded", "Never funded", "same pace as before"], loc="upper right",
               bbox_to_anchor=(0.995, y), frameon=False, fontsize=8.5, ncol=3, labelcolor=INK)


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=200, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    print(f"Saved: {OUT / name}.png / .pdf")


def figure6():
    if not EVENT.is_file():
        print(f"No {EVENT}: rerun script 14 for Figure 6")
        return
    ev = read_tsv(EVENT)
    if "mean_change" not in ev[0]:
        print(f"{EVENT} is from an older script 14: rerun script 14 for Figure 6")
        return
    panels = [("jensen", "Jensen score per year"), ("pubtator", "PubTator score per year (to 2021)"),
              ("r01", "New R01s per year"), ("papers", "New focused papers per year"),
              ("grants", "New NIH grants per year"), ("partners", "New small-scale partners per year"),
              ("go", "New GO terms per year"), ("compounds", "New compounds tested per year"),
              ("structures", "New PDB structures per year")]
    fig, axes = plt.subplots(3, 3, figsize=(14, 11.5), facecolor=SURFACE)
    for ax, (m, title) in zip(axes.flat, panels):
        style(ax)
        e = [r for r in ev if r["measure"] == m]
        if not e:
            ax.set_visible(False)
            continue
        k = [int(r["k"]) for r in e]
        mean = [float(r["mean_change"]) for r in e]
        ax.axhline(0, color=INK2, linewidth=1.2, zorder=1)
        nv = [r for r in e if r.get("never_mean_change", "") != ""]
        if nv:
            kn = [int(r["k"]) for r in nv]
            ax.fill_between(kn, [float(r["never_lo95"]) for r in nv], [float(r["never_hi95"]) for r in nv],
                            color=NEVER, alpha=0.12, linewidth=0, zorder=2)
            ax.plot(kn, [float(r["never_mean_change"]) for r in nv], color=NEVER, linewidth=2.2, marker="o",
                    markersize=3.5, zorder=3, label="Never funded, same years")
        ax.fill_between(k, [float(r["lo95"]) for r in e], [float(r["hi95"]) for r in e], color=FUNDED,
                        alpha=0.15, linewidth=0, zorder=2)
        ax.plot(k, mean, color=FUNDED, linewidth=2.4, marker="o", markersize=4, zorder=4, label="Funded")
        ax.set_xlim(0, 5.15)                                # year 0 at the left edge, no gap
        ax.set_xticks(range(0, 6))
        ax.set_xticklabels(["0 (start)"] + [f"+{x}" for x in range(1, 6)], fontsize=7.5)
        ax.set_title(title, loc="left", fontsize=9.5, color=INK, fontweight="bold", pad=16)
        ax.text(0, 1.02, f"before funding: {float(e[0]['mean_baseline']):.2g} a year on average (= 0)",
                transform=ax.transAxes, fontsize=7.6, color=INK2)
    for ax in axes[2]:
        ax.set_xlabel("years from the start of funding", color=INK, fontsize=8.5)
    for ax in axes[:, 0]:
        ax.set_ylabel("change from the gene's own\npre-funding average", color=INK2, fontsize=8)
    h, lab = axes.flat[0].get_legend_handles_labels()
    fig.legend(h[::-1], lab[::-1], loc="upper left", bbox_to_anchor=(0.005, 0.94), frameon=False, fontsize=9,
               ncol=2, labelcolor=INK)
    fig.suptitle("Figure 6. Each gene against its own past: change from its pre-funding average, year by year",
                 x=0.01, y=1.0, ha="left", fontsize=12, color=INK, fontweight="bold")
    fig.text(0.01, 0.972, "Every gene against its own past: 0 = its own average over the 5 years before funding; year 0 = "
                          "funding starts. Line = average over genes; band = 95% range (genes resampled), same method "
                          "for both lines.\nNever-funded genes borrow the start years of funded genes in their family, "
                          "in the same proportions. Later years rest on fewer funded genes (rounds 2018-2019 only).",
             fontsize=8.3, color=INK2, va="top")
    fig.tight_layout(rect=(0, 0, 1, 0.91), h_pad=2.8, w_pad=2.5)
    save(fig, "fig6_timing")


def main():
    rows = read_tsv(IN)
    padj = padj_lookup()
    random.seed(5)
    nf = sum(r["funded_any"] == "1" for r in rows)
    note = (f"One dot per gene on the 2018 list ({nf} later funded, {len(rows) - nf} never funded). "
            "Yearly averages; axes log(value + small offset) so 0 is shown.")

    # ---- Figure 4
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.4), facecolor=SURFACE)
    for ax, (m, title, c, w) in zip(axes, FIG4):
        scatter_panel(ax, rows, m, title, c, w, padj)
    fig.suptitle("Figure 4. IDG's own criteria before and after the 2018 announcement",
                 x=0.01, y=1.02, ha="left", fontsize=12, color=INK, fontweight="bold")
    fig.text(0.01, 0.955, note, fontsize=8.5, color=INK2)
    legend(fig, 1.03)
    fig.tight_layout(rect=(0, 0, 1, 0.93), w_pad=3)
    save(fig, "fig4_change_idg_criteria")

    # ---- Figure 5
    fig = plt.figure(figsize=(16, 9.6), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 4, hspace=0.45, wspace=0.35)
    spots = [gs[0, 0], gs[0, 1], gs[0, 2], gs[0, 3], gs[1, 0], gs[1, 1]]
    for spot, (m, title, c, w) in zip(spots, FIG5):
        scatter_panel(fig.add_subplot(spot), rows, m, title, c, w, padj)

    ax = fig.add_subplot(gs[1, 2:])
    style(ax)
    ax.grid(axis="x", visible=False)
    x = list(range(len(MILESTONES)))
    wdt = 0.36
    print("\nMilestones (% of all genes in the group whose first came from 2018 on):")
    nlab = {}
    for k, (flag, color, lab) in enumerate((("1", FUNDED, "funded"), ("0", NEVER, "never"))):
        vals, ns = [], []
        for m, _ in MILESTONES:
            g = [r for r in rows if r["funded_any"] == flag]          # all genes; "" = had it by 2017 = no
            vals.append(100 * sum(r[m] == "1" for r in g) / len(g) if g else 0)
            ns.append(len(g))
        xs = [i + (k - 0.5) * wdt for i in x]
        ax.bar(xs, vals, width=wdt, color=color, edgecolor=SURFACE, linewidth=1.5, zorder=2)
        for i, (xx, v, n) in enumerate(zip(xs, vals, ns)):
            ax.text(xx, v + 1, f"{v:.0f}%", ha="center", va="bottom", fontsize=7.5, color=INK)
            nlab.setdefault(i, []).append(str(n))
        print(f"   {lab:7} " + "  ".join(f"{n}: {v:.0f}% (n={c})" for (n, _), v, c in zip(MILESTONES, vals, ns)))
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _, lab in MILESTONES], color=INK, fontsize=7.8)
    ax.set_ylim(0, 40)
    ax.set_yticks([0, 10, 20, 30, 40])
    ax.set_yticklabels(["0", "10", "20", "30", "40%"])
    ax.set_title("Milestones: first ever, from 2018 on", loc="left", fontsize=9.5, color=INK,
                 fontweight="bold", pad=16)
    ax.text(0, 1.015, "% of all genes whose first-ever one came from 2018 on (genes that already had one count as no)", transform=ax.transAxes,
            fontsize=7.6, color=INK2)

    fig.suptitle("Figure 5. Our measures before and after the 2018 announcement",
                 x=0.1, y=0.995, ha="left", fontsize=12, color=INK, fontweight="bold")
    fig.text(0.1, 0.962, note, fontsize=8.5, color=INK2)
    legend(fig, 1.0)
    save(fig, "fig5_change_our_measures")

    figure6()


if __name__ == "__main__":
    main()
