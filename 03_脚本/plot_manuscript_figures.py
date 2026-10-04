#!/usr/bin/env python3
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D

sys.path.insert(0, "/Users/cyibin/.agents/skills/scientific-figures/scripts")
from publication_style import (
    PALETTE,
    FigureStyle,
    apply_publication_style,
    finalize_figure,
    qa_before_save,
)

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "06_真实分析结果" / "state_shift"
OUT = ROOT / "latex" / "manuscript" / "figures"
STAGE_ORDER = ["GS", "IM", "GC"]
STAGE_NAME = {"GS": "Gastritis", "IM": "Metaplasia", "GC": "Cancer"}
STAGE_COLOR = {
    "GS": PALETTE["blue_main"],
    "IM": PALETTE["teal"],
    "GC": PALETTE["red_strong"],
}


def tag(ax, letter, x=-0.14, y=1.04):
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        ha="left",
        va="bottom",
        clip_on=False,
    )


def load():
    patients = pd.read_csv(RES / "patient_qc.csv")
    genes = pd.read_csv(RES / "gene_knockout.csv")
    if len(patients) != 18:
        raise SystemExit(f"expected 18 patients, got {len(patients)}")
    if patients.groupby("stage").size().to_dict() != {"GC": 6, "GS": 6, "IM": 6}:
        raise SystemExit("stage counts are not 6/6/6")
    tested = genes[(genes["weight"] > 0) & (genes["n_im_detected"] >= 4)].copy()
    if len(tested) != 4781:
        raise SystemExit(f"expected 4781 tested genes, got {len(tested)}")
    if int(tested["retained"].sum()) != 416:
        raise SystemExit("retained count is not 416")
    cols = ["im_median_drop", "null_p95", "perm_p", "weight"]
    if not np.isfinite(tested[cols].to_numpy(dtype=float)).all():
        raise SystemExit("non-finite values in tested genes")
    if (tested["im_median_drop"] < 0).any() or (tested["null_p95"] < 0).any():
        raise SystemExit("negative displacement in the tested set")
    return patients, tested


def ecdf(values):
    ordered = np.sort(np.asarray(values, dtype=float))
    if np.any(ordered <= 0):
        raise SystemExit("ECDF has a non-positive displacement")
    frac = np.arange(1, len(ordered) + 1) / len(ordered)
    return ordered, frac


def fig_position(patients):
    gs_min = float(patients.loc[patients.stage == "GS", "position"].min())
    gs_max = float(patients.loc[patients.stage == "GS", "position"].max())
    gc_min = float(patients.loc[patients.stage == "GC", "position"].min())
    fig = plt.figure(figsize=(6.35, 6.7))
    outer = fig.add_gridspec(1, 2, width_ratios=[1.18, 1.0], wspace=0.46)
    left = outer[0].subgridspec(3, 1, hspace=0.28)
    right = outer[1].subgridspec(2, 1, height_ratios=[1.12, 1.0], hspace=0.42)
    axes = []
    for i in range(3):
        ax = fig.add_subplot(left[i], sharex=axes[0] if axes else None)
        axes.append(ax)
    axb = fig.add_subplot(right[0])
    axc = fig.add_subplot(right[1])

    for i, stage in enumerate(STAGE_ORDER):
        ax = axes[i]
        sub = patients.loc[patients.stage == stage].sort_values("position")
        y = np.arange(len(sub))
        color = STAGE_COLOR[stage]
        ax.axvspan(gs_min, gs_max, color=PALETTE["blue_main"], alpha=0.12, lw=0, zorder=0)
        ax.axvline(gc_min, color=PALETTE["red_strong"], lw=0.7, ls=(0, (3, 2)), zorder=1)
        ax.axvline(float(np.median(sub["position"])), color="black", lw=0.85, zorder=2)
        for hp, marker in (("neg", "s"), ("pos", "o")):
            mask = sub["hp"].to_numpy() == hp
            ax.scatter(
                sub.loc[mask, "position"],
                y[mask],
                s=38,
                marker=marker,
                facecolor=color,
                edgecolor="black",
                linewidth=0.45,
                zorder=3,
            )
        ax.set_yticks(y)
        ax.set_yticklabels(sub["gsm"].tolist(), fontsize=7)
        ax.set_ylim(-0.65, 5.65)
        ax.set_title(STAGE_NAME[stage], loc="left", fontsize=9, color=color, pad=1)
        if i < 2:
            ax.tick_params(axis="x", labelbottom=False)
    axes[0].set_xlim(6.4, 20.4)
    axes[2].set_xlabel("Epithelial projection")
    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="0.25", markeredgecolor="black", markersize=5.5, label="H. pylori positive"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="0.25", markeredgecolor="black", markersize=5.5, label="H. pylori negative"),
        Line2D([0], [0], color="black", lw=0.85, label="Group median"),
        Line2D([0], [0], color=PALETTE["red_strong"], lw=0.7, ls=(0, (3, 2)), label="Lowest cancer"),
    ]
    axes[2].legend(
        handles=legend_handles,
        loc="center",
        frameon=False,
        fontsize=6.5,
        borderaxespad=0.1,
        handlelength=1.5,
    )
    tag(axes[0], "A", x=-0.22, y=1.08)

    for stage in STAGE_ORDER:
        sub = patients.loc[patients.stage == stage]
        for hp, marker in (("neg", "s"), ("pos", "o")):
            mask = sub["hp"] == hp
            axb.scatter(
                sub.loc[mask, "position"],
                sub.loc[mask, "n_epithelial"],
                s=36,
                marker=marker,
                facecolor=STAGE_COLOR[stage],
                edgecolor="black",
                linewidth=0.45,
                zorder=3,
            )
    axb.axvspan(gs_min, gs_max, color=PALETTE["blue_main"], alpha=0.12, lw=0, zorder=0)
    axb.axvline(gc_min, color=PALETTE["red_strong"], lw=0.7, ls=(0, (3, 2)), zorder=1)
    axb.set_yscale("log")
    axb.set_xlim(6.4, 20.4)
    axb.set_xlabel("Epithelial projection")
    axb.set_ylabel("Epithelial cells")
    stage_handles = [
        Line2D(
            [0], [0],
            marker="o",
            color="none",
            markerfacecolor=STAGE_COLOR[stage],
            markeredgecolor="black",
            markersize=5.5,
            label=STAGE_NAME[stage],
        )
        for stage in STAGE_ORDER
    ]
    axb.legend(
        handles=stage_handles,
        loc="lower left",
        bbox_to_anchor=(0.18, 1.02),
        ncol=3,
        frameon=False,
        fontsize=6.5,
        borderaxespad=0,
        columnspacing=0.7,
        handletextpad=0.25,
    )
    tag(axb, "B", x=-0.16, y=1.16)

    tested_drop = pd.read_csv(RES / "gene_knockout.csv")
    tested_drop = tested_drop[(tested_drop["weight"] > 0) & (tested_drop["n_im_detected"] >= 4)]
    retained = tested_drop.loc[tested_drop["retained"] == True, "im_median_drop"]  # noqa: E712
    xs, ys = ecdf(tested_drop["im_median_drop"])
    xr, yr = ecdf(retained)
    axc.step(xs, ys, where="post", color="0.35", lw=1.15, label="Tested, 4,781")
    axc.step(xr, yr, where="post", color=PALETTE["red_strong"], lw=1.35, label="Retained, 416")
    axc.set_xscale("log")
    axc.set_xlim(float(xs.min()) * 0.7, 0.35)
    axc.set_ylim(0, 1.04)
    axc.set_xlabel("Median displacement")
    axc.set_ylabel("Cumulative fraction")
    axc.legend(loc="upper left", frameon=False, fontsize=6.5)
    tag(axc, "C", x=-0.16, y=1.04)

    print("position", qa_before_save(fig))
    finalize_figure(fig, OUT / "fig_position.pdf", formats=("pdf", "png"), dpi=300, tight=False, pad=0.12)


def fig_knockout(tested):
    fig = plt.figure(figsize=(6.35, 6.85))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.08, 0.92], hspace=0.40, wspace=0.38)
    axv = fig.add_subplot(grid[0, 0])
    axh = fig.add_subplot(grid[0, 1])
    axd = fig.add_subplot(grid[1, :])

    drop = tested["im_median_drop"].to_numpy(dtype=float)
    logp = -np.log10(tested["perm_p"].to_numpy(dtype=float))
    kept = tested["retained"].to_numpy(dtype=bool)
    axv.scatter(drop[~kept], logp[~kept], s=8, c=PALETTE["neutral"], alpha=0.7, linewidths=0, rasterized=True, zorder=2)
    axv.scatter(drop[kept], logp[kept], s=13, c=PALETTE["red_strong"], alpha=0.9, linewidths=0, rasterized=True, zorder=3)
    axv.axhline(-np.log10(0.05), color="0.2", lw=0.6, ls=(0, (3, 2)), zorder=1)
    axv.set_xlim(-0.004, 0.275)
    axv.set_ylim(0.12, 2.58)
    axv.set_xlabel("Median displacement")
    axv.set_ylabel(r"$-\log_{10}$ permutation $p$")
    texts = []
    for gene in ["HSP90AA1", "HSP90AB1", "RPS8", "RPL13", "LYZ", "CLDN18", "MDM4"]:
        row = tested.loc[tested["gene"] == gene]
        if len(row) != 1:
            raise SystemExit(f"label gene missing: {gene}")
        row = row.iloc[0]
        texts.append(
            axv.text(
                float(row.im_median_drop),
                -np.log10(float(row.perm_p)),
                gene,
                fontsize=7,
                fontstyle="italic",
                zorder=4,
            )
        )
    adjust_text(
        texts,
        ax=axv,
        arrowprops=dict(arrowstyle="-", color="0.35", lw=0.4),
        expand=(1.15, 1.25),
        ensure_inside_axes=True,
    )
    axv.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["red_strong"], markersize=5, label="Retained"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["neutral"], markersize=5, label="Not retained"),
        ],
        loc="lower right",
        frameon=False,
        fontsize=6.5,
    )
    tag(axv, "A", x=-0.18, y=1.04)

    cmap = LinearSegmentedColormap.from_list(
        "null_density",
        ["#F4F7FB", PALETTE["blue_secondary"], PALETTE["blue_main"]],
    )
    bins = axh.hexbin(
        tested["null_p95"],
        tested["im_median_drop"],
        gridsize=26,
        cmap=cmap,
        mincnt=1,
        norm=LogNorm(vmin=1),
        linewidths=0,
        extent=(0, 0.27, 0, 0.27),
        zorder=1,
    )
    axh.plot([0, 0.27], [0, 0.27], color="0.15", lw=0.6, zorder=2)
    kept_df = tested.loc[kept]
    axh.scatter(
        kept_df["null_p95"],
        kept_df["im_median_drop"],
        s=11,
        c=PALETTE["red_strong"],
        linewidths=0,
        rasterized=True,
        zorder=3,
    )
    texts = []
    for gene in ["HSP90AA1", "LYZ", "CLDN18", "RPL13"]:
        row = tested.loc[tested["gene"] == gene].iloc[0]
        texts.append(
            axh.text(
                float(row.null_p95),
                float(row.im_median_drop),
                gene,
                fontsize=7,
                fontstyle="italic",
                zorder=4,
            )
        )
    adjust_text(
        texts,
        ax=axh,
        arrowprops=dict(arrowstyle="-", color="0.35", lw=0.4),
        expand=(1.2, 1.3),
        ensure_inside_axes=True,
    )
    axh.set_xlim(0, 0.27)
    axh.set_ylim(0, 0.27)
    axh.set_xlabel("Null 95th percentile")
    axh.set_ylabel("Median displacement")
    colorbar = fig.colorbar(bins, ax=axh, fraction=0.046, pad=0.03)
    colorbar.set_label("Genes in bin", fontsize=7)
    colorbar.ax.tick_params(labelsize=6)
    tag(axh, "B", x=-0.22, y=1.04)

    top = kept_df.nlargest(12, "im_median_drop").iloc[::-1]
    if len(top) != 12:
        raise SystemExit("top-12 panel is short")
    if float(top["null_p95"].max()) >= float(top["im_median_drop"].max()):
        raise SystemExit("a null percentile exceeds the largest observed displacement")
    y = np.arange(len(top))
    axd.hlines(y, top["null_p95"], top["im_median_drop"], color="0.4", lw=0.9, zorder=2)
    axd.scatter(top["null_p95"], y, s=30, facecolor="white", edgecolor="black", linewidth=0.55, zorder=3)
    axd.scatter(
        top["im_median_drop"],
        y,
        s=30,
        facecolor=PALETTE["red_strong"],
        edgecolor="black",
        linewidth=0.4,
        zorder=4,
    )
    axd.set_yticks(y)
    axd.set_yticklabels(top["gene"].tolist(), fontsize=8, fontstyle="italic")
    axd.set_xlabel("Position units")
    axd.set_xlim(0.09, 0.265)
    axd.set_ylim(-0.7, 11.7)
    axd.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="black", markersize=5.5, label="Null 95th percentile"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["red_strong"], markeredgecolor="black", markersize=5.5, label="Observed median"),
        ],
        loc="lower right",
        frameon=False,
        fontsize=7,
    )
    tag(axd, "C", x=-0.07, y=1.02)

    print("knockout", qa_before_save(fig))
    finalize_figure(fig, OUT / "fig_knockout.pdf", formats=("pdf", "png"), dpi=300, tight=False, pad=0.12)


def fig_axis() -> None:
    genes = pd.read_csv(RES / "gene_knockout.csv")
    axis = genes[genes["note"].fillna("") != "not on the axis"].copy()
    if len(axis) != 7411:
        raise SystemExit(f"expected 7411 axis genes, got {len(axis)}")
    weight = axis["weight"].to_numpy(dtype=float)
    if not np.isfinite(weight).all():
        raise SystemExit("non-finite axis weights")
    length = float(np.linalg.norm(weight))
    if abs(length - 6.907033615435267) > 1e-8:
        raise SystemExit(f"axis length changed: {length}")
    squared = weight ** 2
    share = squared / squared.sum()
    order = np.argsort(share)[::-1]
    cumulative = np.cumsum(share[order])
    if abs(float(cumulative[3]) - 0.10271946986377774) > 1e-8:
        raise SystemExit("four-gene share changed")
    if abs(float(cumulative[153]) - 0.49999572134840253) > 1e-8:
        raise SystemExit("154-gene share changed")

    fig = plt.figure(figsize=(6.35, 4.7))
    grid = fig.add_gridspec(1, 2, width_ratios=[1.02, 1.18], wspace=0.42)
    axa = fig.add_subplot(grid[0])
    axb = fig.add_subplot(grid[1])
    rank = np.arange(1, len(cumulative) + 1)
    axa.plot(rank, cumulative, color=PALETTE["blue_main"], lw=1.35)
    axa.set_xscale("log")
    axa.set_xlim(1, 7411)
    axa.set_ylim(0, 1.04)
    axa.set_xlabel("Genes ranked by squared weight")
    axa.set_ylabel("Cumulative share of axis length")
    marks = ((4, 0.22, 0.10), (154, 0.55, 0.16), (1553, 0.78, -0.02))
    labels = ("4 genes, 10.3%", "154 genes, half", "1,553 genes, 90%")
    for (n, y_text, x_shift), label in zip(marks, labels):
        y = float(cumulative[n - 1])
        axa.scatter([n], [y], s=16, color=PALETTE["red_strong"], zorder=3)
        axa.annotate(
            label,
            xy=(n, y),
            xytext=(n * (10 ** x_shift), y_text),
            fontsize=6.5,
            arrowprops=dict(arrowstyle="-", color="0.4", lw=0.4),
            ha="center",
            va="center",
        )
    tag(axa, "A", x=-0.18, y=1.04)

    signed = np.sign(weight) * squared / length
    axis["signed"] = signed
    top = axis.assign(magnitude=np.abs(signed)).nlargest(16, "magnitude").sort_values("signed")
    if len(top) != 16:
        raise SystemExit("contribution panel is short")
    y = np.arange(len(top))
    colors = [PALETTE["blue_main"] if v < 0 else PALETTE["red_strong"] for v in top["signed"]]
    axb.barh(y, top["signed"], color=colors, edgecolor="black", linewidth=0.35, height=0.72)
    axb.axvline(0, color="black", lw=0.6)
    axb.set_yticks(y)
    axb.set_yticklabels(top["gene"].tolist(), fontsize=7, fontstyle="italic")
    axb.set_xlabel("Contribution to the mean gap")
    axb.set_xlim(-0.36, 0.36)
    tag(axb, "B", x=-0.28, y=1.04)

    print("axis", qa_before_save(fig))
    finalize_figure(fig, OUT / "fig_axis.pdf", formats=("pdf", "png"), dpi=300, tight=False, pad=0.12)


if __name__ == "__main__":
    apply_publication_style(FigureStyle(font_size=8, axes_linewidth=0.7, cjk=False))
    patients, tested = load()
    fig_position(patients)
    fig_knockout(tested)
    fig_axis()
    print("figures ok")
