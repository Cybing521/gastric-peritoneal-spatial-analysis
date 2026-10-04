#!/usr/bin/env python3
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager

font_manager.fontManager.addfont("/System/Library/Fonts/STHeiti Medium.ttc")
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "06_真实分析结果"
EN = ROOT / "latex" / "mdpi" / "figures"
ZH = ROOT / "latex" / "mdpi" / "figures_zh"

PRIMARY = "#0072B2"
CASE = "#D55E00"
NULL = "#4D4D4D"
HEAT_LOW = "#2166AC"
HEAT_MID = "#F7F7F7"
HEAT_HIGH = "#B2182B"

PAIRS = [
    "CCL2–CCR2",
    "C3–C3AR1",
    "SPP1–CD44",
    "TGFB1–TGFBR1",
    "TGFB1–TGFBR2",
    "TGFB1–TGFBR_MIN",
    "APP–CD74",
    "CD99–CD99",
    "SEMA4D–PLXNB2",
]
PRODUCTS = [
    "C3_fib_C3AR1_mac",
    "C3_epi_C3AR1_mac",
    "SPP1_mac_CD44_mac",
    "SPP1_mac_CD44_t",
    "CCL2_epi_CCR2_t",
    "CCL2_fib_CCR2_mac",
    "TGFB1_t_TGFBR2_mac",
]
GENES = ["C3", "SPP1", "CCL2", "TGFB1", "C3AR1", "CD44", "CCR2", "TGFBR2"]
CLASSES = ["fib", "mac", "epi", "t"]


def labels(lang):
    if lang == "zh":
        return {
            "pair": {
                "CCL2–CCR2": "CCL2–CCR2",
                "C3–C3AR1": "C3–C3AR1",
                "SPP1–CD44": "SPP1–CD44",
                "TGFB1–TGFBR1": "TGFB1–TGFBR1",
                "TGFB1–TGFBR2": "TGFB1–TGFBR2",
                "TGFB1–TGFBR_MIN": "TGFB1–TGFBR较小者",
                "APP–CD74": "APP–CD74",
                "CD99–CD99": "CD99–CD99",
                "SEMA4D–PLXNB2": "SEMA4D–PLXNB2",
            },
            "product": {
            "C3_fib_C3AR1_mac": "C3成纤维 × C3AR1",
            "C3_epi_C3AR1_mac": "C3上皮 × C3AR1",
            "SPP1_mac_CD44_mac": "SPP1 × CD44巨噬",
            "SPP1_mac_CD44_t": "SPP1 × CD44 T",
            "CCL2_epi_CCR2_t": "CCL2上皮 × CCR2",
            "CCL2_fib_CCR2_mac": "CCL2成纤维 × CCR2",
            "TGFB1_t_TGFBR2_mac": "TGFB1 T × TGFBR2",
            },
            "cls": {"fib": "成纤维", "mac": "巨噬", "epi": "上皮", "t": "T细胞"},
            "observed": "观察值",
            "null95": "标签打乱第95百分位",
            "contact": "接触比例",
            "n_above": "高于坐标零模型的原发切片数（共9）",
            "excess": "中位超额（评分 − 零模型第95百分位）",
            "radius": "半径（μm）",
            "heat_ex": "中位超额",
            "log_product": "log(1 + 类均值乘积)",
            "primary": "原发肿瘤",
            "pm": "腹膜肿瘤",
            "log2_ratio": "log2（腹膜中位数 / 原发中位数）",
            "diff": "腹膜中位数 − 原发中位数",
            "gse183": "GSE183904",
            "gse308": "GSE308231",
            "one_pm": "1例腹膜",
            "pri_med": "原发中位数",
            "stat3": "巨噬 STAT3",
            "socs3": "巨噬 SOCS3",
            "c3fib": "成纤维 C3",
            "c3ar": "巨噬 C3AR1",
            "c3epi": "上皮 C3",
        }
    return {
        "pair": {
            "CCL2–CCR2": "CCL2–CCR2",
            "C3–C3AR1": "C3–C3AR1",
            "SPP1–CD44": "SPP1–CD44",
            "TGFB1–TGFBR1": "TGFB1–TGFBR1",
            "TGFB1–TGFBR2": "TGFB1–TGFBR2",
            "TGFB1–TGFBR_MIN": "TGFB1–TGFBR min",
            "APP–CD74": "APP–CD74",
            "CD99–CD99": "CD99–CD99",
            "SEMA4D–PLXNB2": "SEMA4D–PLXNB2",
        },
        "product": {
            "C3_fib_C3AR1_mac": "C3 fib × C3AR1 mac",
            "C3_epi_C3AR1_mac": "C3 epi × C3AR1 mac",
            "SPP1_mac_CD44_mac": "SPP1 mac × CD44 mac",
            "SPP1_mac_CD44_t": "SPP1 mac × CD44 T",
            "CCL2_epi_CCR2_t": "CCL2 epi × CCR2 T",
            "CCL2_fib_CCR2_mac": "CCL2 fib × CCR2 mac",
            "TGFB1_t_TGFBR2_mac": "TGFB1 T × TGFBR2 mac",
        },
        "cls": {"fib": "Fibroblast", "mac": "Macrophage", "epi": "Epithelial", "t": "T cell"},
        "observed": "Observed",
        "null95": "Label-shuffle 95th percentile",
        "contact": "Contact fraction",
        "n_above": "Primary sections above the coordinate null (of 9)",
        "excess": "Median excess (score − null 95th percentile)",
        "radius": "Radius (μm)",
        "heat_ex": "Median excess",
        "log_product": "log(1 + class-mean product)",
        "primary": "Primary tumor",
        "pm": "Peritoneal tumor",
        "log2_ratio": "log2 (peritoneal median / primary median)",
        "diff": "Peritoneal median − primary median",
        "gse183": "GSE183904",
        "gse308": "GSE308231",
        "one_pm": "One peritoneal sample",
        "pri_med": "Primary median",
        "stat3": "Macrophage STAT3",
        "socs3": "Macrophage SOCS3",
        "c3fib": "Fibroblast C3",
        "c3ar": "Macrophage C3AR1",
        "c3epi": "Epithelial C3",
    }


def apply_style(lang):
    if lang == "zh":
        plt.rcParams["font.family"] = "Heiti SC" if "Heiti SC" in {f.name for f in font_manager.fontManager.ttflist} else "STHeiti"
        plt.rcParams["axes.unicode_minus"] = False
    else:
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    plt.rcParams.update({
        "font.size": 7,
        "axes.labelsize": 7.5,
        "axes.titlesize": 8,
        "axes.linewidth": 0.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "xtick.major.width": 0.4,
        "ytick.major.width": 0.4,
        "xtick.major.size": 2,
        "ytick.major.size": 2,
        "legend.fontsize": 6.5,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def tag(ax, letter):
    ax.text(-0.08, 1.06, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left", clip_on=False)


def save(fig, folder, stem):
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / f"{stem}.pdf")
    fig.savefig(folder / f"{stem}.png", dpi=300)
    plt.close(fig)


def load():
    contact = pd.read_csv(RES / "colocalization.csv")
    summary = pd.read_csv(RES / "pair_summary.csv")
    primary = pd.read_csv(RES / "界面深化" / "primary_summary.csv")
    scores = pd.read_csv(RES / "界面深化" / "section_scores.csv")
    source = pd.read_csv(RES / "界面深化" / "scrna_source_summary.csv")
    cells = pd.read_csv(RES / "界面深化" / "scrna_cell_source.csv")
    comm = pd.read_csv(RES / "腹膜通讯" / "communication_scores.csv")
    perm = pd.read_csv(RES / "腹膜通讯" / "communication_permutation.csv")
    stat3 = pd.read_csv(RES / "腹膜通讯" / "scrna_stat3_class.csv")
    g308 = pd.read_csv(RES / "公开腹膜" / "sample_class_means.csv")
    g308 = g308[g308["dataset"] == "GSE308231"].copy()
    g163 = pd.read_csv(RES / "公开腹膜" / "gse163558_direction.csv")
    return contact, summary, primary, scores, source, cells, comm, perm, stat3, g308, g163


def heatmap(ax, mat, row_labels, col_labels, cbar_label):
    finite = mat[np.isfinite(mat)]
    if finite.size == 0:
        raise RuntimeError("heatmap has no finite values")
    bound = float(np.max(np.abs(finite)))
    norm = TwoSlopeNorm(vmin=-bound, vcenter=0.0, vmax=bound)
    image = ax.imshow(mat, cmap="RdBu_r", norm=norm, aspect="auto")
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=30, ha="right")
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=6, width=0.4, length=2)
    cbar.set_label(cbar_label, fontsize=6.5)
    cbar.outline.set_linewidth(0.4)
    return image


def strip(ax, frame, value_col, ylabels, lab, legend=True):
    rng = np.random.default_rng(20261005)
    for i, key in enumerate(ylabels):
        sub = frame[frame["pair"] == key]
        pri = sub.loc[sub["role"] == "primary_tumor", value_col].to_numpy(float)
        pm = sub.loc[sub["role"] == "peritoneal_tumor", value_col].to_numpy(float)
        ax.scatter(np.log1p(pri), i + rng.uniform(-0.12, 0.12, len(pri)), s=12, color=PRIMARY, edgecolor="black", linewidth=0.2, zorder=2, label=lab["primary"] if i == 0 else None)
        ax.scatter(np.log1p(pm), np.full(len(pm), i), s=28, color=CASE, edgecolor="black", linewidth=0.3, zorder=3, label=lab["pm"] if i == 0 else None)
        med = float(np.log1p(np.median(pri)))
        ax.plot([med, med], [i - 0.28, i + 0.28], color="black", lw=0.8, zorder=4)
    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels([lab["product"][k] for k in ylabels])
    ax.set_xlabel(lab["log_product"])
    ax.invert_yaxis()
    if legend:
        ax.legend(loc="lower right")


def figure1(lang, contact, summary, primary, scores):
    lab = labels(lang)
    apply_style(lang)
    fig = plt.figure(figsize=(7.2, 8.0))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.05, 1.25, 1.25], hspace=0.55, wspace=0.55)
    ax = fig.add_subplot(gs[0, :])
    x = np.arange(len(contact))
    w = 0.38
    ax.bar(x - w / 2, contact["contact_fraction"], w, color=PRIMARY, edgecolor="black", linewidth=0.3, label=lab["observed"])
    ax.bar(x + w / 2, contact["null_p95"], w, color=NULL, edgecolor="black", linewidth=0.3, label=lab["null95"])
    ax.set_xticks(x)
    ax.set_xticklabels(contact["sample"], rotation=40, ha="right")
    low = float(min(contact["contact_fraction"].min(), contact["null_p95"].min()))
    high = float(max(contact["contact_fraction"].max(), contact["null_p95"].max()))
    ax.set_ylim(low - 0.02, high + 0.03)
    ax.set_ylabel(lab["contact"])
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, borderaxespad=0)
    tag(ax, "A")

    ax = fig.add_subplot(gs[1, 0])
    screen = [p for p in PAIRS if p in set(summary["pair"])]
    sub = summary.set_index("pair").loc[screen]
    y = np.arange(len(screen))
    vals = sub["n_sections_above_null"].to_numpy(float)
    ax.barh(y, vals, color=PRIMARY, edgecolor="black", linewidth=0.3, height=0.62)
    ax.scatter(vals, y, s=8, color=PRIMARY, zorder=3)
    ax.axvline(5, color=CASE, ls="--", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([lab["pair"][p] for p in screen])
    ax.set_xlabel(lab["n_above"])
    ax.set_xlim(0, 9.6)
    ax.invert_yaxis()
    tag(ax, "B")

    ax = fig.add_subplot(gs[1, 1])
    pri = primary.set_index("pair").loc[PAIRS]
    for i, pair in enumerate(PAIRS):
        est = float(pri.loc[pair, "median_excess"])
        lo = float(pri.loc[pair, "median_excess_low"])
        hi = float(pri.loc[pair, "median_excess_high"])
        if not (lo <= est <= hi):
            raise RuntimeError(f"interval order failed for {pair}")
        ax.plot([lo, hi], [i, i], color="black", lw=0.7)
        ax.scatter([est], [i], s=18, color=PRIMARY, edgecolor="black", linewidth=0.3, zorder=3)
    ax.axvline(0, color=CASE, ls="--", lw=0.7)
    ax.set_yticks(range(len(PAIRS)))
    ax.set_yticklabels([lab["pair"][p] for p in PAIRS])
    ax.set_xlabel(lab["excess"])
    ax.invert_yaxis()
    tag(ax, "C")

    ax = fig.add_subplot(gs[2, :])
    use = scores[scores["group"] == "primary"]
    mat = np.full((len(PAIRS), 4), np.nan)
    radii = [100, 150, 200, 300]
    for i, pair in enumerate(PAIRS):
        for j, radius in enumerate(radii):
            vals = use.loc[(use["pair"] == pair) & (use["radius_um"] == radius), "excess"]
            mat[i, j] = float(np.median(vals))
    heatmap(ax, mat, [lab["pair"][p] for p in PAIRS], [str(r) for r in radii], lab["heat_ex"])
    ax.set_xlabel(lab["radius"])
    tag(ax, "D")
    save(fig, EN if lang == "en" else ZH, "fig1_spatial")


def source_matrix(source):
    mat = np.full((len(GENES), len(CLASSES)), np.nan)
    for i, gene in enumerate(GENES):
        for j, cls in enumerate(CLASSES):
            row = source[(source["gene"] == gene) & (source["cell_class"] == cls)]
            if row.empty:
                continue
            pm = float(row["pm_median"].iloc[0])
            pri = float(row["primary_median"].iloc[0])
            if pm > 0 and pri > 0:
                mat[i, j] = np.log2(pm / pri)
    return mat


def figure2(lang, source, comm):
    lab = labels(lang)
    apply_style(lang)
    fig = plt.figure(figsize=(7.2, 7.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.35], hspace=0.42, wspace=0.55)
    ax = fig.add_subplot(gs[0, 0])
    mat = source_matrix(source)
    heatmap(ax, mat, GENES, [lab["cls"][c] for c in CLASSES], lab["log2_ratio"])
    tag(ax, "A")
    ax = fig.add_subplot(gs[0, 1])
    show = [("C3", "fib"), ("C3", "epi"), ("C3", "mac"), ("C3AR1", "mac"), ("SPP1", "mac"), ("CD44", "mac"), ("CCL2", "fib"), ("TGFB1", "t")]
    bars = []
    names = []
    for gene, cls in show:
        row = source[(source["gene"] == gene) & (source["cell_class"] == cls)]
        if row.empty:
            continue
        pm = float(row["pm_median"].iloc[0])
        pri = float(row["primary_median"].iloc[0])
        if pm > 0 and pri > 0:
            bars.append(np.log2(pm / pri))
            names.append(f"{gene} {lab['cls'][cls]}")
    order = np.argsort(bars)
    ax.barh(np.arange(len(order)), np.array(bars)[order], color=[PRIMARY if bars[i] >= 0 else CASE for i in order], edgecolor="black", linewidth=0.25, height=0.7)
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels([names[i] for i in order], fontsize=6.5)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_xlabel(lab["log2_ratio"])
    ax.invert_yaxis()
    tag(ax, "B")
    ax = fig.add_subplot(gs[1, :])
    strip(ax, comm, "score", PRODUCTS, lab)
    tag(ax, "C")
    save(fig, EN if lang == "en" else ZH, "fig2_source")


def figure3(lang, g308, perm183, perm308, g163):
    lab = labels(lang)
    apply_style(lang)
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 4.2), gridspec_kw={"wspace": 0.28})
    g308 = g308.copy()
    g308["pair_role"] = g308["role"]
    long = []
    for pair in PRODUCTS:
        for _, row in g308.iterrows():
            long.append({"pair": pair, "role": row["role"], "score": row[pair]})
    strip(axes[0], pd.DataFrame(long), "score", PRODUCTS, lab, legend=False)
    tag(axes[0], "A")

    ax = axes[1]
    y = np.arange(len(PRODUCTS))
    a = perm183.set_index("pair").loc[PRODUCTS, "diff_median"].to_numpy(float)
    b = perm308.set_index("pair").loc[PRODUCTS, "diff_median"].to_numpy(float)
    ax.scatter(a, y + 0.12, s=22, marker="o", color=PRIMARY, edgecolor="black", linewidth=0.3, label=lab["gse183"], zorder=3)
    ax.scatter(b, y - 0.12, s=22, marker="s", color=CASE, edgecolor="black", linewidth=0.3, label=lab["gse308"], zorder=3)
    for i in range(len(PRODUCTS)):
        ax.plot([a[i], b[i]], [y[i] + 0.12, y[i] - 0.12], color="#BDBDBD", lw=0.6, zorder=1)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlabel(lab["diff"])
    ax.invert_yaxis()
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=1)
    tag(ax, "B")

    ax = axes[2]
    y = np.arange(len(PRODUCTS))
    one = g163.set_index("pair").loc[PRODUCTS]
    ax.scatter(np.log1p(one["primary_median"]), y, s=22, color=PRIMARY, edgecolor="black", linewidth=0.3, label=lab["pri_med"], zorder=3)
    ax.scatter(np.log1p(one["pm_value"]), y, s=28, color=CASE, edgecolor="black", linewidth=0.3, label=lab["one_pm"], zorder=3)
    for i, pair in enumerate(PRODUCTS):
        ax.plot([np.log1p(one.loc[pair, "primary_median"]), np.log1p(one.loc[pair, "pm_value"])], [i, i], color="#BDBDBD", lw=0.6, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlabel(lab["log_product"])
    ax.invert_yaxis()
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=1)
    tag(ax, "C")
    axes[0].legend(
        handles=[
            plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=PRIMARY, markeredgecolor="black", markersize=5, label=lab["primary"]),
            plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=CASE, markeredgecolor="black", markersize=6, label=lab["pm"]),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=1,
    )
    save(fig, EN if lang == "en" else ZH, "fig3_replication")


def figure4(lang, cells, stat3):
    lab = labels(lang)
    apply_style(lang)
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8), gridspec_kw={"hspace": 0.4, "wspace": 0.35})
    pri = cells[cells["role"] == "primary_tumor"]
    panels = [
        (axes[0, 0], "A", "C3__fib", "C3AR1__mac", "n_fib", "n_mac", lab["c3fib"], lab["c3ar"]),
        (axes[0, 1], "B", "C3__epi", "C3AR1__mac", "n_epi", "n_mac", lab["c3epi"], lab["c3ar"]),
    ]
    stats = {}
    for ax, letter, xcol, ycol, nx, ny, xl, yl in panels:
        use = pri[(pri[nx] >= 20) & (pri[ny] >= 20)][[xcol, ycol]].dropna()
        rho, p = spearmanr(use[xcol], use[ycol])
        ax.scatter(use[xcol], use[ycol], s=16, color=PRIMARY, edgecolor="black", linewidth=0.25)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.text(0.04, 0.96, f"Spearman {rho:.2f}\nn = {len(use)}", transform=ax.transAxes, va="top", fontsize=6.5)
        tag(ax, letter)
        stats[f"{xcol}_vs_{ycol}"] = {"rho": float(rho), "n": int(len(use)), "p": float(p)}
    rng = np.random.default_rng(20261005)
    for ax, letter, col, title in (
        (axes[1, 0], "C", "STAT3__mac", lab["stat3"]),
        (axes[1, 1], "D", "SOCS3__mac", lab["socs3"]),
    ):
        for xpos, role, color, name in ((0, "primary_tumor", PRIMARY, lab["primary"]), (1, "peritoneal_tumor", CASE, lab["pm"])):
            vals = stat3.loc[stat3["role"] == role, col].to_numpy(float)
            ax.scatter(np.full(len(vals), xpos) + rng.uniform(-0.08, 0.08, len(vals)), vals, s=16, color=color, edgecolor="black", linewidth=0.25, label=name)
            ax.plot([-0.18 + xpos, 0.18 + xpos], [np.median(vals), np.median(vals)], color="black", lw=0.9)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([lab["primary"], lab["pm"]])
        ax.set_ylabel(title)
        ax.legend(loc="best")
        tag(ax, letter)
    save(fig, EN if lang == "en" else ZH, "fig4_ligand")
    return stats


def main():
    contact, summary, primary, scores, source, cells, comm, perm, stat3, g308, g163 = load()
    perm308 = pd.read_csv(RES / "公开腹膜" / "gse308231_permutation.csv")
    for lang in ("en", "zh"):
        figure1(lang, contact, summary, primary, scores)
        figure2(lang, source, comm)
        figure3(lang, g308, perm, perm308, g163)
        stats = figure4(lang, cells, stat3)
    mac = stat3[stat3["n_mac"] >= 20]
    out = {
        "spearman": stats,
        "stat3_mac_pm": float(mac.loc[mac["role"] == "peritoneal_tumor", "STAT3__mac"].median()),
        "stat3_mac_pri": float(mac.loc[mac["role"] == "primary_tumor", "STAT3__mac"].median()),
        "socs3_mac_pm": float(mac.loc[mac["role"] == "peritoneal_tumor", "SOCS3__mac"].median()),
        "socs3_mac_pri": float(mac.loc[mac["role"] == "primary_tumor", "SOCS3__mac"].median()),
    }
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
