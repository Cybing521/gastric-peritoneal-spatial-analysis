#!/usr/bin/env python3
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
sys.path.insert(0, "/Users/cyibin/.agents/skills/scientific-figures/scripts")
from publication_style import (
    PALETTE,
    FigureStyle,
    apply_publication_style,
    finalize_figure,
    qa_before_save,
)

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "06_真实分析结果" / "界面深化"
OUT = ROOT / "latex" / "analysis" / "figures"
ORDER = [
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
LABELS = [
    "CCL2–CCR2",
    "C3–C3AR1",
    "SPP1–CD44",
    "TGFB1–TGFBR1",
    "TGFB1–TGFBR2",
    "TGFB1–TGFBR复合物",
    "APP–CD74",
    "CD99–CD99",
    "SEMA4D–PLXNB2",
]


def tag(ax, letter):
    ax.text(-0.12, 1.04, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", ha="left", va="bottom")


def main() -> None:
    scores = pd.read_csv(RES / "section_scores.csv")
    colo = pd.read_csv(RES / "colocalization.csv")
    prim = scores[(scores["group"] == "primary") & (scores["radius_um"] == 150)].copy()
    if prim.groupby("pair").size().min() != 9:
        raise SystemExit("primary section count is not 9")
    prim["ratio"] = prim["score"] / prim["null_p95"]
    if not np.isfinite(prim["ratio"]).all():
        raise SystemExit("non-finite score ratio")
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 5.0), gridspec_kw={"width_ratios": [1.15, 1.0], "wspace": 0.55})
    ax = axes[0]
    rng = np.random.default_rng(20261004)
    for i, pair in enumerate(ORDER):
        y = prim.loc[prim["pair"] == pair, "ratio"].to_numpy(float)
        if len(y) != 9:
            raise SystemExit(pair)
        ax.scatter(
            y,
            i + rng.uniform(-0.12, 0.12, size=len(y)),
            s=26,
            color=PALETTE["blue_main"],
            edgecolor="black",
            linewidth=0.4,
            zorder=3,
        )
        ax.plot([np.median(y), np.median(y)], [i - 0.28, i + 0.28], color="black", lw=1.2, zorder=4)
    ax.axvline(1.0, color=PALETTE["red_strong"], lw=0.8, ls=(0, (3, 2)))
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels(LABELS, fontsize=7)
    ax.set_xlabel("加权均值 / 平移空分布第95百分位")
    ax.set_xlim(0.55, 1.28)
    ax.invert_yaxis()
    tag(ax, "A")
    sub = colo[(colo["group"] == "primary") & (colo["radius_um"] == 150)].sort_values("sample")
    if len(sub) != 9:
        raise SystemExit("colocalization rows")
    ax2 = axes[1]
    x = np.arange(len(sub))
    ax2.scatter(x, sub["correlation"], s=32, color=PALETTE["blue_main"], edgecolor="black", linewidth=0.4, zorder=3)
    ax2.scatter(x, sub["null_p95"], s=32, facecolor="white", edgecolor=PALETTE["red_strong"], linewidth=1.0, zorder=3)
    ax2.set_xticks(x)
    ax2.set_xticklabels(sub["sample"], rotation=40, ha="right", fontsize=6.5)
    ax2.set_ylabel("成纤维比例与邻域巨噬比例的相关")
    tag(ax2, "B")
    print(qa_before_save(fig))
    finalize_figure(fig, OUT / "fig_interface.pdf", formats=("pdf", "png"), dpi=300, tight=False, pad=0.15)


if __name__ == "__main__":
    apply_publication_style(FigureStyle(font_size=8, axes_linewidth=0.7, cjk=True))
    main()
