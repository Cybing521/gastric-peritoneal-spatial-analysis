#!/usr/bin/env python3
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, "/Users/cyibin/.agents/skills/scientific-figures/scripts")
from publication_style import PALETTE, FigureStyle, apply_publication_style, finalize_figure, qa_before_save

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "06_真实分析结果" / "公开腹膜" / "sample_class_means.csv"
OUT = ROOT / "latex" / "analysis" / "figures"
ORDER = ["C3_fib_C3AR1_mac", "C3_epi_C3AR1_mac", "SPP1_mac_CD44_mac"]
LABELS = ["C3成纤维 × C3AR1巨噬", "C3上皮 × C3AR1巨噬", "SPP1巨噬 × CD44巨噬"]


def main() -> None:
    frame = pd.read_csv(RES)
    sub = frame[frame["dataset"] == "GSE308231"]
    if len(sub) != 6:
        raise SystemExit(len(sub))
    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    rng = np.random.default_rng(20261005)
    for i, col in enumerate(ORDER):
        for role, color in (("primary_tumor", PALETTE["blue_main"]), ("peritoneal_tumor", PALETTE["red_strong"])):
            y = sub.loc[sub["role"] == role, col].to_numpy(float)
            ax.scatter(
                np.log1p(y),
                i + rng.uniform(-0.08, 0.08, size=len(y)),
                s=42,
                color=color,
                edgecolor="black",
                linewidth=0.4,
                zorder=3,
            )
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels(LABELS, fontsize=8)
    ax.set_xlabel("log(1 + 类均值乘积)")
    ax.invert_yaxis()
    print(qa_before_save(fig))
    finalize_figure(fig, OUT / "fig_public_pm.pdf", formats=("pdf", "png"), dpi=300, tight=False, pad=0.25)


if __name__ == "__main__":
    apply_publication_style(FigureStyle(font_size=8, axes_linewidth=0.7, cjk=True))
    main()
