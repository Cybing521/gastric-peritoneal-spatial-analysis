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
RES = ROOT / "06_真实分析结果" / "腹膜通讯"
OUT = ROOT / "latex" / "analysis" / "figures"
ORDER = [
    "C3_fib_C3AR1_mac",
    "C3_epi_C3AR1_mac",
    "SPP1_mac_CD44_mac",
    "SPP1_mac_CD44_t",
    "CCL2_epi_CCR2_t",
    "CCL2_fib_CCR2_mac",
    "TGFB1_t_TGFBR2_mac",
]
LABELS = [
    "C3成纤维 × C3AR1巨噬",
    "C3上皮 × C3AR1巨噬",
    "SPP1巨噬 × CD44巨噬",
    "SPP1巨噬 × CD44 T细胞",
    "CCL2上皮 × CCR2 T细胞",
    "CCL2成纤维 × CCR2巨噬",
    "TGFB1 T细胞 × TGFBR2巨噬",
]


def main() -> None:
    scores = pd.read_csv(RES / "communication_scores.csv")
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    rng = np.random.default_rng(20261005)
    for i, pair in enumerate(ORDER):
        sub = scores[scores["pair"] == pair]
        y_pri = sub.loc[sub["role"] == "primary_tumor", "score"].to_numpy(float)
        y_pm = sub.loc[sub["role"] == "peritoneal_tumor", "score"].to_numpy(float)
        if len(y_pm) != 3:
            raise SystemExit(pair)
        ax.scatter(
            np.log1p(y_pri),
            i + rng.uniform(-0.12, 0.12, size=len(y_pri)),
            s=16,
            color=PALETTE["blue_main"],
            edgecolor="black",
            linewidth=0.3,
            zorder=2,
        )
        ax.scatter(
            np.log1p(y_pm),
            np.full(len(y_pm), i),
            s=36,
            color=PALETTE["red_strong"],
            edgecolor="black",
            linewidth=0.4,
            zorder=3,
        )
        ax.plot(
            [np.log1p(np.median(y_pri)), np.log1p(np.median(y_pri))],
            [i - 0.28, i + 0.28],
            color="black",
            lw=1.1,
            zorder=4,
        )
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels(LABELS, fontsize=7)
    ax.set_xlabel("log(1 + 类均值乘积)")
    ax.invert_yaxis()
    print(qa_before_save(fig))
    finalize_figure(fig, OUT / "fig_communication.pdf", formats=("pdf", "png"), dpi=300, tight=False, pad=0.2)


if __name__ == "__main__":
    apply_publication_style(FigureStyle(font_size=8, axes_linewidth=0.7, cjk=True))
    main()
