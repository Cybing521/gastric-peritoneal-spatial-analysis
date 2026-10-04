#!/usr/bin/env python3
"""Report figures from verified result tables. case_id: compare_bar (ImmunoStruct plot_bars)."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path.home() / ".agents/skills/scientific-figures/scripts"))
from house_flowchart import apply_flowchart_style, draw_vertical_workflow
from publication_style import (
    PALETTE,
    FigureStyle,
    apply_publication_style,
    create_subplots,
    finalize_figure,
    make_grouped_bar,
    qa_before_save,
)

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "06_真实分析结果"
OUT = ROOT / "latex" / "analysis" / "figures"
OUT.mkdir(parents=True, exist_ok=True)


def save(fig, stem: str) -> None:
    qa_before_save(fig, require_cjk_font=True)
    finalize_figure(fig, OUT / stem, formats=("pdf", "png"), dpi=300)


def fig_contact() -> None:
    apply_publication_style(FigureStyle(font_size=12, axes_linewidth=1.4, cjk=True))
    colo = pd.read_csv(RES / "colocalization.csv")
    fig, axes = create_subplots(1, 1, figsize=(12.8, 4.4))
    ax = axes[0]
    make_grouped_bar(
        ax,
        colo["sample"].tolist(),
        [colo["contact_fraction"].tolist(), colo["null_p95"].tolist()],
        ["观察值", "标签打乱第95百分位"],
        ylabel="接触比例",
        colors=[PALETTE["blue_main"], PALETTE["red_strong"]],
        edgecolor="black",
        linewidth=1.1,
    )
    ax.set_ylim(0.75, 1.02)
    ax.tick_params(axis="x", labelrotation=35)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.16), ncol=2, frameon=False)
    save(fig, "fig_contact")


def fig_retention() -> None:
    apply_publication_style(FigureStyle(font_size=12, axes_linewidth=1.4, cjk=True))
    summary = pd.read_csv(RES / "pair_summary.csv")
    wanted = [
        "APP–CD74",
        "CD99–CD99",
        "SEMA4D–PLXNB2",
        "TGFB1–TGFBR2",
        "TGFB1–TGFBR1",
        "C3–C3AR1",
        "SPP1–CD44",
        "CCL2–CCR2",
    ]
    sub = summary.set_index("pair").loc[wanted].reset_index()
    sub = sub.iloc[::-1]
    labels = [p.replace("–", "-") for p in sub["pair"]]
    vals = sub["n_sections_above_null"].to_numpy(float)
    fig, axes = create_subplots(1, 1, figsize=(8.4, 5.0))
    ax = axes[0]
    y = np.arange(len(labels))
    ax.barh(
        y,
        vals,
        color=PALETTE["blue_main"],
        edgecolor="black",
        linewidth=1.1,
        height=0.62,
    )
    ax.axvline(5, color=PALETTE["red_strong"], linestyle="--", linewidth=1.3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(0, 9.4)
    ax.set_xlabel("超过坐标打乱空分布的原发灶切片数（共9例）")
    save(fig, "fig_retention")


def fig_peritoneal() -> None:
    apply_publication_style(FigureStyle(font_size=12, axes_linewidth=1.4, cjk=True))
    direction = pd.read_csv(RES / "peritoneal_direction.csv")
    genes = ["C3", "SPP1", "CCL2", "CCR2", "TGFBR2", "TGFBR1", "TGFB1", "CD44"]
    sub = direction.drop_duplicates("gene").set_index("gene").loc[genes].iloc[::-1]
    fig, axes = create_subplots(1, 1, figsize=(8.2, 4.8))
    ax = axes[0]
    vals = sub["log2fc_median"].to_numpy(float)
    y = np.arange(len(sub))
    colors = [PALETTE["blue_main"] if v >= 0 else PALETTE["red_strong"] for v in vals]
    ax.barh(y, vals, color=colors, edgecolor="black", linewidth=1.1, height=0.62)
    ax.axvline(0, color="black", linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(sub.index.tolist())
    ax.set_xlabel("腹膜肿瘤相对原发肿瘤的 log2 中位数倍数")
    save(fig, "fig_peritoneal")


def fig_flow() -> None:
    apply_flowchart_style(lang="zh", base_size=9)
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = [
        "PingFang SC",
        "Heiti SC",
        "Songti SC",
        "Noto Sans CJK SC",
        "Arial Unicode MS",
    ]
    stages = [
        ("1  质控与细胞评分", "组织内点位\n基因数≥200，线粒体比例≤25%\n成纤维、巨噬评分各取前25%"),
        ("2  邻域接触", "半径150 μm\n对照：保持两类数目，打乱标签200次"),
        ("3  配体受体虚拟阻断", "界面中位数 = 受体 × 成纤维邻居配体\n对照：打乱坐标200次"),
        ("4  腹膜方向核对", "3份腹膜肿瘤，26份原发肿瘤\n比较样本均值，不参与空间排序"),
    ]
    notes = [
        "GSE251950\n9例原发灶进入排序",
        "判定：观察值高于\n第95百分位",
        "保留：≥5/9例，\n且留一后仍≥5/8",
        "GSE183904\nGC6-PM不参与排序",
    ]
    for suffix in ("pdf", "png"):
        draw_vertical_workflow(
            stages,
            notes,
            title="",
            fig_w=9.8,
            fig_h=7.6,
            out=OUT / f"fig_flow.{suffix}",
            dpi=300,
        )
    plt.close("all")


if __name__ == "__main__":
    fig_contact()
    fig_retention()
    fig_peritoneal()
    fig_flow()
    print("wrote", OUT)
