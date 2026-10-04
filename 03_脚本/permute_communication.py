#!/usr/bin/env python3
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "06_真实分析结果" / "界面深化" / "scrna_cell_source.csv"
OUT = ROOT / "06_真实分析结果" / "腹膜通讯"
PAIRS = [
    ("C3_fib_C3AR1_mac", "C3__fib", "C3AR1__mac"),
    ("C3_epi_C3AR1_mac", "C3__epi", "C3AR1__mac"),
    ("SPP1_mac_CD44_mac", "SPP1__mac", "CD44__mac"),
    ("SPP1_mac_CD44_t", "SPP1__mac", "CD44__t"),
    ("CCL2_epi_CCR2_t", "CCL2__epi", "CCR2__t"),
    ("CCL2_fib_CCR2_mac", "CCL2__fib", "CCR2__mac"),
    ("TGFB1_t_TGFBR2_mac", "TGFB1__t", "TGFBR2__mac"),
]


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, idx in enumerate(order[::-1], start=1):
        i = n - rank
        val = p[idx] * n / (i + 1)
        prev = min(prev, val)
        q[idx] = prev
    return np.clip(q, 0, 1)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    wide = pd.read_csv(SRC)
    use = wide[wide["role"].isin(["peritoneal_tumor", "primary_tumor"])].copy()
    rows = []
    long_rows = []
    for name, lig, rec in PAIRS:
        score = use[lig] * use[rec]
        frame = use[["sample", "gsm", "role"]].copy()
        frame["pair"] = name
        frame["score"] = score
        finite = frame[np.isfinite(frame["score"])].reset_index(drop=True)
        pm = finite[finite["role"] == "peritoneal_tumor"]
        pr = finite[finite["role"] == "primary_tumor"]
        if len(pm) != 3:
            raise SystemExit(f"{name} peritoneal n={len(pm)}")
        obs = float(pm["score"].median() - pr["score"].median())
        values = finite["score"].to_numpy(float)
        k = 3
        n_comb = 0
        n_ge = 0
        for idx in combinations(range(len(values)), k):
            n_comb += 1
            chosen = values[list(idx)]
            rest = np.delete(values, idx)
            stat = float(np.median(chosen) - np.median(rest))
            if stat >= obs - 1e-12:
                n_ge += 1
        paired_pm = finite[(finite["sample"] == 38)]
        paired_pr = finite[(finite["sample"] == 36)]
        paired = float("nan")
        if len(paired_pm) == 1 and len(paired_pr) == 1:
            paired = float(paired_pm["score"].iloc[0] - paired_pr["score"].iloc[0])
        rows.append(
            {
                "pair": name,
                "pm_n": int(len(pm)),
                "primary_n": int(len(pr)),
                "n_comb": int(n_comb),
                "pm_median": float(pm["score"].median()),
                "primary_median": float(pr["score"].median()),
                "diff_median": obs,
                "n_pm_above_primary_median": int((pm["score"] > pr["score"].median()).sum()),
                "perm_p": n_ge / n_comb,
                "paired_38_minus_36": paired,
            }
        )
        long_rows.append(finite)
    table = pd.DataFrame(rows)
    table["bh_q"] = bh(table["perm_p"])
    table.to_csv(OUT / "communication_permutation.csv", index=False)
    pd.concat(long_rows, ignore_index=True).to_csv(OUT / "communication_scores.csv", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
