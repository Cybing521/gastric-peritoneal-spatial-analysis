#!/usr/bin/env python3
from __future__ import annotations

import gzip
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

ROOT = Path("/root/autodl-tmp/gc_spatial_niche")
DATA = ROOT / "data" / "public_pm"
OUT = ROOT / "results" / "public_pm"
OUT.mkdir(parents=True, exist_ok=True)

PROGRAMS = {
    "fib": ["DCN", "LUM", "COL1A1", "FAP"],
    "mac": ["CD68", "CSF1R", "C1QA", "C1QB"],
    "epi": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "t": ["CD3D", "CD3E", "TRAC"],
}
EXTRA = ["C3", "SPP1", "CCL2", "TGFB1", "CCR2", "C3AR1", "CD44", "TGFBR1", "TGFBR2", "STAT3", "SOCS3"]
PAIRS = [
    ("C3_fib_C3AR1_mac", "C3__fib", "C3AR1__mac"),
    ("C3_epi_C3AR1_mac", "C3__epi", "C3AR1__mac"),
    ("SPP1_mac_CD44_mac", "SPP1__mac", "CD44__mac"),
    ("SPP1_mac_CD44_t", "SPP1__mac", "CD44__t"),
    ("CCL2_epi_CCR2_t", "CCL2__epi", "CCR2__t"),
    ("CCL2_fib_CCR2_mac", "CCL2__fib", "CCR2__mac"),
    ("TGFB1_t_TGFBR2_mac", "TGFB1__t", "TGFBR2__mac"),
]
MIN_CLASS = 20

SAMPLES = [
    ("GSE308231", "GSM9240670", "Ca_1", "primary_tumor"),
    ("GSE308231", "GSM9240671", "Ca_2", "primary_tumor"),
    ("GSE308231", "GSM9240672", "Ca_3", "primary_tumor"),
    ("GSE308231", "GSM9240673", "F_14", "peritoneal_tumor"),
    ("GSE308231", "GSM9240674", "F_15", "peritoneal_tumor"),
    ("GSE308231", "GSM9240675", "F_16", "peritoneal_tumor"),
    ("GSE163558", "GSM5004180", "PT1", "primary_tumor"),
    ("GSE163558", "GSM5004181", "PT2", "primary_tumor"),
    ("GSE163558", "GSM5004182", "PT3", "primary_tumor"),
    ("GSE163558", "GSM5004187", "P1", "peritoneal_tumor"),
    ("GSE228598", "GSM7133741", "P1", "peritoneal_fluid"),
    ("GSE228598", "GSM7133742", "P2", "peritoneal_fluid"),
    ("GSE228598", "GSM7133743", "P3", "peritoneal_fluid"),
    ("GSE228598", "GSM7133744", "P4", "peritoneal_fluid"),
    ("GSE228598", "GSM7133745", "P5", "peritoneal_fluid"),
    ("GSE228598", "GSM7133746", "P6", "peritoneal_fluid"),
    ("GSE228598", "GSM7133747", "P7", "peritoneal_fluid"),
    ("GSE228598", "GSM7133748", "P8", "peritoneal_fluid"),
    ("GSE228598", "GSM7133749", "P9", "peritoneal_fluid"),
    ("GSE228598", "GSM7133750", "P10", "peritoneal_fluid"),
    ("GSE228598", "GSM7133751", "P11", "peritoneal_fluid"),
    ("GSE228598", "GSM7133752", "P12", "peritoneal_fluid"),
    ("GSE228598", "GSM7133753", "P13", "peritoneal_fluid"),
    ("GSE228598", "GSM7133754", "P14", "peritoneal_fluid"),
    ("GSE228598", "GSM7133755", "P16", "peritoneal_fluid"),
    ("GSE228598", "GSM7133756", "P17", "peritoneal_fluid"),
    ("GSE228598", "GSM7133757", "P18", "peritoneal_fluid"),
    ("GSE228598", "GSM7133758", "P19", "peritoneal_fluid"),
    ("GSE228598", "GSM7133759", "P20", "peritoneal_fluid"),
    ("GSE228598", "GSM7133760", "P21", "peritoneal_fluid"),
    ("GSE228598", "GSM7133761", "P22", "peritoneal_fluid"),
    ("GSE228598", "GSM7133762", "P23", "peritoneal_fluid"),
    ("GSE228598", "GSM7133763", "P24", "peritoneal_fluid"),
    ("GSE228598", "GSM7133764", "P25", "peritoneal_fluid"),
    ("GSE228598", "GSM7133765", "P26", "peritoneal_fluid"),
    ("GSE228598", "GSM7133766", "P27", "peritoneal_fluid"),
    ("GSE228598", "GSM7133767", "P28", "peritoneal_fluid"),
    ("GSE228598", "GSM7133768", "P29", "peritoneal_fluid"),
]


def log(msg: str) -> None:
    print(msg, flush=True)


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


def read_table(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        return pd.read_csv(handle, sep="\t", header=None)


def sample_dir(dataset: str) -> Path:
    return DATA / dataset


def load_counts(dataset: str, gsm: str, code: str):
    folder = sample_dir(dataset)
    prefix = f"{gsm}_{code}"
    matrix_path = folder / f"{prefix}_matrix.mtx.gz"
    feature_path = folder / f"{prefix}_features.tsv.gz"
    if not feature_path.exists():
        feature_path = folder / f"{prefix}_genes.tsv.gz"
    with gzip.open(matrix_path, "rb") as handle:
        matrix = mmread(handle).tocsr()
    features = read_table(feature_path)
    barcodes = read_table(folder / f"{prefix}_barcodes.tsv.gz")
    if matrix.shape[0] != len(features) and matrix.shape[1] == len(features):
        matrix = matrix.T.tocsr()
    if matrix.shape[0] != len(features) or matrix.shape[1] != len(barcodes):
        raise SystemExit(f"{prefix} shape {matrix.shape} genes {len(features)} cells {len(barcodes)}")
    symbols = features.iloc[:, 1].astype(str) if features.shape[1] > 1 else features.iloc[:, 0].astype(str)
    return matrix, symbols.to_numpy()


def class_row(matrix, symbols) -> dict:
    want = set(EXTRA)
    for genes in PROGRAMS.values():
        want.update(genes)
    positions = {}
    for i, gene in enumerate(symbols):
        if gene in want:
            positions.setdefault(gene, []).append(i)
    keep_rows = sorted({i for rows in positions.values() for i in rows})
    sub = matrix[keep_rows]
    row_gene = [symbols[i] for i in keep_rows]
    grouped = {}
    for gene, rows in positions.items():
        idx = [keep_rows.index(i) for i in rows]
        grouped[gene] = np.asarray(sub[idx].sum(axis=0)).ravel()
    lib = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    raw_probe = matrix.data[:10000] if matrix.data.size else np.array([0.0])
    integer_like = float(np.mean(np.abs(raw_probe - np.rint(raw_probe)))) < 1e-6
    peak = float(matrix.data.max()) if matrix.data.size else 0.0
    if integer_like and peak > 30:
        lib[lib == 0] = np.nan
        for gene in grouped:
            grouped[gene] = np.log1p(grouped[gene] / lib * 1e4)
    n = matrix.shape[1]
    scores = []
    labels = list(PROGRAMS)
    for label in labels:
        mats = [grouped[g] for g in PROGRAMS[label] if g in grouped]
        scores.append(np.vstack(mats).mean(axis=0) if mats else np.full(n, np.nan))
    stack = np.vstack(scores)
    winner = np.argmax(stack, axis=0)
    best = stack[winner, np.arange(n)]
    second = np.partition(stack, -2, axis=0)[-2]
    assigned = (best > 0) & (best > second) & np.isfinite(best)
    out = {"n_cells": int(n), "n_genes_used": int(len(grouped))}
    for i, label in enumerate(labels):
        mask = assigned & (winner == i)
        out[f"n_{label}"] = int(mask.sum())
        for gene in EXTRA:
            if gene not in grouped or int(mask.sum()) < MIN_CLASS:
                out[f"{gene}__{label}"] = float("nan")
            else:
                out[f"{gene}__{label}"] = float(np.nanmean(grouped[gene][mask]))
    out["n_unlabeled"] = int((~assigned).sum())
    return out


def add_products(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    for name, lig, rec in PAIRS:
        out[name] = out[lig] * out[rec]
    return out


def perm_table(frame: pd.DataFrame, pm_role: str, pr_role: str) -> pd.DataFrame:
    rows = []
    for name, _, _ in PAIRS:
        sub = frame[np.isfinite(frame[name])].reset_index(drop=True)
        pm = sub[sub["role"] == pm_role]
        pr = sub[sub["role"] == pr_role]
        if len(pm) == 0 or len(pr) == 0:
            continue
        obs = float(pm[name].median() - pr[name].median())
        values = sub[name].to_numpy(float)
        n_ge = 0
        n_comb = 0
        for idx in combinations(range(len(values)), len(pm)):
            n_comb += 1
            chosen = values[list(idx)]
            rest = np.delete(values, list(idx))
            if float(np.median(chosen) - np.median(rest)) >= obs - 1e-12:
                n_ge += 1
        rows.append(
            {
                "pair": name,
                "pm_n": int(len(pm)),
                "primary_n": int(len(pr)),
                "n_comb": int(n_comb),
                "pm_median": float(pm[name].median()),
                "primary_median": float(pr[name].median()),
                "diff_median": obs,
                "n_pm_above_primary_median": int((pm[name] > pr[name].median()).sum()),
                "perm_p": n_ge / n_comb,
            }
        )
    table = pd.DataFrame(rows)
    if len(table):
        table["bh_q"] = bh(table["perm_p"])
    return table


def direction_one(frame: pd.DataFrame) -> pd.DataFrame:
    pm = frame[frame["role"] == "peritoneal_tumor"]
    pr = frame[frame["role"] == "primary_tumor"]
    rows = []
    for name, _, _ in PAIRS:
        rows.append(
            {
                "pair": name,
                "pm_value": float(pm[name].iloc[0]) if len(pm) == 1 else float("nan"),
                "primary_median": float(pr[name].median()),
                "primary_n": int(pr[name].notna().sum()),
                "above_primary_median": bool(pm[name].iloc[0] > pr[name].median()) if len(pm) == 1 and np.isfinite(pm[name].iloc[0]) else False,
            }
        )
    return pd.DataFrame(rows)


def fluid_summary(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    cols = [name for name, _, _ in PAIRS] + [f"{gene}__{label}" for gene in ["C3", "C3AR1", "SPP1", "CD44", "CCL2", "CCR2", "STAT3", "SOCS3"] for label in PROGRAMS]
    for col in cols:
        vals = frame[col].to_numpy(float)
        finite = vals[np.isfinite(vals)]
        rows.append(
            {
                "measure": col,
                "n": int(len(finite)),
                "median": float(np.median(finite)) if len(finite) else float("nan"),
                "min": float(np.min(finite)) if len(finite) else float("nan"),
                "max": float(np.max(finite)) if len(finite) else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    rows = []
    for dataset, gsm, code, role in SAMPLES:
        log(f"{dataset} {gsm}_{code}")
        matrix, symbols = load_counts(dataset, gsm, code)
        row = {"dataset": dataset, "gsm": gsm, "sample": code, "role": role}
        row.update(class_row(matrix, symbols))
        rows.append(row)
        del matrix
    wide = add_products(pd.DataFrame(rows))
    wide.to_csv(OUT / "sample_class_means.csv", index=False)
    solid = wide[wide["dataset"] == "GSE308231"]
    perm_table(solid, "peritoneal_tumor", "primary_tumor").to_csv(OUT / "gse308231_permutation.csv", index=False)
    one = wide[wide["dataset"] == "GSE163558"]
    direction_one(one).to_csv(OUT / "gse163558_direction.csv", index=False)
    fluid = wide[wide["dataset"] == "GSE228598"]
    fluid_summary(fluid).to_csv(OUT / "gse228598_fluid_summary.csv", index=False)
    log(f"public pm samples {len(wide)}")
    log("public pm done")


if __name__ == "__main__":
    main()
