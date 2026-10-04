#!/usr/bin/env python3
"""Gastric epithelial gastritis-to-cancer axis and additive virtual knockout.

Rules are locked in 01_方案/状态转换_需求锁.md. Do not retune after results.
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io, sparse

ROOT = Path("/root/autodl-tmp/gc_state_shift")
DATA = ROOT / "data"
RESULTS = ROOT / "results"

SEED = 20261003
N_PERM = 200
MIN_GENES = 200
MIN_COUNTS = 1000
MIN_LOG10_RATIO = 0.7
MAX_MT = 0.20
MAX_HB = 0.05
MIN_EPI_CELLS = 50
DETECT_FRAC = 0.10
MIN_POLE_PATIENTS = 4
MIN_IM_PATIENTS = 4
LOO_PASS_MIN = 5
P_CUT = 0.05

EPI_GENES = ["EPCAM", "KRT8", "KRT18", "KRT19"]
IMM_GENES = ["PTPRC", "CD3D", "CD79A", "MS4A1"]
FIB_GENES = ["COL1A1", "DCN", "LUM"]
HB_GENES = ["HBA1", "HBA2", "HBB", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"]
PRESPEC = ["HNF4G", "CDX2", "MUC2", "CLDN18"]

# library order matches GEO sample1-18. No personal names.
SAMPLES = [
    ("sample1", "GSM7966226", "GC", "neg"),
    ("sample2", "GSM7966227", "GC", "neg"),
    ("sample3", "GSM7966228", "GC", "neg"),
    ("sample4", "GSM7966229", "GC", "pos"),
    ("sample5", "GSM7966230", "GC", "pos"),
    ("sample6", "GSM7966231", "GC", "pos"),
    ("sample7", "GSM7966232", "GS", "neg"),
    ("sample8", "GSM7966233", "GS", "neg"),
    ("sample9", "GSM7966234", "GS", "neg"),
    ("sample10", "GSM7966235", "GS", "pos"),
    ("sample11", "GSM7966236", "GS", "pos"),
    ("sample12", "GSM7966237", "GS", "pos"),
    ("sample13", "GSM7966238", "IM", "neg"),
    ("sample14", "GSM7966239", "IM", "neg"),
    ("sample15", "GSM7966240", "IM", "neg"),
    ("sample16", "GSM7966241", "IM", "pos"),
    ("sample17", "GSM7966242", "IM", "pos"),
    ("sample18", "GSM7966243", "IM", "pos"),
]


BLOCK = 6_794_880


def load_adata() -> tuple[ad.AnnData, pd.DataFrame]:
    """Load the 10x matrix without materializing 122 million whitelist barcodes.

    Each library is a contiguous block of 6,794,880 barcodes, suffix -1 through -18.
    Cells below the locked count and gene thresholds cannot pass QC, so they are
    dropped before AnnData is built. Mitochondrial and hemoglobin gates still run
    on the remaining cells.
    """
    features = pd.read_csv(DATA / "features.tsv.gz", sep="\t", header=None)
    symbols = features.iloc[:, 1].astype(str) if features.shape[1] > 1 else features.iloc[:, 0].astype(str)
    print("reading matrix", flush=True)
    matrix = io.mmread(DATA / "matrix.mtx.gz").tocsr()
    if matrix.shape[1] != 18 * BLOCK and matrix.shape[0] == 18 * BLOCK:
        matrix = matrix.T.tocsr()
    if matrix.shape[1] != 18 * BLOCK:
        raise SystemExit(f"unexpected matrix shape {matrix.shape}")
    if matrix.shape[0] != len(symbols):
        raise SystemExit(f"features {len(symbols)} != genes {matrix.shape[0]}")
    n_counts = np.asarray(matrix.sum(axis=0)).ravel()
    n_genes = np.bincount(matrix.indices, minlength=matrix.shape[1])
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.log10(n_genes) / np.log10(n_counts)
    pre = (n_genes >= MIN_GENES) & (n_counts >= MIN_COUNTS) & (ratio >= MIN_LOG10_RATIO)
    detected = n_counts > 0
    rows = []
    for i, (library, gsm, stage, hp) in enumerate(SAMPLES):
        sl = slice(i * BLOCK, (i + 1) * BLOCK)
        rows.append(
            {
                "gsm": gsm,
                "stage": stage,
                "hp": hp,
                "n_detected": int(detected[sl].sum()),
                "n_prefilter": int(pre[sl].sum()),
            }
        )
    n_raw = pd.DataFrame(rows)
    idx = np.flatnonzero(pre)
    print(f"prefilter cells {len(idx)}", flush=True)
    sub = matrix[:, idx].T.tocsr()
    del matrix
    adata = ad.AnnData(X=sub)
    adata.var_names = pd.Index(symbols.to_numpy())
    adata.var_names_make_unique()
    adata.obs["library_index"] = (idx // BLOCK) + 1
    meta = pd.DataFrame(SAMPLES, columns=["library", "gsm", "stage", "hp"])
    meta["library_index"] = np.arange(1, 19)
    adata.obs = adata.obs.join(meta.set_index("library_index"), on="library_index")
    if adata.obs["gsm"].isna().any():
        raise SystemExit("some cells did not match a library")
    return adata, n_raw


def qc(adata: ad.AnnData) -> ad.AnnData:
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    hb = [g for g in HB_GENES if g in adata.var_names]
    adata.var["hb"] = adata.var_names.isin(hb)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "hb"], log1p=False, inplace=True)
    n_counts = np.asarray(adata.obs["total_counts"])
    n_genes = np.asarray(adata.obs["n_genes_by_counts"])
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.log10(n_genes) / np.log10(n_counts)
    keep = (
        (n_genes >= MIN_GENES)
        & (n_counts >= MIN_COUNTS)
        & (ratio >= MIN_LOG10_RATIO)
        & (np.asarray(adata.obs["pct_counts_mt"]) <= MAX_MT * 100)
        & (np.asarray(adata.obs["pct_counts_hb"]) <= MAX_HB * 100)
    )
    before = adata.n_obs
    adata = adata[keep].copy()
    print(f"qc {before} -> {adata.n_obs}", flush=True)
    return adata, before


def detect_fraction(block: sparse.spmatrix) -> np.ndarray:
    return np.asarray((block > 0).mean(axis=0)).ravel()


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    adata, n_raw = load_adata()
    adata, _n_prefilter = qc(adata)
    n_qc = (
        adata.obs.groupby(["gsm", "stage", "hp"], observed=True)
        .size()
        .rename("n_qc")
        .reset_index()
    )
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    for key, genes in [("epi", EPI_GENES), ("imm", IMM_GENES), ("fib", FIB_GENES)]:
        use = [g for g in genes if g in adata.var_names]
        if len(use) < 2:
            raise SystemExit(f"missing score genes for {key}: {genes}")
        sc.tl.score_genes(adata, use, score_name=f"{key}_score")
    epi = (
        (adata.obs["epi_score"] > 0)
        & (adata.obs["epi_score"] > adata.obs["imm_score"])
        & (adata.obs["epi_score"] > adata.obs["fib_score"])
    )
    adata = adata[epi].copy()
    print(f"epithelial cells {adata.n_obs}", flush=True)

    patients = []
    profiles = []
    detects = []
    genes = adata.var_names.to_numpy()
    for gsm, sub in adata.obs.groupby("gsm", observed=True):
        stage = sub["stage"].iloc[0]
        hp = sub["hp"].iloc[0]
        n_epi = len(sub)
        row = {"gsm": gsm, "stage": stage, "hp": hp, "n_epithelial": int(n_epi)}
        patients.append(row)
        if n_epi < MIN_EPI_CELLS:
            continue
        block = adata[sub.index].X
        if sparse.issparse(block):
            block = block.tocsr()
            mean = np.asarray(block.mean(axis=0)).ravel()
            det = detect_fraction(block)
        else:
            mean = np.asarray(block.mean(axis=0)).ravel()
            det = np.mean(block > 0, axis=0)
        profiles.append((gsm, mean))
        detects.append((gsm, det))

    patient_df = pd.DataFrame(patients)
    patient_df = patient_df.merge(
        n_raw[["gsm", "n_detected", "n_prefilter"]], on="gsm", how="left"
    )
    patient_df = patient_df.merge(n_qc[["gsm", "n_qc"]], on="gsm", how="left")
    order = [s[1] for s in SAMPLES]
    patient_df["gsm"] = pd.Categorical(patient_df["gsm"], order, ordered=True)
    patient_df = patient_df.sort_values("gsm")

    usable = patient_df[patient_df["n_epithelial"] >= MIN_EPI_CELLS]
    counts = usable["stage"].value_counts().to_dict()
    print("usable patients", counts, flush=True)
    if counts.get("GS", 0) != 6 or counts.get("GC", 0) != 6 or counts.get("IM", 0) != 6:
        patient_df.to_csv(RESULTS / "patient_qc.csv", index=False)
        summary = {
            "stopped": True,
            "reason": "a stage does not have 6 patients with at least 50 epithelial cells",
            "n_usable": {k: int(v) for k, v in counts.items()},
        }
        (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2))
        print("stopped", summary, flush=True)
        return

    prof = pd.DataFrame(
        {gsm: vec for gsm, vec in profiles},
        index=genes,
    ).T
    det = pd.DataFrame({gsm: vec for gsm, vec in detects}, index=genes).T
    prof = prof.loc[usable["gsm"].astype(str)]
    det = det.loc[usable["gsm"].astype(str)]
    meta = usable.set_index(usable["gsm"].astype(str))

    gs = meta.index[meta["stage"] == "GS"]
    gc = meta.index[meta["stage"] == "GC"]
    im = meta.index[meta["stage"] == "IM"]

    axis_ok = (det.loc[gs] >= DETECT_FRAC).sum(axis=0) >= MIN_POLE_PATIENTS
    axis_ok &= (det.loc[gc] >= DETECT_FRAC).sum(axis=0) >= MIN_POLE_PATIENTS
    banned = set(EPI_GENES + IMM_GENES + FIB_GENES + HB_GENES)
    banned |= {g for g in genes if str(g).upper().startswith("MT-")}
    keep_genes = [g for g in genes if axis_ok[g] and g not in banned]
    print(f"axis genes {len(keep_genes)}", flush=True)

    X = prof[keep_genes].to_numpy(dtype=np.float64)
    gene_index = {g: i for i, g in enumerate(keep_genes)}
    stage = meta.loc[prof.index, "stage"].to_numpy()
    gs_rows = np.where(stage == "GS")[0]
    gc_rows = np.where(stage == "GC")[0]
    im_rows = np.where(stage == "IM")[0]
    w = X[gc_rows].mean(axis=0) - X[gs_rows].mean(axis=0)
    norm = np.linalg.norm(w)
    if not np.isfinite(norm) or norm == 0:
        raise SystemExit("axis norm is zero")
    position = X @ w / norm
    patient_df["position"] = patient_df["gsm"].astype(str).map(
        dict(zip(prof.index, position))
    )
    patient_df.to_csv(RESULTS / "patient_qc.csv", index=False)

    pos_med = {
        st: float(np.median(position[stage == st])) for st in ["GS", "IM", "GC"]
    }
    im_between = pos_med["GS"] < pos_med["IM"] < pos_med["GC"] or pos_med["GC"] < pos_med["IM"] < pos_med["GS"]

    X_im = X[im_rows]
    observed = np.median(X_im * w / norm, axis=0)

    rng = np.random.default_rng(SEED)
    labels = np.array(["GS"] * 6 + ["GC"] * 6)
    pole_rows = np.concatenate([gs_rows, gc_rows])
    null = np.empty((N_PERM, X.shape[1]), dtype=np.float64)
    loo_null_p95 = np.empty((6, X.shape[1]), dtype=np.float64)
    for b in range(N_PERM):
        shuffled = rng.permutation(labels)
        fake_gc = pole_rows[shuffled == "GC"]
        fake_gs = pole_rows[shuffled == "GS"]
        w_b = X[fake_gc].mean(axis=0) - X[fake_gs].mean(axis=0)
        n_b = np.linalg.norm(w_b)
        if n_b == 0:
            null[b] = np.nan
            continue
        drop_b = X_im * w_b / n_b
        null[b] = np.median(drop_b, axis=0)
        for j in range(6):
            loo_null_p95[j]  # placeholder to keep shape; filled after the loop more cheaply
        if (b + 1) % 50 == 0:
            print(f"perm {b + 1}", flush=True)
    # LOO null percentiles need the per-subset medians. Recompute from stored weights
    # by a second pass that only keeps the 6 subset medians' 95th percentile.
    # Store is too big if we keep all. Second pass is 200 * 6 medians, cheap.
    loo_null = np.empty((N_PERM, 6, X.shape[1]), dtype=np.float64)
    rng2 = np.random.default_rng(SEED)
    for b in range(N_PERM):
        shuffled = rng2.permutation(labels)
        fake_gc = pole_rows[shuffled == "GC"]
        fake_gs = pole_rows[shuffled == "GS"]
        w_b = X[fake_gc].mean(axis=0) - X[fake_gs].mean(axis=0)
        n_b = np.linalg.norm(w_b)
        drop_b = X_im * w_b / n_b if n_b else np.full_like(X_im, np.nan)
        for j in range(6):
            loo_null[b, j] = np.median(np.delete(drop_b, j, axis=0), axis=0)
    # The block above duplicates the permutation. null[] already used the first stream.
    # loo must use the same shuffles as null. The second rng restarts at SEED, so the
    # shuffles match. Good.

    finite = np.isfinite(observed) & np.isfinite(null).all(axis=0)
    perm_p = np.full(observed.shape, np.nan)
    perm_p[finite] = (1.0 + np.sum(null[:, finite] >= observed[finite], axis=0)) / (N_PERM + 1)
    null_p95 = np.nanquantile(null, 0.95, axis=0)

    im_detect_n = (det.loc[im, keep_genes] >= DETECT_FRAC).sum(axis=0).to_numpy()
    loo_pass = np.zeros(len(keep_genes), dtype=int)
    loo_obs = np.empty((6, len(keep_genes)))
    for j in range(6):
        loo_obs[j] = np.median(np.delete(X_im * w / norm, j, axis=0), axis=0)
        p95_j = np.nanquantile(loo_null[:, j, :], 0.95, axis=0)
        loo_pass += (loo_obs[j] > p95_j).astype(int)

    retained = (
        (w > 0)
        & (im_detect_n >= MIN_IM_PATIENTS)
        & (perm_p <= P_CUT)
        & (loo_pass >= LOO_PASS_MIN)
    )
    summary_df = pd.DataFrame(
        {
            "gene": keep_genes,
            "weight": w,
            "im_median_drop": observed,
            "null_p95": null_p95,
            "perm_p": perm_p,
            "n_im_detected": im_detect_n.astype(int),
            "n_loo_pass": loo_pass,
            "retained": retained,
        }
    )
    # Prespec genes absent from the axis still get a row.
    have = set(summary_df["gene"])
    extra = []
    for gene in PRESPEC:
        if gene in have:
            continue
        extra.append(
            {
                "gene": gene,
                "weight": np.nan,
                "im_median_drop": np.nan,
                "null_p95": np.nan,
                "perm_p": np.nan,
                "n_im_detected": 0,
                "n_loo_pass": 0,
                "retained": False,
                "note": "not on the axis",
            }
        )
    if extra:
        summary_df["note"] = ""
        summary_df = pd.concat([summary_df, pd.DataFrame(extra)], ignore_index=True)
    else:
        summary_df["note"] = ""
    summary_df.loc[summary_df["gene"].isin(PRESPEC) & (summary_df["note"] == ""), "note"] = "prespec"
    summary_df = summary_df.sort_values(
        ["retained", "im_median_drop"], ascending=[False, False]
    )
    summary_df.to_csv(RESULTS / "gene_knockout.csv", index=False)

    prespec_rows = summary_df[summary_df["gene"].isin(PRESPEC)][
        ["gene", "weight", "im_median_drop", "perm_p", "n_loo_pass", "retained", "note"]
    ]
    hp = meta.loc[im, "hp"]
    retained_genes = summary_df.loc[summary_df["retained"] == True, "gene"].tolist()  # noqa: E712
    hp_rows = []
    for gene in retained_genes:
        if gene not in gene_index:
            continue
        j = gene_index[gene]
        drop_i = X_im[:, j] * w[j] / norm
        hp_rows.append(
            {
                "gene": gene,
                "n_hp_pos_drop_gt0": int(np.sum(drop_i[hp.to_numpy() == "pos"] > 0)),
                "n_hp_neg_drop_gt0": int(np.sum(drop_i[hp.to_numpy() == "neg"] > 0)),
            }
        )
    pd.DataFrame(hp_rows).to_csv(RESULTS / "retained_hp.csv", index=False)

    summary = {
        "n_cells_detected": int(n_raw["n_detected"].sum()),
        "n_epithelial": int(patient_df["n_epithelial"].sum()),
        "n_axis_genes": len(keep_genes),
        "n_tested_positive_weight": int(np.sum((w > 0) & (im_detect_n >= MIN_IM_PATIENTS))),
        "n_retained": int(np.sum(retained)),
        "position_median": pos_med,
        "im_between_gs_and_gc": bool(im_between),
        "prespec": prespec_rows.to_dict(orient="records"),
        "n_perm": N_PERM,
        "seed": SEED,
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(json.dumps(summary, indent=2, default=float), flush=True)


if __name__ == "__main__":
    main()
