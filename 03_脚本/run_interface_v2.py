#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.optimize import nnls
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

ROOT = Path("/root/autodl-tmp/gc_spatial_niche")
DATA = ROOT / "data"
VISIUM = DATA / "visium"
OUT = ROOT / "results" / "interface_v2"
OUT.mkdir(parents=True, exist_ok=True)

RADIUS_PRIMARY = 150.0
RADII = (100.0, 150.0, 200.0, 300.0)
QUANTILES = (0.60, 0.75, 0.90)
SPOT_DIAMETER_UM = 55.0
N_PERM = 200
N_BOOT = 2000
SEED = 20261004
MIN_GENES = 200
MAX_MT_FRAC = 0.25
SHIFT_DIST = (200.0, 300.0, 400.0, 500.0)
SHIFT_MATCH_UM = 60.0
MIN_CLASS_CELLS = 20

PROGRAMS = {
    "fib": ["DCN", "LUM", "COL1A1", "FAP"],
    "mac": ["CD68", "CSF1R", "C1QA", "C1QB"],
    "epi": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "t": ["CD3D", "CD3E", "TRAC"],
}
STAT3_GENES = ["STAT3", "SOCS3"]
PAIRS = [
    ("CCL2", "CCR2", "prespec"),
    ("C3", "C3AR1", "prespec"),
    ("SPP1", "CD44", "prespec"),
    ("TGFB1", "TGFBR1", "prespec"),
    ("TGFB1", "TGFBR2", "prespec"),
    ("APP", "CD74", "borderline"),
    ("CD99", "CD99", "borderline"),
    ("SEMA4D", "PLXNB2", "borderline"),
    ("TGFB1", "TGFBR_MIN", "complex"),
]
SOURCE_GENES = ["C3", "SPP1", "CCL2", "TGFB1", "CCR2", "C3AR1", "CD44", "TGFBR1", "TGFBR2"]
SAMPLES = [
    ("GC1", "GSM7990473", "primary"),
    ("GC2", "GSM7990474", "primary"),
    ("GC3", "GSM7990475", "primary"),
    ("GC4", "GSM7990476", "primary"),
    ("GC5", "GSM7990477", "primary"),
    ("GC6", "GSM7990478", "primary"),
    ("GC6-PM", "GSM7990479", "metastasis_site_unspecified"),
    ("GC7", "GSM7990480", "primary"),
    ("GC8", "GSM7990481", "primary"),
    ("GC9", "GSM7990482", "primary"),
]


def log(msg: str) -> None:
    print(msg, flush=True)


def find_sample_dir(gsm: str) -> Path:
    hits = [p for p in VISIUM.iterdir() if p.is_dir() and gsm in p.name]
    if len(hits) != 1:
        raise FileNotFoundError(f"{gsm}: {hits}")
    return hits[0]


def find_file(root: Path, patterns: list[str]) -> Path:
    for pattern in patterns:
        hits = sorted(root.rglob(pattern))
        if hits:
            return hits[0]
    raise FileNotFoundError(f"{patterns} under {root}")


def load_section(sample: str, gsm: str) -> ad.AnnData:
    root = find_sample_dir(gsm)
    h5 = find_file(root, ["*filtered_feature_bc_matrix.h5"])
    adata = sc.read_10x_h5(h5)
    if (
        "gene_symbols" in adata.var.columns
        and pd.Index(adata.var_names.astype(str)).str.startswith("ENSG").mean() > 0.5
    ):
        adata.var_names = adata.var["gene_symbols"].astype(str)
    adata.var_names_make_unique()
    pos = find_file(root, ["*tissue_positions.csv", "*tissue_positions_list.csv"])
    sf = json.loads(find_file(root, ["*scalefactors_json.json"]).read_text())
    table = pd.read_csv(pos, header=None)
    if str(table.iloc[0, 0]) in {"barcode", "barcodes"}:
        table = pd.read_csv(pos)
        table.columns = [c.strip() for c in table.columns]
    else:
        table.columns = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
        ]
    table["barcode"] = table["barcode"].astype(str)
    table = table.set_index("barcode")
    common = adata.obs_names.intersection(table.index)
    adata = adata[common].copy()
    table = table.loc[adata.obs_names]
    adata = adata[(table["in_tissue"].astype(int) == 1).to_numpy()].copy()
    table = table.loc[adata.obs_names]
    mpp = SPOT_DIAMETER_UM / float(sf["spot_diameter_fullres"])
    adata.obsm["spatial_um"] = np.column_stack(
        [
            table["pxl_col_in_fullres"].to_numpy(float) * mpp,
            table["pxl_row_in_fullres"].to_numpy(float) * mpp,
        ]
    )
    adata.uns["sample"] = sample
    return adata


def qc_filter(adata: ad.AnnData) -> ad.AnnData:
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    keep = (adata.obs["n_genes_by_counts"] >= MIN_GENES) & (
        adata.obs["pct_counts_mt"] <= MAX_MT_FRAC * 100
    )
    return adata[keep].copy()


def gene_vector(adata: ad.AnnData, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.full(adata.n_obs, np.nan, dtype=np.float64)
    x = adata[:, gene].X
    if sparse.issparse(x):
        x = x.toarray()
    return np.asarray(x, dtype=np.float64).ravel()


def proportions(adata: ad.AnnData) -> np.ndarray:
    names = []
    blocks = []
    for label, genes in PROGRAMS.items():
        have = [g for g in genes if g in adata.var_names]
        if len(have) < 2:
            raise RuntimeError(f"{adata.uns['sample']} missing {label}: {genes}")
        names.append(label)
        blocks.append(have)
    order = []
    for genes in blocks:
        for gene in genes:
            if gene not in order:
                order.append(gene)
    design = np.zeros((len(order), len(blocks)), dtype=np.float64)
    index = {gene: i for i, gene in enumerate(order)}
    for j, genes in enumerate(blocks):
        for gene in genes:
            design[index[gene], j] = 1.0
    expr = np.zeros((adata.n_obs, len(order)), dtype=np.float64)
    for j, gene in enumerate(order):
        expr[:, j] = gene_vector(adata, gene)
    out = np.zeros((adata.n_obs, len(blocks)), dtype=np.float64)
    for i in range(adata.n_obs):
        coef, _ = nnls(design, expr[i])
        total = float(coef.sum())
        if total > 0:
            out[i] = coef / total
    return out


def adjacency(coords: np.ndarray, radius: float) -> sparse.csr_matrix:
    tree = cKDTree(coords)
    pairs = tree.query_pairs(radius, output_type="ndarray")
    n = coords.shape[0]
    if len(pairs) == 0:
        return sparse.csr_matrix((n, n), dtype=np.float64)
    rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
    cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
    data = np.ones(len(rows), dtype=np.float64)
    return sparse.csr_matrix((data, (rows, cols)), shape=(n, n))


def make_shifts(coords: np.ndarray, rng: np.random.Generator) -> list[tuple[np.ndarray, np.ndarray]]:
    tree = cKDTree(coords)
    shifts = []
    for _ in range(N_PERM):
        dist = float(rng.choice(SHIFT_DIST))
        ang = float(rng.uniform(0, 2 * np.pi))
        delta = np.array([np.cos(ang) * dist, np.sin(ang) * dist])
        dist_nn, index = tree.query(coords - delta, k=1)
        ok = dist_nn <= SHIFT_MATCH_UM
        shifts.append((index.astype(np.int32), ok))
    return shifts


def shifted(values: np.ndarray, index: np.ndarray, ok: np.ndarray) -> np.ndarray:
    out = np.full(values.shape[0], np.nan, dtype=np.float64)
    out[ok] = values[index[ok]]
    return out


def neighbor_mean(values: np.ndarray, adj: sparse.csr_matrix) -> np.ndarray:
    known = np.isfinite(values).astype(np.float64)
    filled = np.where(np.isfinite(values), values, 0.0)
    counts = np.asarray(adj @ known).ravel()
    total = np.asarray(adj @ filled).ravel()
    out = np.full(values.shape[0], np.nan, dtype=np.float64)
    ok = counts > 0
    out[ok] = total[ok] / counts[ok]
    return out


def weighted_mean(signal: np.ndarray, weight: np.ndarray) -> float:
    ok = np.isfinite(signal) & np.isfinite(weight) & (weight > 0)
    if int(ok.sum()) == 0 or float(weight[ok].sum()) == 0:
        return float("nan")
    return float(np.sum(weight[ok] * signal[ok]) / np.sum(weight[ok]))


def zero_weight_fraction(signal: np.ndarray, weight: np.ndarray) -> float:
    ok = np.isfinite(signal) & np.isfinite(weight) & (weight > 0)
    if int(ok.sum()) == 0 or float(weight[ok].sum()) == 0:
        return float("nan")
    return float(np.sum(weight[ok] * (signal[ok] == 0)) / np.sum(weight[ok]))


def null_summary(observed: float, nulls: np.ndarray) -> tuple[float, float, float]:
    finite = nulls[np.isfinite(nulls)]
    if not np.isfinite(observed) or len(finite) < N_PERM:
        return float("nan"), float("nan"), float("nan")
    p95 = float(np.quantile(finite, 0.95))
    p = float((1 + np.sum(finite >= observed)) / (len(finite) + 1))
    return p95, p, float(observed - p95)


def receptor_vector(adata: ad.AnnData, receptor: str) -> np.ndarray:
    if receptor == "TGFBR_MIN":
        a = gene_vector(adata, "TGFBR1")
        b = gene_vector(adata, "TGFBR2")
        return np.minimum(a, b)
    return gene_vector(adata, receptor)


def bootstrap_median(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float, float]:
    finite = values[np.isfinite(values)]
    if len(finite) == 0:
        return float("nan"), float("nan"), float("nan")
    draw = rng.integers(0, len(finite), size=(N_BOOT, len(finite)))
    med = np.median(finite[draw], axis=1)
    low, high = np.quantile(med, [0.025, 0.975])
    return float(np.median(finite)), float(low), float(high)


def analyze_sections() -> None:
    rng = np.random.default_rng(SEED)
    score_rows = []
    colo_rows = []
    binary_rows = []
    assoc_rows = []
    for sample, gsm, group in SAMPLES:
        log(f"load {sample}")
        adata = qc_filter(load_section(sample, gsm))
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        props = proportions(adata)
        fib = props[:, 0]
        mac = props[:, 1]
        coords = np.asarray(adata.obsm["spatial_um"], dtype=np.float64)
        weight = fib * mac
        shifts = make_shifts(coords, rng)
        match = float(np.mean([ok.mean() for _, ok in shifts]))
        log(f"  spots {adata.n_obs} shift-match {match:.3f}")
        ligand_cache = {}
        receptor_cache = {}
        for ligand, receptor, _panel in PAIRS:
            if ligand not in ligand_cache:
                ligand_cache[ligand] = gene_vector(adata, ligand)
            if receptor not in receptor_cache:
                receptor_cache[receptor] = receptor_vector(adata, receptor)
        adj_primary = None
        for radius in RADII:
            adj = adjacency(coords, radius)
            if radius == RADIUS_PRIMARY:
                adj_primary = adj
            mac_neighbor = neighbor_mean(mac, adj)
            both = np.isfinite(fib) & np.isfinite(mac_neighbor)
            observed_corr = float(np.corrcoef(fib[both], mac_neighbor[both])[0, 1]) if both.sum() > 5 else float("nan")
            null_corr = np.empty(N_PERM, dtype=np.float64)
            for b, (index, ok) in enumerate(shifts):
                moved = neighbor_mean(shifted(mac, index, ok), adj)
                mask = np.isfinite(fib) & np.isfinite(moved)
                null_corr[b] = float(np.corrcoef(fib[mask], moved[mask])[0, 1]) if mask.sum() > 5 else float("nan")
            p95, p, excess = null_summary(observed_corr, null_corr)
            colo_rows.append(
                {
                    "sample": sample,
                    "gsm": gsm,
                    "group": group,
                    "radius_um": radius,
                    "correlation": observed_corr,
                    "null_p95": p95,
                    "perm_p": p,
                    "excess": excess,
                    "shift_match": match,
                    "n_spots": int(adata.n_obs),
                }
            )
            for ligand, receptor, panel in PAIRS:
                lig = ligand_cache[ligand]
                rec = receptor_cache[receptor]
                signal = neighbor_mean(lig, adj) * rec
                observed = weighted_mean(signal, weight)
                nulls = np.empty(N_PERM, dtype=np.float64)
                for b, (index, ok) in enumerate(shifts):
                    moved = neighbor_mean(shifted(lig, index, ok), adj) * rec
                    nulls[b] = weighted_mean(moved, weight)
                p95, p, excess = null_summary(observed, nulls)
                score_rows.append(
                    {
                        "sample": sample,
                        "gsm": gsm,
                        "group": group,
                        "radius_um": radius,
                        "pair": f"{ligand}–{receptor}",
                        "panel": panel,
                        "score": observed,
                        "null_p95": p95,
                        "perm_p": p,
                        "excess": excess,
                        "zero_weight_fraction": zero_weight_fraction(signal, weight),
                        "weight_sum": float(np.nansum(weight)),
                        "shift_match": match,
                        "n_spots": int(adata.n_obs),
                    }
                )
                if radius == RADIUS_PRIMARY:
                    for q in QUANTILES:
                        mask = (fib >= np.quantile(fib, q)) & (mac >= np.quantile(mac, q))
                        usable = mask & np.isfinite(signal)
                        binary_rows.append(
                            {
                                "sample": sample,
                                "gsm": gsm,
                                "group": group,
                                "pair": f"{ligand}–{receptor}",
                                "quantile": q,
                                "n_spots": int(usable.sum()),
                                "median_signal": float(np.median(signal[usable])) if usable.sum() else float("nan"),
                                "mean_signal": float(np.mean(signal[usable])) if usable.sum() else float("nan"),
                            }
                        )
        stat3 = np.nanmean(np.column_stack([gene_vector(adata, g) for g in STAT3_GENES]), axis=1)
        focus = weight >= np.nanmedian(weight)
        for ligand, receptor, panel in PAIRS:
            if panel != "prespec" or ligand not in {"C3", "SPP1", "CCL2"}:
                continue
            signal = neighbor_mean(ligand_cache[ligand], adj_primary) * receptor_cache[receptor]
            mask = focus & np.isfinite(signal) & np.isfinite(stat3)
            if mask.sum() < 10 or np.nanstd(signal[mask]) == 0 or np.nanstd(stat3[mask]) == 0:
                rho, p = float("nan"), float("nan")
            else:
                rho, p = spearmanr(signal[mask], stat3[mask])
            assoc_rows.append(
                {
                    "sample": sample,
                    "gsm": gsm,
                    "group": group,
                    "pair": f"{ligand}–{receptor}",
                    "spearman": float(rho),
                    "spot_p": float(p),
                    "n_spots": int(mask.sum()),
                }
            )
    scores = pd.DataFrame(score_rows)
    colo = pd.DataFrame(colo_rows)
    binary = pd.DataFrame(binary_rows)
    assoc = pd.DataFrame(assoc_rows)
    scores.to_csv(OUT / "section_scores.csv", index=False)
    colo.to_csv(OUT / "colocalization.csv", index=False)
    binary.to_csv(OUT / "binary_sensitivity.csv", index=False)
    assoc.to_csv(OUT / "stat3_association.csv", index=False)
    boot_rng = np.random.default_rng(SEED)
    summary_rows = []
    primary = scores[(scores["group"] == "primary") & (scores["radius_um"] == RADIUS_PRIMARY)]
    for pair, sub in primary.groupby("pair"):
        med, low, high = bootstrap_median(sub["score"].to_numpy(float), boot_rng)
        emed, elow, ehigh = bootstrap_median(sub["excess"].to_numpy(float), boot_rng)
        summary_rows.append(
            {
                "pair": pair,
                "radius_um": RADIUS_PRIMARY,
                "n_sections": int(np.isfinite(sub["score"]).sum()),
                "median_score": med,
                "median_score_low": low,
                "median_score_high": high,
                "median_excess": emed,
                "median_excess_low": elow,
                "median_excess_high": ehigh,
                "n_sections_p_le_0.05": int(np.nansum(sub["perm_p"] <= 0.05)),
                "median_zero_weight_fraction": float(np.nanmedian(sub["zero_weight_fraction"])),
            }
        )
    pd.DataFrame(summary_rows).to_csv(OUT / "primary_summary.csv", index=False)
    log(f"spatial rows {len(scores)}")


def class_means(path: Path) -> dict:
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    values = frame.to_numpy(dtype=np.float64)
    finite = np.isfinite(values)
    integer_like = np.nanmean(np.abs(values[finite] - np.rint(values[finite]))) < 1e-6
    if integer_like and np.nanmax(values) > 30:
        lib = values.sum(axis=0)
        lib[lib == 0] = np.nan
        values = np.log1p(values / lib * 1e4)
    look = SOURCE_GENES + [g for genes in PROGRAMS.values() for g in genes]
    have = {gene: values[frame.index.get_loc(gene)] for gene in look if gene in frame.index}
    n = values.shape[1]
    scores = {}
    for label, genes in PROGRAMS.items():
        mats = [have[g] for g in genes if g in have]
        scores[label] = np.vstack(mats).mean(axis=0) if mats else np.full(n, np.nan)
    stack = np.vstack([scores[k] for k in PROGRAMS])
    winner = np.argmax(stack, axis=0)
    labels = np.array(list(PROGRAMS))
    best = stack[winner, np.arange(n)]
    second = np.partition(stack, -2, axis=0)[-2]
    assigned = (best > 0) & (best > second) & np.isfinite(best)
    out = {"n_cells": int(n)}
    for i, label in enumerate(labels):
        mask = assigned & (winner == i)
        out[f"n_{label}"] = int(mask.sum())
        for gene in SOURCE_GENES:
            if gene not in have or mask.sum() < MIN_CLASS_CELLS:
                out[f"{gene}__{label}"] = float("nan")
            else:
                out[f"{gene}__{label}"] = float(np.nanmean(have[gene][mask]))
    out["n_unlabeled"] = int((~assigned).sum())
    return out


def analyze_scrna() -> None:
    manifest = pd.read_csv(DATA / "scrna" / "manifest.tsv", sep="\t")
    rows = []
    for rec in manifest.itertuples(index=False):
        path = DATA / "scrna" / f"{rec.gsm}_sample{rec.sample}.csv.gz"
        if not path.exists():
            log(f"missing {path.name}")
            continue
        log(f"scrna {path.name}")
        row = {"sample": int(rec.sample), "gsm": rec.gsm, "role": rec.role}
        row.update(class_means(path))
        rows.append(row)
    wide = pd.DataFrame(rows)
    wide.to_csv(OUT / "scrna_cell_source.csv", index=False)
    pm = wide[wide["role"] == "peritoneal_tumor"]
    pt = wide[wide["role"] == "primary_tumor"]
    summary = []
    for gene in SOURCE_GENES:
        for label in PROGRAMS:
            col = f"{gene}__{label}"
            if col not in wide.columns:
                continue
            pm_vals = pm[col].to_numpy(float)
            pt_vals = pt[col].to_numpy(float)
            summary.append(
                {
                    "gene": gene,
                    "cell_class": label,
                    "pm_n": int(np.isfinite(pm_vals).sum()),
                    "primary_n": int(np.isfinite(pt_vals).sum()),
                    "pm_median": float(np.nanmedian(pm_vals)) if np.isfinite(pm_vals).any() else float("nan"),
                    "primary_median": float(np.nanmedian(pt_vals)) if np.isfinite(pt_vals).any() else float("nan"),
                }
            )
    pd.DataFrame(summary).to_csv(OUT / "scrna_source_summary.csv", index=False)
    log(f"scrna samples {len(wide)}")


if __name__ == "__main__":
    analyze_sections()
    analyze_scrna()
    log("interface v2 done")
