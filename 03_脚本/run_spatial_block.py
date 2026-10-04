#!/usr/bin/env python3
"""Primary gastric-cancer Visium: fibroblast-macrophage interface LR scores.

Endpoint, locked before looking at results:
  For each ligand-receptor pair, the section score is the median, over
  macrophage-high spots that have at least one fibroblast-high neighbor
  within 150 um, of (receptor at the spot) * (mean ligand on those
  fibroblast-high neighbors). Expression is log1p of counts per 10k.
  Under this additive receiver model, removing the pair drops the score
  by exactly that quantity. The null shuffles spot coordinates (200
  times, seed 20261002) and takes the 95th percentile.

A primary section supports spatial contact when the fraction of
macrophage-high spots with a fibroblast-high neighbor exceeds the 95th
percentile of a label shuffle that preserves both counts. A pair is
above the spatial null in a section when its score exceeds that
section's null 95th percentile. Across the 9 primary sections, a pair
is retained when this happens in at least 5 sections, and it remains
so in leave-one-section-out (at least 5 of the remaining 8).

GC6-PM is scored with the same locked rules and is not used in the
ranking. CellChat heteromeric complexes (underscore in the gene symbol)
are not scored.
"""

from __future__ import annotations

import json
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.spatial import cKDTree

ROOT = Path("/root/autodl-tmp/gc_spatial_niche")
DATA = ROOT / "data"
VISIUM = DATA / "visium"
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

RADIUS_UM = 150.0
SPOT_DIAMETER_UM = 55.0
HIGH_Q = 0.75
N_PERM = 200
SEED = 20261002
MIN_GENES = 200
MAX_MT_FRAC = 0.25
MIN_INTERFACE = 10
DETECT_FRAC = 0.05
MIN_SECTIONS_DETECT = 6
MIN_SECTIONS_ABOVE = 5

FIB_GENES = ["DCN", "LUM", "COL1A1", "FAP"]
MAC_GENES = ["CD68", "CSF1R", "C1QA", "C1QB"]
PRESPEC = [
    ("CCL2", "CCR2", "prespec"),
    ("C3", "C3AR1", "prespec"),
    ("SPP1", "CD44", "prespec"),
    ("TGFB1", "TGFBR1", "prespec"),
    ("TGFB1", "TGFBR2", "prespec"),
]

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


def unpack_visium() -> None:
    tar_path = DATA / "GSE251950_RAW.tar"
    marker = DATA / "visium_unpacked.ok"
    if marker.exists():
        log("visium already unpacked")
        return
    log(f"unpack {tar_path}")
    VISIUM.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path) as outer:
        outer.extractall(DATA / "visium_tar", filter="data")
    nested = list((DATA / "visium_tar").rglob("*.tar.gz")) + list((DATA / "visium_tar").rglob("*.tgz"))
    if not nested:
        nested = [tar_path]
        src_root = DATA / "visium_tar"
    else:
        src_root = None
    for item in nested:
        dest = VISIUM / item.name.replace(".tar.gz", "").replace(".tgz", "")
        dest.mkdir(parents=True, exist_ok=True)
        with tarfile.open(item) as inner:
            inner.extractall(dest, filter="data")
        log(f"  unpacked {item.name}")
    if src_root is not None:
        pass
    marker.write_text("ok\n")


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
    sf_path = find_file(root, ["*scalefactors_json.json"])
    sf = json.loads(sf_path.read_text())
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
    in_tissue = table["in_tissue"].astype(int) == 1
    adata = adata[in_tissue.to_numpy()].copy()
    table = table.loc[adata.obs_names]
    mpp = SPOT_DIAMETER_UM / float(sf["spot_diameter_fullres"])
    xy = np.column_stack(
        [
            table["pxl_col_in_fullres"].to_numpy(float) * mpp,
            table["pxl_row_in_fullres"].to_numpy(float) * mpp,
        ]
    )
    adata.obsm["spatial_um"] = xy
    adata.uns["sample"] = sample
    adata.uns["microns_per_pixel"] = mpp
    return adata


def qc_filter(adata: ad.AnnData) -> tuple[ad.AnnData, dict]:
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    n_before = adata.n_obs
    keep = (adata.obs["n_genes_by_counts"] >= MIN_GENES) & (
        adata.obs["pct_counts_mt"] <= MAX_MT_FRAC * 100
    )
    info = {
        "n_before": int(n_before),
        "n_after": int(keep.sum()),
        "min_genes": MIN_GENES,
        "max_pct_mt": MAX_MT_FRAC * 100,
        "median_counts": float(adata.obs["total_counts"].median()),
        "median_pct_mt": float(adata.obs["pct_counts_mt"].median()),
    }
    return adata[keep].copy(), info


def normalize(adata: ad.AnnData) -> None:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)


def present_genes(adata: ad.AnnData, genes: list[str]) -> list[str]:
    have = set(adata.var_names)
    return [g for g in genes if g in have]


def high_mask(adata: ad.AnnData, genes: list[str], key: str) -> np.ndarray:
    use = present_genes(adata, genes)
    if len(use) < 2:
        raise RuntimeError(f"{adata.uns['sample']} missing genes for {key}: {genes}")
    sc.tl.score_genes(adata, use, score_name=key, use_raw=False)
    thr = float(np.quantile(adata.obs[key].to_numpy(), HIGH_Q))
    adata.uns[f"{key}_threshold"] = thr
    adata.uns[f"{key}_genes"] = use
    return adata.obs[key].to_numpy() >= thr


def radius_adjacency(coords: np.ndarray, radius: float) -> sparse.csr_matrix:
    tree = cKDTree(coords)
    pairs = tree.query_pairs(radius, output_type="ndarray")
    n = coords.shape[0]
    if len(pairs) == 0:
        return sparse.csr_matrix((n, n), dtype=np.float32)
    rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
    cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
    data = np.ones(len(rows), dtype=np.float32)
    return sparse.csr_matrix((data, (rows, cols)), shape=(n, n))


def expr_matrix(adata: ad.AnnData, genes: list[str]) -> np.ndarray:
    idx = adata.var_names.get_indexer(genes)
    if np.any(idx < 0):
        missing = [g for g, i in zip(genes, idx) if i < 0]
        raise KeyError(missing)
    x = adata[:, genes].X
    if sparse.issparse(x):
        x = x.toarray()
    return np.asarray(x, dtype=np.float32)


def pair_signal(
    ligand: np.ndarray,
    receptor: np.ndarray,
    fib: np.ndarray,
    mac: np.ndarray,
    adj: sparse.csr_matrix,
) -> tuple[np.ndarray, int]:
    fib_f = fib.astype(np.float32)
    neigh_counts = np.asarray(adj @ fib_f).ravel()
    weighted = ligand * fib_f[:, None]
    neigh_sum = np.asarray(adj @ weighted)
    mean_l = np.zeros_like(neigh_sum, dtype=np.float32)
    ok = neigh_counts > 0
    mean_l[ok] = neigh_sum[ok] / neigh_counts[ok, None]
    raw = mean_l * receptor
    interface = mac & ok
    n_int = int(interface.sum())
    if n_int < MIN_INTERFACE:
        return np.full(ligand.shape[1], np.nan, dtype=np.float64), n_int
    return np.median(raw[interface], axis=0).astype(np.float64), n_int


def contact_fraction(fib: np.ndarray, mac: np.ndarray, adj: sparse.csr_matrix) -> tuple[float, int, int]:
    neigh_counts = np.asarray(adj @ fib.astype(np.float32)).ravel()
    n_mac = int(mac.sum())
    if n_mac == 0:
        return float("nan"), 0, 0
    n_contact = int((mac & (neigh_counts > 0)).sum())
    return n_contact / n_mac, n_contact, n_mac


def load_lr_table(gene_names: set[str]) -> pd.DataFrame:
    path = DATA / "cellchat_lr.csv"
    table = pd.read_csv(path)
    table["ligand"] = table["ligand"].astype(str)
    table["receptor"] = table["receptor"].astype(str)
    table = table[~table["ligand"].str.contains("_") & ~table["receptor"].str.contains("_")]
    table = table[table["ligand"].isin(gene_names) & table["receptor"].isin(gene_names)]
    table = table.drop_duplicates(["ligand", "receptor"])
    pres = pd.DataFrame(PRESPEC, columns=["ligand", "receptor", "source"])
    table["source"] = "cellchat_single_gene"
    keys = set(zip(table["ligand"], table["receptor"]))
    extra = pres[~pres.apply(lambda r: (r["ligand"], r["receptor"]) in keys, axis=1)]
    extra = extra[extra["ligand"].isin(gene_names) & extra["receptor"].isin(gene_names)]
    table = pd.concat([table, extra.drop(columns=[])], ignore_index=True)
    # mark prespec
    pres_keys = set(zip(pres["ligand"], pres["receptor"]))
    table["panel"] = [
        "prespec" if (a, b) in pres_keys else "cellchat"
        for a, b in zip(table["ligand"], table["receptor"])
    ]
    table = table.drop_duplicates(["ligand", "receptor"]).reset_index(drop=True)
    table["pair"] = table["ligand"] + "–" + table["receptor"]
    return table


def gene_universe(sections: list[ad.AnnData]) -> set[str]:
    names = None
    for adata in sections:
        cur = set(adata.var_names)
        names = cur if names is None else names & cur
    return names or set()


def analyze_spatial() -> None:
    log("loading sections")
    loaded = []
    qc_rows = []
    for sample, gsm, group in SAMPLES:
        log(f"  load {sample}")
        adata = load_section(sample, gsm)
        adata, info = qc_filter(adata)
        normalize(adata)
        info.update({"sample": sample, "gsm": gsm, "group": group})
        fib = high_mask(adata, FIB_GENES, "fib_score")
        mac = high_mask(adata, MAC_GENES, "mac_score")
        info["n_fib_high"] = int(fib.sum())
        info["n_mac_high"] = int(mac.sum())
        info["fib_genes"] = ",".join(adata.uns["fib_score_genes"])
        info["mac_genes"] = ",".join(adata.uns["mac_score_genes"])
        qc_rows.append(info)
        loaded.append((sample, gsm, group, adata, fib, mac))
        log(f"    spots {info['n_after']} fib {info['n_fib_high']} mac {info['n_mac_high']}")
    pd.DataFrame(qc_rows).to_csv(RESULTS / "section_qc.csv", index=False)

    universe = gene_universe([a for *_, a, _, _ in loaded])
    pairs = load_lr_table(universe)
    pairs.to_csv(RESULTS / "lr_panel.csv", index=False)
    log(f"LR pairs with both genes in all sections: {len(pairs)}")
    ligands = pairs["ligand"].tolist()
    receptors = pairs["receptor"].tolist()
    # unique gene order for matrices
    lig_genes = list(dict.fromkeys(ligands))
    rec_genes = list(dict.fromkeys(receptors))
    lig_index = {g: i for i, g in enumerate(lig_genes)}
    rec_index = {g: i for i, g in enumerate(rec_genes)}
    lig_col = np.array([lig_index[g] for g in ligands])
    rec_col = np.array([rec_index[g] for g in receptors])

    rng = np.random.default_rng(SEED)
    colo_rows = []
    score_rows = []

    for sample, gsm, group, adata, fib, mac in loaded:
        coords = np.asarray(adata.obsm["spatial_um"], dtype=np.float64)
        adj = radius_adjacency(coords, RADIUS_UM)
        frac, n_contact, n_mac = contact_fraction(fib, mac, adj)
        null_frac = np.empty(N_PERM, dtype=np.float64)
        labels_f = fib.copy()
        labels_m = mac.copy()
        for b in range(N_PERM):
            perm = rng.permutation(adata.n_obs)
            null_frac[b], _, _ = contact_fraction(labels_f[perm], labels_m[perm], adj)
        null_cut = float(np.quantile(null_frac, 0.95))
        colo_rows.append(
            {
                "sample": sample,
                "gsm": gsm,
                "group": group,
                "contact_fraction": frac,
                "n_contact": n_contact,
                "n_mac_high": n_mac,
                "null_p95": null_cut,
                "above_null": bool(frac > null_cut),
            }
        )
        log(f"  contact {sample}: {frac:.3f} vs null95 {null_cut:.3f}")

        L_all = expr_matrix(adata, lig_genes)
        R_all = expr_matrix(adata, rec_genes)
        L = L_all[:, lig_col]
        R = R_all[:, rec_col]
        observed, n_int = pair_signal(L, R, fib, mac, adj)
        null_mat = np.empty((N_PERM, len(pairs)), dtype=np.float64)
        for b in range(N_PERM):
            shuffled = coords[rng.permutation(adata.n_obs)]
            adj_b = radius_adjacency(shuffled, RADIUS_UM)
            null_mat[b], _ = pair_signal(L, R, fib, mac, adj_b)
        null_p95 = np.nanquantile(null_mat, 0.95, axis=0)
        perm_p = np.full(len(pairs), np.nan)
        finite = np.isfinite(observed)
        if finite.any():
            perm_p[finite] = (1.0 + np.sum(null_mat[:, finite] >= observed[None, finite], axis=0)) / (N_PERM + 1)
        # detection on the unshuffled interface
        fib_f = fib.astype(bool)
        neigh_counts = np.asarray(adj @ fib.astype(np.float32)).ravel()
        interface = mac & (neigh_counts > 0)
        if interface.sum() == 0 or fib_f.sum() == 0:
            lig_det = np.zeros(len(pairs))
            rec_det = np.zeros(len(pairs))
        else:
            lig_det = (L[fib_f] > 0).mean(axis=0)
            rec_det = (R[interface] > 0).mean(axis=0)
        for i, row in pairs.iterrows():
            score_rows.append(
                {
                    "sample": sample,
                    "group": group,
                    "pair": row["pair"],
                    "ligand": row["ligand"],
                    "receptor": row["receptor"],
                    "panel": row["panel"],
                    "score": observed[i],
                    "null_p95": null_p95[i],
                    "above_null": bool(observed[i] > null_p95[i]) if np.isfinite(observed[i]) else False,
                    "perm_p": float(perm_p[i]) if np.isfinite(perm_p[i]) else np.nan,
                    "n_interface": n_int,
                    "ligand_detect_in_fib": float(lig_det[i]),
                    "receptor_detect_in_interface": float(rec_det[i]),
                    "detected": bool(lig_det[i] >= DETECT_FRAC and rec_det[i] >= DETECT_FRAC),
                }
            )
        log(f"  scored {sample} interface {n_int}")

    colo = pd.DataFrame(colo_rows)
    colo.to_csv(RESULTS / "colocalization.csv", index=False)
    scores = pd.DataFrame(score_rows)
    scores.to_csv(RESULTS / "pair_section_scores.csv", index=False)
    # null matrices are large; keep the p95 already stored. Write a compact check only.
    log(f"wrote section scores: {scores.shape}")
    summarize(colo, scores)


def summarize(colo: pd.DataFrame, scores: pd.DataFrame) -> None:
    primary_colo = colo[colo["group"] == "primary"]
    n_contact_ok = int(primary_colo["above_null"].sum())
    premise_ok = n_contact_ok >= MIN_SECTIONS_ABOVE

    prim = scores[scores["group"] == "primary"].copy()
    rows = []
    for pair, sub in prim.groupby("pair", sort=False):
        detected_n = int(sub["detected"].sum())
        above_n = int(sub["above_null"].sum())
        panel = sub["panel"].iloc[0]
        ligand = sub["ligand"].iloc[0]
        receptor = sub["receptor"].iloc[0]
        eligible = detected_n >= MIN_SECTIONS_DETECT or panel == "prespec"
        loo_ok = True
        loo_fail = []
        if eligible:
            for left in sub["sample"]:
                held = sub[sub["sample"] != left]
                if int(held["above_null"].sum()) < MIN_SECTIONS_ABOVE:
                    loo_ok = False
                    loo_fail.append(left)
        rows.append(
            {
                "pair": pair,
                "ligand": ligand,
                "receptor": receptor,
                "panel": panel,
                "n_sections_detected": detected_n,
                "n_sections_above_null": above_n,
                "median_score": float(sub["score"].median()),
                "median_perm_p": float(sub["perm_p"].median()),
                "eligible": eligible,
                "retained": bool(eligible and above_n >= MIN_SECTIONS_ABOVE and loo_ok),
                "loo_stable": bool(loo_ok),
                "loo_failed_when_dropping": ",".join(loo_fail),
            }
        )
    summary = pd.DataFrame(rows).sort_values(
        ["retained", "n_sections_above_null", "median_score"], ascending=[False, False, False]
    )
    summary.to_csv(RESULTS / "pair_summary.csv", index=False)

    meta = {
        "radius_um": RADIUS_UM,
        "n_perm": N_PERM,
        "seed": SEED,
        "high_quantile": HIGH_Q,
        "primary_sections_contact_above_null": n_contact_ok,
        "primary_sections": int(primary_colo.shape[0]),
        "spatial_contact_premise": premise_ok,
        "n_pairs_tested": int(summary.shape[0]),
        "n_retained": int(summary["retained"].sum()),
        "ccl2_ccr2": summary.loc[summary["pair"] == "CCL2–CCR2"].to_dict(orient="records"),
    }
    (RESULTS / "summary.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False))
    log(json.dumps(meta, indent=2, ensure_ascii=False))


def read_sample_means(path: Path, genes: list[str]) -> pd.Series:
    df = pd.read_csv(path, index_col=0)
    idx_as_genes = pd.Index(df.index.astype(str)).str.match(r"^[A-Za-z][A-Za-z0-9\.\-]*$").mean() > 0.5
    if not idx_as_genes:
        df = df.T
    df.index = df.index.astype(str)
    values = df.to_numpy(dtype=np.float64)
    finite = np.isfinite(values)
    integer_like = np.nanmean(np.abs(values[finite] - np.rint(values[finite]))) < 1e-6
    if integer_like and np.nanmax(values) > 30:
        lib = values.sum(axis=0)
        lib[lib == 0] = np.nan
        values = np.log1p(values / lib * 1e4)
    means = []
    for gene in genes:
        hits = np.flatnonzero(df.index == gene)
        if len(hits) == 0:
            means.append(np.nan)
        else:
            means.append(float(np.nanmean(values[hits[0]])))
    return pd.Series(means, index=genes)


def analyze_peritoneal(summary: pd.DataFrame) -> None:
    manifest = pd.read_csv(DATA / "scrna" / "manifest.tsv", sep="\t")
    focus = summary[summary["eligible"]].copy()
    if focus.empty:
        focus = summary.head(0)
    genes = sorted(set(focus["ligand"]).union(focus["receptor"]).union(FIB_GENES).union(MAC_GENES))
    # always include prespec genes
    for a, b, _ in PRESPEC:
        if a not in genes:
            genes.append(a)
        if b not in genes:
            genes.append(b)
    rows = []
    for rec in manifest.itertuples(index=False):
        path = DATA / "scrna" / f"{rec.gsm}_sample{rec.sample}.csv.gz"
        if not path.exists():
            log(f"missing {path.name}")
            continue
        log(f"  pseudobulk {path.name}")
        means = read_sample_means(path, genes)
        row = {"sample": int(rec.sample), "gsm": rec.gsm, "role": rec.role}
        row.update(means.to_dict())
        rows.append(row)
    wide = pd.DataFrame(rows)
    wide.to_csv(RESULTS / "scrna_sample_means.csv", index=False)

    pm = wide[wide["role"] == "peritoneal_tumor"]
    pt = wide[wide["role"] == "primary_tumor"]
    out = []
    pair_genes = focus[["pair", "ligand", "receptor", "panel", "retained"]].drop_duplicates()
    if pair_genes.empty:
        pair_genes = pd.DataFrame(PRESPEC, columns=["ligand", "receptor", "panel"])
        pair_genes["pair"] = pair_genes["ligand"] + "–" + pair_genes["receptor"]
        pair_genes["retained"] = False
    for rec in pair_genes.itertuples(index=False):
        for gene, role_g in [(rec.ligand, "ligand"), (rec.receptor, "receptor")]:
            if gene not in wide.columns:
                continue
            pm_vals = pm[gene].to_numpy(float)
            pt_vals = pt[gene].to_numpy(float)
            pm_med = float(np.nanmedian(pm_vals))
            pt_med = float(np.nanmedian(pt_vals))
            paired = {}
            # sample 38 peritoneal tumor vs sample 36 primary tumor, same patient in the series map
            if {36, 38}.issubset(set(wide["sample"])):
                v38 = float(wide.loc[wide["sample"] == 38, gene].iloc[0])
                v36 = float(wide.loc[wide["sample"] == 36, gene].iloc[0])
                paired = {"paired_sample38": v38, "paired_sample36": v36, "paired_log2fc": float(np.log2((v38 + 1e-3) / (v36 + 1e-3)))}
            n_pm_above_pt_median = int(np.nansum(pm_vals > pt_med))
            out.append(
                {
                    "pair": rec.pair,
                    "panel": rec.panel,
                    "retained_in_spatial": bool(rec.retained),
                    "gene": gene,
                    "gene_role": role_g,
                    "pm_n": int(np.isfinite(pm_vals).sum()),
                    "primary_n": int(np.isfinite(pt_vals).sum()),
                    "pm_median": pm_med,
                    "primary_median": pt_med,
                    "log2fc_median": float(np.log2((pm_med + 1e-3) / (pt_med + 1e-3))),
                    "n_pm_above_primary_median": n_pm_above_pt_median,
                    **paired,
                }
            )
    direc = pd.DataFrame(out)
    direc.to_csv(RESULTS / "peritoneal_direction.csv", index=False)
    log(f"peritoneal rows {len(direc)}")


def plot_results() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig_dir = RESULTS / "figures"
    fig_dir.mkdir(exist_ok=True)
    colo = pd.read_csv(RESULTS / "colocalization.csv")
    summary = pd.read_csv(RESULTS / "pair_summary.csv")
    scores = pd.read_csv(RESULTS / "pair_section_scores.csv")

    order = colo.sort_values("sample")
    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(order))
    ax.bar(x, order["contact_fraction"], color="#4C78A8", label="Observed")
    ax.plot(x, order["null_p95"], color="#E45756", marker="o", linestyle="none", label="Label-shuffle null 95th")
    ax.set_xticks(x)
    ax.set_xticklabels(order["sample"], rotation=45, ha="right")
    ax.set_ylabel("Fraction of macrophage-high spots\nwith a fibroblast-high neighbor")
    ax.set_title("150 μm neighborhood, primary Visium plus GC6-PM")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "colocalization.png", dpi=160)
    plt.close(fig)

    show = summary[summary["eligible"]].sort_values("median_score", ascending=False).head(20)
    if show.empty:
        show = summary.sort_values("median_score", ascending=False).head(15)
    prim = scores[scores["group"] == "primary"]
    med_null = prim.groupby("pair")["null_p95"].median()
    fig, ax = plt.subplots(figsize=(8, 6))
    y = np.arange(len(show))
    ax.barh(y, show["median_score"], color="#72B7B2", label="Median section score")
    ax.plot(show["pair"].map(med_null), y, color="#E45756", marker="o", linestyle="none", label="Median null 95th")
    ax.set_yticks(y)
    ax.set_yticklabels(show["pair"])
    ax.invert_yaxis()
    ax.set_xlabel("Interface signal (median across 9 primary sections)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "pair_scores.png", dpi=160)
    plt.close(fig)


def main() -> None:
    unpack_visium()
    analyze_spatial()
    summary = pd.read_csv(RESULTS / "pair_summary.csv")
    analyze_peritoneal(summary)
    plot_results()
    log("done")


if __name__ == "__main__":
    main()
