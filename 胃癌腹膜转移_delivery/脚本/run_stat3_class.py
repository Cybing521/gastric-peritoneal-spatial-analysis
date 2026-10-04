#!/usr/bin/env python3
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/root/autodl-tmp/gc_spatial_niche")
DATA = ROOT / "data" / "scrna"
OUT = ROOT / "results" / "interface_v2"
PROGRAMS = {
    "fib": ["DCN", "LUM", "COL1A1", "FAP"],
    "mac": ["CD68", "CSF1R", "C1QA", "C1QB"],
    "epi": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "t": ["CD3D", "CD3E", "TRAC"],
}
EXTRA = ["STAT3", "SOCS3"]
MIN_CLASS = 20


def log(msg: str) -> None:
    print(msg, flush=True)


def class_extra(path: Path) -> dict:
    want = set(EXTRA)
    for genes in PROGRAMS.values():
        want.update(genes)
    kept = {}
    lib = None
    with gzip.open(path, "rt") as handle:
        handle.readline()
        for line in handle:
            gene, _, rest = line.partition(",")
            gene = gene.strip().strip('"')
            arr = np.fromstring(rest, sep=",")
            if lib is None:
                lib = np.zeros(len(arr), dtype=np.float64)
            if len(arr) != len(lib):
                raise SystemExit(f"{path.name} {gene}")
            lib += np.nan_to_num(arr)
            if gene in want and gene not in kept:
                kept[gene] = arr
    if not kept:
        raise SystemExit(path.name)
    n = len(lib)
    values = {g: kept[g] for g in kept}
    integer_like = np.mean([np.mean(np.abs(v - np.rint(v))) for v in values.values()]) < 1e-6
    if integer_like and max(float(np.nanmax(v)) for v in values.values()) > 30:
        lib[lib == 0] = np.nan
        for g in values:
            values[g] = np.log1p(values[g] / lib * 1e4)
    scores = []
    labels = list(PROGRAMS)
    for label in labels:
        mats = [values[g] for g in PROGRAMS[label] if g in values]
        scores.append(np.vstack(mats).mean(axis=0) if mats else np.full(n, np.nan))
    stack = np.vstack(scores)
    winner = np.argmax(stack, axis=0)
    best = stack[winner, np.arange(n)]
    second = np.partition(stack, -2, axis=0)[-2]
    assigned = (best > 0) & (best > second) & np.isfinite(best)
    out = {}
    for i, label in enumerate(labels):
        mask = assigned & (winner == i)
        out[f"n_{label}"] = int(mask.sum())
        for gene in EXTRA:
            if gene not in values or mask.sum() < MIN_CLASS:
                out[f"{gene}__{label}"] = float("nan")
            else:
                out[f"{gene}__{label}"] = float(np.nanmean(values[gene][mask]))
    return out


def main() -> None:
    manifest = pd.read_csv(DATA / "manifest.tsv", sep="\t")
    rows = []
    for rec in manifest.itertuples(index=False):
        path = DATA / f"{rec.gsm}_sample{rec.sample}.csv.gz"
        log(path.name)
        row = {"sample": int(rec.sample), "gsm": rec.gsm, "role": rec.role}
        row.update(class_extra(path))
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "scrna_stat3_class.csv", index=False)
    log("stat3 class done")


if __name__ == "__main__":
    main()
