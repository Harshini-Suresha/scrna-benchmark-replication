"""Cross-protocol similarity (paper Fig. 3c,f), merged into site/data/results.json.

Per protocol (human matrices in config matrix_dir), one chunked pass:
  - per-cell detected genes -> QC mask (n_genes >= min_genes; mito gate
    skipped: ENSEMBL IDs, same artefact as the main pipeline)
  - per-gene total counts over QC-pass cells (+ grand total for CP10K)
  - dense sampled submatrix of <=100 QC-pass cells (seed 0) for curves
Then across protocols (no RNG beyond the fixed sampling):
  - common genes (intersection over all protocols)
  - mean CP10K per common gene -> log1p -> Pearson correlation matrix
  - cumulative detected genes over 1..100 sampled cells, mean of 3 shuffles
Stores d["similarity"] = {genes_common, order, correlation, cumulative}.
Usage: python pipeline/similarity.py [results.json]
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
import yaml

RES = sys.argv[1] if len(sys.argv) > 1 else "site/data/results.json"
rng = np.random.default_rng(0)
N_CUM = 100


def stream(path, min_genes):
    sep = "\t" if (".tsv" in path or ".txt" in path) else ","
    # pass 1: per-cell detected genes (cheap) -> QC mask + fixed sample
    ncells, ng = 0, None
    for ch in pd.read_csv(path, index_col=0, sep=sep, chunksize=2000):
        v = ch.values
        if ng is None:
            ncells = v.shape[1]
            ng = np.zeros(ncells)
        ng += (v > 0).sum(0)
    qc = np.where(ng >= min_genes)[0]
    take = np.sort(rng.choice(qc, min(N_CUM, len(qc)), replace=False))
    # pass 2: per-gene sums over QC cells + sampled presence submatrix
    gsums, total, gnames, sblocks = None, 0.0, [], []
    for ch in pd.read_csv(path, index_col=0, sep=sep, chunksize=2000):
        v = ch.values.astype("float64")
        if gsums is None:
            gsums = np.zeros(v.shape[0] * 0)  # grown below
            gsums = np.zeros(0)
        q = v[:, qc]
        gsums = np.concatenate([gsums, q.sum(1)])
        total += q.sum()
        gnames += ch.index.astype(str).tolist()
        sblocks.append((v[:, take] > 0).astype(np.int8))
    sample = np.vstack(sblocks)
    gnames = np.array(gnames)
    return {
        "cells": ncells,
        "n_qc": int(len(qc)),
        "gsums": gsums,
        "total": float(total),
        "genes": gnames,
        "sample": sample,
    }


def main():
    d = json.load(open(RES))
    cfg = d["config"]
    min_genes = cfg.get("min_genes", 200)
    files = {
        os.path.basename(f.rstrip("/")).split(".")[0]: f
        for f in glob.glob(cfg["matrix_dir"] + "/*")
        if not f.endswith((".md", ".js"))
    }
    order = [p["protocol"] for p in d["protocols"]]
    per = {}
    for name in order:
        print("similarity: streaming", name, flush=True)
        per[name] = stream(files[name], min_genes)
    common = set(per[order[0]]["genes"])
    for name in order[1:]:
        common &= set(per[name]["genes"])
    common = sorted(common)
    print(f"common genes: {len(common)}", flush=True)
    means = []
    for name in order:
        st = per[name]
        pos = {g: i for i, g in enumerate(st["genes"])}
        idx = np.array([pos[g] for g in common])
        m = st["gsums"][idx] / max(1, st["n_qc"])
        m = m / max(1e-9, (st["total"] / max(1, st["n_qc"]))) * 1e4
        means.append(np.log1p(m))
    corr = np.corrcoef(np.array(means))
    cum = {}
    for name in order:
        S = per[name]["sample"]
        curves = []
        for _ in range(3):
            o = rng.permutation(S.shape[1])
            curves.append((S[:, o].cumsum(axis=1) > 0).sum(axis=0))
        cum[name] = [int(round(float(v))) for v in np.mean(curves, axis=0)]
        print(f"  {name}: {cum[name][-1]} genes by 100 cells", flush=True)
    d["similarity"] = {
        "genes_common": len(common),
        "qc_note": f"means over cells with >={min_genes} genes (mito gate skipped: ENSEMBL IDs)",
        "order": order,
        "correlation": [[round(float(v), 3) for v in row] for row in corr],
        "cumulative": cum,
    }
    json.dump(d, open(RES, "w"))
    bn = RES; var = "window.DEMO" if "demo" in bn else ("window.RESULTS_MOUSE" if "mouse" in bn else "window.RESULTS")
    open(os.path.splitext(RES)[0] + ".js", "w").write(
        var + " = " + json.dumps(d, separators=(",", ":")) + ";"
    )
    print("updated", RES)


main()
