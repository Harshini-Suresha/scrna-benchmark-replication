"""Per-cell QC distributions + library-dominance stats, merged into site/data/results.json.

Reads data/processed/* in gene-blocks (memory-safe), accumulates per-cell
n_genes / n_umis PRE-QC and per-gene totals, then stores compact percentiles
and top-gene shares per protocol. Run after pipeline/run_analysis.py:
    python pipeline/qc_distributions.py [results.json] [--dtype float32]
Accumulator dtype flag (float64 default; float32 ≈2x faster on large
matrices with negligible percentile differences). Does not touch
embeddings, clusters, or t-SNE output.
"""
import glob
import json
import os
import sys

VIOLIN_N = 1000
viol_rng = __import__('numpy').random.default_rng(0)

import numpy as np
import pandas as pd

RAW = sys.argv[1:]
DTYPE = "float64"
ARGS = []
skip_next = False
for i, a in enumerate(RAW):
    if skip_next:
        skip_next = False
        continue
    if a == "--dtype":
        if i + 1 < len(RAW) and not RAW[i + 1].startswith("--"):
            DTYPE = RAW[i + 1].strip()
            skip_next = True
        continue
    if a.startswith("--dtype="):
        DTYPE = a.split("=", 1)[1].strip()
        continue
    ARGS.append(a)
if DTYPE not in ("float32", "float64"):
    sys.exit("--dtype must be float32 or float64")
print(f"accumulator dtype: {DTYPE}")

RES = ARGS[0] if ARGS else "site/data/results.json"
QS = [0, 10, 25, 50, 75, 90, 100]


def stats(path):
    sep = "\t" if (".tsv" in path or ".txt" in path) else ","
    ncells, umis, ngenes, gtot, gcnt, gnames = 0, None, None, None, None, []
    for ch in pd.read_csv(path, index_col=0, sep=sep, chunksize=2000):
        v = ch.values.astype(DTYPE)
        if umis is None:
            ncells = v.shape[1]
            umis = np.zeros(ncells)
            ngenes = np.zeros(ncells)
            gtot = np.zeros(0)
        umis += v.sum(0)
        ngenes += (v > 0).sum(0)
        gtot = np.concatenate([gtot, v.sum(1)])
        gcnt = np.concatenate([gcnt, (v > 0).sum(1)]) if gcnt is not None else (v > 0).sum(1)
        gnames += ch.index.astype(str).tolist()
    order = np.argsort(gtot)[::-1]
    tot = gtot.sum()
    qm = np.where(ngenes >= 200)[0]
    take = np.sort(viol_rng.choice(qm, min(VIOLIN_N, len(qm)), replace=False))
    return {
        "input_cells": int(ncells),
        "genes_total_1cell": int((gtot > 0).sum()),
        "genes_total_3cells": int((gcnt >= 3).sum()),
        "genes_q": [float(x) for x in np.percentile(ngenes, QS)],
        "umis_q": [float(x) for x in np.percentile(umis, QS)],
        "top50_share": float(gtot[order[:50]].sum() / tot),
        "top200_share": float(gtot[order[:200]].sum() / tot),
        "top5_genes": [str(gnames[i]) for i in order[:5]],
        "cell_sample": {"genes": ngenes[take].astype(int).tolist(),
                        "umis": umis[take].astype(int).tolist()},
    }


def main():
    d = json.load(open(RES))
    cfg_dir = d["config"]["matrix_dir"]
    by_p = {p["protocol"]: p for p in d["protocols"]}
    files = {
        os.path.basename(f.rstrip("/")).split(".")[0]: f
        for f in glob.glob(cfg_dir + "/*")
        if not f.endswith((".md", ".js"))
    }
    for name, proto in by_p.items():
        if name not in files:
            print(f"{name}: matrix not found, skipped")
            continue
        print("QC-dist", name)
        st = stats(files[name])
        proto["qc_dist"] = st
        print(" ", {k: v for k, v in st.items() if k != "top5_genes"})
    json.dump(d, open(RES, "w"))
    bn = os.path.basename(RES); var = "window.DEMO" if "demo" in bn else ("window.RESULTS_MOUSE" if "mouse" in bn else "window.RESULTS")
    open(os.path.splitext(RES)[0] + ".js", "w").write(
        var + " = " + json.dumps(d, separators=(",", ":")) + ";"
    )
    print("updated", RES)


main()
