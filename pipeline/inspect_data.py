"""Diagnose data/processed/*: shape, orientation, counts vs normalised, and whether cell names match data/metadata.csv."""
import glob, os, pandas as pd, numpy as np
meta = pd.read_csv("data/metadata.csv", index_col=0) if os.path.exists("data/metadata.csv") else None
print("metadata:", "none" if meta is None else f"{len(meta)} cells, columns {list(meta.columns)}")
for f in sorted(glob.glob("data/processed/*")):
    if f.endswith(".md"): continue
    print("\n==", f)
    if os.path.isdir(f): print("  10x folder:", os.listdir(f)); continue
    sep = "\t" if (".tsv" in f or ".txt" in f) else ","
    head = pd.read_csv(f, index_col=0, sep=sep, nrows=20)
    names, n = [], 0
    for ch in pd.read_csv(f, usecols=[0], sep=sep, chunksize=100000):
        n += len(ch); names += ch.iloc[:, 0].astype(str).tolist()[:3000]
    num = head.select_dtypes("number")
    print(f"  rows={n} cols={head.shape[1]} numeric_cols={num.shape[1]}/{head.shape[1]}")
    print("  row names :", names[:4], "\n  col names :", list(head.columns[:4]))
    print("  integer counts:", bool(np.all(num.values % 1 == 0)), "(False => already normalised; QC/normalisation assumptions break)")
    print("  likely orientation:", "genes_by_cells" if n > head.shape[1] else "cells_by_genes (or very few genes)")
    if meta is not None:
        for lab, ids in (("rows", names), ("cols", list(head.columns))):
            print(f"  metadata overlap with {lab}: {np.mean([i in meta.index for i in ids]):.0%}")
