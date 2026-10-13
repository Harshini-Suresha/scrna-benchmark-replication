"""Scanpy parallel track: standard-methods reanalysis merged into a results JSON.

Per protocol (narrow column-sampled reads, memory-safe):
  - QC (min_genes) -> up to N_CELLS cells
  - AnnData -> normalize_total(1e4) -> log1p
  - highly_variable_genes (Seurat flavor, n_top_genes = n_hvg)
  - PCA (n_pcs) -> neighbors (knn) -> leiden (resolution)
  - UMAP + rank_genes_groups (t-test) top markers per Leiden cluster
  - agreement with GEO reference labels: ARI / NMI / silhouette (>=2 classes)
  - marker-set Jaccard vs the main pipeline's union markers (same file)

Stores d["scanpy"] = {per_protocol: {...}, note}. Needs scanpy
(+ leidenalg via its deps). Usage:
    python pipeline/scanpy_compare.py [results.json]
Does not touch existing blocks.
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

RES = sys.argv[1] if len(sys.argv) > 1 else "site/data/results.json"
N_CANDIDATES = 2000
N_CELLS = 1500
N_UMAP = 300


def header_barcodes(path, sep):
    hdr = pd.read_csv(path, sep=sep, nrows=0).columns.astype(str).values
    nfields = pd.read_csv(path, sep=sep, nrows=1, header=None).shape[1]
    if nfields == len(hdr) + 1:
        return hdr
    return hdr[1:]


def main():
    import scipy.sparse as sp
    import scanpy as sc
    from sklearn.metrics import (adjusted_rand_score as ari,
                                 normalized_mutual_info_score as nmi,
                                 silhouette_score)
    sc.settings.verbosity = 0
    d = json.load(open(RES))
    cfg = d["config"]
    min_genes = cfg.get("min_genes", 200)
    seed = int(cfg.get("seed", 0))
    lrng = np.random.default_rng(seed)
    files = {os.path.basename(f.rstrip("/")).split(".")[0]: f
             for f in glob.glob(cfg["matrix_dir"] + "/*")
             if not f.endswith((".md", ".js"))}
    order = [p["protocol"] for p in d["protocols"]]
    by_p = {p["protocol"]: p for p in d["protocols"]}
    meta = None
    try:
        mp = cfg.get("metadata", "")
        if mp and os.path.exists(mp):
            meta = pd.read_csv(mp, index_col=0, low_memory=False)
            meta.index = meta.index.astype(str)
    except Exception as e:
        print("metadata skipped:", e, flush=True)
    ccol = cfg.get("celltype_col", "celltype")

    out = {}
    for name in order:
        print("scanpy:", name, flush=True)
        path = files[name]
        sep = "\t" if (".tsv" in path or ".txt" in path) else ","
        all_bc = header_barcodes(path, sep)
        cand = np.sort(lrng.choice(all_bc, min(N_CANDIDATES, len(all_bc)),
                                   replace=False).astype(str))
        dt = {b: np.float32 for b in cand}
        genes, blocks = [], []
        for ch in pd.read_csv(path, index_col=0, sep=sep, chunksize=2000,
                               usecols=list(cand), dtype=dt):
            genes += ch.index.astype(str).tolist()
            blocks.append(ch.values.astype(np.float32))
        genes = np.array(genes)
        V = np.vstack(blocks)  # genes x candidates
        ng = np.asarray((V > 0).sum(0)).ravel()
        qc = np.where(ng >= min_genes)[0]
        if len(qc) < 50:
            print(f"  {name}: only {len(qc)} QC cells, skipped", flush=True)
            continue
        take = np.sort(lrng.choice(qc, min(N_CELLS, len(qc)), replace=False))
        M = sp.csr_matrix(V[:, take].T)  # cells x genes
        ad = sc.AnnData(M, var=pd.DataFrame(index=genes),
                        obs=pd.DataFrame(index=cand[take]))
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
        sc.pp.highly_variable_genes(
            ad, n_top_genes=min(int(cfg.get("n_hvg", 2000)), ad.shape[1] - 1),
            flavor="seurat", subset=True)
        sc.pp.pca(ad, n_comps=min(int(cfg.get("n_pcs", 30)),
                                   ad.shape[1] - 1),
                  random_state=seed)
        sc.pp.neighbors(ad, n_neighbors=min(int(cfg.get("knn", 15)),
                                            ad.shape[0] - 1),
                        random_state=seed)
        sc.tl.leiden(ad, resolution=float(cfg.get("resolution", 1.0)),
                     random_state=seed, key_added="leiden")
        sc.tl.umap(ad, random_state=seed)
        lab = np.array(ad.obs["leiden"].astype(int).tolist())
        row = {"n_cells": int(ad.shape[0]),
               "n_leiden": int(lab.max() + 1)}
        # labels
        cover, aval = 0.0, None
        if meta is not None and ccol in meta.columns:
            hits = meta.reindex(ad.obs.index)[ccol].astype(str).values
            ok = ~pd.isna(hits)
            cover = float(ok.mean())
            if ok.sum() > 0.5 * len(hits):
                aval = hits[ok].astype(str)
                row.update(ari=float(ari(aval, lab[ok])),
                           nmi=float(nmi(aval, lab[ok])))
                if pd.Series(aval).nunique() >= 2:
                    row["silhouette"] = float(silhouette_score(
                        ad.obsm["X_pca"][ok],
                        aval, sample_size=min(2000, ok.sum()),
                        random_state=seed))
        row["label_coverage"] = cover
        # markers + overlap with main-pipeline union set
        sc.tl.rank_genes_groups(ad, groupby="leiden", method="t-test")
        rnames = ad.uns["rank_genes_groups"]["names"]
        union = set()
        for grp in rnames.dtype.names:
            union |= set(np.atleast_1d(rnames[grp][:5]).astype(str).tolist())
        row["n_union_markers"] = int(len(union))
        main_union = set((by_p[name].get("markers") or {}).get("genes", []))
        if main_union and union:
            row["marker_jaccard"] = round(
                len(main_union & union) / len(main_union | union), 3)
        # UMAP wall subset
        u = np.sort(lrng.choice(ad.shape[0], min(N_UMAP, ad.shape[0]),
                                replace=False))
        uu = ad.obsm["X_umap"][u]
        uspec = (meta.reindex(ad.obs.index[u])[ccol].astype(str).values
                 if meta is not None and ccol in meta.columns
                 else np.array(["NA"] * len(u)))
        uspec = ["NA" if pd.isna(v) else str(v) for v in uspec]
        row["umap"] = {"x": [round(float(v), 2) for v in uu[:, 0]],
                       "y": [round(float(v), 2) for v in uu[:, 1]],
                       "cluster": [int(v) for v in lab[u]],
                       "species": uspec}
        out[name] = row
        print(" ", {k: v for k, v in row.items()
                    if k not in ("umap",)}, flush=True)
    d["scanpy"] = {"per_protocol": out,
                   "note": "scanpy 2nd track: normalize_total/log1p, "
                           "Seurat-flavor HVG, PCA, kNN, Leiden, UMAP, "
                           "t-test markers; seeds from config"}
    json.dump(d, open(RES, "w"))
    var = "window.DEMO" if "demo" in os.path.basename(RES) \
        else ("window.RESULTS_MOUSE" if "mouse" in os.path.basename(RES)
              else "window.RESULTS")
    open(os.path.splitext(RES)[0] + ".js", "w").write(
        var + " = " + json.dumps(d, separators=(",", ":")) + ";")
    print("updated", RES, "scanpy protocols:", len(out), flush=True)


if __name__ == "__main__":
    main()
