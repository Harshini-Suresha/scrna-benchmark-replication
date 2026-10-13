"""Paper-figure panels (Mereu et al. 2020, Figs. 1/2/3/5 analogues), merged into a results JSON.

What it computes, per species results file:
  - composition (Fig. 1): per-protocol Species fractions from metadata.csv
  - mapping (Fig. 1): per-protocol exonic/intronic/intergenic/unmapped read
    fractions from the GEO *_metada.tsv.gz supplements
  - depth (Fig. 2d-e): median genes detected at thinned UMI depths
    [1k, 2.5k, 5k, 10k, 20k] per protocol, overall + Human/Mouse split.
    Thinning = binomial resampling of each cell's counts (seeded).
    Approximation of the paper's zUMIs 20K-read downsampling.
  - dropout (Fig. 2f): mean zero-fraction per expression bin per protocol
    (empirical analogue of the paper's SCDE dropout probabilities)
  - joint (Fig. 3a-b): joint HVG-PCA over sampled cells across protocols
    (common genes only): PC1/PC2 coords + protocol/species/ngenes
  - mixability (Fig. 5g-h): kBET-style chi-square neighbourhood test on the
    joint PC space: per-protocol rejection rate (low = well mixed)

I/O strategy: matrices are genes (rows) x cells (cols) with up to 80k
columns, so full parses are prohibitive. Per protocol we read the header
for barcodes, then stream only ~800 sampled cell columns in chunks
(name-based usecols + index_col=0, verified exact). Steps needing
metadata/raw supplements are skipped gracefully when those files are
absent (e.g. the synthetic demo). Usage:
    python pipeline/paper_figures.py [results.json]
Does not touch embeddings, clusters, or t-SNE output.
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

RES = sys.argv[1] if len(sys.argv) > 1 else "site/data/results.json"

ACC2PROTO = {"GSE133535": "10X_2x5K_250K", "GSE133536": "10X_8x10K_25K",
             "GSE141469": "10X_V3", "GSE133537": "C1HT_medium",
             "GSE133538": "C1HT_small", "GSE133539": "CELseq2",
             "GSE133540": "Dropseq", "GSE133541": "ICELL8",
             "GSE133542": "MARSseq", "GSE133543": "QUARTZseq2",
             "GSE133544": "gmSCRBseq", "GSE133545": "SMARTseq2",
             "GSE133546": "SingleNuclei", "GSE133547": "ddSEQ",
             "GSE133548": "inDrop"}
DEPTHS = [1000, 2500, 5000, 10000, 20000]
N_CANDIDATES = 800
N_DEPTH_CELLS = 60
N_JOINT_CELLS = 100
N_HVG = 2000
N_PCS = 30
KNN_MIX = 25


def header_barcodes(path, sep):
    """Cell barcodes from the header, handling both layouts.

    Real GEO matrices: header row holds ONLY cell barcodes (each data row
    has one extra leading gene-ID field). Synthetic CSVs: header holds a
    gene column + cell barcodes. Distinguished by field counts.
    """
    hdr = pd.read_csv(path, sep=sep, nrows=0).columns.astype(str).values
    nfields = pd.read_csv(path, sep=sep, nrows=1, header=None).shape[1]
    if nfields == len(hdr) + 1:
        return hdr  # gene-index layout: every header field is a cell
    return hdr[1:]  # first header field names the gene column


def read_columns(path, sep, barcodes):
    """Streamed genes x len(barcodes) float32 blocks (name selection = exact)."""
    dt = {b: np.float32 for b in barcodes}
    genes, blocks = [], []
    for ch in pd.read_csv(path, index_col=0, sep=sep, chunksize=2000,
                           usecols=list(barcodes), dtype=dt):
        genes += ch.index.astype(str).tolist()
        blocks.append(ch.values.astype(np.float32))
    return np.array(genes), np.vstack(blocks)


def main():
    import scipy.sparse as sp
    d = json.load(open(RES))
    cfg = d["config"]
    min_genes = cfg.get("min_genes", 200)
    seed = int(cfg.get("seed", 0))
    lrng = np.random.default_rng(seed)
    files = {os.path.basename(f.rstrip("/")).split(".")[0]: f
             for f in glob.glob(cfg["matrix_dir"] + "/*")
             if not f.endswith((".md", ".js"))}
    order = [p["protocol"] for p in d["protocols"]]
    species_tok = "mouse" if "mouse" in cfg.get("matrix_dir", "") else "human"
    paper = d.get("paper", {})

    # ---- metadata: composition (Fig. 1) ----
    meta = None
    try:
        mp = cfg.get("metadata", "")
        if mp and os.path.exists(mp):
            meta = pd.read_csv(mp, index_col=0, low_memory=False)
            meta.index = meta.index.astype(str)
    except Exception as e:
        print("metadata skipped:", e, flush=True)
    ccol = cfg.get("celltype_col", "celltype")
    if meta is not None and ccol in meta.columns:
        comp = {}
        for name in order:
            ms = meta[meta["protocol"] == name] if "protocol" in meta.columns \
                else meta
            vc = ms[ccol].astype(str).value_counts(normalize=True).round(4)
            comp[name] = vc.to_dict()
        paper["composition"] = comp
        print("composition: stored for", len(comp), "protocols", flush=True)

    # ---- mapping stats (Fig. 1) ----
    mapping = {}
    for acc, proto in ACC2PROTO.items():
        if proto not in order:
            continue
        fs = sorted(glob.glob(f"data/raw/{acc}/*{species_tok}*metada*"))
        if not fs:
            continue
        try:
            m = pd.read_csv(fs[0], sep="\t", index_col=0, low_memory=True)
            cols = {c: m[c].astype(float).values
                    for c in ("nExonReads", "nIntronReads", "nIntergenicReads",
                              "nUnmappedReads") if c in m.columns}
            if len(cols) < 2:
                continue
            tot = sum(cols.values()) + 1e-9
            lbl = {"nExonReads": "exonic", "nIntronReads": "intronic",
                   "nIntergenicReads": "intergenic",
                   "nUnmappedReads": "unmapped"}
            mapping[proto] = {}
            for k, v in cols.items():
                mu = float(np.mean(v / tot))
                # Non-finite means (empty/all-zero columns) become null so the
                # JSON record stays strictly parseable; the P2 panel skips them.
                mapping[proto][lbl[k]] = round(mu, 4) if np.isfinite(mu) else None
        except Exception as e:
            print(f"mapping {proto} skipped:", e, flush=True)
    if mapping:
        paper["mapping"] = mapping
        print("mapping: stored for", len(mapping), "protocols", flush=True)

    # ---- phase 1: narrow sample per protocol -> depth + dropout ----
    cache = {}
    depths, dropout = {}, {}
    for name in order:
        print("phase1:", name, flush=True)
        path = files[name]
        sep = "\t" if (".tsv" in path or ".txt" in path) else ","
        all_bc = header_barcodes(path, sep)
        cand = np.sort(lrng.choice(all_bc, min(N_CANDIDATES, len(all_bc)),
                                   replace=False).astype(str))
        genes, V = read_columns(path, sep, cand)  # genes x candidates
        ng = np.asarray((V > 0).sum(0)).ravel()
        qcpos = np.where(ng >= min_genes)[0]
        if not len(qcpos):
            continue
        spec = np.array(["NA"] * len(cand))
        if meta is not None and ccol in meta.columns:
            hits = meta.reindex(cand)[ccol].astype(str).values
            spec = np.where(pd.isna(hits), "NA", hits)
        cache[name] = {"genes": genes, "cand": cand, "qcpos": qcpos,
                       "spec": spec, "ng": ng}
        # depth curves on <=N_DEPTH_CELLS QC cells
        take = np.sort(lrng.choice(qcpos, min(N_DEPTH_CELLS, len(qcpos)),
                                   replace=False))
        S = V[:, take].astype("float64")
        totals = S.sum(0)
        tkspec = spec[take]
        med = {}
        for D in DEPTHS:
            ok = np.where(totals >= D)[0]
            if not len(ok):
                continue
            det = np.empty(len(ok))
            for j, ci in enumerate(ok):
                th = lrng.binomial(S[:, ci].astype(np.int64),
                                   min(1.0, D / totals[ci]))
                det[j] = (th > 0).sum()
            row = {"n": int(len(ok)), "all": round(float(np.median(det)), 1)}
            for lab in ("Human", "Mouse"):
                m = tkspec[ok] == lab
                if m.sum() >= 5:
                    row[lab] = round(float(np.median(det[m])), 1)
            med[str(D)] = row
        depths[name] = med
        # dropout: zero fraction vs mean expression over sampled cells
        mu = S.mean(1)
        zf = (S == 0).mean(1)
        keep = mu > 0
        mu, zf = np.log1p(mu[keep]), zf[keep]
        if len(mu) > 100:
            qs = np.quantile(mu, np.linspace(0, 1, 11))
            qs[-1] += 1e-9
            bi = np.clip(np.digitize(mu, qs) - 1, 0, 9)
            xs = [round(float(mu[bi == b].mean()), 3) if (bi == b).any()
                  else None for b in range(10)]
            ys = [round(float(zf[bi == b].mean()), 3) if (bi == b).any()
                  else None for b in range(10)]
            dropout[name] = {"x": xs, "y": ys}
    paper["depths"] = {"depths": DEPTHS, "per_protocol": depths,
                       "note": "median genes at thinned UMI depth (binomial, "
                               f"seed {seed}, <={N_DEPTH_CELLS} QC cells/protocol)"}
    if dropout:
        paper["dropout"] = dropout

    # ---- phase 2: joint HVG-PCA + mixability (Figs. 3a-b, 5g-h) ----
    try:
        from sklearn.decomposition import TruncatedSVD
        from sklearn.neighbors import NearestNeighbors
        from scipy.stats import chi2 as chi2dist
        have = [n for n in order if n in cache]
        common = set(cache[have[0]]["genes"].tolist())
        for n in have[1:]:
            common &= set(cache[n]["genes"].tolist())
        common = sorted(common)
        print(f"joint: {len(common)} common genes", flush=True)
        if len(common) >= 1000:
            mats, labs, ngs, protos = [], [], [], []
            for n in have:
                c = cache[n]
                path = files[n]
                sep = "\t" if (".tsv" in path or ".txt" in path) else ","
                pool = c["qcpos"]
                jt = np.sort(lrng.choice(pool, min(N_JOINT_CELLS, len(pool)),
                                         replace=False))
                jbc = c["cand"][jt]
                genes2, V2 = read_columns(path, sep, jbc)
                pos = {g: i for i, g in enumerate(genes2)}
                Jp = V2[[pos[g] for g in common]]
                mats.append(sp.csr_matrix(Jp))
                labs += [c["spec"][j] for j in jt]
                ngs += c["ng"][jt].astype(int).tolist()
                protos += [n] * len(jt)
                print(f"joint: {n} {len(jt)} cells", flush=True)
            X = sp.hstack(mats, format="csr").T.tocsr()  # cells x genes
            tot = np.asarray(X.sum(1)).ravel()
            X = X.multiply((1e4 / np.maximum(tot, 1))[:, None]).log1p().tocsr()
            mu = np.asarray(X.mean(0)).ravel()
            vr = np.asarray(X.multiply(X).mean(0)).ravel() - mu ** 2
            hv = np.argsort(vr)[::-1][:min(N_HVG, X.shape[1])]
            svd = TruncatedSVD(min(N_PCS, X.shape[1] - 1), random_state=seed)
            P = svd.fit_transform(X[:, hv])
            var = [round(float(v), 4)
                   for v in svd.explained_variance_ratio_[:2]]
            # kBET-style chi-square mixing test
            up = sorted(set(protos))
            gf = np.array([protos.count(u) for u in up]) / len(protos)
            nbrs = NearestNeighbors(
                n_neighbors=min(KNN_MIX, len(protos) - 1)).fit(P)
            rej = {u: [0, 0] for u in up}
            rej_all = 0
            parr = np.array(protos)
            for i, nb in enumerate(nbrs.kneighbors(P, return_distance=False)):
                obs = np.array([(parr[nb] == u).sum() for u in up])
                exp = len(nb) * gf
                msk = exp > 0
                stat = float(((obs[msk] - exp[msk]) ** 2 / exp[msk]).sum())
                df = int(msk.sum() - 1)
                r = bool(df > 0 and stat > chi2dist.ppf(0.95, df))
                rej[protos[i]][int(r)] += 1
                rej_all += r
            mix = {u: round(v[1] / max(1, v[0] + v[1]), 3)
                   for u, v in rej.items()}
            paper["joint"] = {
                "n_cells": int(X.shape[0]), "n_genes_common": len(common),
                "hvg": int(len(hv)), "pc_var12": var,
                "cells": [{"p": p, "s": s, "g": int(g),
                           "x": round(float(a), 3), "y": round(float(b), 3)}
                          for p, s, g, a, b in
                          zip(protos, labs, ngs, P[:, 0], P[:, 1])],
                "note": f"joint log1p-CP10K HVG PCA, <={N_JOINT_CELLS} "
                        f"QC cells/protocol, seed {seed}"}
            paper["mixability"] = {
                "per_protocol_rejection": mix,
                "overall_rejection": round(rej_all / len(protos), 3),
                "note": "kBET-style chi2 neighbourhood test "
                        f"(k={KNN_MIX}, alpha=0.05); low = well mixed"}
            print("joint: stored", X.shape[0], "cells; mixability overall",
                  round(rej_all / len(protos), 3), flush=True)
    except Exception as e:
        print("joint/mixability skipped:", repr(e), flush=True)

    d["paper"] = paper
    json.dump(d, open(RES, "w"))
    var = "window.DEMO" if "demo" in os.path.basename(RES) \
        else ("window.RESULTS_MOUSE" if "mouse" in os.path.basename(RES)
              else "window.RESULTS")
    open(os.path.splitext(RES)[0] + ".js", "w").write(
        var + " = " + json.dumps(d, separators=(",", ":")) + ";")
    print("updated", RES, "paper keys:", sorted(paper.keys()), flush=True)


if __name__ == "__main__":
    main()
