"""Per-protocol QC -> normalise -> HVG -> PCA -> kNN graph -> Louvain; benchmark metrics -> JSON.
Dependencies: numpy, pandas, scipy, scikit-learn, networkx, pyyaml (+ umap-learn if available, else UMAP is skipped)."""
import glob, json, os, sys, datetime
import numpy as np, pandas as pd, scipy.sparse as sp, yaml, networkx as nx
from scipy.io import mmread
from sklearn.decomposition import TruncatedSVD
from sklearn.manifold import TSNE
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics import adjusted_rand_score as ari, normalized_mutual_info_score as nmi, silhouette_score
try:
    from umap import UMAP
    HAS_UMAP = True
except Exception as e:
    UMAP, HAS_UMAP = None, False
    print("umap-learn unavailable, UMAP will be skipped:", e)

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config.yaml"))
rng = np.random.default_rng(cfg["seed"])

def pick(d, *names):
    for n in names:
        g = glob.glob(os.path.join(d, n + "*"))
        if g: return g[0]
    raise FileNotFoundError(f"{names} not in {d}")

def load(path):
    """Returns cells x genes CSR matrix, gene names, cell names. Text files are read in chunks (memory-safe)."""
    if os.path.isdir(path):  # 10x-style folder
        X = sp.csr_matrix(mmread(pick(path, "matrix.mtx")).T)
        f = pd.read_csv(pick(path, "features.tsv", "genes.tsv"), sep="\t", header=None)
        genes = f[1 if f.shape[1] > 1 else 0].astype(str).values
        cells = pd.read_csv(pick(path, "barcodes.tsv"), header=None)[0].astype(str).values
        return X.astype("float32"), genes, cells
    sep = "\t" if (".tsv" in path or ".txt" in path) else ","
    blocks, rows = [], []
    for ch in pd.read_csv(path, index_col=0, sep=sep, chunksize=1000):
        cols = ch.columns.astype(str).values; rows += ch.index.astype(str).tolist()
        blocks.append(sp.csr_matrix(ch.values.astype("float32")))
    M, rows, o = sp.vstack(blocks).tocsr(), np.array(rows), cfg["orientation"]
    genes_rows = o == "genes_by_cells" or (o == "auto" and M.shape[0] > M.shape[1])
    print(f"  {os.path.basename(path)}: {M.shape} read as {'genes x cells' if genes_rows else 'cells x genes'}")
    return (M.T.tocsr(), rows, cols) if genes_rows else (M, cols, rows)

def qc(X, genes, cells):
    ng = np.asarray((X > 0).sum(1)).ravel(); um = np.asarray(X.sum(1)).ravel()
    mt = np.array([g.upper().startswith(tuple(cfg["mito_prefixes"])) for g in genes])
    mito = np.asarray(X[:, mt].sum(1)).ravel() / np.maximum(um, 1) * 100 if mt.any() else np.zeros(len(um))
    keep = (ng >= cfg["min_genes"]) & (mito <= cfg["max_mito_pct"])
    X, cells, ng, um, mito = X[keep], cells[keep], ng[keep], um[keep], mito[keep]
    gk = np.asarray((X > 0).sum(0)).ravel() >= cfg["min_cells"]
    return X[:, gk], genes[gk], cells, ng, um, mito

def embed(X):
    L = np.log1p(X.multiply(1e4 / np.maximum(np.asarray(X.sum(1)), 1)).toarray())
    Z = L[:, np.argsort(L.var(0))[::-1][:cfg["n_hvg"]]]
    Z = np.clip((Z - Z.mean(0)) / (Z.std(0) + 1e-6), -10, 10)
    svd = TruncatedSVD(min(cfg["n_pcs"], Z.shape[1] - 1), random_state=cfg["seed"])
    return svd.fit_transform(Z), svd.explained_variance_ratio_

def cluster(P):
    A = kneighbors_graph(P, min(cfg["knn"], len(P) - 1), include_self=False)
    G = nx.from_scipy_sparse_array(A.maximum(A.T))
    lab = np.zeros(len(P), int)
    for i, c in enumerate(nx.community.louvain_communities(G, resolution=cfg["resolution"], seed=cfg["seed"])):
        lab[list(c)] = i
    return lab

meta = pd.read_csv(cfg["metadata"], index_col=0) if os.path.exists(cfg["metadata"]) else pd.DataFrame()
meta.index = meta.index.astype(str)
if meta.empty: print("No metadata file: running without reference labels (ARI/NMI/silhouette skipped).")
files = {os.path.basename(f.rstrip("/")).split(".")[0]: f for f in glob.glob(cfg["matrix_dir"] + "/*") if not f.endswith(".md")}
names = cfg["protocols"] or sorted(files)

data, errors = {}, []
for p in names:
    try:
        print("Loading", p); X, genes, cells, ng, um, mito = qc(*load(files[p]))
        if X.shape[0] < 50 or X.shape[1] < 100:
            raise ValueError(f"only {X.shape[0]} cells / {X.shape[1]} genes pass QC; check orientation, thresholds, file choice")
        data[p] = (X, genes, cells, ng, um, mito)
    except Exception as e:
        print("  SKIPPED:", repr(e)); errors.append({"protocol": p, "error": repr(e)[:200]})
if not data: sys.exit("No usable protocols. Run: python pipeline/inspect_data.py")
n_use = min([cfg["max_cells"]] + ([len(v[2]) for v in data.values()] if cfg["equalize_cells"] else []))

old_out = {}
try:
    old_out = json.load(open(cfg["output"]))
except Exception:
    pass
def _env():
    import platform
    e = {"python": platform.python_version(), "packages": {}, "umap_available": HAS_UMAP}
    for mod, dist in [("numpy", "numpy"), ("pandas", "pandas"), ("scipy", "scipy"), ("sklearn", "scikit-learn"),
                      ("networkx", "networkx"), ("yaml", "PyYAML"), ("umap", "umap-learn")]:
        try:
            m = __import__(mod)
            v = getattr(m, "__version__", None) or __import__("importlib.metadata").metadata.version(dist)
            e["packages"][dist] = str(v)
        except Exception as ex:
            e["packages"][dist] = "unknown (%s)" % ex
    return e
try:
    _environment = _env()
except Exception as ex:
    _environment = {"error": repr(ex)}
    if isinstance(old_out.get("environment"), dict):
        _environment = old_out["environment"]
out = {"generated": datetime.datetime.now(datetime.timezone.utc).isoformat(), "config": cfg, "environment": _environment, "cells_per_protocol": int(n_use), "protocols": [], "discrepancies": old_out.get("discrepancies", []), "similarity": old_out.get("similarity", {}), "paper": old_out.get("paper", {}), "errors": errors}
for p, (X, genes, cells, ng, um, mito) in data.items():
    row = {"protocol": p, "cells_pass_qc": int(len(cells)), "median_genes": float(np.median(ng)),
           "median_umis": float(np.median(um)), "median_mito_pct": float(np.median(mito))}
    idx = np.sort(rng.choice(len(cells), min(n_use, len(cells)), replace=False))
    Xs, cs = X[idx], cells[idx]
    P, pvar = embed(Xs); full = cluster(P); row["n_clusters"] = int(full.max() + 1)
    Ln = np.log1p((Xs.multiply(1e4 / np.maximum(np.asarray(Xs.sum(1)), 1))).toarray().astype("float32"))
    hv = np.argsort(Ln.var(0))[::-1][:cfg["n_hvg"]]
    G, M = genes[hv], Ln[:, hv]
    top, union = {}, []
    for c in sorted(set(full.tolist())):
        m = full == c
        mi, mo = M[m].mean(0), M[~m].mean(0)
        fi = (M[m] > 0).mean(0)
        sc = (mi - mo) * fi
        tk = np.argsort(sc)[::-1][:5]
        top[c] = tk
        union += [t for t in tk if t not in union]
    union = sorted(union)
    pos = {t: i for i, t in enumerate(union)}
    row["markers"] = {"genes": [str(G[t]) for t in union],
                      "top": {str(c): [pos[t] for t in tk] for c, tk in top.items()},
                      "clusters": {str(c): {"means": [round(float(M[full == c][:, union].mean(0)[i]), 3) for i in range(len(union))],
                                            "fracs": [round(float((M[full == c][:, union] > 0).mean(0)[i]), 3) for i in range(len(union))]}
                                   for c in sorted(set(full.tolist()))}}
    print(f"  {p}: markers for {len(top)} clusters ({len(union)} union genes)")
    lab = meta.reindex(cs)[cfg["celltype_col"]].values if cfg["celltype_col"] in meta else np.full(len(cs), np.nan)
    ok = pd.notna(lab)
    row["label_coverage"] = float(ok.mean()); print(f"  {p}: {ok.mean():.0%} of cells have reference labels")
    if ok.sum() > 0.5 * len(cs):
        l = lab[ok].astype(str)
        row.update(ari=float(ari(l, full[ok])), nmi=float(nmi(l, full[ok])))
        if pd.Series(l).nunique() >= 2:
            row["silhouette"] = float(silhouette_score(P[ok], l, sample_size=min(2000, ok.sum()), random_state=cfg["seed"]))
        row["composition"] = pd.Series(l).value_counts(normalize=True).round(4).to_dict()
    st = []
    for _ in range(cfg["stability_repeats"]):
        sub = np.sort(rng.choice(len(cs), int(0.8 * len(cs)), replace=False))
        st.append(ari(full[sub], cluster(embed(Xs[sub])[0])))
    row["stability_mean"], row["stability_sd"] = float(np.mean(st)), float(np.std(st))
    t = np.sort(rng.choice(len(cs), min(cfg["tsne_cells"], len(cs)), replace=False))
    xy = TSNE(2, perplexity=min(30, (len(t) - 1) / 3), init="pca", random_state=cfg["seed"]).fit_transform(P[t])
    row["tsne"] = {"x": xy[:, 0].round(2).tolist(), "y": xy[:, 1].round(2).tolist(),
                   "cluster": full[t].tolist(), "label": [str(v) if pd.notna(v) else "NA" for v in lab[t]]}
    row["pca"] = {"x": P[t, 0].round(2).tolist(), "y": P[t, 1].round(2).tolist(),
                  "cluster": full[t].tolist(),
                  "var": [round(float(v), 4) for v in pvar]}
    row["cellstat"] = {"cluster": full[t].tolist(), "genes": ng[idx[t]].astype(int).tolist(),
                       "umis": um[idx[t]].astype(int).tolist()}
    try:
        uu = UMAP(n_components=2, random_state=cfg["seed"]).fit_transform(P) if HAS_UMAP else None
    except Exception as e:
        print(f"  {p}: UMAP skipped ({e})"); uu = None
    row["umap"] = {"x": uu[t, 0].round(2).tolist(), "y": uu[t, 1].round(2).tolist(),
                   "cluster": full[t].tolist()} if uu is not None else None
    out["protocols"].append(row); print({k: v for k, v in row.items() if k not in ("tsne", "composition", "pca", "cellstat", "umap", "markers")})

os.makedirs(os.path.dirname(cfg["output"]), exist_ok=True)
json.dump(out, open(cfg["output"], "w"))
# Offline fallback for file:// viewing: script-tag loadable copy next to the JSON.
try:
    js_path = os.path.splitext(cfg["output"])[0] + ".js"
    bn = os.path.basename(js_path)
    var = "window.DEMO" if "demo" in bn else ("window.RESULTS_MOUSE" if "mouse" in bn else "window.RESULTS")
    open(js_path, "w").write(var + " = " + json.dumps(out, separators=(",", ":")) + ";")
    print("wrote", js_path)
except Exception as e:
    print("  (js fallback skipped:", e, ")")
