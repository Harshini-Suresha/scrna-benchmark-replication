"""Tiny synthetic dataset to smoke-test the pipeline. Output is NEVER real data."""
import os, numpy as np, pandas as pd
rng = np.random.default_rng(1); os.makedirs("data/synthetic/m", exist_ok=True)
G, types, meta = 1200, ["T", "B", "Mono", "HEK"], []
base = rng.gamma(0.5, 1, (len(types), G)); base[np.arange(4)[:, None], rng.integers(0, G, (4, 120))] *= 8
genes = [f"MT-{i}" if i < 13 else f"G{i}" for i in range(G)]
for name, depth in [("synA", 4000), ("synB", 1500), ("synC", 700)]:
    ct = rng.integers(0, 4, 400); mu = base[ct] / base[ct].sum(1, keepdims=True) * depth
    cells = [f"{name}_c{i}" for i in range(400)]
    pd.DataFrame(rng.poisson(mu).T, index=genes, columns=cells).to_csv(f"data/synthetic/m/{name}.csv.gz")
    meta += [(c, types[t], name) for c, t in zip(cells, ct)]
pd.DataFrame(meta, columns=["cell", "celltype", "protocol"]).set_index("cell").to_csv("data/synthetic/metadata.csv")
