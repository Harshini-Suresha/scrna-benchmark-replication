# Paper figures vs this repo (from the bioRxiv v1 legends; verify against the published figures)
| Fig | Content in the paper | Status |
|---|---|---|
| 1 | Reference sample design, species composition, mapping stats | Covered (P1 species composition from GEO labels; P2 exon/intron/intergenic/unmapped mapping stats; Paper tab) |
| 2 | UMAPs of references; genes detected (20K); genes vs depth; dropout vs expression | Covered at species level (P3 thinned-depth genes curves overall + Human/Mouse; P4 empirical dropout-vs-expression; Paper tab). Reference UMAPs and per-cell-type (HEK293T/monocyte/B) splits not reproduced — no matchSCore2/Table S2 labels |
| 3 | HVG PCA; Pearson correlation; marker heatmaps; cumulative genes | Covered (P5 joint HVG-PCA by protocol/species/genes; R11 correlation; R10 markers; R12 cumulative; Paper tab) |
| 4 | t-SNE clusters per protocol; cluster accuracy and ASW | Approximate (ARI/NMI/silhouette vs species-mix labels, t-SNE in Explorer). Paper accuracy is vs matchSCore2 cell types |
| 5 | Integration UMAPs; signature scores; separation vs mixing | Partial (P6 pre- + P8 post-integration kBET-style mixability around Harmony batch correction P7/P9/P11/P13 + pre/post pair P10/P12; Paper tab). No Seurat/scMerge integration, no projection probabilities |
| 6 | Overall benchmarking score | Approximate (R5 arithmetic-mean composite + R9 rank agreement; paper used harmonic mean of 6 metrics) |

# Replication-only figures (no paper counterpart; all full-depth unless stated)
| Fig | Content | Tab |
|---|---|---|
| R1a/b | Median genes / UMIs per cell, log-scale bars | Quality |
| R2 | Genes-vs-UMIs bubble scatter (log–log, area ∝ QC-pass cells) | Quality |
| R3 | QC retention: input → pass bars | Quality |
| R4a | Clustering stability, mean ± SD over subsamples | Clustering |
| R4b | Louvain community counts per protocol | Clustering |
| R5 | Benchmark composite dot plot (5 metrics + mean) | Benchmark |
| R6/R7 | Per-cell gene / count distributions (min–max whiskers, IQR box, median, P10/P90) | Quality |
| R8 | Count concentration in top-50 / top-200 genes | Quality |
| R9 | Paper-vs-replication rank bump chart | Benchmark |
| R10 | Cluster marker dot plot (union top markers × clusters) | Markers |
| R11 | Cross-protocol Pearson correlation heatmap (common genes) | Similarity |
| R12 | Cumulative detected genes over 1–100 cells | Similarity |
| R13a/b | Per-cell gene/count violins (KDE, ≤1,000 QC-pass cells) | Quality |
| R14/R15 | Genes/counts per Louvain cluster (Tukey boxes) | Clustering |
| R16a–c | Same-cell t-SNE/PCA/UMAP gallery | Clustering |
| R16d | PCA scree (per-PC + cumulative variance) | Clustering |
| R17 | Communities vs median genes, size ∝ stability | Clustering |
| R18 | Median genes human vs mouse (log–log, Pearson r) | Generalizability |
| R19 | Human vs mouse rank bump (median genes) | Generalizability |
| R20 | Stability human (disc) vs mouse (square) | Generalizability |
| R21 | Depth-adjusted complexity (log₁₀ residual vs depth trend) | Signatures |
| R22 | Cell heterogeneity (CV of genes/counts, cellstat) | Signatures |
| R23 | Cluster evenness (k vs effective k = exp entropy) | Signatures |
| R24 | Accuracy vs precision (ARI vs stability, area ∝ k) | Signatures |
| R25 | Species capture (stacked Human/Mouse/Doublet/Dog) | Signatures |
| R26 | Marker definition (union count vs strength, area ∝ k) | Signatures |
| R27 | Median counts human vs mouse (log–log, Pearson r) | Generalizability |
| R28 | Stability gap, mouse − human (diverging, ARI) | Generalizability |
| R29 | ARI human vs mouse (area ∝ mean k) | Generalizability |
| R30 | Per-species composite-rank bump | Generalizability |
| R31 | Genes per 1,000 counts (efficiency) | Signatures |
| R32 | Capture skew (mean/median, genes + counts) | Signatures |
| R33 | Rarefaction increment, 50→100 cells | Signatures |
| R34 | k vs stability | Signatures |
| R35 | Silhouette vs ARI | Signatures |
| R36 | NMI vs ARI concordance | Signatures |
| R37 | Stability SD over 5 repeats | Signatures |
| R38 | Largest-cluster share (dominance) | Signatures |
| R39 | Total PCA variance retained | Signatures |
| R40 | Union markers per cluster (economy) | Signatures |
| R41 | Mean top-marker strength ranking | Signatures |
| R42 | Marker-set Jaccard heatmap (by composite rank) | Signatures |
| R43 | Mean off-diagonal correlation (hubness) | Signatures |
| R44 | Rarefaction endpoint (genes @100 cells) | Signatures |
| R45 | Early-saturation ratio (genes@10/genes@100) | Signatures |
| R46 | Doublet share (composition) | Signatures |
| R47 | Dog spike-in share (composition) | Signatures |
| R48 | Silhouette ranking | Signatures |
| R49 | Rank instability (max−min component rank) | Signatures |
| R50 | PC1 variance share | Signatures |
| R51 | t-SNE wall: all protocols, cluster colours | Explorer |
| R52 | PCA wall: all protocols + variance shares | Explorer |
| R53a | Louvain k vs scanpy Leiden k (paired bars) | Clustering |
| R53b | Custom vs scanpy ARI vs labels (paired bars) | Clustering |
| R54 | Scanpy Leiden UMAP wall: all protocols | Explorer |
| R55 | Custom vs scanpy NMI parity (+ diagonal) | Clustering |
| R56 | Custom vs scanpy silhouette parity (+ diagonal) | Clustering |
| R57 | Marker-set Jaccard custom vs scanpy, ranked | Clustering |
| R58 | Union marker counts custom vs scanpy (paired bars) | Clustering |
| R59 | Scanpy UMAP wall: species colours | Explorer |

# Paper-figure analogues (Paper tab, `pipeline/paper_figures.py`)
| Fig | Paper panel | Analogue |
|---|---|---|
| P1 | Fig. 1 species mix | Species shares from GEO labels (stratified by species, not cell type) |
| P2 | Fig. 1 mapping stats | Exonic/intronic/intergenic/unmapped shares from QC supplements |
| P3 | Fig. 2d–e genes at 20K reads | Median genes after binomial thinning to equal UMI depths |
| P4 | Fig. 2f dropout vs expression | Empirical dropout fraction vs mean log1p CP10K (deciles) |
| P5 | Fig. 3a–b HVG PCA | Joint HVG-PCA over common genes, coloured by protocol/species/genes |
| P6 | Fig. 5g–h integration/mixability | kBET-style χ² neighbourhood rejection rate on joint PCs |
| P7 | Fig. 5a–d integrated embedding | Joint embedding after Harmony batch correction (protocol = batch, harmonypy) |
| P8 | Fig. 5g–h post-integration mixing | Same kBET-style test on Harmony-corrected PCs (pre vs post comparison) |
| P9 | Harmony wall per protocol | Integrated coordinates split one panel per protocol, shared axes, species colours |
| P10 | Pre vs post integration | Joint PCA vs Harmony embedding side by side, protocol colours |
| P11 | Harmony wall: protocol colours | Integrated panels coloured by protocol |
| P12 | Rejection pre → post per protocol | Paired bars sorted by improvement |
| P13 | Harmony wall: genes per cell | Integrated panels shaded by library complexity |
