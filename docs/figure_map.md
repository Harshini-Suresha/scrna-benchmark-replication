# Paper figures vs this repo (from the bioRxiv v1 legends; verify against the published figures)
| Fig | Content in the paper | Status |
|---|---|---|
| 1 | Reference sample design, species composition, mapping stats | Covered (P1 species composition from GEO labels; P2 exon/intron/intergenic/unmapped mapping stats; Paper tab) |
| 2 | UMAPs of references; genes detected (20K); genes vs depth; dropout vs expression | Covered at species level (P3 thinned-depth genes curves overall + Human/Mouse; P4 empirical dropout-vs-expression; Paper tab). Reference UMAPs and per-cell-type (HEK293T/monocyte/B) splits not reproduced — no matchSCore2/Table S2 labels |
| 3 | HVG PCA; Pearson correlation; marker heatmaps; cumulative genes | Covered (P5 joint HVG-PCA by protocol/species/genes; R11 correlation; R10 markers; R12 cumulative; Paper tab) |
| 4 | t-SNE clusters per protocol; cluster accuracy and ASW | Approximate (ARI/NMI/silhouette vs species-mix labels, t-SNE in Explorer). Paper accuracy is vs matchSCore2 cell types |
| 5 | Integration UMAPs; signature scores; separation vs mixing | Partial (P6 pre-integration kBET-style mixability baseline; Paper tab). No Seurat/Harmony/scMerge integration, no projection probabilities |
| 6 | Overall benchmarking score | Approximate (R5 arithmetic-mean composite + R9 rank agreement; paper used harmonic mean of 6 metrics) |

# Replication-only figures (no paper counterpart; all full-depth unless stated)
| Fig | Content | Tab |
|---|---|---|
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

# Paper-figure analogues (Paper tab, `pipeline/paper_figures.py`)
| Fig | Paper panel | Analogue |
|---|---|---|
| P1 | Fig. 1 species mix | Species shares from GEO labels (stratified by species, not cell type) |
| P2 | Fig. 1 mapping stats | Exonic/intronic/intergenic/unmapped shares from QC supplements |
| P3 | Fig. 2d–e genes at 20K reads | Median genes after binomial thinning to equal UMI depths |
| P4 | Fig. 2f dropout vs expression | Empirical dropout fraction vs mean log1p CP10K (deciles) |
| P5 | Fig. 3a–b HVG PCA | Joint HVG-PCA over common genes, coloured by protocol/species/genes |
| P6 | Fig. 5g–h integration/mixability | kBET-style χ² neighbourhood rejection rate on joint PCs |
