# Methodology

Replication target: Mereu, Lafzi *et al.*, *Nat Biotechnol* 38:747–755 (2020), via the open [bioRxiv preprint 630087v1](https://www.biorxiv.org/content/10.1101/630087v1), the published abstract and figure legends, and the [authors' R code](https://github.com/ati-lz/HCA_Benchmarking). **The Online Methods and supplementary parameter settings were not accessible at the time of writing; no equivalence with the authors' implementation is claimed.** The complete execution record is `site/methodology.html`.

## 1. Objective and estimands

Objective: compare single-cell/single-nucleus RNA-seq protocols on a shared reference sample with respect to (i) molecular yield, (ii) quality distributions, (iii) partition structure, and (iv) cross-condition stability of protocol ordering.

Per-protocol estimands: median detected genes cell⁻¹; median counts cell⁻¹; QC retention; Louvain partition count k; subsample stability (mean ARI); and, where external labels exist, agreement with reference labels (ARI, NMI, silhouette). Cross-protocol estimands: correlation of mean expression, cumulative gene discovery, composite rank, human–mouse rank concordance, and signature diagnostics (R21–R50).

## 2. Reference study design (paper)

A single reference suspension was prepared in one batch, cryopreserved, and distributed: human peripheral-blood mononuclear cells (60%), mouse colon (30%), HEK293T (human, 6%), NIH3T3 (mouse, 3%), and MDCK (dog, 1%). Thirteen sc/snRNA-seq protocols were applied (~3,000 cells each; 5 plate-based, 7 microfluidic, 1 single-nucleus), plus a large Chromium V2 reference (30,807 human + 19,749 mouse cells) used only to define cell types and markers — never to score protocols. The paper standardises to 20,000 reads cell⁻¹ (lowest first quartile, `zUMIs` downsampling) with fractional depths and equal cell numbers.

## 3. Replication scope and deviations

| Paper step | Paper | This replication | Fidelity |
|---|---|---|---|
| Primary processing | Joint-genome mapping, per-species counts | Starts from GEO count matrices | Not reproduced |
| Cell-type labels | Chromium reference + matchSCore2 (p ≥ 0.5) | External `metadata.csv` (barcode index, `celltype` column); species labels in current outputs | Partial |
| Standardisation | 20K reads cell⁻¹ + equal cells | Equal cells only (`cells_per_protocol` = min QC-pass n, capped at `max_cells`) | Partial |
| Fig. 2 yield | Per cell type at 20K reads | Full-depth overall medians | Partial; depth-confounded |
| Fig. 3 structure | HVG PCA, correlation, markers, cumulative genes | Correlation (R11), markers (R10), cumulative genes (R12); no HVG-PCA panel | Partial |
| Fig. 4 clustering | Cluster accuracy + average silhouette width | ARI/NMI/silhouette where labels exist + subsample stability (added) | Approximate |
| Fig. 5 integration | Separation vs mixing; projection | Not implemented | Not implemented |
| Fig. 6 overall | Harmonic mean of 6 scores | Arithmetic-mean composite over 5 available full-depth metrics | Approximate |

"Cluster accuracy" is undefined in the accessible text; ARI/NMI are reported as chance-corrected stand-ins. Subsample stability is an addition, not a paper metric. Current `results.json` carries species-composition labels (Human/Mouse/Doublet/Dog); ARI against species labels is therefore low by construction and must not be read as cell-type accuracy.

## 4. Data provenance

Source: GEO SuperSeries GSE133549 (no files of its own); count matrices reside in 15 per-protocol SubSeries plus GSE141469 (10X V3, postdates the paper, unranked). Retrieved 2026-10-13 (`data/raw/`, 15 SubSeries directories); processed matrices written the same day (`data/processed/`, `data/processed_mouse/`). Each SubSeries ships gzipped TSVs `{ACC}_{PROTO}_human_exp_mat.tsv.gz` with a truncated `{ACC}_{PROTO}_human_metada.tsv.gz` companion, and the mouse pair — except GSE141469, which ships species-separated matrices (`…_human_exp_mat` / `…_mouse_exp_mat`, correctly spelled `…_metadata`) and therefore yields single-species 10X_V3 tracks by construction (see §11). Identifiers are versioned ENSEMBL IDs (e.g. `ENSG00000000003.14`), not symbols; the header row is corner-less (barcodes only), handled by `pandas.read_csv(index_col=0)` in chunks. Counts are integers (UMI or non-UMI reads depending on chemistry — see §6). Retrieval: `pipeline/download_geo.sh`; accession→protocol mapping in `pipeline/prepare_data.py`; orientation fixed to `genes_by_cells`.

## 5. Environment and configuration

Python ≥ 3.9, virtualenv via `bash run.sh setup` (`requirements.txt`: numpy, pandas, scipy, scikit-learn, networkx, PyYAML, umap-learn optional). No scanpy/Seurat stack. All stochastic operators consume `config.yaml: seed` (default 0): subsampling (`numpy.random.default_rng`), SVD (`TruncatedSVD(random_state)`), Louvain (`seed`), t-SNE (`random_state`). Key parameters: `min_genes: 200`, `max_mito_pct: 20` (inoperative — §6), `min_cells: 3`, `equalize_cells: true`, `max_cells: 5000`, `n_hvg: 2000`, `n_pcs: 30`, `knn: 15`, `resolution: 1.0`, `stability_repeats: 5`, `tsne_cells: 1500`. Every run embeds its software environment in `results.json: environment` (Python + package versions + UMAP flag, via `run_analysis.py: _env`); existing `site/data/results*.json` were backfilled from the project `.venv` (Python 3.9.6, numpy 2.0.2, pandas 2.3.3, scipy 1.13.1, scikit-learn 1.6.1, networkx 3.2.1, PyYAML 6.0.3, umap-learn 0.5.12).

## 6. Pipeline (in order)

1. **Load** (`run_analysis.py:load`): stream text matrices in 1,000-row gene blocks to CSR cells×genes (float32); 10x folders via `scipy.io.mmread(...).T`. Peak memory ≈ one block.
2. **QC**: per cell `n_genes = #{X>0}`, `umis = ΣX`, `mito% = 100·ΣX_mito/max(umis,1)`; retain `n_genes ≥ 200 ∧ mito% ≤ 20`; drop genes in <3 cells. The mito prefix match (`MT-`/`MT.`, case-insensitive) never fires on ENSEMBL IDs (MT-ND1 = ENSG00000198888 is present but invisible); retention therefore indexes the `min_genes` rule only. Non-UMI chemistries (Smart-seq2, C1HT, ICELL8) report reads; "UMI" figures mix molecular estimands across chemistries.
3. **Equalise**: seeded subsample without replacement to `n_use = min(max_cells, min QC-pass counts)` per protocol (human 1,451; mouse 1,154 in the current run).
4. **Normalise → HVG → PCA**: CP10K `X′ = X·10⁴/ΣX`, `L = ln(1+X′)`, top-2,000 variance genes, z-clip ±10, `TruncatedSVD(≤30)`.
5. **Graph → cluster**: `kneighbors_graph(k ≤ 15)` symmetrised by element-wise maximum, `louvain_communities(resolution=1.0, seed)`.
6. **Markers**: on log1p CP10K HVGs, per-community score `(mean_in − mean_out) × fraction expressing`; top 5 retained; union stored with per-cluster means/fractions (R10). Deterministic given the equalised set.
7. **Labels (where present)**: barcode reindex; ARI/NMI/silhouette (silhouette on ≤2,000 cells) require >50% labelled cells.
8. **Stability (always)**: 5× seeded 80% subsamples, full re-embedding + re-clustering, ARI(full labels, sub labels) restricted to the subsample; mean ± SD. Precision, not accuracy.
9. **Embeddings**: t-SNE (≤1,500 cells, perplexity min(30,(n−1)/3), init PCA, seed 0, 2 dp) + PCA scatter with variance ratios + UMAP iff available + `cellstat` (cluster, genes, UMIs).
10. **Post-passes**: `qc_distributions.py` (pre-QC percentiles, concentration, ≤1,000-cell samples; R6–R8, R13; `--dtype float32` halves peak memory with identical rounded outputs), `similarity.py` (common-gene Pearson matrix; cumulative union over 1–100 cells, mean of 3 shuffles; R11–R12), `map_symbols.py` (ENSEMBL→symbol where resolvable).

## 7. Metric definitions

Median genes/UMIs: medians over QC-pass cells at full depth (paper: 20K-equalised per cell type). Retention: QC-pass/input. ARI ∈ [−1,1] (≈0 random); NMI ∈ [0,1]; silhouette ∈ [−1,1] (≈ paper ASW). Stability: mean ± SD of subsample ARI. Composite: per-metric min–max normalisation over processed protocols, arithmetic mean (the paper's harmonic mean degenerates under few-protocol min–max). Concentration: top-50/200 count shares. Heterogeneity (R22): CV = SD/mean across cells. Evenness (R23): effective k = exp(Shannon entropy) of cluster sizes; evenness = H/log k. Residuals (R21): OLS residual of log₁₀ genes on log₁₀ counts.

## 8. Statistical reporting standard

Every dashboard legend states the display, estimator, sample size (n), transform, ordering, and error representation. Axes carry units or explicit log₁₀ scales. Error bars are ±1 SD over 5 seeded repeats unless stated. Orderings are by the displayed statistic (descending) unless rank-based. Seeds, equalised n, and full-depth status are stated per view. Claims are scoped to the processed set; min–max composites rescale when protocols are added.

## 9. Generalizability: more tissues, not more species

Human + mouse suffices for protocol benchmarking — ordering signal derives from chemistry, not species count. Generalise by varying tissue, complexity, and depth within the same species (heterogeneous vs homogeneous tissue; high- vs low-RNA-content cells; full depth vs 20K downsampling), not by adding species, which confounds annotation completeness, mito-prefix rules, orthology, and dissociation bias. Add species only for species-mixing doublet or cross-species mapping questions.

Dashboard tab `Generalizability` (both JSONs; 15 joint protocols): R18 human–mouse median-genes scatter (log–log, Pearson r ≈ 0.84); R19 rank-bump (parallel = preserved); R20 paired stability (bar = cross-tissue instability); R27 median-counts concordance (log–log, Pearson r); R28 stability gap mouse − human (diverging, ARI units); R29 human vs mouse ARI (area ∝ mean k); R30 per-species composite-rank bump (5-metric mean recomputed per species).

## 10. Signatures (R21–R26 and R31–R50, client-side, `results.json` only)

- R21 depth-adjusted complexity: OLS residual of log₁₀ genes on log₁₀ counts; + = exceeds depth expectation.
- R22 heterogeneity: CV of genes/counts across cells (`cellstat`); high = variable capture or mixed states.
- R23 evenness: nominal vs effective k; deficit = dominance by few clusters.
- R24 accuracy vs precision: ARI vs stability, area ∝ k; high-precision/low-accuracy is expected when partitions are orthogonal to available labels.
- R25 species capture: 100% stacked Human/Mouse/Doublet/Dog shares from `composition`, ordered by Human share.
- R26 marker definition: union marker count vs mean top-marker strength (mean×frac), area ∝ k.
- R31 count efficiency: median genes per 1,000 counts, ordered (non-UMI reads noted).
- R32 capture skew: mean/median ratio of genes and counts across cells (`cellstat`); >1 = right-skewed.
- R33 rarefaction increment: union genes gained 50 → 100 sampled cells (mean of 3 shuffles).
- R34 resolution vs precision: k vs subsample stability (y ∈ [0.5, 1.0]).
- R35 separation vs agreement: silhouette vs ARI (protocols with both).
- R36 agreement concordance: NMI vs ARI (x ∈ [0, 1]).
- R37 stability uncertainty: SD of subsample ARI over 5 repeats, ordered.
- R38 dominance: largest-cluster share of t-SNE-subset cells, ordered.
- R39 PCA variance retained: sum of variance ratios ×100 (%), ordered.
- R40 marker economy: union markers per cluster, ordered.
- R41 marker strength ranking: mean top-marker strength, ordered.
- R42 marker-set overlap: Jaccard heatmap of union marker sets (rows/cols by composite rank).
- R43 hubness: mean off-diagonal Pearson r per protocol (from R11), ordered.
- R44 rarefaction endpoint: union genes at 100 sampled cells, ordered.
- R45 early saturation: genes@10/genes@100, ordered.
- R46 doublet rate: Doublet share from `composition`, ordered.
- R47 spike-in recovery: Dog (MDCK) share from `composition`, ordered.
- R48 silhouette ranking: mean silhouette width, ordered (protocols with ≥2 label classes).
- R49 rank instability: max − min component rank across the 5 composite metrics, ordered.
- R50 PC1 dominance: PC1 variance share (%), ordered.

Fifty numbered figures in total (R1–R50); every legend states display, estimator, n, transform, ordering, and error representation.

## 10b. Paper analogues (P1–P6, `pipeline/paper_figures.py` → `results.json: paper`)

Direct per-figure analogues computed on replication data; species-mix labels stand in for the paper's matchSCore2 cell types throughout, and depths are binomial-thinned UMIs approximating the 20K-read zUMIs standard: P1 species composition (Fig. 1); P2 read mapping by genomic location (Fig. 1); P3 depth-equalised genes detected (Fig. 2d–e); P4 dropout probability vs expression (Fig. 2f); P5 joint HVG-PCA across protocols (Fig. 3a–b); P6 kBET-style protocol mixability (Fig. 5g–h). With R1–R50 the dashboard holds 56 figures.

## 10c. Data exposure (Datasets tab)

Beyond figures, the Datasets tab serves the underlying numbers: a download centre (complete-record JSON plus client-side full-metrics, markers, composition, correlation, cumulative-genes, and QC-distribution CSVs for the active species track) and full tables (per-protocol metrics, label agreement, species composition, pre-QC percentiles, embedded software environment). The Overview CSV button exports the same full-metrics table.

## 11. Threats to validity

Construct: full-depth medians conflate depth with efficiency; UMI/read mixing; mito gate void; Louvain/30-PC/1.0 stack differs from Seurat/8-PC/0.6; variance HVGs differ from paper; species labels stand in for cell types. Structural: GSE141469 ships species-separated matrices, so 10X_V3 compositions are single-species by construction (Human 1.0 / Mouse 1.0) and its R25/R29/R30 points are not comparable on species shares — verified against `data/metadata*.csv` protocol×label crosstabs. Internal: single seed (0); 5 stability repeats; min–max composites are set-relative. External: one batch, two tissues, no 20K standardisation — rankings are conditional, not universal. Mitigations: R21 residualisation, R18–R20 cross-tissue checks, explicit per-legend scoping.

## 12. Reproduce / extend

`bash run.sh setup && source .venv/bin/activate`; `bash run.sh test` (synthetic, asserts ARI > 0.5); `download → prepare → inspect → analyse → qc_distributions → map_symbols → similarity → serve`. Mouse track via `config.mouse.yaml`. Software versions are embedded automatically in `results.json: environment`; commit `site/data/results*.json`, record GEO retrieval dates, and curate `discrepancies` entries (`docs/discrepancies.md`).
