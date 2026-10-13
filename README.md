# Benchmarking scRNA-seq protocols: partial replication + dashboard
Partial replication of Mereu, Lafzi et al., *Nat Biotechnol* 38:747–755 (2020), [doi](https://doi.org/10.1038/s41587-020-0469-4).
Data: [GSE133549](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133549), a SuperSeries of 15 per-protocol SubSeries.
Authors' R code: [ati-lz/HCA_Benchmarking](https://github.com/ati-lz/HCA_Benchmarking). See `docs/figure_map.md` for scope.

**New here? Read [docs/METHODOLOGY.md](docs/METHODOLOGY.md)** (plain-language method, diagrams) or open `site/methodology.html`.

## Run it (Python ≥3.9)
```bash
bash run.sh setup && source .venv/bin/activate
bash run.sh test                       # synthetic smoke test; preview: bash run.sh serve → /?demo
bash pipeline/download_geo.sh list     # which SubSeries have files
bash run.sh download GSE133540 GSE133535 GSE133539   # start with 2-3 protocols
bash run.sh prepare GSE133540 GSE133535 GSE133539
bash run.sh inspect                    # CHECK shapes, orientation, integer counts, label overlap
# add data/metadata.csv (index = cell barcode, column `celltype`) if labels are available
bash run.sh analyse                    # writes site/data/results.json (+ offline .js fallback)
python pipeline/qc_distributions.py  # pre-QC distributions for Quality Figs. R6-R8
python pipeline/map_symbols.py       # ENSEMBL→symbol for marker genes (needs network)
python pipeline/similarity.py        # Fig. 3c,f: correlation heatmap + cumulative curves
bash run.sh serve
```
Currently processed: **all 15 human + all 15 mouse SubSeries** (Human/Mouse toggle in the header) (see `site/methodology.html` §2.2 for the matrix table).
Dashboard tabs: Overview · Quality (R1–R3, R6–R8) · Clustering (R4) · Markers (R10 dot plot) · Benchmark (R5 composite + R9 ranks) · Protocols · Explorer (t-SNE) · Glossary · Datasets · Discrepancies.
If `inspect` shows problems, fix `config.yaml` (orientation, thresholds) or the file choice, then rerun. Without labels the run still works
but ARI/NMI/silhouette are skipped; the dashboard shows label coverage per protocol.

## Publish
Push to `main`, then Settings → Pages → Source: GitHub Actions. Commit `site/data/results.json` (and add real entries to the Discrepancies tab).
## Method notes
Equalised cell numbers · mito/gene QC · log-normalised counts · top-variance genes · PCA · kNN graph · Louvain · cluster markers · stability over 80% subsamples.
