# Project report: PBMC immune-cell profiling

## Abstract
This project uses single-cell RNA-seq to explore cellular heterogeneity in a public, healthy-donor PBMC sample. The workflow includes quality control, normalization, dimensionality reduction, clustering and assessment of marker genes. **Results are pending execution of the analysis**; numerical claims and cluster identities should be added only after reviewing generated output.

## Introduction
Bulk RNA sequencing measures the average expression across many cells. Single-cell RNA sequencing instead provides gene-expression profiles for individual cells, making it possible to examine differences among immune populations. This project examines a healthy blood dataset as a foundation for later work with cancer-associated immune cells.

## Methods
The 10x Genomics PBMC 3k dataset was selected. See `analysis.py` for the exact computational protocol and recorded QC thresholds. Data preprocessing uses Scanpy, followed by PCA, a neighbor graph, UMAP, Leiden clustering, and Wilcoxon-based marker ranking. Candidate labels are derived from multiple published immune marker genes and require review.

## Results — complete after running
- Original cells and genes: refer to `results/run_summary.json`.
- Cells retained after QC: refer to `results/run_summary.json`.
- Number and sizes of clusters: refer to `results/tables/cluster_sizes.csv`.
- QC histograms: `results/figures/01_qc_before.png`, `02_qc_after.png`.
- UMAP: `results/figures/03_umap_clusters.png`.
- Marker evidence: `results/figures/04_marker_dotplot.png` and `results/tables/ranked_cluster_markers.csv`.
- Manual annotation decision: record supporting **and conflicting** markers for each cluster before claiming a cell identity.

## Discussion — complete after interpreting results
Explain which populations appear well supported by multiple markers, where overlap occurs, and how QC or clustering choices might influence patterns. Avoid causal claims or cancer interpretations: this sample is healthy blood and from a single donor.

## Limitations
The dataset is from a healthy donor; no tumor comparisons, disease associations, biological replicates or clinical conclusions are justified. Cluster marker testing should be treated as exploratory, not formal between-patient inference. The preliminary cell-type suggestions are not externally validated.

## Reproducibility
Use the project requirements and `analysis.py`; preserve the `results/run_summary.json` generated locally. The project's tracked files do not contain synthetic results.
