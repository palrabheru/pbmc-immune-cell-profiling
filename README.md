# PBMC Immune Cell Profiling | Single-Cell RNA-seq

A reproducible Python/Scanpy project exploring immune populations in a public single-cell gene-expression dataset. My goal is to understand **how clusters are formed and how marker genes support possible cell identities**, rather than assigning labels just because a UMAP looks separated.

## Research question
Can we distinguish major immune-cell populations from gene-expression measurements of individual peripheral blood cells?

## Data
- 10x Genomics **PBMC 3k**, loaded with `scanpy.datasets.pbmc3k()`.
- Healthy donor blood sample. **Not a cancer dataset**; it does not establish cancer biomarkers or patient outcomes.
- [Dataset documentation](https://scanpy.readthedocs.io/en/stable/generated/scanpy.datasets.pbmc3k.html).

## Analysis workflow
1. Filter cells with fewer than 200 detected genes and genes expressed in fewer than three cells.
2. Calculate QC metrics, including mitochondrial transcript percentage; inspect before/after plots. Remove cells with >=2,500 genes or >=5% mitochondrial RNA (example thresholds, not universally optimal).
3. Normalize counts per cell, log-transform and select 2,000 highly variable genes.
4. Scale expression, run PCA, build a 10-nearest-neighbor graph, embed cells with UMAP and cluster using Leiden at resolution 0.5.
5. Rank cluster-associated genes with Wilcoxon tests.
6. Compare candidate cell identities using panels of T, B, NK, monocyte, dendritic and platelet markers. Labels are explicitly marked **review**, not definitive.

## How to run (Windows PowerShell)
Python 3.11 recommended:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe analysis.py
```

The public dataset will download the first time. Outputs are written into `results/`:

| Output | Description |
|---|---|
| `figures/01_qc_before.png` / `02_qc_after.png` | Quality-control histograms |
| `figures/03_umap_clusters.png` | Leiden clusters |
| `figures/04_marker_dotplot.png` | Multiple immune marker genes across clusters |
| `figures/05_marker_umaps.png` | Selected gene expression on UMAP |
| `tables/ranked_cluster_markers.csv` | Cluster marker tests |
| `tables/cluster_sizes.csv` | Cluster cell counts |
| `tables/candidate_cluster_annotations.csv` | Preliminary marker-based labels, for manual review |
| `tables/marker_reference.csv` | Marker gene panel |
| `processed_pbmc3k.h5ad` | Processed AnnData object |
| `run_summary.json` | QC settings and run summary |

To make figures appear on GitHub, copy *selected* figures into the tracked `figures/` directory (see `figures/README.md`). Don't commit the large processed `.h5ad` file.

## Findings and limitations
**Current status: implementation complete, real-data execution and interpretation pending.** This README intentionally does not invent cell counts, identified cell types, p-values, or final figures. See `REPORT.md` for the report framework to fill after running.

Cluster labels are suggestions based on expression of marker genes, **not validated identities**. Scores are not probabilities and should not be interpreted as calibrated confidence. Marker overlap, QC thresholds, clustering resolution and donor-specific effects may affect conclusions.

## Future extensions
- Validate annotations with multiple markers and known PBMC references.
- Perform QC sensitivity analysis and resolution comparisons.
- Investigate differential expression between **biologically justified** groups.
- Transfer the workflow to actual cancer tumor-microenvironment datasets.

## References
- [Scanpy clustering tutorial](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering-2017.html)
- [Scanpy PBMC3k dataset](https://scanpy.readthedocs.io/en/stable/generated/scanpy.datasets.pbmc3k.html)
- 10x Genomics PBMC 3k public data.
