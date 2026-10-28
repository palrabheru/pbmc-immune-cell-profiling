# PBMC Immune Cell Profiling with Single-Cell RNA-seq

This project explores immune-cell populations using single-cell RNA sequencing and Python/Scanpy.

## Research question
Can expression profiles from individual blood cells reveal distinct immune cell populations?

## Data
[10x Genomics PBMC 3k](https://scanpy.readthedocs.io/en/stable/generated/scanpy.datasets.pbmc3k.html) contains around 2,700 cells from a healthy donor. This is **not cancer data**.

## Methods
1. Filter cells and genes; inspect gene counts and mitochondrial RNA.
2. Normalize and log-transform gene expression.
3. Select variable genes, run PCA, calculate nearest neighbors, UMAP and Leiden clustering.
4. Rank marker genes and visualize CD3D, MS4A1, NKG7, LYZ, FCGR3A and PPBP.

## Run with Python 3.11
```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python analysis.py
```

Outputs are written to `results/`, including quality-control plots, cluster UMAP, marker UMAP, cluster sizes, ranked genes, and a processed h5ad file. The first run downloads the public data.

## Status and interpretation
**Analysis code uploaded; execution and biological validation still needed.** Clusters are computational groups, not confirmed immune-cell labels. Cell identities require careful review of multiple markers. QC cutoffs are starting choices rather than universal rules.

## Next steps
Run the workflow, interpret the figures, label candidate T/B/NK/myeloid populations and test sensitivity to QC settings.

Dataset credit: 10x Genomics. Background: [Scanpy clustering tutorial](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering-2017.html).
