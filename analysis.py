"""PBMC3k immune-cell profiling. Run with Python 3.11: python analysis.py"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
FIG = OUT / "figures"
TAB = OUT / "tables"
for directory in (FIG, TAB):
    directory.mkdir(parents=True, exist_ok=True)
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=110, facecolor="white")
SEED = 42

# This public dataset comes from a healthy donor, not a cancer study.
adata = sc.datasets.pbmc3k()
adata.var_names_make_unique()
raw_cells, raw_genes = adata.shape
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)

def qc_plot(filename):
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.5))
    for ax, key in zip(axes, ["n_genes_by_counts", "total_counts", "pct_counts_mt"]):
        ax.hist(adata.obs[key].to_numpy(), bins=40)
        ax.set(title=key, ylabel="Cells")
    fig.tight_layout()
    fig.savefig(FIG / filename, dpi=170)
    plt.close(fig)

qc_plot("01_qc_before.png")
adata = adata[(adata.obs["n_genes_by_counts"] < 2500) & (adata.obs["pct_counts_mt"] < 5)].copy()
qc_plot("02_qc_after.png")
if adata.n_obs < 100:
    raise ValueError("QC filtering left too few cells. Inspect your input and thresholds.")
filtered_cells, filtered_genes = adata.shape

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
adata = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata, random_state=SEED)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden", flavor="igraph",
             n_iterations=2, directed=False)

sc.pl.umap(adata, color="leiden", show=False, title="PBMC3k | Leiden clusters")
plt.savefig(FIG / "03_umap_clusters.png", bbox_inches="tight", dpi=180)
plt.close("all")

marker_panel = {
    "T cells": ["CD3D", "CD3E", "IL7R"],
    "CD8 / cytotoxic T": ["CD8A", "CD8B", "GZMK"],
    "B cells": ["MS4A1", "CD79A", "CD79B"],
    "NK cells": ["NKG7", "GNLY", "KLRD1"],
    "Classical monocytes": ["LYZ", "S100A8", "S100A9"],
    "FCGR3A+ monocytes": ["FCGR3A", "LST1", "MS4A7"],
    "Dendritic cells": ["FCER1A", "CST3", "CLEC10A"],
    "Platelets": ["PPBP", "PF4", "GP9"],
}
markers = list(dict.fromkeys(gene for genes in marker_panel.values()
                             for gene in genes if gene in adata.raw.var_names))
if markers:
    sc.pl.dotplot(adata, markers, groupby="leiden", use_raw=True, show=False)
    plt.savefig(FIG / "04_marker_dotplot.png", dpi=180, bbox_inches="tight")
    plt.close("all")
    subset = [gene for gene in ["CD3D", "MS4A1", "NKG7", "LYZ", "FCGR3A", "PPBP"]
              if gene in adata.raw.var_names]
    if subset:
        sc.pl.umap(adata, color=subset, use_raw=True, show=False, ncols=3)
        plt.savefig(FIG / "05_marker_umaps.png", dpi=180, bbox_inches="tight")
        plt.close("all")

sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon",
                        use_raw=True, pts=True)
markers_df = sc.get.rank_genes_groups_df(adata, group=None)
markers_df.to_csv(TAB / "ranked_cluster_markers.csv", index=False)
sizes = adata.obs["leiden"].value_counts().sort_index()
sizes.rename_axis("cluster").rename("n_cells").to_csv(TAB / "cluster_sizes.csv")

# Interpret clusters with multiple markers. These are *candidate* annotations;
# mixed profiles or close scores are kept as uncertain.
raw = adata.raw.to_adata()
rows = []
for cluster in sorted(adata.obs["leiden"].unique(), key=int):
    mask = (adata.obs["leiden"].astype(str).to_numpy() == str(cluster))
    scores = {}
    for cell_type, genes in marker_panel.items():
        present = [gene for gene in genes if gene in raw.var_names]
        if len(present) < 2:
            continue
        sub = raw[mask, present].X
        scores[cell_type] = float(np.asarray(sub.mean(axis=0)).ravel().mean())
    ordered = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    winner = ordered[0][0] if ordered else "Unresolved"
    # Scores are useful for review but aren't calibrated confidence probabilities.
    rows.append({"cluster": cluster, "cells": int(mask.sum()),
                 "candidate_label": winner + " (review)" if ordered else "Unresolved",
                 "top_marker_score": ordered[0][1] if ordered else np.nan,
                 "second_marker_score": ordered[1][1] if len(ordered) > 1 else np.nan})
pd.DataFrame(rows).to_csv(TAB / "candidate_cluster_annotations.csv", index=False)
pd.DataFrame({"marker_category": [k for k,v in marker_panel.items() for _ in v],
              "gene": [g for v in marker_panel.values() for g in v]}).to_csv(
    TAB / "marker_reference.csv", index=False)

adata.write_h5ad(OUT / "processed_pbmc3k.h5ad", compression="gzip")
metadata = {"dataset": "10x Genomics PBMC 3k (healthy blood)",
            "cells_loaded": raw_cells, "genes_loaded": raw_genes,
            "cells_after_qc": filtered_cells, "genes_after_qc": filtered_genes,
            "n_clusters": int(adata.obs["leiden"].nunique()),
            "qc": {"min_genes": 200, "min_cells_per_gene": 3,
                   "max_genes_per_cell_exclusive": 2500,
                   "max_mito_percent_exclusive": 5},
            "umap_seed": SEED,
            "note": "Candidate annotations require manual inspection; not a cancer cohort"}
(OUT / "run_summary.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
print(json.dumps(metadata, indent=2))
print("Figures:", FIG, "| Tables:", TAB)
