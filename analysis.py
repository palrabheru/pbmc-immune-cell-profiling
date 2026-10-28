"""PBMC 3k immune-cell exploration in Scanpy. Run: python analysis.py"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

out = Path(__file__).parent / "results"
out.mkdir(exist_ok=True)
adata = sc.datasets.pbmc3k()
adata.var_names_make_unique()
print("Original shape:", adata.shape)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)

fig, axes = plt.subplots(1, 3, figsize=(12, 3))
for ax, key in zip(axes, ["n_genes_by_counts", "total_counts", "pct_counts_mt"]):
    ax.hist(adata.obs[key], bins=40)
    ax.set_title(key)
fig.tight_layout()
fig.savefig(out / "qc.png", dpi=150)
plt.close(fig)

adata = adata[(adata.obs.n_genes_by_counts < 2500) & (adata.obs.pct_counts_mt < 5)].copy()
print("After QC:", adata.shape)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata, random_state=42)
sc.tl.leiden(adata, resolution=0.5, flavor="igraph", directed=False, n_iterations=2, key_added="leiden")

sc.pl.umap(adata, color="leiden", show=False)
plt.savefig(out / "umap_clusters.png", dpi=180, bbox_inches="tight")
plt.close("all")

markers = [g for g in ["CD3D", "MS4A1", "NKG7", "LYZ", "FCGR3A", "PPBP"] if g in adata.raw.var_names]
sc.pl.umap(adata, color=markers, use_raw=True, ncols=3, show=False)
plt.savefig(out / "umap_markers.png", dpi=180, bbox_inches="tight")
plt.close("all")

sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)
sc.get.rank_genes_groups_df(adata, group=None).to_csv(out / "cluster_markers.csv", index=False)
adata.obs["leiden"].value_counts().sort_index().rename_axis("cluster").rename("cells").to_csv(out / "cluster_sizes.csv")
adata.write_h5ad(out / "processed_pbmc3k.h5ad", compression="gzip")
print("Finished. Results saved in:", out)
