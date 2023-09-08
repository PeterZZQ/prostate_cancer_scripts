# In[]
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc
from anndata import AnnData
import scipy as sci

# In[]
# ---------------------------------------------------------------- #
#
# Preprocessing
#
# ---------------------------------------------------------------- #

adata = sc.read_h5ad("Cell_Ranger_output/adata_seurat.h5ad")
# sc.pp.filter_cells(adata, min_genes = 200)
sc.pp.filter_genes(adata, min_cells = 3)
# keep the original count before normalization and log transform
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

adata.obs["annotation"] = "other"
adata.obs.loc[adata.obs["seurat_clusters"] == "0", "annotation"] = "Luminal 1"
adata.obs.loc[adata.obs["seurat_clusters"] == "1", "annotation"] = "Luminal 1"
adata.obs.loc[adata.obs["seurat_clusters"] == "2", "annotation"] = "Hillock epithelia (Basal)"
adata.obs.loc[adata.obs["seurat_clusters"] == "3", "annotation"] = "Club epithelia (Luminal)"
adata.obs.loc[adata.obs["seurat_clusters"] == "4", "annotation"] = "Macrophage (Myeloid, Total immune)"
adata.obs.loc[adata.obs["seurat_clusters"] == "5", "annotation"] = "Monocytes (Myeloid, Total immune)"
adata.obs.loc[adata.obs["seurat_clusters"] == "6", "annotation"] = "Mesenchymal"
adata.obs.loc[adata.obs["seurat_clusters"] == "7", "annotation"] = "Spink1+ (Luminal)"
adata.obs.loc[adata.obs["seurat_clusters"] == "8", "annotation"] = "Basal"
adata.obs.loc[adata.obs["seurat_clusters"] == "9", "annotation"] = "Lymphoid (Total immune)"
adata.obs.loc[adata.obs["seurat_clusters"] == "10", "annotation"] = "SV, Spink1+"
adata.obs.loc[adata.obs["seurat_clusters"] == "11", "annotation"] = "Endothelial"
adata.obs.loc[adata.obs["seurat_clusters"] == "12", "annotation"] = "Luminal 2"
adata.obs.loc[adata.obs["seurat_clusters"] == "13", "annotation"] = "Luminal 3"

# In[]
adata_ti = adata[(adata.obs["annotation"] == "Basal") |
                 (adata.obs["annotation"] == "Hillock epithelia (Basal)") |
                 (adata.obs["annotation"] == "Luminal 1")|
                 (adata.obs["annotation"] == "Spink1+ (Luminal)"), :]

# adata_ti = adata.copy()
sc.pp.highly_variable_genes(adata_ti, n_top_genes = 5000)
adata_ti = adata_ti[:, adata_ti.var.highly_variable]
sc.tl.pca(adata_ti)
sc.pp.neighbors(adata_ti, n_neighbors = 30, n_pcs = 30)
sc.tl.leiden(adata_ti, resolution = 0.05)
sc.tl.umap(adata_ti, min_dist = 0.1)
sc.pl.umap(adata_ti, color = ["annotation", "sample", "leiden"])

adata_ti.write_h5ad("Cell_Ranger_output/adata_ti.h5ad")
# In[]
# # Diffusion pseudotime
# adata_ti.obs["dpt_pseudotime"] = 0
# for sample in np.unique(adata_ti.obs["sample"].values):
#     adata_sample = adata_ti[adata_ti.obs["sample"] == sample, :]

#     sc.pp.neighbors(adata_sample, n_neighbors = 30, n_pcs = 30)
#     sc.tl.diffmap(adata_sample)
#     sc.pl.diffmap(adata_sample, color = "annotation")
#     sc.pp.neighbors(adata_sample, n_neighbors=10, use_rep='X_diffmap')

#     adata_sample.uns["iroot"] = np.flatnonzero(adata_sample.obs['leiden']  == "3")[0]
#     sc.tl.dpt(adata_sample)
#     adata_ti.obs.loc[adata_ti.obs["sample"] == sample, "dpt_pseudotime"] = adata_sample.obs["dpt_pseudotime"]

# In[]
# use slingshot result
adata_ti = sc.read_h5ad("Cell_Ranger_output/adata_ti.h5ad")
adata_ti_sample = adata_ti[adata_ti.obs["sample"] == "M1417-PP12", :]
sc.pl.umap(adata_ti_sample, color = ["pseudotime", "pseudotime_joint", "annotation"], save = "_PT_M1417-PP12.png")

adata_ti_sample = adata_ti[adata_ti.obs["sample"] == "M1416-PP18", :]
sc.pl.umap(adata_ti_sample, color = ["pseudotime", "pseudotime_joint", "annotation"], save = "_PT_M1416-PP18.png")

adata_ti_sample = adata_ti[adata_ti.obs["sample"] == "M1436-PPC12", :]
sc.pl.umap(adata_ti_sample, color = ["pseudotime", "pseudotime_joint", "annotation"], save = "_PT_M1436-PPC12.png")

adata_ti_sample = adata_ti[adata_ti.obs["sample"] == "M1437-PPC18", :]
sc.pl.umap(adata_ti_sample, color = ["pseudotime", "pseudotime_joint", "annotation"], save = "_PT_M1437-PPC18.png")

# In[]
# Plot the change of marker gene expression along the pseudotime axis 
adata_ti_raw = adata[(adata.obs["annotation"] == "Basal") |
                     (adata.obs["annotation"] == "Hillock epithelia (Basal)") |
                     (adata.obs["annotation"] == "Luminal 1")|
                     (adata.obs["annotation"] == "Spink1+ (Luminal)"), :]
adata_ti_raw.obs["pseudotime"] = adata_ti.obs["pseudotime"]
adata_ti_raw.obs["pseudotime_joint"] = adata_ti.obs["pseudotime_joint"]

plt.rcParams["font.size"] = 20

cmap = plt.cm.get_cmap("tab10")
for marker in ["Trp63", "Krt5", "Krt14", "Spink1"]:
    fig = plt.figure(figsize = (30, 10))
    ax = fig.subplots(nrows = 2, ncols = 2)
    for i, sample in enumerate(["M1417-PP12", "M1416-PP18", "M1436-PPC12", "M1437-PPC18"]):
        adata_ti_sample = adata_ti_raw[adata_ti_raw.obs["sample"] == sample, :]
        marker_expr = adata_ti_sample[:, marker].X.toarray().squeeze()
        pseudo_time = adata_ti_sample.obs["pseudotime_joint"].values.squeeze()
        for j, celltype in enumerate(np.unique(adata_ti_sample.obs["annotation"].values)):
            ax[i//2, i%2].scatter(pseudo_time[adata_ti_sample.obs["annotation"] == celltype], marker_expr[adata_ti_sample.obs["annotation"] == celltype], 
                                  label = celltype, s = 6, color = cmap(j))
        ax[i//2, i%2].set_xlabel("Pseudotime (joint)")
        ax[i//2, i%2].set_ylabel("Expression")
        ax[i//2, i%2].set_title(sample)
        ax[i//2, i%2].legend(loc='upper left', prop={'size': 20}, frameon = False, ncol = 1, bbox_to_anchor=(1.04, 1), markerscale = 3)
    plt.suptitle("Marker gene: " + marker, fontsize = 30)
    plt.tight_layout()
    fig.savefig(f"figures/expr_change_{marker}_joint.png", bbox_inches = "tight")


# In[]
import seaborn as sns

# Check Basal cell density
pseudo_time = adata_ti_raw.obs[["pseudotime", "sample"]]
pseudo_time["pseudotime"] = pseudo_time["pseudotime"]/np.max(pseudo_time["pseudotime"])
    
basal_cell_idx = (adata_ti_raw.obs["annotation"] == "Basal") | (adata_ti_raw.obs["annotation"] == "Hillock epithelia (Basal)")
pseudo_time_basal = pseudo_time[basal_cell_idx]

g = sns.displot(data = pseudo_time_basal[(pseudo_time_basal["sample"] == "M1416-PP18") | (pseudo_time_basal["sample"] == "M1437-PPC18")], x = "pseudotime", hue = "sample", kind = "kde", common_norm = False)
g = sns.displot(data = pseudo_time_basal[(pseudo_time_basal["sample"] == "M1417-PP12") | (pseudo_time_basal["sample"] == "M1436-PPC12")], x = "pseudotime", hue = "sample", kind = "kde", common_norm = False)

# In[]
pseudo_time = adata_ti_raw.obs[["pseudotime_joint", "sample"]]

basal_cell_idx = (adata_ti_raw.obs["annotation"] == "Basal") | (adata_ti_raw.obs["annotation"] == "Hillock epithelia (Basal)")
pseudo_time_basal = pseudo_time[basal_cell_idx]
pseudo_time_basal.ca

g = sns.displot(data = pseudo_time_basal[(pseudo_time_basal["sample"] == "M1416-PP18") | (pseudo_time_basal["sample"] == "M1437-PPC18")], x = "pseudotime_joint", hue = "sample", kind = "kde", common_norm = False)
g = sns.displot(data = pseudo_time_basal[(pseudo_time_basal["sample"] == "M1417-PP12") | (pseudo_time_basal["sample"] == "M1436-PPC12")], x = "pseudotime_joint", hue = "sample", kind = "kde", common_norm = False)


# %%
