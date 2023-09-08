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
# ---------------------------------------------------------------- #
#
# Plot differentially expressed gene heatmap
#
# ---------------------------------------------------------------- #
# read in the enriched and depleted genes of each cluster
DEs_enriched = []
DEs_depleted = []
DEs = []
for idx in range(14):
    DE_enriched = pd.read_csv(f"integration/DE/DE_cluster{idx}_enriched.csv").index.values[:10]
    DE_depleted = pd.read_csv(f"integration/DE/DE_cluster{idx}_depleted.csv").index.values[:10]
    DEs_enriched.append(DE_enriched)
    DEs_depleted.append(DE_depleted)
    DEs.append(np.concatenate([DE_enriched, DE_depleted], axis = 0))

# group by seurat clusters
adata_clusters = []
for idx in range(0,14):
    adata_clusters.append(adata[adata.obs["seurat_clusters"] == str(idx),:])
adata_merge = AnnData.concatenate(*adata_clusters, join = "inner")
counts = adata_merge[:,np.concatenate(DEs_enriched, axis = 0)].X.toarray()

# z-score
counts = (counts - np.mean(counts, axis = 1, keepdims = True))/np.sqrt(np.mean((counts - np.mean(counts, axis = 1, keepdims = True))**2, axis = 1, keepdims = True) + 1e-9)


# plot heatmap
plt.rcParams["font.size"] = 10
lut = plt.cm.get_cmap("tab20b", 14)
yticks = []
yticklabels = []
count = 0
for idx in range(0,14):

    yticks.append(count + int(np.sum(adata_merge.obs["seurat_clusters"] == str(idx))/2))
    count += int(np.sum(adata.obs["seurat_clusters"] == str(idx)))
    yticklabels.append(idx)
yticks = np.array(yticks)
yticklabels = np.array(yticklabels)

fig = plt.figure(figsize = (70, 40))
ax = fig.add_subplot()
ax = sns.heatmap(counts, ax = ax)
ax.set_yticks(yticks)
ax.set_xticks(np.arange(np.concatenate(DEs_enriched, axis = 0).shape[0]))
ax.tick_params(axis='y', which='major', pad=30, length=0) # extra padding to leave room for the row colors
ax.set_yticklabels(yticklabels, rotation=45)
ax.set_xticklabels(np.concatenate(DEs_enriched, axis = 0), rotation = 45)

count = 0
for idx in range(0,14):
    ax.add_patch(plt.Rectangle(xy=(-0.007, count), width=0.007, height=int(np.sum(adata.obs["seurat_clusters"] == str(idx))), 
                         color=lut(idx), lw=0, transform=ax.get_yaxis_transform(), clip_on=False))
    count += int(np.sum(adata.obs["seurat_clusters"] == str(idx)))

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=40)
# fig.savefig("integration/DE/heatmap.pdf", bbox_inches = "tight")
fig.savefig("integration/DE/heatmap_zscore.png", bbox_inches = "tight", dpi = 100)


# In[]
# read in the enriched and depleted genes of each cluster
DEs_enriched = []
DEs_depleted = []
DEs = []
celltypes = ["Luminal 1", "Hillock epithelia (Basal)", "Club epithelia (Luminal)", "Macrophage (Myeloid, Total immune)", "Monocytes (Myeloid, Total immune)", "Mesenchymal", 
            "Spink1+ (Luminal)", "Basal", "Lymphoid (Total immune)", "SV, Spink1+", "Endothelial", "Luminal 2", "Luminal 3"]
for idx in celltypes:
    DE_enriched = pd.read_csv(f"integration/DE/DE_{idx}_enriched.csv").index.values[:10]
    DE_depleted = pd.read_csv(f"integration/DE/DE_{idx}_depleted.csv").index.values[:10]
    DEs_enriched.append(DE_enriched)
    DEs_depleted.append(DE_depleted)
    DEs.append(np.concatenate([DE_enriched, DE_depleted], axis = 0))

# group by seurat clusters
adata_clusters = []
for idx in celltypes:
    adata_clusters.append(adata[adata.obs["annotation"] == idx,:])
adata_merge = AnnData.concatenate(*adata_clusters, join = "inner")
counts = adata_merge[:,np.concatenate(DEs_enriched, axis = 0)].X.toarray()

# z-score
# counts = (counts - np.mean(counts, axis = 1, keepdims = True))/np.sqrt(np.mean((counts - np.mean(counts, axis = 1, keepdims = True))**2, axis = 1, keepdims = True) + 1e-9)


# plot heatmap
plt.rcParams["font.size"] = 10
lut = plt.cm.get_cmap("tab20b", 14)
yticks = []
yticklabels = []
count = 0
for idx in celltypes:

    yticks.append(count + int(np.sum(adata_merge.obs["annotation"] == idx)/2))
    count += int(np.sum(adata.obs["annotation"] == idx))
    yticklabels.append(idx)
yticks = np.array(yticks)
yticklabels = np.array(yticklabels)

fig = plt.figure(figsize = (70, 60))
ax = fig.add_subplot()
ax = sns.heatmap(counts, ax = ax)
ax.set_yticks(yticks)
ax.set_xticks(np.arange(np.concatenate(DEs_enriched, axis = 0).shape[0]))
ax.tick_params(axis='y', which='major', pad=30, length=0) # extra padding to leave room for the row colors
ax.set_yticklabels(yticklabels, rotation=45)
ax.set_xticklabels(np.concatenate(DEs_enriched, axis = 0), rotation = 45)

count = 0
for i, idx in enumerate(celltypes):
    ax.add_patch(plt.Rectangle(xy=(-0.007, count), width=0.007, height=int(np.sum(adata.obs["annotation"] == idx)), 
                         color=lut(i), lw=0, transform=ax.get_yaxis_transform(), clip_on=False))
    count += int(np.sum(adata.obs["annotation"] == idx))

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=40)
# fig.savefig("integration/DE/heatmap.pdf", bbox_inches = "tight")
fig.savefig("integration/DE/heatmap_annotation.png", bbox_inches = "tight", dpi = 100)


# In[]
# ---------------------------------------------------------------- #
#
# Plot cell composition pie charts
#
# ---------------------------------------------------------------- #

# immune cells count under 12wk
adata_ctrl_12wk = adata[(adata.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata.obs["age"] == "12wk"),:]
adata_knockout_12wk = adata[(adata.obs["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata.obs["age"] == "12wk"),:]
n_immune_ctrl_12wk = np.sum((adata_ctrl_12wk.obs["seurat_clusters"] == "4") | (adata_ctrl_12wk.obs["seurat_clusters"] == "5") | (adata_ctrl_12wk.obs["seurat_clusters"] == "9"))
n_other_ctrl_12wk = adata_ctrl_12wk.shape[0] - n_immune_ctrl_12wk
n_immune_knockout_12wk = np.sum((adata_knockout_12wk.obs["seurat_clusters"] == "4") | (adata_knockout_12wk.obs["seurat_clusters"] == "5") | (adata_knockout_12wk.obs["seurat_clusters"] == "9"))
n_other_knockout_12wk = adata_knockout_12wk.shape[0] - n_immune_knockout_12wk

plt.rcParams["font.size"] = 15
fig = plt.figure(figsize = (12,6))
axs = fig.subplots(nrows = 1, ncols = 2)
axs[0].pie([n_immune_ctrl_12wk, n_other_ctrl_12wk], labels = ["Immune", "Other"], autopct='%.0f%%')
axs[1].pie([n_immune_knockout_12wk, n_other_knockout_12wk], labels = ["Immune", "Other"], autopct='%.0f%%')
axs[0].set_title("PbCre(+/-),Pten(-/-),P53(-/-)")
axs[1].set_title("PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)")
fig.savefig("immune_percentage_12wk.png", bbox_inches = "tight")


# immune cells count under 18wk
adata_ctrl_18wk = adata[(adata.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata.obs["age"] == "18wk"),:]
adata_knockout_18wk = adata[(adata.obs["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata.obs["age"] == "18wk"),:]
n_immune_ctrl_18wk = np.sum((adata_ctrl_18wk.obs["seurat_clusters"] == "4") | (adata_ctrl_18wk.obs["seurat_clusters"] == "5") | (adata_ctrl_18wk.obs["seurat_clusters"] == "9"))
n_other_ctrl_18wk = adata_ctrl_18wk.shape[0] - n_immune_ctrl_18wk
n_immune_knockout_18wk = np.sum((adata_knockout_18wk.obs["seurat_clusters"] == "4") | (adata_knockout_18wk.obs["seurat_clusters"] == "5") | (adata_knockout_18wk.obs["seurat_clusters"] == "9"))
n_other_knockout_18wk = adata_knockout_18wk.shape[0] - n_immune_knockout_18wk

plt.rcParams["font.size"] = 15
fig = plt.figure(figsize = (12,6))
axs = fig.subplots(nrows = 1, ncols = 2)
axs[0].pie([n_immune_ctrl_18wk, n_other_ctrl_18wk], labels = ["Immune", "Other"], autopct='%.0f%%')
axs[1].pie([n_immune_knockout_18wk, n_other_knockout_18wk], labels = ["Immune", "Other"], autopct='%.0f%%')
axs[0].set_title("PbCre(+/-),Pten(-/-),P53(-/-)")
axs[1].set_title("PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)")
fig.savefig("immune_percentage_18wk.png", bbox_inches = "tight")


# %%
