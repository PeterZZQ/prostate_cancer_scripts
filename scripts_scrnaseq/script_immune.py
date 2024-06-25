# In[]
import sys
sys.path.append("..")
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc
import anndata
import scipy as sci
import seaborn as sns

import utils
import warnings
warnings.filterwarnings("ignore")


# Volcano plots
import matplotlib.pyplot as plt
import seaborn as sns


# In[]
# ---------------------------------------------------------------- #
#
# 1. Data loading
#
# ---------------------------------------------------------------- #
raw_dir = "../dataset/data_scrnaseq/data/qc_data/"
seurat_dir = "../dataset/data_scrnaseq/seurat_integration/"
vis_intact_dir = "../results_scrnaseq/annot_immune_subtype/"

# raw data
adata_intact_raw = sc.read_h5ad(raw_dir + "adata_qc_intact.h5ad")
# log-normalized data
adata_intact_seurat = sc.read_h5ad(seurat_dir + "adata_intact_seurat.h5ad")
adata_intact_seurat.layers["raw"] = adata_intact_raw.X.copy()

# # load the interest genes
# interest_genes = pd.read_csv("../dataset/markers_info/interest_genes.csv", index_col = 0).index.values
# interest_genes = [x for x in interest_genes if (x != "Cxcl11") & (x != "Tnfa") & (x != "Tgfb") & (x != "Il12") & (x != "Pd1") & (x != "Tim3") & (x != "Ccl21")]
# adata_intact_seurat.var["interest_gene"] = False
# adata_intact_seurat.var.loc[interest_genes, "interest_gene"] = True

# load markers of macrophage and T cell subtypes 
markers = {}
# Ptprc is immune cell marker (can be removed)
markers["immune"] = ["Ptprc"]
# Cd3e is the lymphoid marker
markers["lymphoid"] = ["Cd3e"]#, "Ms4a1", "Klrb1c"]
# lymphoid cell include T cell, B cell, and NK cell
# T cells can be further separated, now just cd4 and cd8
markers["CD4_T"] = ["Cd4"]
markers["CD8_T"] = ["Cd8a"]
markers["cytotoxicity T"] = ["Gzmb", "Prf1"]
markers["exhausted T"] = ["Tigit", "Havcr2"]
markers["NK"] = ["Klrd1"]
markers["B"] = ["Ms4a1"]
# macrophage marker
markers["macrophage"] = ["Itgam", "Adgre1"]
# M1: not exist in expr "Cd64", "Cd16", "Cd32", chatgpt: "Nos2", "Tnf", "Il1b", "Il6", "Cxcl10", "Hla-dr" (not)
# provided list (refined):  "Nos2", "Cd80", "Cd86"
markers["macrophage_m1"] = ["Cd80", "Cd86", "Nos2", "Tnf", "Il1b", "Il6", "Cxcl10"]
# M2: not exist in expr "Cd206", chatgpt: "Arg1", "Mrc1", "Retnla", "Ccl18" (not), "Il10", "TGFβ" (not), "Chi3l3" (not)
# provided list (refined):  "Arg1", "Mrc1", "Cd163"
markers["macrophage_m2"] = ["Cd163", "Arg1", "Mrc1", "Retnla", "Il10"]
# Select immune cell, including mainly macrophage and T cells (Lymphoid)
adata_macrophage = adata_intact_seurat[adata_intact_seurat.obs["annot"] == "Macrophage", :]
adata_lymphoid = adata_intact_seurat[adata_intact_seurat.obs["annot"] == "Lymphoid", :]


# load all immune cell
# seurat embedding
cca_pca_immune = pd.read_csv(seurat_dir + "immune_subset/cca_pca_immune.csv", sep = "\t")
cca_umap_immune = pd.read_csv(seurat_dir + "immune_subset/cca_umap_immune.csv", sep = "\t")
meta_immune_seurat = pd.read_csv(seurat_dir + "immune_subset/meta_immune.csv", sep = "\t")
adata_immune_raw = adata_intact_raw[meta_immune_seurat.index.values,:]

adata_immune_raw.obsm["X_seurat_pca"] = cca_pca_immune.loc[adata_immune_raw.obs.index,:].values
adata_immune_raw.obsm["X_seurat_umap"] = cca_umap_immune.loc[adata_immune_raw.obs.index,:].values
adata_immune_raw.obs["seurat_cluster"] = meta_immune_seurat.loc[adata_immune_raw.obs.index,"seurat_clusters"].astype("category")

# In[]
# ---------------------------------------------------------------- #
#
# 2. Calculate umap on lymphoid and macrophage separately
#
# ---------------------------------------------------------------- #
# 
print(f"number of lymphoid cells: {adata_lymphoid.shape[0]}")
print(f"number of macrophage cells: {adata_macrophage.shape[0]}")

# Could over-smooth the visualization, no further gene filtering
# for ct in markers.keys():
#     for gene in markers[ct]:
#         if gene in adata_lymphoid.var.index.values:
#             adata_lymphoid.var.loc[gene, "highly_variable"] = True
#             adata_macrophage.var.loc[gene, "highly_variable"] = True
#         else:
#             print(gene + " not in the expression data")
            

# adata_lymphoid = adata_lymphoid[:, adata_lymphoid.var["highly_variable"]]
# adata_macrophage = adata_macrophage[:, adata_macrophage.var["highly_variable"]]

# # optional, unit variance and clip the count
# sc.pp.scale(adata_lymphoid, max_value=10)
# sc.pp.scale(adata_macrophage, max_value=10)

sc.tl.pca(adata_lymphoid, n_comps = 100)
sc.pp.neighbors(adata_lymphoid, n_neighbors = 15, n_pcs = 100)
sc.tl.umap(adata_lymphoid, min_dist = 0.1)

sc.tl.pca(adata_macrophage, n_comps = 100)
sc.pp.neighbors(adata_macrophage, n_neighbors = 15, n_pcs = 100)
sc.tl.umap(adata_macrophage, min_dist = 0.1)

# In[]
# TODO: the figure looks ugly, improve later.
# markers of lymphoid (removed) subtypes
marker = markers["CD4_T"] + markers["CD8_T"] + markers["cytotoxicity T"] + markers["exhausted T"] + markers["NK"] + markers["B"]
fig = plt.figure(figsize = (len(marker) * 7, 10))
axs = fig.subplots(nrows = 2, ncols = len(marker))

for idx, gene in enumerate(marker):
    # NOTE: here we use the standardized expression score (max 1), since we only use the relative level to separate subtypes
    adata_lymphoid.obs[gene + "_score"] = adata_lymphoid[:,gene].X.toarray()/(np.max(adata_lymphoid[:,gene].X.toarray()) + 1e-7)
    # PP
    sc.pl.umap(adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"], color = gene + "_score", color_map = utils.SUPER_MAGMA, ax = axs[0, idx], show = False, s = 60, alpha = 1)
    # PPC
    sc.pl.umap(adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"], color = gene + "_score", color_map = utils.SUPER_MAGMA, ax = axs[1, idx], show = False, s = 60, alpha = 1)
    axs[0, idx].set_xlim([np.min(adata_lymphoid.obsm["X_umap"][:,0]) - 0.5,np.max(adata_lymphoid.obsm["X_umap"][:,0]) + 0.5])
    axs[0, idx].set_ylim([np.min(adata_lymphoid.obsm["X_umap"][:,1]) - 0.5,np.max(adata_lymphoid.obsm["X_umap"][:,1]) + 0.5])
    axs[1, idx].set_xlim([np.min(adata_lymphoid.obsm["X_umap"][:,0]) - 0.5,np.max(adata_lymphoid.obsm["X_umap"][:,0]) + 0.5])
    axs[1, idx].set_ylim([np.min(adata_lymphoid.obsm["X_umap"][:,1]) - 0.5,np.max(adata_lymphoid.obsm["X_umap"][:,1]) + 0.5])
# fig.suptitle(ct, fontsize = 30)
fig.savefig(vis_intact_dir + f"markers_lymphoid.pdf", bbox_inches = "tight")



# markers of Macrophage and its subtypes
marker = markers["macrophage_m1"] + markers["macrophage_m2"] #+ markers["macrophage"]
fig = plt.figure(figsize = (len(marker) * 7, 10))
axs = fig.subplots(nrows = 2, ncols = len(marker))

for idx, gene in enumerate(marker):
    # NOTE: here we use the standardized expression score (max 1), since we only use the relative level to separate subtypes
    adata_macrophage.obs[gene + "_score"] = adata_macrophage[:,gene].X.toarray()/(np.max(adata_macrophage[:,gene].X.toarray()) + 1e-7)
    # PP
    sc.pl.umap(adata_macrophage[adata_macrophage.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"], color = gene + "_score", color_map = utils.SUPER_MAGMA, ax = axs[0, idx], show = False, s = 60, alpha = 1)
    # PPC
    sc.pl.umap(adata_macrophage[adata_macrophage.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"], color = gene + "_score", color_map = utils.SUPER_MAGMA, ax = axs[1, idx], show = False, s = 60, alpha = 1)
    axs[0, idx].set_xlim([np.min(adata_macrophage.obsm["X_umap"][:,0]) - 0.5,np.max(adata_macrophage.obsm["X_umap"][:,0]) + 0.5])
    axs[0, idx].set_ylim([np.min(adata_macrophage.obsm["X_umap"][:,1]) - 0.5,np.max(adata_macrophage.obsm["X_umap"][:,1]) + 0.5])
    axs[1, idx].set_xlim([np.min(adata_macrophage.obsm["X_umap"][:,0]) - 0.5,np.max(adata_macrophage.obsm["X_umap"][:,0]) + 0.5])
    axs[1, idx].set_ylim([np.min(adata_macrophage.obsm["X_umap"][:,1]) - 0.5,np.max(adata_macrophage.obsm["X_umap"][:,1]) + 0.5])
# fig.suptitle(ct, fontsize = 30)
fig.savefig(vis_intact_dir + f"markers_macrophage.pdf", bbox_inches = "tight")


# In[]
# ---------------------------------------------------------------- #
#
# 3. separate T/macrophage cell subtypes
#  * count the number of cells expression the T cell subtype markers (percentage) under two different genotypes
#  * boxplot showing the nonzero marker gene expression (remove zeros) level under two different genotypes
#
# ---------------------------------------------------------------- #

# 1. Lymphoid cell subtype comparison
# extract expression
adata_lymphoid_total = adata_lymphoid.copy()
lymphoid_sub_df = pd.DataFrame(columns = ["Subtype", "GT", "Percentage", "Mean expr (nonzero)"])
for age in ["12wk", "18wk"]:
    adata_lymphoid = adata_lymphoid_total[adata_lymphoid_total.obs["age"] == age]
    x_pp_cd4T = adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)", "Cd4"].X.toarray().squeeze()
    x_ppc_cd4T = adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", "Cd4"].X.toarray().squeeze()
    x_pp_cd8T = adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)", "Cd8a"].X.toarray().squeeze()
    x_ppc_cd8T = adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", "Cd8a"].X.toarray().squeeze()
    x_pp_cytoT = adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)", ["Gzmb", "Prf1"]].X.toarray().squeeze()
    x_ppc_cytoT = adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ["Gzmb", "Prf1"]].X.toarray().squeeze()
    x_pp_exhaustT = adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)", ["Tigit", "Havcr2"]].X.toarray().squeeze()
    x_ppc_exhaustT = adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ["Tigit", "Havcr2"]].X.toarray().squeeze()

    # calculate percentage
    pct_pp_cd4T = np.sum((x_pp_cd4T > 0))/x_pp_cd4T.shape[0]
    pct_ppc_cd4T = np.sum((x_ppc_cd4T > 0))/x_ppc_cd4T.shape[0]
    pct_pp_cd8T = np.sum((x_pp_cd8T > 0))/x_pp_cd8T.shape[0]
    pct_ppc_cd8T = np.sum((x_ppc_cd8T > 0))/x_ppc_cd8T.shape[0]
    pct_pp_cytoT = np.sum((x_pp_cytoT[:,0] > 0) & (x_pp_cytoT[:,1] > 0))/x_pp_cytoT.shape[0]
    pct_ppc_cytoT = np.sum((x_ppc_cytoT[:,0] > 0) & (x_ppc_cytoT[:,1] > 0))/x_ppc_cytoT.shape[0]
    pct_pp_exhaustT = np.sum((x_pp_exhaustT[:,0] > 0) & (x_pp_exhaustT[:,1] > 0))/x_pp_exhaustT.shape[0]
    pct_ppc_exhaustT = np.sum((x_ppc_exhaustT[:,0] > 0) & (x_ppc_exhaustT[:,1] > 0))/x_ppc_exhaustT.shape[0]

    expr_pp_cd4T = np.mean(x_pp_cd4T[x_pp_cd4T > 0])
    expr_ppc_cd4T = np.mean(x_ppc_cd4T[x_ppc_cd4T > 0])
    expr_pp_cd8T = np.mean(x_pp_cd8T[x_pp_cd8T > 0])
    expr_ppc_cd8T = np.mean(x_ppc_cd8T[x_ppc_cd8T > 0])
    expr_pp_cytoT = np.mean(x_pp_cytoT[(x_pp_cytoT[:,0] > 0) & (x_pp_cytoT[:,1] > 0)])
    expr_ppc_cytoT = np.mean(x_ppc_cytoT[(x_ppc_cytoT[:,0] > 0) & (x_ppc_cytoT[:,1] > 0)])
    expr_pp_exhaustT = np.mean(x_pp_exhaustT[(x_pp_exhaustT[:,0] > 0) & (x_pp_exhaustT[:,1] > 0)])
    expr_ppc_exhaustT = np.mean(x_ppc_exhaustT[(x_ppc_exhaustT[:,0] > 0) & (x_ppc_exhaustT[:,1] > 0)])

    lymphoid_sub_df = lymphoid_sub_df._append({"Subtype": "CD4-T", "GT": "PP " + age, "Percentage": pct_pp_cd4T, "Mean expr (nonzero)": expr_pp_cd4T}, ignore_index = True)
    lymphoid_sub_df = lymphoid_sub_df._append({"Subtype": "CD4-T", "GT": "PPC " + age, "Percentage": pct_ppc_cd4T, "Mean expr (nonzero)": expr_ppc_cd4T}, ignore_index = True)
    lymphoid_sub_df = lymphoid_sub_df._append({"Subtype": "CD8-T", "GT": "PP " + age, "Percentage": pct_pp_cd8T, "Mean expr (nonzero)": expr_pp_cd8T}, ignore_index = True)
    lymphoid_sub_df = lymphoid_sub_df._append({"Subtype": "CD8-T", "GT": "PPC " + age, "Percentage": pct_ppc_cd8T, "Mean expr (nonzero)": expr_ppc_cd8T}, ignore_index = True)
    lymphoid_sub_df = lymphoid_sub_df._append({"Subtype": "CytoT", "GT": "PP " + age, "Percentage": pct_pp_cytoT, "Mean expr (nonzero)": expr_pp_cytoT}, ignore_index = True)
    lymphoid_sub_df = lymphoid_sub_df._append({"Subtype": "CytoT", "GT": "PPC " + age, "Percentage": pct_ppc_cytoT, "Mean expr (nonzero)": expr_ppc_cytoT}, ignore_index = True)
    lymphoid_sub_df = lymphoid_sub_df._append({"Subtype": "ExhaustT", "GT": "PP " + age, "Percentage": pct_pp_exhaustT, "Mean expr (nonzero)": expr_pp_exhaustT}, ignore_index = True)
    lymphoid_sub_df = lymphoid_sub_df._append({"Subtype": "ExhaustT", "GT": "PPC " + age, "Percentage": pct_ppc_exhaustT, "Mean expr (nonzero)": expr_ppc_exhaustT}, ignore_index = True)
lymphoid_sub_df.to_csv(vis_intact_dir + "lymphoid_subtypes.csv")

# In[]
sns.set_theme()
fig = plt.figure(figsize = (20, 7))
ax = fig.subplots(nrows = 1, ncols = 2)
ax[0] = sns.barplot(data = lymphoid_sub_df, x = "Subtype", y = "Percentage", hue = "GT", ax = ax[0], width=0.8)
sns.move_legend(ax[0], "upper left", bbox_to_anchor=(1.04, 1), frameon = False, fontsize = 15, title = None)
ax[1] = sns.barplot(data = lymphoid_sub_df, x = "Subtype", y = "Mean expr (nonzero)", hue = "GT", ax = ax[1], width=0.8)
sns.move_legend(ax[1], "upper left", bbox_to_anchor=(1.04, 1), frameon = False, fontsize = 15, title = None)
ax[0].set_xlabel("Lymphoid Subtype", fontsize = 15)
ax[0].set_ylabel("Proportion", fontsize = 15)
ax[1].set_xlabel("Lymphoid Subtype", fontsize = 15)
ax[1].set_ylabel("Mean expr", fontsize = 15)
for i in ax[0].containers:
    ax[0].bar_label(i, fmt = "%.2f", fontsize = 15)
for i in ax[1].containers:
    ax[1].bar_label(i, fmt = "%.2f", fontsize = 15)
plt.tight_layout()
# sns.reset_orig()
fig.savefig(vis_intact_dir + "lymphoid_subtypes_sep.png", bbox_inches = "tight", dpi = 150)

# In[]
# 2. Macrophage subtypes
adata_macrophage_total = adata_macrophage.copy()
macrophage_sub_df = pd.DataFrame(columns = ["Subtype", "GT", "Percentage", "Mean expr (nonzero)"])
for age in ["12wk", "18wk"]:
    adata_macrophage = adata_macrophage_total[adata_macrophage_total.obs["age"] == age]
    x_pp_m1 = adata_macrophage[adata_macrophage.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)", ["Cd80", "Cd86"]].X.toarray().squeeze() # NOTE: expression of Nos2 is too few
    x_ppc_m1 = adata_macrophage[adata_macrophage.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ["Cd80", "Cd86"]].X.toarray().squeeze()
    x_pp_m2 = adata_macrophage[adata_macrophage.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)", ["Arg1", "Mrc1", "Cd163"]].X.toarray().squeeze()
    x_ppc_m2 = adata_macrophage[adata_macrophage.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ["Arg1", "Mrc1", "Cd163"]].X.toarray().squeeze()

    pct_pp_m1 = np.sum((x_pp_m1[:,0] > 0) & (x_pp_m1[:,1] > 0))/x_pp_m1.shape[0]
    pct_ppc_m1 = np.sum((x_ppc_m1[:,0] > 0) & (x_ppc_m1[:,1] > 0))/x_ppc_m1.shape[0]
    pct_pp_m2 = np.sum((x_pp_m2[:,0] > 0) & (x_pp_m2[:,1] > 0) & (x_pp_m2[:,2] > 0))/x_pp_m2.shape[0]
    pct_ppc_m2 = np.sum((x_ppc_m2[:,0] > 0) & (x_ppc_m2[:,1] > 0) & (x_ppc_m2[:,2] > 0))/x_ppc_m2.shape[0]

    expr_pp_m1 = np.mean(x_pp_m1[(x_pp_m1[:,0] > 0) & (x_pp_m1[:,1] > 0)])
    expr_ppc_m1 = np.mean(x_ppc_m1[(x_ppc_m1[:,0] > 0) & (x_ppc_m1[:,1] > 0)])
    expr_pp_m2 = np.mean(x_pp_m2[(x_pp_m2[:,0] > 0) & (x_pp_m2[:,1] > 0) & (x_pp_m2[:,2] > 0)])
    expr_ppc_m2 = np.mean(x_ppc_m2[(x_ppc_m2[:,0] > 0) & (x_ppc_m2[:,1] > 0) & (x_ppc_m2[:,2] > 0)])

    macrophage_sub_df = macrophage_sub_df._append({"Subtype": "M1 Macro", "GT": "PP " + age, "Percentage": pct_pp_m1, "Mean expr (nonzero)": expr_pp_m1}, ignore_index = True)
    macrophage_sub_df = macrophage_sub_df._append({"Subtype": "M1 Macro", "GT": "PPC " + age, "Percentage": pct_ppc_m1, "Mean expr (nonzero)": expr_ppc_m1}, ignore_index = True)
    macrophage_sub_df = macrophage_sub_df._append({"Subtype": "M2 Macro", "GT": "PP " + age, "Percentage": pct_pp_m2, "Mean expr (nonzero)": expr_pp_m2}, ignore_index = True)
    macrophage_sub_df = macrophage_sub_df._append({"Subtype": "M2 Macro", "GT": "PPC " + age, "Percentage": pct_ppc_m2, "Mean expr (nonzero)": expr_ppc_m2}, ignore_index = True)

macrophage_sub_df.to_csv(vis_intact_dir + "macrophage_subtypes.csv")

# In[]
sns.set_theme()
fig = plt.figure(figsize = (20, 7))
ax = fig.subplots(nrows = 1, ncols = 2)
ax[0] = sns.barplot(data = macrophage_sub_df, x = "Subtype", y = "Percentage", hue = "GT", ax = ax[0], width=0.4)
sns.move_legend(ax[0], "upper left", bbox_to_anchor=(1.04, 1), frameon = False, fontsize = 15, title = None)
ax[1] = sns.barplot(data = macrophage_sub_df, x = "Subtype", y = "Mean expr (nonzero)", hue = "GT", ax = ax[1], width=0.4)
sns.move_legend(ax[1], "upper left", bbox_to_anchor=(1.04, 1), frameon = False, fontsize = 15, title = None)
ax[0].set_xlabel("Macrophage Subtype", fontsize = 15)
ax[0].set_ylabel("Proportion", fontsize = 15)
ax[1].set_xlabel("Macrophage Subtype", fontsize = 15)
ax[1].set_ylabel("Mean expr", fontsize = 15)
for i in ax[0].containers:
    ax[0].bar_label(i, fmt = "%.3f", fontsize = 15)
for i in ax[1].containers:
    ax[1].bar_label(i, fmt = "%.3f", fontsize = 15)
plt.tight_layout()
# sns.reset_orig()
fig.savefig(vis_intact_dir + "macropahge_subtypes_sep.png", bbox_inches = "tight", dpi = 150)



# %%
