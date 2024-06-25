# In[]
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm

# Define the colormap from grey to white to purple
colors = [(0.7, 0.7, 0.7), (1, 1, 1), (1, 0, 0)]  # Grey to White to Purple
cmap_name = 'grey_to_white_to_purple'
cmap_expr = LinearSegmentedColormap.from_list(cmap_name, colors)

result_dir = "results_seurat_visium/"
# In[]
# Read VISIUM data
adata_pp12 = sc.read_visium(path = "spatial_data/Visium_python/225_PP12/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_pp18 = sc.read_visium(path = "spatial_data/Visium_python/1687_PP18/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_ppc12 = sc.read_visium(path = "spatial_data/Visium_python/1161_PPC/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_ppc18 = sc.read_visium(path = "spatial_data/Visium_python/1660_PPC18/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_pp12.obsm["spatial"] = adata_pp12.obsm["spatial"].astype(np.float64)
adata_pp18.obsm["spatial"] = adata_pp18.obsm["spatial"].astype(np.float64)
adata_ppc12.obsm["spatial"] = adata_ppc12.obsm["spatial"].astype(np.float64)
adata_ppc18.obsm["spatial"] = adata_ppc18.obsm["spatial"].astype(np.float64)

# Load Seurat transferred labels, and select the label with maximum value for each spot
transfer_label_pp12 = pd.read_csv("spatial_data/Visium/transfer_labels_225_pp12.csv", index_col = 0, sep = " ")
transfer_label_pp18 = pd.read_csv("spatial_data/Visium/transfer_labels_1687_pp18.csv", index_col = 0, sep = " ")
transfer_label_ppc12 = pd.read_csv("spatial_data/Visium/transfer_labels_1161_ppc.csv", index_col = 0, sep = " ")
transfer_label_ppc18 = pd.read_csv("spatial_data/Visium/transfer_labels_1660_ppc18.csv", index_col = 0, sep = " ")
adata_pp12.obs["predict_label"] = transfer_label_pp12["predict.label"]
adata_pp18.obs["predict_label"] = transfer_label_pp18["predict.label"]
adata_ppc12.obs["predict_label"] = transfer_label_ppc12["predict.label"]
adata_ppc18.obs["predict_label"] = transfer_label_ppc18["predict.label"]
adata_pp12.obs.loc[adata_pp12.obs["predict_label"] == "SV, Spink1+", "predict_label"] = "SV"
adata_pp18.obs.loc[adata_pp18.obs["predict_label"] == "SV, Spink1+", "predict_label"] = "SV"
adata_ppc12.obs.loc[adata_ppc12.obs["predict_label"] == "SV, Spink1+", "predict_label"] = "SV"
adata_ppc18.obs.loc[adata_ppc18.obs["predict_label"] == "SV, Spink1+", "predict_label"] = "SV"
adata_pp12.obs.loc[adata_pp12.obs["predict_label"] == "Luminal 1", "predict_label"] = "Luminal"
adata_pp18.obs.loc[adata_pp18.obs["predict_label"] == "Luminal 1", "predict_label"] = "Luminal"
adata_ppc12.obs.loc[adata_ppc12.obs["predict_label"] == "Luminal 1", "predict_label"] = "Luminal"
adata_ppc18.obs.loc[adata_ppc18.obs["predict_label"] == "Luminal 1", "predict_label"] = "Luminal"

category = ['Basal', 'Club epithelia (Luminal)', 'Endothelial', 'Hillock epithelia (Basal)', 'Luminal', 'Lymphoid (Total immune)', 
            'Macrophage (Myeloid, Total immune)', 'Mesenchymal', 'Monocytes (Myeloid, Total immune)', 'SV', 'Spink1+ (Luminal)']

adata_pp12.obs["predict_label"] = pd.Categorical(adata_pp12.obs["predict_label"], categories=category)
adata_pp18.obs["predict_label"] = pd.Categorical(adata_pp18.obs["predict_label"], categories=category)
adata_ppc12.obs["predict_label"] = pd.Categorical(adata_ppc12.obs["predict_label"], categories=category)
adata_ppc18.obs["predict_label"] = pd.Categorical(adata_ppc18.obs["predict_label"], categories=category)
# Make unique var names
adata_pp12.var_names_make_unique()
adata_pp18.var_names_make_unique()
adata_ppc12.var_names_make_unique()
adata_ppc18.var_names_make_unique()

# In[]
# preprocessing
sc.pp.normalize_total(adata_pp12, inplace=True, target_sum = 100000)
sc.pp.log1p(adata_pp12)
adata_pp12.raw = adata_pp12
# sc.pp.highly_variable_genes(adata_pp12, flavor="seurat", n_top_genes=2000)
sc.pp.normalize_total(adata_pp18, inplace=True, target_sum = 100000)
sc.pp.log1p(adata_pp18)
adata_pp18.raw = adata_pp18
sc.pp.normalize_total(adata_ppc12, inplace=True, target_sum = 100000)
sc.pp.log1p(adata_ppc12)
adata_ppc12.raw = adata_ppc12
sc.pp.normalize_total(adata_ppc18, inplace=True, target_sum = 100000)
sc.pp.log1p(adata_ppc18)
adata_ppc18.raw = adata_ppc18

# # In[]
# # Plot the transferred labels
# fig = plt.figure(figsize = (12,7))
# ax = fig.add_subplot()
# sc.pl.spatial(adata_pp12, img_key = "hires", color = "predict_label", ax = ax)
# fig.savefig(result_dir + f"seurat_transfer_pp12.pdf", bbox_inches = "tight")
# fig.savefig(result_dir + f"seurat_transfer_pp12.png", bbox_inches = "tight")
# fig = plt.figure(figsize = (12,7))
# ax = fig.add_subplot()
# sc.pl.spatial(adata_pp18, img_key = "hires", color = "predict_label", ax = ax)
# fig.savefig(result_dir + f"seurat_transfer_pp18.pdf", bbox_inches = "tight")
# fig.savefig(result_dir + f"seurat_transfer_pp18.png", bbox_inches = "tight")
# fig = plt.figure(figsize = (12,7))
# ax = fig.add_subplot()
# sc.pl.spatial(adata_ppc12, img_key = "hires", color = "predict_label", ax = ax)
# fig.savefig(result_dir + f"seurat_transfer_ppc12.pdf", bbox_inches = "tight")
# fig.savefig(result_dir + f"seurat_transfer_ppc12.png", bbox_inches = "tight")
# fig = plt.figure(figsize = (12,7))
# ax = fig.add_subplot()
# sc.pl.spatial(adata_ppc18, img_key = "hires", color = "predict_label", ax = ax)
# fig.savefig(result_dir + f"seurat_transfer_ppc18.pdf", bbox_inches = "tight")
# fig.savefig(result_dir + f"seurat_transfer_ppc18.png", bbox_inches = "tight")


# In[]
# ------------------------------------------------------------------------------------------
# 
# Spatial distribution of the expression of interest genes
#
# ------------------------------------------------------------------------------------------
# 
# Extract the interest genes
interest_genes = pd.read_csv("spatial_data/interest_genes.csv", index_col = 0)
gene_names =[x for x in interest_genes["Genes"] if (x != "Cxcl11") & (x != "Tnfa") & (x != "Tgfb") & (x != "Il12") & (x != "Pd1") & (x != "Tim3")]

# SUPER_MAGMA = LinearSegmentedColormap.from_list('super_magma', colors=['#e0e0e0', '#dedede', '#fff68f', '#ffec8b', '#ffc125', '#ee7600', '#ee5c42', '#cd3278', '#c71585', '#68228b'], N=500)
SUPER_MAGMA = LinearSegmentedColormap.from_list('green_to_red', [(0, 0.7, 0), (1, 1, 1), (1, 0, 0), (0.5, 0, 0.5)], N=500)

# plot the spatial expression of important genes
for gene in gene_names:
    fig = plt.figure(figsize = (20, 10))
    ax = fig.subplots(nrows = 2, ncols = 2)

    vmin = min(np.min(adata_pp12[:, gene].X), np.min(adata_ppc12[:, gene].X), np.min(adata_pp18[:, gene].X), np.min(adata_ppc18[:, gene].X))
    vmax = max(np.max(adata_pp12[:, gene].X), np.max(adata_ppc12[:, gene].X), np.max(adata_pp18[:, gene].X), np.max(adata_ppc18[:, gene].X))

    sc.pl.spatial(adata_pp12, img_key="hires", color = gene, save = None, ax = ax[0,0], show = False, vmin = vmin, vmax = vmax, color_map = SUPER_MAGMA)
    sc.pl.spatial(adata_pp18, img_key="hires", color = gene, save = None, ax = ax[1,0], show = False, vmin = vmin, vmax = vmax, color_map = SUPER_MAGMA)
    sc.pl.spatial(adata_ppc12, img_key="hires", color = gene, save = None, ax = ax[0,1], show = False, vmin = vmin, vmax = vmax, color_map = SUPER_MAGMA)
    sc.pl.spatial(adata_ppc18, img_key="hires", color = gene, save = None, ax = ax[1,1], show = False, vmin = vmin, vmax = vmax, color_map = SUPER_MAGMA)
    
    fig.suptitle(gene, fontsize = 20)
    ax[0,0].set_title("PP12")
    ax[0,1].set_title("PP18")
    ax[1,0].set_title("PPC12")
    ax[1,1].set_title("PPC18")
    fig.tight_layout()

    fig.savefig(f"results_seurat_scrnaseq/figure_markers/important_genes/visium_{gene}.png", dpi = 200, bbox_inches = "tight")

    

# # In[]
# # interest ligand-receptor pairs
# lr_cxc = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "CXC")
# lr_ccl = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "CCL")
# lr_il = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "IL")
# lr_tgfb = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "TGFB")
# lr_ifn = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "IFN")
# lr_tnf = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "TNF family")
# lr_immune = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "immune_checkpoint")
# lr_cd = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "CD")
# lr_mhc = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "MHC")
# lr_EGF_FGF_IGF = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "EGF_FGF_IGF")

# lr_cxc["sheet_name"] = "CXC"
# lr_ccl["sheet_name"] = "CCL"
# lr_il["sheet_name"] = "IL"
# lr_tgfb["sheet_name"] = "TGFB"
# lr_ifn["sheet_name"] = "IFN"
# lr_tnf["sheet_name"] = "TNF family"
# lr_immune["sheet_name"] = "immune_checkpoint"
# lr_cd["sheet_name"] = "CD"
# lr_mhc["sheet_name"] = "MHC"
# lr_EGF_FGF_IGF["sheet_name"] = "EGF_FGF_IGF"

# lr = pd.concat([lr_cxc, lr_ccl, lr_il, lr_tgfb, lr_ifn, lr_tnf, lr_immune, lr_cd, lr_mhc, lr_EGF_FGF_IGF], axis = 0, ignore_index = True)
# lr.index = np.arange(lr.shape[0])
# lr.to_csv("spatial_data/ligand_receptor_reformat.csv")    

# In[]
# ------------------------------------------------------------------------------------------
# 
# Boxplot of interest genes across genotypes
#
# ------------------------------------------------------------------------------------------
#
# Trp53, Pten, Ackr3: perturbed genes
# Cxcr4: interacts
# Cxcl12, Mif: potential ligand
# Spink1: affected
# Mki67: proliferation

# markers for other cell types
# Ptprc: immune cells
# Cd4: CD4+ T cells
# Cd8a: CD8+ T cells
# Foxp3: Treg T cells
# Ctla4: T cells
# Itgax: Myeloid cells (DCs)
# Il2, Il4, Il6: Interleukin
# Cxcl2, Cxcl12: Chemokine Ligand
# Tmprss4: Transmembrane Serine Protease
# Syp, Chga, Sox2, Mycn: NE 
# Dcn, FN1: Fibroblasts
# Nkx3-1, Pbsn: L1
# LGR5: Stem cells
# Ccl21: Endothelial cells
# Cxcr5, Ms4a1: B-cells 
# Klrd1, Klrb1c: NK cell
# Gzmb, Prf1, Tnfa, Ifng, Tgfb, Il6: Cytotoxic

interest_genes = pd.read_csv("../dataset/markers_info/interest_genes.csv", index_col = 0).index.values
functions = pd.read_csv("../dataset/markers_info/interest_genes.csv", index_col = 0)["Note"].values
interest_genes = [x for x in interest_genes if (x != "Cxcl11") & (x != "Tnfa") & (x != "Tgfb") & (x != "Il12") & (x != "Pd1") & (x != "Tim3") & (x != "Ccl21")]
full_names = [x + "_" + y for x,y in zip(interest_genes, functions) if (x != "Cxcl11") & (x != "Tnfa") & (x != "Tgfb") & (x != "Il12") & (x != "Pd1") & (x != "Tim3") & (x != "Ccl21")]


# In[]
# ------------------------------------------------------------------------------------------
#
# Visualize CellChat cci result, CellChat requires cell type annotation for each 
# spot in visium data (could be inaccurate)
#
# ------------------------------------------------------------------------------------------
# 
cci_dir = "results_cellchat_visium/"
cci_pp18 = pd.read_csv(cci_dir + "cci_pp18.csv", sep = "\t")
clusters = ["Luminal 1", "Mesenchymal"]
cci_clusts_pp18 = {}
cci_clusts_top_pp18 = {}
for source in clusters:
    for target in clusters:
        cci_clust = cci_pp18[(cci_pp18["source"] == source) & (cci_pp18["target"] == target)]
        cci_clust = cci_clust.sort_values(by = ["prob"], ascending = False)
        cci_clusts_pp18[source + "-" + target] = cci_clust
        cci_clust_top = cci_clust.iloc[:20,:]
        cci_clusts_top_pp18[source + "-" + target] = cci_clust_top


cci_ppc18 = pd.read_csv(cci_dir + "cci_ppc18.csv", sep = "\t")
clusters = ["Luminal 1", "Mesenchymal"]
cci_clusts_ppc18 = {}
cci_clusts_top_ppc18 = {}
for source in clusters:
    for target in clusters:
        cci_clust = cci_ppc18[(cci_ppc18["source"] == source) & (cci_ppc18["target"] == target)]
        cci_clust = cci_clust.sort_values(by = ["prob"], ascending = False)
        cci_clusts_ppc18[source + "-" + target] = cci_clust
        cci_clust_top = cci_clust.iloc[:20,:]
        cci_clusts_top_ppc18[source + "-" + target] = cci_clust_top


scores_dict = {}

columns = np.union1d(cci_clusts_top_pp18["Luminal 1-Luminal 1"]["interaction_name_2"].values, cci_clusts_top_ppc18["Luminal 1-Luminal 1"]["interaction_name_2"].values)
scores = pd.DataFrame(data = 0, columns = columns, index = ["pp18", "ppc18"])
for col in columns:
    if col in cci_clusts_pp18["Luminal 1-Luminal 1"]["interaction_name_2"].values:
        scores.loc["pp18", col] = np.sum(cci_clusts_pp18["Luminal 1-Luminal 1"].loc[cci_clusts_pp18["Luminal 1-Luminal 1"]["interaction_name_2"] == col, "prob"])
    if col in cci_clusts_ppc18["Luminal 1-Luminal 1"]["interaction_name_2"].values:
        scores.loc["ppc18", col] = np.sum(cci_clusts_ppc18["Luminal 1-Luminal 1"].loc[cci_clusts_ppc18["Luminal 1-Luminal 1"]["interaction_name_2"] == col, "prob"])
scores_dict["Luminal 1-Luminal 1"] = scores

columns = np.union1d(cci_clusts_top_pp18["Luminal 1-Mesenchymal"]["interaction_name_2"].values, cci_clusts_top_ppc18["Luminal 1-Mesenchymal"]["interaction_name_2"].values)
scores = pd.DataFrame(data = 0, columns = columns, index = ["pp18", "ppc18"])
for col in columns:
    if col in cci_clusts_pp18["Luminal 1-Mesenchymal"]["interaction_name_2"].values:
        scores.loc["pp18", col] = np.sum(cci_clusts_pp18["Luminal 1-Mesenchymal"].loc[cci_clusts_pp18["Luminal 1-Mesenchymal"]["interaction_name_2"] == col, "prob"])
    if col in cci_clusts_ppc18["Luminal 1-Mesenchymal"]["interaction_name_2"].values:
        scores.loc["ppc18", col] = np.sum(cci_clusts_ppc18["Luminal 1-Mesenchymal"].loc[cci_clusts_ppc18["Luminal 1-Mesenchymal"]["interaction_name_2"] == col, "prob"])
scores_dict["Luminal 1-Mesenchymal"] = scores

columns = np.union1d(cci_clusts_top_pp18["Mesenchymal-Mesenchymal"]["interaction_name_2"].values, cci_clusts_top_ppc18["Mesenchymal-Mesenchymal"]["interaction_name_2"].values)
scores = pd.DataFrame(data = 0, columns = columns, index = ["pp18", "ppc18"])
for col in columns:
    if col in cci_clusts_pp18["Mesenchymal-Mesenchymal"]["interaction_name_2"].values:
        scores.loc["pp18", col] = np.sum(cci_clusts_pp18["Mesenchymal-Mesenchymal"].loc[cci_clusts_pp18["Mesenchymal-Mesenchymal"]["interaction_name_2"] == col, "prob"])
    if col in cci_clusts_ppc18["Mesenchymal-Mesenchymal"]["interaction_name_2"].values:
        scores.loc["ppc18", col] = np.sum(cci_clusts_ppc18["Mesenchymal-Mesenchymal"].loc[cci_clusts_ppc18["Mesenchymal-Mesenchymal"]["interaction_name_2"] == col, "prob"])
scores_dict["Mesenchymal-Mesenchymal"] = scores

columns = np.union1d(cci_clusts_top_pp18["Mesenchymal-Luminal 1"]["interaction_name_2"].values, cci_clusts_top_ppc18["Mesenchymal-Luminal 1"]["interaction_name_2"].values)
scores = pd.DataFrame(data = 0, columns = columns, index = ["pp18", "ppc18"])
for col in columns:
    if col in cci_clusts_pp18["Mesenchymal-Luminal 1"]["interaction_name_2"].values:
        scores.loc["pp18", col] = np.sum(cci_clusts_pp18["Mesenchymal-Luminal 1"].loc[cci_clusts_pp18["Mesenchymal-Luminal 1"]["interaction_name_2"] == col, "prob"])
    if col in cci_clusts_ppc18["Mesenchymal-Luminal 1"]["interaction_name_2"].values:
        scores.loc["ppc18", col] = np.sum(cci_clusts_ppc18["Mesenchymal-Luminal 1"].loc[cci_clusts_ppc18["Mesenchymal-Luminal 1"]["interaction_name_2"] == col, "prob"])
scores_dict["Mesenchymal-Luminal 1"] = scores

import seaborn as sns
fig = plt.figure(figsize = (30, 7))
ax = fig.subplots(nrows = 2, ncols = 2)
sns.heatmap(scores_dict["Luminal 1-Luminal 1"], ax = ax[0,0], annot=True)
sns.heatmap(scores_dict["Luminal 1-Mesenchymal"], ax = ax[0,1], annot=True)
sns.heatmap(scores_dict["Mesenchymal-Luminal 1"], ax = ax[1,0], annot=True)
sns.heatmap(scores_dict["Mesenchymal-Mesenchymal"], ax = ax[1,1], annot=True)
ax[0,0].set_title("Luminal 1-Luminal 1")
ax[0,1].set_title("Luminal 1-Mesenchymal")
ax[1,0].set_title("Mesenchymal-Luminal 1")
ax[1,1].set_title("Mesenchymal-Mesenchymal")
plt.tight_layout()

fig.savefig(cci_dir + "cci_compare_18.pdf", bbox_inches = "tight")

# In[]
cci_pp12 = pd.read_csv(cci_dir + "cci_pp12.csv", sep = "\t")
clusters = ["Luminal 1", "Mesenchymal"]
cci_clusts_pp12 = {}
cci_clusts_top_pp12 = {}
for source in clusters:
    for target in clusters:
        cci_clust = cci_pp12[(cci_pp12["source"] == source) & (cci_pp12["target"] == target)]
        cci_clust = cci_clust.sort_values(by = ["prob"], ascending = False)
        cci_clusts_pp12[source + "-" + target] = cci_clust
        cci_clust_top = cci_clust.iloc[:20,:]
        cci_clusts_top_pp12[source + "-" + target] = cci_clust_top


cci_ppc12 = pd.read_csv(cci_dir + "cci_ppc12.csv", sep = "\t")
clusters = ["Luminal 1", "Mesenchymal"]
cci_clusts_ppc12 = {}
cci_clusts_top_ppc12 = {}
for source in clusters:
    for target in clusters:
        cci_clust = cci_ppc12[(cci_ppc12["source"] == source) & (cci_ppc12["target"] == target)]
        cci_clust = cci_clust.sort_values(by = ["prob"], ascending = False)
        cci_clusts_ppc12[source + "-" + target] = cci_clust
        cci_clust_top = cci_clust.iloc[:20,:]
        cci_clusts_top_ppc12[source + "-" + target] = cci_clust_top

scores_dict = {}

columns = np.union1d(cci_clusts_top_pp12["Luminal 1-Luminal 1"]["interaction_name_2"].values, cci_clusts_top_ppc12["Luminal 1-Luminal 1"]["interaction_name_2"].values)
scores = pd.DataFrame(data = 0, columns = columns, index = ["pp12", "ppc12"])
for col in columns:
    if col in cci_clusts_pp12["Luminal 1-Luminal 1"]["interaction_name_2"].values:
        scores.loc["pp12", col] = np.sum(cci_clusts_pp12["Luminal 1-Luminal 1"].loc[cci_clusts_pp12["Luminal 1-Luminal 1"]["interaction_name_2"] == col, "prob"])
    if col in cci_clusts_ppc12["Luminal 1-Luminal 1"]["interaction_name_2"].values:
        scores.loc["ppc12", col] = np.sum(cci_clusts_ppc12["Luminal 1-Luminal 1"].loc[cci_clusts_ppc12["Luminal 1-Luminal 1"]["interaction_name_2"] == col, "prob"])
scores_dict["Luminal 1-Luminal 1"] = scores

columns = np.union1d(cci_clusts_top_pp12["Luminal 1-Mesenchymal"]["interaction_name_2"].values, cci_clusts_top_ppc12["Luminal 1-Mesenchymal"]["interaction_name_2"].values)
scores = pd.DataFrame(data = 0, columns = columns, index = ["pp12", "ppc12"])
for col in columns:
    if col in cci_clusts_pp12["Luminal 1-Mesenchymal"]["interaction_name_2"].values:
        scores.loc["pp12", col] = np.sum(cci_clusts_pp12["Luminal 1-Mesenchymal"].loc[cci_clusts_pp12["Luminal 1-Mesenchymal"]["interaction_name_2"] == col, "prob"])
    if col in cci_clusts_ppc12["Luminal 1-Mesenchymal"]["interaction_name_2"].values:
        scores.loc["ppc12", col] = np.sum(cci_clusts_ppc12["Luminal 1-Mesenchymal"].loc[cci_clusts_ppc12["Luminal 1-Mesenchymal"]["interaction_name_2"] == col, "prob"])
scores_dict["Luminal 1-Mesenchymal"] = scores

columns = np.union1d(cci_clusts_top_pp12["Mesenchymal-Mesenchymal"]["interaction_name_2"].values, cci_clusts_top_ppc12["Mesenchymal-Mesenchymal"]["interaction_name_2"].values)
scores = pd.DataFrame(data = 0, columns = columns, index = ["pp12", "ppc12"])
for col in columns:
    if col in cci_clusts_pp12["Mesenchymal-Mesenchymal"]["interaction_name_2"].values:
        scores.loc["pp12", col] = np.sum(cci_clusts_pp12["Mesenchymal-Mesenchymal"].loc[cci_clusts_pp12["Mesenchymal-Mesenchymal"]["interaction_name_2"] == col, "prob"])
    if col in cci_clusts_ppc12["Mesenchymal-Mesenchymal"]["interaction_name_2"].values:
        scores.loc["ppc12", col] = np.sum(cci_clusts_ppc12["Mesenchymal-Mesenchymal"].loc[cci_clusts_ppc12["Mesenchymal-Mesenchymal"]["interaction_name_2"] == col, "prob"])
scores_dict["Mesenchymal-Mesenchymal"] = scores

columns = np.union1d(cci_clusts_top_pp12["Mesenchymal-Luminal 1"]["interaction_name_2"].values, cci_clusts_top_ppc12["Mesenchymal-Luminal 1"]["interaction_name_2"].values)
scores = pd.DataFrame(data = 0, columns = columns, index = ["pp12", "ppc12"])
for col in columns:
    if col in cci_clusts_pp12["Mesenchymal-Luminal 1"]["interaction_name_2"].values:
        scores.loc["pp12", col] = np.sum(cci_clusts_pp12["Mesenchymal-Luminal 1"].loc[cci_clusts_pp12["Mesenchymal-Luminal 1"]["interaction_name_2"] == col, "prob"])
    if col in cci_clusts_ppc12["Mesenchymal-Luminal 1"]["interaction_name_2"].values:
        scores.loc["ppc12", col] = np.sum(cci_clusts_ppc12["Mesenchymal-Luminal 1"].loc[cci_clusts_ppc12["Mesenchymal-Luminal 1"]["interaction_name_2"] == col, "prob"])
scores_dict["Mesenchymal-Luminal 1"] = scores

import seaborn as sns
fig = plt.figure(figsize = (30, 7))
ax = fig.subplots(nrows = 2, ncols = 2)
sns.heatmap(scores_dict["Luminal 1-Luminal 1"], ax = ax[0,0], annot=True)
sns.heatmap(scores_dict["Luminal 1-Mesenchymal"], ax = ax[0,1], annot=True)
sns.heatmap(scores_dict["Mesenchymal-Luminal 1"], ax = ax[1,0], annot=True)
sns.heatmap(scores_dict["Mesenchymal-Mesenchymal"], ax = ax[1,1], annot=True)
ax[0,0].set_title("Luminal 1-Luminal 1")
ax[0,1].set_title("Luminal 1-Mesenchymal")
ax[1,0].set_title("Mesenchymal-Luminal 1")
ax[1,1].set_title("Mesenchymal-Mesenchymal")
plt.tight_layout()
fig.savefig(cci_dir + "cci_compare_12.pdf", bbox_inches = "tight")




# %%
