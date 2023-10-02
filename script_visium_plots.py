# In[]
import numpy as np
import scanpy as sc
import pandas as pd
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt

def calc_dist_thr(loc, K):
    D = pairwise_distances(loc)
    knn_index = np.argpartition(D, kth = K-1, axis = 1)[:, (K-1)]
    kth_dist = np.take_along_axis(D, knn_index[:,None], axis = 1)
    # maximum knn distance as dist cut-off, change with K selection
    dis_thr = np.max(kth_dist)
    return dis_thr

def save_commot_results(adata, lr_index_path, cci_path):
    LR_pairs = []
    cci_predict = []
    for LR_pair in adata.obsp.keys():
        LR_pairs.append(LR_pair)
        cci_predict.append(adata.obsp[LR_pair].toarray()[None,:,:])
    LR_pairs = np.array(LR_pairs)
    cci_predict = np.concatenate(cci_predict, axis = 0)
    print(cci_predict.shape)
    print(LR_pairs)
    np.savetxt(lr_index_path, LR_pairs, fmt = "%s")
    np.save(cci_path, cci_predict)

result_dir = "results_seurat_visium/"
# In[]
adata_pp12 = sc.read_visium(path = "spatial_data/Visium_python/225_PP12/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_pp18 = sc.read_visium(path = "spatial_data/Visium_python/1687_PP18/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_ppc12 = sc.read_visium(path = "spatial_data/Visium_python/1161_PPC/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_ppc18 = sc.read_visium(path = "spatial_data/Visium_python/1660_PPC18/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
transfer_label_pp12 = pd.read_csv("spatial_data/Visium/transfer_labels_225_pp12.csv", index_col = 0, sep = " ")
transfer_label_pp18 = pd.read_csv("spatial_data/Visium/transfer_labels_1687_pp18.csv", index_col = 0, sep = " ")
transfer_label_ppc12 = pd.read_csv("spatial_data/Visium/transfer_labels_1161_ppc.csv", index_col = 0, sep = " ")
transfer_label_ppc18 = pd.read_csv("spatial_data/Visium/transfer_labels_1660_ppc18.csv", index_col = 0, sep = " ")

adata_pp12.obs["predict_label"] = transfer_label_pp12["predict.label"]
adata_pp18.obs["predict_label"] = transfer_label_pp18["predict.label"]
adata_ppc12.obs["predict_label"] = transfer_label_ppc12["predict.label"]
adata_ppc18.obs["predict_label"] = transfer_label_ppc18["predict.label"]

adata_pp12.obsm["spatial"] = adata_pp12.obsm["spatial"].astype(np.float64)
adata_pp18.obsm["spatial"] = adata_pp18.obsm["spatial"].astype(np.float64)
adata_ppc12.obsm["spatial"] = adata_ppc12.obsm["spatial"].astype(np.float64)
adata_ppc18.obsm["spatial"] = adata_ppc18.obsm["spatial"].astype(np.float64)

adata_pp12.var_names_make_unique()
adata_pp18.var_names_make_unique()
adata_ppc12.var_names_make_unique()
adata_ppc18.var_names_make_unique()

adata_pp12.obs.loc[adata_pp12.obs["predict_label"] == "SV, Spink1+", "predict_label"] = "SV"
adata_pp18.obs.loc[adata_pp18.obs["predict_label"] == "SV, Spink1+", "predict_label"] = "SV"
adata_ppc12.obs.loc[adata_ppc12.obs["predict_label"] == "SV, Spink1+", "predict_label"] = "SV"
adata_ppc18.obs.loc[adata_ppc18.obs["predict_label"] == "SV, Spink1+", "predict_label"] = "SV"

category = ['Basal', 'Club epithelia (Luminal)', 'Endothelial',
            'Hillock epithelia (Basal)', 'Luminal 1',
            'Lymphoid (Total immune)', 'Macrophage (Myeloid, Total immune)',
            'Mesenchymal', 'Monocytes (Myeloid, Total immune)', 'SV',
            'Spink1+ (Luminal)']

adata_pp12.obs["predict_label"] = pd.Categorical(adata_pp12.obs["predict_label"], categories=category)
adata_pp18.obs["predict_label"] = pd.Categorical(adata_pp18.obs["predict_label"], categories=category)
adata_ppc12.obs["predict_label"] = pd.Categorical(adata_ppc12.obs["predict_label"], categories=category)
adata_ppc18.obs["predict_label"] = pd.Categorical(adata_ppc18.obs["predict_label"], categories=category)


# In[]
fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_pp12, img_key = "hires", color = "predict_label", ax = ax)
fig.savefig(result_dir + f"seurat_transfer_pp12.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"seurat_transfer_pp12.png", bbox_inches = "tight")
fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_pp18, img_key = "hires", color = "predict_label", ax = ax)
fig.savefig(result_dir + f"seurat_transfer_pp18.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"seurat_transfer_pp18.png", bbox_inches = "tight")
fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_ppc12, img_key = "hires", color = "predict_label", ax = ax)
fig.savefig(result_dir + f"seurat_transfer_ppc12.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"seurat_transfer_ppc12.png", bbox_inches = "tight")
fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_ppc18, img_key = "hires", color = "predict_label", ax = ax)
fig.savefig(result_dir + f"seurat_transfer_ppc18.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"seurat_transfer_ppc18.png", bbox_inches = "tight")


# In[]
interest_genes = pd.read_csv("spatial_data/interest_genes.csv", index_col = 0)
gene_names =[x for x in interest_genes["Genes"] if (x != "Cxcl11") & (x != "Tnfa") & (x != "Tgfb") & (x != "Il12")]

sc.pl.spatial(adata_pp12, img_key="hires", color = gene_names, save = "visium_interest_genes.pdf")


# In[]
# interest ligand-receptor pairs
lr_cxc = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "CXC")
lr_ccl = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "CCL")
lr_il = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "IL")
lr_tgfb = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "TGFB")
lr_ifn = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "IFN")
lr_tnf = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "TNF family")
lr_immune = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "immune_checkpoint")
lr_cd = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "CD")
lr_mhc = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "MHC")
lr_EGF_FGF_IGF = pd.read_excel("spatial_data/ligand-receptor-pairs.xlsx", sheet_name = "EGF_FGF_IGF")

lr_cxc["sheet_name"] = "CXC"
lr_ccl["sheet_name"] = "CCL"
lr_il["sheet_name"] = "IL"
lr_tgfb["sheet_name"] = "TGFB"
lr_ifn["sheet_name"] = "IFN"
lr_tnf["sheet_name"] = "TNF family"
lr_immune["sheet_name"] = "immune_checkpoint"
lr_cd["sheet_name"] = "CD"
lr_mhc["sheet_name"] = "MHC"
lr_EGF_FGF_IGF["sheet_name"] = "EGF_FGF_IGF"

lr = pd.concat([lr_cxc, lr_ccl, lr_il, lr_tgfb, lr_ifn, lr_tnf, lr_immune, lr_cd, lr_mhc, lr_EGF_FGF_IGF], axis = 0, ignore_index = True)
lr.index = np.arange(lr.shape[0])
lr.to_csv("spatial_data/ligand_receptor_reformat.csv")

# In[]
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
