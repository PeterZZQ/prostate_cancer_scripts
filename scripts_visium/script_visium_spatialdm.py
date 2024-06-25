# In[]
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import spatialdm as sdm
import spatialdm.plottings as pl
import warnings
warnings.filterwarnings("ignore")

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
# Load Seurat transferred labels
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

category = ['Basal', 'Club epithelia (Luminal)', 'Endothelial',
            'Hillock epithelia (Basal)', 'Luminal',
            'Lymphoid (Total immune)', 'Macrophage (Myeloid, Total immune)',
            'Mesenchymal', 'Monocytes (Myeloid, Total immune)', 'SV',
            'Spink1+ (Luminal)']

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
adata_pp12.raw = adata_pp12
sc.pp.normalize_total(adata_pp12, inplace=True)
sc.pp.log1p(adata_pp12)
# sc.pp.highly_variable_genes(adata_pp12, flavor="seurat", n_top_genes=2000)
adata_pp18.raw = adata_pp18
sc.pp.normalize_total(adata_pp18, inplace=True)
sc.pp.log1p(adata_pp18)
adata_ppc12.raw = adata_ppc12
sc.pp.normalize_total(adata_ppc12, inplace=True)
sc.pp.log1p(adata_ppc12)
adata_ppc18.raw = adata_ppc18
sc.pp.normalize_total(adata_ppc18, inplace=True)
sc.pp.log1p(adata_ppc18)


# In[]
lr_df = pd.read_excel("markers_info/ligand receptor genes of interest.xlsx") 
lr_genes = lr_df["Gene"].values.squeeze()
lr_genes = [x for x in lr_genes if (x != "H2-Q8") & (x != "H2-T3") & (x != "H2-T9") & (x != "H2-M10.5") \
            & (x != "H2-T10") & (x != "H2-Bl") & (x != "H2-T18") & (x != "H2-Q9") & (x != "H2-L")\
                & (x != "Gm10499") & (x != "H2-T-ps") & (x != "H2-BI") & (x != "H2-D") & (x != "Raet1a") & (x != "H60a") & (x != "Kir3dl1") & (x != "Klra")]
# Plot the spatial expression pattern of interest genes
sc.pl.spatial(adata_pp12, img_key="hires", color = lr_genes, vmin = 0, vmax = 6,  save = "pathway_genes_pp12.pdf")
sc.pl.spatial(adata_pp18, img_key="hires", color = lr_genes, vmin = 0, vmax = 6, save = "pathway_genes_pp18.pdf")
sc.pl.spatial(adata_ppc12, img_key="hires", color = lr_genes, vmin = 0, vmax = 6, save = "pathway_genes_ppc12.pdf")
sc.pl.spatial(adata_ppc18, img_key="hires", color = lr_genes, vmin = 0, vmax = 6, save = "pathway_genes_ppc18.pdf")

# In[]
# calculate correlation
adata = adata_ppc18.copy()
result_dir = "results_spatialdm_visium/ppc18/"
# weight_matrix by rbf kernel, l need to be selected carefully
sdm.weight_matrix(adata, l=200, cutoff=0.2, single_cell=False) 
# check weight scales
plt.scatter(adata.obsm['spatial'][:,0], adata.obsm['spatial'][:,1], c=adata.obsp['weight'].A[50])
# find overlapping LRs from CellChatDB
sdm.extract_lr(adata, 'mouse', min_cell=3) 
# optional: update the LR pairs to include only the MHC-I, IFN-II, and IL6 pathways
adata.uns["ligand"] = adata.uns["ligand"][adata.uns["geneInter"]["pathway_name"].isin(["MHC-I", "IFN-II", "IL6"])]    
adata.uns["receptor"] = adata.uns["receptor"][adata.uns["geneInter"]["pathway_name"].isin(["MHC-I", "IFN-II", "IL6"])]    
adata.uns["geneInter"] = adata.uns["geneInter"][adata.uns["geneInter"]["pathway_name"].isin(["MHC-I", "IFN-II", "IL6"])]    
pathway_lrs = adata.uns["geneInter"][adata.uns["geneInter"]["pathway_name"].isin(["MHC-I", "IFN-II", "IL6"])]    

# global Moran selection, method z-score, permutation or both
sdm.spatialdm_global(adata, n_perm=1000, specified_ind=None, method='both', nproc=1)     
# select significant pairs method z-score, permutation, 
# MHC-I, IFN-II, IL6 pathway related gene set is small, no need for pathway correction
sdm.sig_pairs(adata, method='permutation', fdr=True, threshold=0.1)
# print the number of selected pairs, locals are carried out on the selected
print(adata.uns['global_res'].selected.sum())
# adata.uns['global_res'].sort_values(by='fdr').head()
# local spot selection, method z-score, permutation or both, spot number more than 1000
sdm.spatialdm_local(adata, n_perm=1000, method='both', specified_ind=None, nproc=1)     
# significant local spots, method z-score, permutation or both, 
# selected pairs is small, no need for FDR correction
sdm.sig_spots(adata, method='permutation', fdr=False, threshold=0.1)     

# adata.uns["ligand"]["Ligand1"] = adata.uns["ligand"]["Ligand1"].astype("str")
# adata.uns["ligand"]["Ligand2"] = adata.uns["ligand"]["Ligand2"].astype("str")
# adata.uns["receptor"]["Receptor1"] = adata.uns["receptor"]["Receptor1"].astype("str")
# adata.uns["receptor"]["Receptor2"] = adata.uns["receptor"]["Receptor2"].astype("str")
# adata.uns["receptor"]["Receptor3"] = adata.uns["receptor"]["Receptor3"].astype("str")
# adata.uns["geneInter"]["agonist"] = adata.uns["geneInter"]["agonist"].astype("str")
# adata.uns["geneInter"]["antagonist"] = adata.uns["geneInter"]["antagonist"].astype("str")
# adata.uns["geneInter"]["co_A_receptor"] = adata.uns["geneInter"]["co_A_receptor"].astype("str")
# adata.uns["geneInter"]["co_I_receptor"] = adata.uns["geneInter"]["co_I_receptor"].astype("str")
# adata.uns["global_res"]["Ligand1"] = adata.uns["global_res"]["Ligand1"].astype("str")
# adata.uns["global_res"]["Ligand2"] = adata.uns["global_res"]["Ligand2"].astype("str")
# adata.uns["global_res"]["Receptor1"] = adata.uns["global_res"]["Receptor1"].astype("str")
# adata.uns["global_res"]["Receptor2"] = adata.uns["global_res"]["Receptor2"].astype("str")
# adata.uns["global_res"]["Receptor3"] = adata.uns["global_res"]["Receptor3"].astype("str")
# adata.write_h5ad("spatial_data/Visium_python/spatialdm_pp12.h5ad")
# In[]
global_scores = adata.uns['global_res'].sort_values(by='fdr')
global_scores["pathway"] = pathway_lrs.loc[global_scores.index.values, "pathway_name"]
global_scores.to_csv(result_dir + "coexpr_MHC-I_IFN-II_IL6.csv")
pl.plot_pairs(adata, [x for x in adata.uns['global_res'].index[adata.uns['global_res'].selected].values], marker='s', pdf = result_dir + f"coexpr_MHC-I_IFN-II_IL6_local")

# In[]
from functools import reduce
coexpr_pp12_pval = pd.read_csv("results_spatialdm_visium/pp12/coexpr_MHC-I_IFN-II_IL6.csv", sep = ",", index_col = 0)
coexpr_pp18_pval = pd.read_csv("results_spatialdm_visium/pp18/coexpr_MHC-I_IFN-II_IL6.csv", sep = ",", index_col = 0)
coexpr_ppc12_pval = pd.read_csv("results_spatialdm_visium/ppc12/coexpr_MHC-I_IFN-II_IL6.csv", sep = ",", index_col = 0)
coexpr_ppc18_pval = pd.read_csv("results_spatialdm_visium/ppc18/coexpr_MHC-I_IFN-II_IL6.csv", sep = ",", index_col = 0)

intersect = reduce(np.intersect1d, [coexpr_pp12_pval.index.values, coexpr_ppc12_pval.index.values, coexpr_pp18_pval.index.values, coexpr_ppc18_pval.index.values])
pway = coexpr_pp12_pval.loc[intersect, "pathway"]

coexpr_pp12_pval = coexpr_pp12_pval.loc[intersect,"fdr"]
coexpr_pp18_pval = coexpr_pp18_pval.loc[intersect,"fdr"]
coexpr_ppc12_pval = coexpr_ppc12_pval.loc[intersect,"fdr"]
coexpr_ppc18_pval = coexpr_ppc18_pval.loc[intersect,"fdr"]

coexpr = 1-np.concatenate([coexpr_pp12_pval.values[:,None], coexpr_ppc12_pval.values[:,None], coexpr_pp18_pval.values[:,None], coexpr_ppc18_pval.values[:,None]], axis = 1)
coexpr = pd.DataFrame(data = coexpr, index = [x + " (" + y + ")" for x,y in zip(intersect, pway)], columns = ["pp12", "ppc12", "pp18", "ppc18"])

cmap = sns.clustermap(coexpr, col_cluster = False, row_cluster = True, yticklabels=True)
cmap.figure.savefig("results_spatialdm_visium/compare_coexpr_MHC-I_IFN-II_IL6.pdf", bbox_inches = "tight")

# In[]
# calculate correlation
for adata in [adata_pp12, adata_pp18, adata_ppc12, adata_ppc18]:
    result_dir = "results_spatialdm_visium/pp12/"
    # weight_matrix by rbf kernel, l need to be selected carefully
    sdm.weight_matrix(adata, l=200, cutoff=0.2, single_cell=False) 
    # check weight scales
    plt.scatter(adata.obsm['spatial'][:,0], adata.obsm['spatial'][:,1], c=adata.obsp['weight'].A[50])
    # find overlapping LRs from CellChatDB
    sdm.extract_lr(adata, 'mouse', min_cell=3) 
    pathway_lrs = adata.uns["geneInter"][adata.uns["geneInter"]["pathway_name"].isin(["MHC-I", "IFN-II", "IL6"])]    

    # global Moran selection, method z-score, permutation or both
    sdm.spatialdm_global(adata, n_perm=1000, specified_ind=None, method='z-score', nproc=1)     
    # select significant pairs method z-score, permutation, 
    # MHC-I, IFN-II, IL6 pathway related gene set is small, no need for pathway correction
    sdm.sig_pairs(adata, method='z-score', fdr=True, threshold=0.1)
    # print the number of selected pairs, locals are carried out on the selected
    print(adata.uns['global_res'].selected.sum())
    # local spot selection, method z-score, permutation or both, spot number more than 1000
    sdm.spatialdm_local(adata, n_perm=1000, method='z-score', specified_ind=None, nproc=1)     
    # significant local spots, method z-score, permutation or both, 
    # selected pairs is small, no need for FDR correction
    sdm.sig_spots(adata, method='z-score', fdr=False, threshold=0.1)     

# In[]
from spatialdm.diff_utils import *
concat=concat_obj([adata_pp12, adata_ppc12, adata_pp18, adata_ppc18], ["pp12", "ppc12", "pp18", "ppc18"], 'mouse', 'z-score', fdr=False)
sns.clustermap(1-concat.uns['p_df'], yticklabels=True)

# conduct differential test across conditions, .uns["p_val"] stores the differential p-value, 
# and .uns["diff_fdr"] stores the fdr-corrected p-value
differential_test(concat, ["pp12", "ppc12", "pp18", "ppc18"], np.array([0,1,0,1]))
group_differential_pairs(concat, 'CTRL', 'CXCR7-KO')

# plot p_df, add pathway names on p_df
pathway = concat.uns["geneInter"].loc[concat.uns["p_df"].index.values, "pathway_name"].values
concat.uns["p_df"].index = [x + " (" + y + ")" for x,y in zip(concat.uns["p_df"].index.values, pathway)]
pathway = concat.uns["geneInter"].loc[concat.uns["tf_df"].index.values, "pathway_name"].values
concat.uns["tf_df"].index = [x + " (" + y + ")" for x,y in zip(concat.uns["tf_df"].index.values, pathway)]

def differential_dendrogram(sample):
    _range = np.arange(1, sample.uns['n_sub'])
    ax = sns.clustermap(1-sample.uns['p_df'].loc[(sample.uns['p_val']<0.1) & (sample.uns['tf_df'].sum(1).isin(_range)),
                                     sample.uns['subset']], figsize = (10,55), yticklabels = True)
    return ax



ax = differential_dendrogram(concat)
ax.figure.savefig("results_spatialdm_visium/compare_coexpr.pdf", bbox_inches = "tight")
pl.differential_volcano(concat, legend=['CTRL', 'CXCR7-KO'])


# %%
