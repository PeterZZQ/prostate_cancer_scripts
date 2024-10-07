# In[]
import pandas as pd
import numpy as np 
import scanpy as sc
import anndata
import scvi
import sys, os
sys.path.append("../")
import utils
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from sklearn.metrics import confusion_matrix
import seaborn as sns
from matplotlib.colors import LogNorm
# In[]
# -----------------------------------------------------------------------------------------
# NOTE: SECTION 1: Loading dataset
# load reference and query dataset
# 1. Reference: scRNA-seq dataset, include 4 samples/batches.
# 2. Query: merfish dataset, include 2 samples/batches: Merfish-PP15, Merfish-PPC15
#
# NOTE: For label transfer, we use annot instead of annot.sub, annot.sub provide detailed 
# annotation of immune subtypes (including Lymphoid and Macrophage), but the accuracy is low
# For annot, the high-level immune cell annotation is highly accurate 
# (Lymphoid, Macrophage, myeloid are separate clusters)
#
# ----------------------------------------------------------------------------------------
use_proseg = True
dataset = "merfish0626"

# first batch of merfish dataset
if dataset == "merfish1":
    # load merfish query dataset, region 0: PP15, region 1: PPC15
    # 1. load region 0
    if use_proseg:
        adata_merfish_pp = sc.read_h5ad("../dataset/data_merfish/region_0/proseg_results/adata_qc.h5ad")
    else:
        adata_merfish_pp = sc.read_h5ad("../dataset/data_merfish/region_0/adata_qc.h5ad")
    adata_merfish_pp.obs["annot"] = pd.Categorical(["Unknown"] * adata_merfish_pp.shape[0], categories = ["Unknown"])
    adata_merfish_pp.obs["genotype"] = pd.Categorical(["PbCre(+/-),Pten(-/-),P53(-/-)"] * adata_merfish_pp.shape[0], categories = ["PbCre(+/-),Pten(-/-),P53(-/-)", "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"])
    adata_merfish_pp.obs["sample"] = pd.Categorical([f"{dataset}-PP"] * adata_merfish_pp.shape[0], categories = [f"{dataset}-PP", f"{dataset}-PPC"])
    # 2. load region 1
    if use_proseg:
        adata_merfish_ppc = sc.read_h5ad("../dataset/data_merfish/region_1/proseg_results/adata_qc.h5ad")
    else:
        adata_merfish_ppc = sc.read_h5ad("../dataset/data_merfish/region_1/adata_qc.h5ad")
    adata_merfish_ppc.obs["annot"] = pd.Categorical(["Unknown"] * adata_merfish_ppc.shape[0], categories = ["Unknown"])
    adata_merfish_ppc.obs["genotype"] = pd.Categorical(["PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"] * adata_merfish_ppc.shape[0], categories = ["PbCre(+/-),Pten(-/-),P53(-/-)", "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"])
    adata_merfish_ppc.obs["sample"] = pd.Categorical([f"{dataset}-PPC"] * adata_merfish_ppc.shape[0], categories = [f"{dataset}-PP", f"{dataset}-PPC"])

else:
    # merfish0617 query dataset, region 0: PPC15, region 1: PP15
    # merfish0626 query dataset, region 0: PPC19, region 1: PP19
    # 1. load region 0
    if use_proseg:
        adata_merfish_ppc = sc.read_h5ad(f"../dataset/{dataset}/region_0/proseg_results/adata_qc.h5ad")
    else:
        adata_merfish_ppc = sc.read_h5ad(f"../dataset/{dataset}/region_0/adata_qc.h5ad")
    adata_merfish_ppc.obs["annot"] = pd.Categorical(["Unknown"] * adata_merfish_ppc.shape[0], categories = ["Unknown"])
    adata_merfish_ppc.obs["genotype"] = pd.Categorical(["PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"] * adata_merfish_ppc.shape[0], categories = ["PbCre(+/-),Pten(-/-),P53(-/-)", "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"])
    adata_merfish_ppc.obs["sample"] = pd.Categorical([f"{dataset}-PPC"] * adata_merfish_ppc.shape[0], categories = [f"{dataset}-PP", f"{dataset}-PPC"])
    # 2. load region 1
    if use_proseg:
        if dataset == "merfish0626":
            adata_merfish_pp = sc.read_h5ad(f"../dataset/{dataset}/region_1/proseg_results/adata_qc_f.h5ad")
        else:
            adata_merfish_pp = sc.read_h5ad(f"../dataset/{dataset}/region_1/proseg_results/adata_qc.h5ad")
    else:
        adata_merfish_pp = sc.read_h5ad(f"../dataset/{dataset}/region_1/adata_qc.h5ad")
    adata_merfish_pp.obs["annot"] = pd.Categorical(["Unknown"] * adata_merfish_pp.shape[0], categories = ["Unknown"])
    adata_merfish_pp.obs["genotype"] = pd.Categorical(["PbCre(+/-),Pten(-/-),P53(-/-)"] * adata_merfish_pp.shape[0], categories = ["PbCre(+/-),Pten(-/-),P53(-/-)", "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"])
    adata_merfish_pp.obs["sample"] = pd.Categorical([f"{dataset}-PP"] * adata_merfish_pp.shape[0], categories = [f"{dataset}-PP", f"{dataset}-PPC"])

# NOTE: load reference annotation dataset, ref1 is the generated dataset, ref2 is the public available dataset
ref = "combined"

# NOTE: reference 1, self-generated
# load reference scRNA-seq dataset, NOTE: use raw data, include 4 samples/batches
adata_ref1 =  sc.read_h5ad("../dataset/data_scrnaseq/data/qc_data/adata_qc_intact.h5ad")
# The smallest cell cluster in reference: Lymphoid, include 349 cells
# lt_correct is the corrected version of the label 
adata_ref1.obs["annot"] = sc.read_h5ad("../dataset/data_scrnaseq/seurat_integration/adata_intact_seurat_lt_correct.h5ad").obs.loc[adata_ref1.obs.index,"annot (lt_correct)"].values

# NOTE: reference 2: public healthy data
# NOTE: the public dataset is already qc-ed
adata_ref2 = sc.read_h5ad("../dataset/reference_scrnaseq/reference_raw.h5ad")
# update the batch and annotation labels for LT
adata_ref2.obs["sample"] = adata_ref2.obs["batchID"].values
print("number of batches in reference: {:d}".format(len(np.unique(adata_ref2.obs["batchID"]))))
# fulltypemerged separate immune cells compared to the IntType
adata_ref2.obs["annot"] = adata_ref2.obs["FullTypeMerged"].values.astype(object)
# remove the predicted doublets
adata_ref2 = adata_ref2[~adata_ref2.obs["annot"].isin(["PredDoublet_Epi_Imm", "PredDoublet_Str_Epi", "PredDoublet_Str_Imm"]),:]

# process the annotations
combined_annos = [x + "-" + y for x,y in zip(adata_ref2.obs["IntType"].values, adata_ref2.obs["FullTypeMerged"].values)]
for x in np.unique(combined_annos):
    print(x)
# combine the SV cells, not of interest, better to remove SV cells when calculating the cell type proportion
adata_ref2.obs.loc[adata_ref2.obs["annot"].isin(["Epi_SV_Basal", "Epi_SV_Ionocyte", "Epi_SV_Luminal"]), "annot"] = "SV"
# combine endothelium
adata_ref2.obs.loc[adata_ref2.obs["annot"].isin(["Str_Endothelium_Lymphatic", "Str_Endothelium_Vascular"]), "annot"] = "Endothelial"
# combine luminal
adata_ref2.obs.loc[adata_ref2.obs["annot"].isin(["Epi_Luminal", "Epi_Luminal_2Psca", "Epi_Luminal_3Foxi1"]), "annot"] = "Luminal"
# combine lymphoid, NOTE: the lymphoid subtype annotation is not very accurate according to marker gene checking
adata_ref2.obs.loc[adata_ref2.obs["annot"].isin(["Imm_B", "Imm_NK", "Imm_Tcell"]), "annot"] = "Lymphoid"
# change name
adata_ref2.obs.loc[adata_ref2.obs["annot"].isin(["Str_Mesenchymal"]), "annot"] = "Mesenchymal"
adata_ref2.obs.loc[adata_ref2.obs["annot"].isin(["Epi_Basal"]), "annot"] = "Basal"
adata_ref2.obs.loc[adata_ref2.obs["annot"].isin(["Imm_Macrophage"]), "annot"] = "Macrophage"
adata_ref2.obs.loc[adata_ref2.obs["annot"].isin(["Imm_DC"]), "annot"] = "DC"
adata_ref2.obs["annot"] = adata_ref2.obs["annot"].astype("category")
print(np.unique(adata_ref2.obs["annot"], return_counts = True))


if ref == "ref1":
    adata_ref = adata_ref1.copy()
    adata_ref.obs["modality"] = "Ref"
elif ref == "ref2":
    adata_ref = adata_ref2.copy()
    adata_ref.obs["modality"] = "Ref"
elif ref == "combined":    
    # combined luminal spink1+ into luminal, club epithelial is a separate cluster, not combined
    adata_ref1.obs["annot"] = adata_ref1.obs["annot"].values.astype(object)
    adata_ref1.obs.loc[adata_ref1.obs["annot"].isin(["Luminal (Spink1+)"]), "annot"] = "Luminal"
    adata_ref = anndata.concat([adata_ref1, adata_ref2], axis = 0, join = "inner", label = "modality", keys = ["Ref1", "Ref2"])
    adata_ref.obs["annot"] = adata_ref.obs["annot"].astype("category")
# find overlapping features
gene_overlap = np.intersect1d(adata_merfish_pp.var.index.values, adata_ref.var.index.values)
adata_ref = adata_ref[:, gene_overlap]
adata_merfish_pp = adata_merfish_pp[:, gene_overlap]
adata_merfish_ppc = adata_merfish_ppc[:, gene_overlap]
adata_ref.layers["counts"] = adata_ref.X.copy()
adata_merfish_pp.layers["counts"] = adata_merfish_pp.X.copy()
adata_merfish_ppc.layers["counts"] = adata_merfish_ppc.X.copy()

# result directory
# number of samples, maybe try different values for ref2, the population of ref2 is large, could run 500
nsamples_scanvi = 300
if use_proseg:
    if ref == "ref1":
        res_dir = f"../results_merfish/results_{dataset}_proseg/annot_{nsamples_scanvi}_ref1/"
    elif ref == "ref2":
        res_dir = f"../results_merfish/results_{dataset}_proseg/annot_{nsamples_scanvi}_ref2/"
    elif ref == "combined":
        res_dir = f"../results_merfish/results_{dataset}_proseg/annot_{nsamples_scanvi}_combined/"
else:
    if ref == "ref1":
        res_dir = f"../results_merfish/results_{dataset}/annot_{nsamples_scanvi}_ref1/"
    elif ref == "ref2":
        res_dir = f"../results_merfish/results_{dataset}/annot_{nsamples_scanvi}_ref2/"
    elif ref == "combined":
        res_dir = f"../results_merfish/results_{dataset}/annot_{nsamples_scanvi}_combined/"
if not os.path.exists(res_dir):
    os.makedirs(res_dir + "plots_marker/")

# In[] 
# -----------------------------------------------------------------------------------------
# NOTE: SECTION 2: Running scANVI
# Label transfer -- scANVI
# obtain label annotation
#
# -----------------------------------------------------------------------------------------
# set the seed before running the model, for reproductivity
scvi.settings.seed = 0

# NOTE: Step 1, Train scVI on reference scRNA-seq dataset
# Train scVI on the reference dataset 
scvi.model.SCVI.setup_anndata(adata_ref, batch_key = "sample", layer = "counts", labels_key = "annot")
scvi_ref = scvi.model.SCVI(adata_ref, use_layer_norm = "both", use_batch_norm = "none", 
                           encode_covariates = True, dropout_rate = 0.2, n_layers = 2)
scvi_ref.train()
# Save the trained model on the reference scRNA-seq data
scvi_ref.save(res_dir + "scvi_reference", overwrite=True)

# NOTE: Step 2, Train scANVI based on scVI model
# Reference
# use scvi initialized with the reference dataset for better result
scvi_model = scvi.model.SCVI.load(res_dir + "scvi_reference", adata = adata_ref)
# use the same label, layer, batch_label as in scvi
scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_model, labels_key = "annot", unlabeled_category="Unknown")
# default for n_samples_per_label is None, here we use the example value 100
scanvi_model.train(max_epochs = 20, n_samples_per_label = nsamples_scanvi)
scanvi_model.save(res_dir + "scanvi_reference", overwrite=True)
SCANVI_LATENT_KEY = "X_scANVI"

# Joint training: transfer labels jointly to pp and ppc
# combine pp and ppc, no need for label and keys
adata_merfish = anndata.concat([adata_merfish_pp, adata_merfish_ppc], join = "inner")
# Query classification for combined merfish data
scvi.model.SCANVI.prepare_query_anndata(adata_merfish, scanvi_model)
scanvi_query = scvi.model.SCANVI.load_query_data(adata_merfish, scanvi_model)
scanvi_query.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0}, check_val_every_n_epoch=10,)
SCANVI_PREDICTIONS_KEY = "predictions_scanvi"
adata_merfish.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation()
adata_merfish.obs[SCANVI_PREDICTIONS_KEY] = scanvi_query.predict()
scanvi_query.save(res_dir + "scanvi", overwrite = True)
adata_merfish.write_h5ad(res_dir + "scanvi/adata_merfish.h5ad")

# In[]
# -----------------------------------------------------------------------------------------
#
# NOTE: SECTION 3: Downstream analysis, Load annotation result and marker information for downstream analysis
#
# -----------------------------------------------------------------------------------------
adata_merfish = sc.read_h5ad(res_dir + "scanvi/adata_merfish.h5ad")
SCANVI_LATENT_KEY = "X_scANVI"
SCANVI_PREDICTIONS_KEY = "predictions_scanvi"
if use_proseg:
    adata_merfish.obsm["spatial"] = adata_merfish.obsm["X_spatial"]
# separate
adata_merfish_pp = adata_merfish[adata_merfish.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"]
adata_merfish_ppc = adata_merfish[adata_merfish.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"]

# Read in the gene panel information
gene_panel_info = pd.read_csv("../dataset/data_merfish/MERSCOPE_gene_panel_info.csv", index_col = 0)
markers = {}
# picked only the one that are included in spatial data
markers["luminal"] = ["Ar", "Krt8"]
markers["basal"] = ["Trp63", "Krt5"]
markers["club_epithelia"] = ["Agr2", "Krt7"]
markers["endothelial"] = ["Ackr1", "Cldn5"]
markers["mesenchymal"] = ["Fgf10", "Wnt10a", "Wnt2"]
# immune markers
markers["lymphoid"] = ["Cd3e", "Ms4a1", "Klrb1c"]
# markers["myeloid"] = ["Ptprc", "Itgam"]
markers["monocytes"] = ["Ptprc", "Itgam", "Cd14", "S100a8", "S100a9"] 
markers["macrophage"] = ["Ptprc", "Itgam", "Adgre1"]

# macrophage sub markers in MerFISH dataset
# Il6: M1, inflammation
# Cd68, Il23a, Cd80, Cd86, Nos2, Cxcl9, Cxcl10: M1
# Ccl22, Ccl24: M2 cytokines
# Msr1: M2 TAM
# Cd163: M2 macrophage
# Mrc1: M2 macrophage, Reg-TAM
markers["macrophage_m1"] = ["Il6", "Il23a", "Cd80", "Cd86", "Cxcl9", "Cxcl10", "Cd68", "Nos2"]
markers["macrophage_m2"] = ["Ccl22", "Ccl24", "Msr1", "Cd163", "Mrc1"]
# useless, just for sanity check
# markers["sv"] = ["Pax2"]

# NOTE: Sanity check, umap of data with predict label
# # log-normalize the cells
# adata_merfish = anndata.concat([adata_merfish_pp, adata_merfish_ppc])
# sc.pp.normalize_total(adata_merfish, 1e4)
# sc.pp.log1p(adata_merfish)
# # calculate the umap embedding for the query dataset
# # sc.tl.pca(adata_merfish, n_comps = 100)
# sc.pp.neighbors(adata_merfish, n_neighbors = 50)
# sc.tl.umap(adata_merfish, min_dist = 0.4)
# Visualize using UMAP embedding, batch effect is not strong here, not useful
# fig = utils.plot_by_batch(adata_merfish.obsm["X_umap"], annos = np.array([x for x in adata_merfish.obs[SCANVI_PREDICTIONS_KEY].values]), 
#                           batches = np.array([x for x in adata_merfish.obs["sample"].values]), figsize = (12, 6), ncols = 2, s = 1, alpha = 0.3, markerscale = 10)
# fig.savefig(res_dir + "umap_scanvi_merfish.png", bbox_inches = "tight", dpi = 150)
# fig = utils.plot_embeds(adata_merfish.obsm["X_umap"], annos = ct_df, figsize = (12, 6), s = 1, alpha = 0.7, markerscale = 10, colormap = cmap2colors, ncols = 3)
# fig.savefig(res_dir + "umap_scanvi_merfish_ct.png", bbox_inches = "tight", dpi = 150)

# In[]
# -----------------------------------------------------------------------------------------
#
# Visualize the integrated latent embedding
#
# -----------------------------------------------------------------------------------------
plt.ioff()
# NOTE: Sanity check, visualize using joint scANVI embedding, see if the query dataset distribution truly is aligned to reference
# scANVI can perform bad when there is mismatch
# JOINT
adata_merfish_pp.obs["modality"] = "Merfish-pp"
adata_merfish_ppc.obs["modality"] = "Merfish-ppc"
adata_ref_merfish = anndata.concat([adata_ref, adata_merfish_pp, adata_merfish_ppc], join = "inner")
adata_ref_merfish.obs.loc[adata_merfish.obs.index,"annot"] = [x for x in adata_merfish.obs[SCANVI_PREDICTIONS_KEY].values]
scanvi_query = scvi.model.SCANVI.load(res_dir + "scanvi/", adata = adata_ref_merfish)
# obtain the integrated latent embedding of ref, pp and ppc
adata_ref_merfish.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation(adata_ref_merfish)
sc.pp.neighbors(adata_ref_merfish, use_rep=SCANVI_LATENT_KEY)
sc.tl.umap(adata_ref_merfish)

fig = utils.plot_by_batch(adata_ref_merfish.obsm["X_umap"], annos = np.array([x for x in adata_ref_merfish.obs["annot"].values]),
                          batches = np.array([x for x in adata_ref_merfish.obs["modality"].values]), figsize = (12, 6),
                          ncols = 2, s = 1, alpha = 0.3, markerscale = 10)
fig.savefig(res_dir + "scanvi_embed.png", bbox_inches = "tight", dpi = 150)

# update the umap of adata_merfish
# adata_merfish.obsm["X_scanvi_umap"] = np.concatenate([adata_ref_pp[adata_ref_pp.obs["modality"] == "Merfish"].obsm["X_umap"], 
#                                                       adata_ref_ppc[adata_ref_ppc.obs["modality"] == "Merfish"].obsm["X_umap"]], axis = 0)
# Joint, update the umap of merfish data with integrated space
adata_merfish.obsm["X_scanvi_umap"] = adata_ref_merfish[adata_merfish.obs.index.values,:].obsm["X_umap"]


# In[]
# -----------------------------------------------------------------------------------------
#
# [Only for combined] Combine annotations from different sources
# NOTE: the purpose of this step is to find the consensus of annotations across different resources
#
# -----------------------------------------------------------------------------------------
# NOTE: [Cannot run for multiple times!!] need to make sure res_dir points to combined
if ref == "combined":
    adata_merfish_ref1 = sc.read_h5ad(f"../results_merfish/results_{dataset}_proseg/annot_{nsamples_scanvi}_ref1/scanvi/adata_merfish.h5ad")
    # TODO: probably need to retrain the model since lymphoid subtype annotation is not accurate in pub dataset.
    adata_merfish_ref2 = sc.read_h5ad(f"../results_merfish/results_{dataset}_proseg/annot_{nsamples_scanvi}_ref2/scanvi/adata_merfish.h5ad")
    adata_merfish_combined = sc.read_h5ad(f"../results_merfish/results_{dataset}_proseg/annot_{nsamples_scanvi}_combined/scanvi/adata_merfish.h5ad")
    adata_merfish.obs["annot_ref1"] = adata_merfish_ref1.obs[SCANVI_PREDICTIONS_KEY].astype(object)
    adata_merfish.obs["annot_ref2"] = adata_merfish_ref2.obs[SCANVI_PREDICTIONS_KEY].astype(object)

    # unify the annotation
    # adata_merfish.obs.loc[adata_merfish.obs["annot_ref2"].isin(["Imm_B", "Imm_NK", "Imm_Tcell"]), "annot_ref2"] = "Lymphoid"
    # adata_merfish.obs.loc[adata_merfish.obs["annot_ref2"].isin(["Imm_Macrophage"]), "annot_ref2"] = "Macrophage"
    # adata_merfish.obs.loc[adata_merfish.obs["annot_ref2"].isin(["Imm_DC"]), "annot_ref2"] = "DC"    
    adata_merfish.obs.loc[adata_merfish.obs["annot_ref1"].isin(["Luminal (Spink1+)"]), "annot_ref1"] = "Luminal"
    adata_merfish.obs["annot_comb"] = adata_merfish_combined.obs[SCANVI_PREDICTIONS_KEY]

    adata_merfish.obs["annot_ref1"] = adata_merfish.obs["annot_ref1"].astype("category")
    adata_merfish.obs["annot_ref2"] = adata_merfish.obs["annot_ref2"].astype("category")
    adata_merfish.obs["annot_comb"] = adata_merfish.obs["annot_comb"].astype("category")
    adata_merfish.obs["annot_ref1"] = adata_merfish.obs["annot_ref1"].cat.set_categories(adata_merfish.obs["annot_comb"].cat.categories)
    adata_merfish.obs["annot_ref2"] = adata_merfish.obs["annot_ref2"].cat.set_categories(adata_merfish.obs["annot_comb"].cat.categories)

    del adata_merfish_ref1, adata_merfish_ref2, adata_merfish_combined

    # compare annotations 
    fig = utils.plot_embeds(adata_merfish.obsm["X_scanvi_umap"], annos = adata_merfish.obs[["annot_ref1", "annot_ref2", "annot_comb"]],
                             figsize = (12, 6), ncols = 3, s = 1, alpha = 0.3, markerscale = 10)
    fig.savefig(res_dir + "compare_annot.png", bbox_inches = "tight", dpi = 150)
    
    # calculate confusion matrix
    # 1. combine with ref1
    conf_mtx1 = confusion_matrix(y_true = np.array([x for x in adata_merfish.obs["annot_comb"]]), y_pred = np.array([x for x in adata_merfish.obs["annot_ref1"]]), labels = adata_merfish.obs["annot_comb"].cat.categories)
    conf_mtx1 = pd.DataFrame(data = conf_mtx1.astype(int), index = adata_merfish.obs["annot_comb"].cat.categories, columns = adata_merfish.obs["annot_comb"].cat.categories)

    # 1. combine with ref2
    conf_mtx2 = confusion_matrix(y_true = np.array([x for x in adata_merfish.obs["annot_comb"]]), y_pred = np.array([x for x in adata_merfish.obs["annot_ref2"]]), labels = adata_merfish.obs["annot_comb"].cat.categories)
    conf_mtx2 = pd.DataFrame(data = conf_mtx2.astype(int), index = adata_merfish.obs["annot_comb"].cat.categories, columns = adata_merfish.obs["annot_comb"].cat.categories)

    fig = plt.figure(figsize = (20 * 2, 16))
    ax = fig.subplots(nrows = 1, ncols = 2)
    sns.heatmap(conf_mtx1 + 1e-4, norm=LogNorm(), annot = True, ax = ax[0])
    ax[0].set_xlabel("Ref1", fontsize = 20)
    ax[0].set_ylabel("Combined", fontsize = 20)
    ax[0].set_title("Combined & Ref1", fontsize = 25)

    sns.heatmap(conf_mtx2 + 1e-4, norm=LogNorm(), annot = True, ax = ax[1])
    ax[1].set_xlabel("Ref2", fontsize = 20)
    ax[1].set_ylabel("Combined", fontsize = 20)
    ax[1].set_title("Combined & Ref2", fontsize = 25)

    fig.savefig(res_dir + "confusion_matrix.png", bbox_inches = "tight", dpi = 150)

    # filter inconsistent annotations, except for Club epithelia, Monoctyes, Str_Glial, Str_SmoothMuscle that do not exist in all clusters 
    adata_merfish.obs["annot_comb"] = adata_merfish.obs["annot_comb"].astype(object)
    for cell_id in adata_merfish.obs.index.values:
        if not adata_merfish.obs.loc[cell_id, "annot_comb"] == adata_merfish.obs.loc[cell_id, "annot_ref1"] == adata_merfish.obs.loc[cell_id, "annot_ref2"]:
            if (adata_merfish.obs.loc[cell_id, "annot_comb"] != "Club epithelia") & (adata_merfish.obs.loc[cell_id, "annot_comb"] != "Monocytes")\
                & (adata_merfish.obs.loc[cell_id, "annot_comb"] != "Str_Glial") & (adata_merfish.obs.loc[cell_id, "annot_comb"] != "Str_SmoothMuscle"):
                adata_merfish.obs.loc[cell_id, "annot_comb"] = "Unknown"
    
    adata_merfish.obs["annot_comb"] = adata_merfish.obs["annot_comb"].astype("category") 

    # filtering and remove unknown cells
    adata_merfish = adata_merfish[adata_merfish.obs["annot_comb"] != "Unknown",:]
    fig = utils.plot_embeds(adata_merfish.obsm["X_scanvi_umap"], annos = adata_merfish.obs[["annot_comb"]],
                             figsize = (12, 6), s = 1, alpha = 0.3, markerscale = 10)
    
    fig.savefig(res_dir + "annot_consensus.png", bbox_inches = "tight", dpi = 150)
    # update the prediction key to be filtered version
    SCANVI_PREDICTIONS_KEY = "annot_comb"

    # save the filtered adata
    adata_merfish.write_h5ad(res_dir + "scanvi/adata_merfish_combfilter.h5ad")


# In[]
# -----------------------------------------------------------------------------------------
#
# NOTE: Visualize the spatial location of the predicted cell types
#
# -----------------------------------------------------------------------------------------
plt.ioff()
# cts = ["Lymphoid", "Macrophage", "Monocytes", "Luminal (Spink1+)", "Luminal", "Endothelial", "Basal", "Club epithelia", "Mesenchymal"]
adata_merfish_pp = adata_merfish[adata_merfish.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"]
adata_merfish_ppc = adata_merfish[adata_merfish.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"]

cts = [x for x in adata_merfish.obs[SCANVI_PREDICTIONS_KEY].cat.categories if x != "SV"]
ct_df = pd.DataFrame(columns = cts, index = adata_merfish.obs.index, data = "0. Other")
for ct in cts:
    ct_df.loc[adata_merfish.obs[SCANVI_PREDICTIONS_KEY] == ct, ct] = f"1. {ct}"
    ct_df[ct] = ct_df[ct].astype("category")
cmap2colors = ListedColormap(colors = ["#e0e0e0", "#68228b"])

fig = utils.plot_embeds(adata_merfish_pp.obsm["spatial"], annos = adata_merfish_pp.obs[[SCANVI_PREDICTIONS_KEY]], figsize = (25, 13), s = 1, alpha = 0.4, markerscale = 15)
fig.savefig(res_dir + "spatial_scanvi_merfish_pp.png", bbox_inches = "tight", dpi = 150)
fig = utils.plot_embeds(adata_merfish_pp.obsm["spatial"], annos = ct_df.loc[adata_merfish_pp.obs.index.values,:], figsize = (10, 6), s = 1, alpha = 0.2, markerscale = 15, colormap = cmap2colors, ncols = 3)
fig.savefig(res_dir + "spatial_scanvi_merfish_pp_ct.png", bbox_inches = "tight", dpi = 150)
fig = utils.plot_embeds(adata_merfish_ppc.obsm["spatial"], annos = adata_merfish_ppc.obs[[SCANVI_PREDICTIONS_KEY]], figsize = (25, 13), s = 1, alpha = 0.4, markerscale = 15)
fig.savefig(res_dir + "spatial_scanvi_merfish_ppc.png", bbox_inches = "tight", dpi = 150)
fig = utils.plot_embeds(adata_merfish_ppc.obsm["spatial"], annos = ct_df.loc[adata_merfish_ppc.obs.index.values,:], figsize = (10, 6), s = 1, alpha = 0.2, markerscale = 15, colormap = cmap2colors, ncols = 3)
fig.savefig(res_dir + "spatial_scanvi_merfish_ppc_ct.png", bbox_inches = "tight", dpi = 150)


# In[]
# -----------------------------------------------------------------------------------------
#
# NOTE: Plot marker gene expression
#
# -----------------------------------------------------------------------------------------
# plt.ioff()

samples = adata_merfish.obs["sample"].cat.categories

# NOTE: Sanity check: markers for the correctness of coarse level annotations
for ct, marker in markers.items():
# for marker in [["Il6", "Ifng"]]:
    # ct = "il6_ifng"    
    nrows = len(marker)
    # adata_intact is already log-normalized
    # spatial level
    gene_expr = pd.DataFrame(data = np.log1p(adata_merfish[:, marker].X.toarray()), columns = marker, index = adata_merfish.obs.index.values.squeeze())
    vmax = [np.max(gene_expr[x].squeeze()) for x in marker]
    for sample in samples:
        fig = utils.plot_embeds(adata_merfish[adata_merfish.obs["sample"] == sample, :].obsm["X_scanvi_umap"], 
                                annos = adata_merfish[adata_merfish.obs["sample"] == sample, :].obs[[SCANVI_PREDICTIONS_KEY]], figsize = (13,6), s = 2, markerscale = 10)
        fig.suptitle(sample, fontsize = 35)
        fig.savefig(res_dir + f"plots_marker/{sample}_umap.png", bbox_inches = "tight", dpi = 150)

        fig = utils.plot_embeds_continuous(adata_merfish[adata_merfish.obs["sample"] == sample, :].obsm["spatial"], 
                                                    annos = gene_expr.loc[adata_merfish.obs["sample"] == sample, :],
                                                    colormap = utils.SUPER_MAGMA, vmax = vmax, figsize = (10,6), ncols = 2, s = 2)
        fig.suptitle(sample, fontsize = 35)
        fig.savefig(res_dir + f"plots_marker/{ct}_marker_{sample}_spatial.png", bbox_inches = "tight", dpi = 150)

    # umap level
    for sample in samples:
        fig = utils.plot_embeds_continuous(adata_merfish[adata_merfish.obs["sample"] == sample, :].obsm["X_scanvi_umap"], 
                                                    annos = gene_expr.loc[adata_merfish.obs["sample"] == sample, :],
                                                    colormap = utils.SUPER_MAGMA, vmax = vmax, figsize = (10,6), ncols = 2, s = 2)
        fig.suptitle(sample, fontsize = 35)
        fig.savefig(res_dir + f"plots_marker/{ct}_marker_{sample}_umap.png", bbox_inches = "tight", dpi = 150)

# In[]
# -----------------------------------------------------------------------------------------
#
# NOTE: Analysis on cell type composition
#
# -----------------------------------------------------------------------------------------
# Count cell type composition, coarse level
# SCANVI_PREDICTIONS_KEY = "annot_comb"
annos_pp = np.array([x for x in adata_merfish_pp.obs[SCANVI_PREDICTIONS_KEY]])
annos_ppc = np.array([x for x in adata_merfish_ppc.obs[SCANVI_PREDICTIONS_KEY]])
# NOTE: remove SV cells 
annos_pp = annos_pp[annos_pp != "SV"]
annos_ppc = annos_ppc[annos_ppc != "SV"]

ct_pp, ct_count_pp = np.unique(annos_pp, return_counts = True)
ct_pct_pp = ct_count_pp/np.sum(ct_count_pp)
ct_ppc, ct_count_ppc = np.unique(annos_ppc, return_counts = True)
ct_pct_ppc = ct_count_ppc/np.sum(ct_count_ppc)

unique_ct = np.union1d(ct_pp, ct_ppc)
count_ct = pd.DataFrame(data = 0.0, index = unique_ct, columns = [f"PP ({dataset})", f"PPC ({dataset})"])
count_ct.loc[ct_pp, f"PP ({dataset})"] = ct_count_pp
count_ct.loc[ct_pp, f"PPC ({dataset})"] = ct_count_ppc

pct_ct = pd.DataFrame(data = 0.0, index = unique_ct, columns = [f"PP ({dataset})", f"PPC ({dataset})"])
pct_ct.loc[ct_pp, f"PP ({dataset})"] = ct_pct_pp
pct_ct.loc[ct_pp, f"PPC ({dataset})"] = ct_pct_ppc

plt.rcParams["font.size"] = 10
fig = plt.figure(figsize = (20,10))
axs = fig.subplots(nrows = 1, ncols = 2)
axs[0].pie(pct_ct[f"PP ({dataset})"].values.squeeze(), labels = pct_ct.index.values.squeeze(), autopct='%.1f%%', pctdistance = 0.8, labeldistance = 1.1)
axs[1].pie(pct_ct[f"PPC ({dataset})"].values.squeeze(), labels = pct_ct.index.values.squeeze(), autopct='%.1f%%', pctdistance = 0.8, labeldistance = 1.1)
axs[0].set_title(f"PP ({dataset})")
axs[1].set_title(f"PPC ({dataset})")
fig.suptitle("Cell percentage", fontsize = 30)
fig.savefig(res_dir + "percent_cell.png", bbox_inches = "tight")

# %%
