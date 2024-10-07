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

# In[]
# -----------------------------------------------------------------------------------------
#
# load reference and query dataset
# 1. Reference: scRNA-seq dataset public ver
# 2. Query: self-generated scRNA-seq dataset
#
# NOTE: For label transfer, we use annot instead of annot.sub, annot.sub provide detailed 
# annotation of immune subtypes (including Lymphoid and Macrophage), but the accuracy is low
# For annot, the high-level immune cell annotation is highly accurate 
# (Lymphoid, Macrophage, myeloid are separate clusters)
#
# -----------------------------------------------------------------------------------------

# load the reference dataset
# NOTE: the public dataset is already qc-ed
adata_ref = sc.read_h5ad("../dataset/reference_scrnaseq/reference_raw.h5ad")
# update the batch and annotation labels for LT
adata_ref.obs["sample"] = adata_ref.obs["batchID"].values
print("number of batches in reference: {:d}".format(len(np.unique(adata_ref.obs["batchID"]))))
for batch in np.unique(adata_ref.obs["batchID"]):
    print(np.sum(adata_ref.obs["batchID"] == batch))

# fulltypemerged separate immune cells compared to the IntType
adata_ref.obs["annot"] = adata_ref.obs["FullTypeMerged"].values.astype(object)
# remove the predicted doublets
adata_ref = adata_ref[~adata_ref.obs["annot"].isin(["PredDoublet_Epi_Imm", "PredDoublet_Str_Epi", "PredDoublet_Str_Imm"]),:]
# process the annotations
combined_annos = [x + "-" + y for x,y in zip(adata_ref.obs["IntType"].values, adata_ref.obs["FullTypeMerged"].values)]
# for x in np.unique(combined_annos):
#     print(x)
# combine the SV cells, not of interest
adata_ref.obs.loc[adata_ref.obs["annot"].isin(["Epi_SV_Basal", "Epi_SV_Ionocyte", "Epi_SV_Luminal"]), "annot"] = "SV"
# combine endothelium
adata_ref.obs.loc[adata_ref.obs["annot"].isin(["Str_Endothelium_Lymphatic", "Str_Endothelium_Vascular"]), "annot"] = "Endothelial"
# combine luminal
adata_ref.obs.loc[adata_ref.obs["annot"].isin(["Epi_Luminal", "Epi_Luminal_2Psca", "Epi_Luminal_3Foxi1"]), "annot"] = "Luminal"
# change name
adata_ref.obs.loc[adata_ref.obs["annot"].isin(["Str_Mesenchymal"]), "annot"] = "Mesenchymal"
adata_ref.obs.loc[adata_ref.obs["annot"].isin(["Epi_Basal"]), "annot"] = "Basal"
adata_ref.obs["annot"] = adata_ref.obs["annot"].astype("category")
# NOTE: keep the original immune name, unless when combining with ref1
# adata_ref.obs.loc[adata_ref.obs["annot"].isin(["Imm_B", "Imm_NK", "Imm_Tcell"])] = "Lymphoid"
# adata_ref.obs.loc[adata_ref.obs["annot"].isin(["Imm_Macrophage", "Immune-Imm_DC"])] = "Macrophage"
print(np.unique(adata_ref.obs["annot"], return_counts = True))


# load query scRNA-seq dataset, NOTE: use raw data, include 4 samples/batches
adata_query =  sc.read_h5ad("../dataset/data_scrnaseq/data/qc_data/adata_qc_intact.h5ad")
# load the original annotation for comparison purpose
adata_ref.obs["annot_orig"] = adata_ref.obs["annot"].values
adata_query.obs["annot_orig"] = sc.read_h5ad("../dataset/data_scrnaseq/seurat_integration/adata_intact_seurat.h5ad").obs.loc[adata_query.obs.index,"annot"].values

# find overlapping features
gene_overlap = np.intersect1d(adata_query.var.index.values, adata_ref.var.index.values)
adata_ref = adata_ref[:, gene_overlap]
adata_query = adata_query[:, gene_overlap]
adata_ref.layers["counts"] = adata_ref.X.copy()
adata_query.layers["counts"] = adata_query.X.copy()

# result directory
# number of samples, maybe try different values for ref2, the population of ref2 is large, could run 500
nsamples_scanvi = 300
res_dir = f"../results_scrnaseq/annot_scanvi/annot_{nsamples_scanvi}/"
if not os.path.exists(res_dir):
    os.makedirs(res_dir)
    os.makedirs(res_dir + "plots_marker/")


# In[] 
# -----------------------------------------------------------------------------------------
#
# Label transfer -- scANVI
# obtain label annotation
#
# -----------------------------------------------------------------------------------------
# set the seed before running the model, for reproductivity
scvi.settings.seed = 0

# NOTE: Step 1, Train scVI on reference scRNA-seq dataset
scvi.model.SCVI.setup_anndata(adata_ref, batch_key = "sample", layer = "counts", labels_key = "annot")
scvi_ref = scvi.model.SCVI(adata_ref, use_layer_norm = "both", use_batch_norm = "none", 
                           encode_covariates = True, dropout_rate = 0.2, n_layers = 2)
scvi_ref.train()
# Save the trained model on the reference scRNA-seq data
scvi_ref.save(res_dir + "scvi_reference", overwrite=True)

# # Sanity check, visualize the trained reference result (scVI)
# SCVI_LATENT_KEY = "X_scVI"
# adata_ref.obsm[SCVI_LATENT_KEY] = scvi_ref.get_latent_representation()
# sc.pp.neighbors(adata_ref, use_rep = SCVI_LATENT_KEY)
# sc.tl.umap(adata_ref)
# sc.pl.umap(adata_ref, color=["sample", "annot"], frameon=False, ncols=1)
# utils.plot_latent(adata_ref.obsm["X_umap"], annos = np.array([x for x in adata_ref.obs["annot"].values]), batches = np.array([x for x in adata_ref.obs["sample"].values]), mode = "separate", figsize = (10, 25))

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

# # Sanity check, visualize the trained reference result (scANVI)
# adata_ref.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation()
# sc.pp.neighbors(adata_ref, use_rep=SCANVI_LATENT_KEY)
# sc.tl.umap(adata_ref)
# sc.pl.umap(adata_ref, color=["sample", "annot"], frameon=False, ncols=1)
# utils.plot_latent(adata_ref.obsm["X_umap"], annos = np.array([x for x in adata_ref.obs["annot"].values]), batches = np.array([x for x in adata_ref.obs["sample"].values]), mode = "separate", figsize = (10, 25))

# query annotation
# scanvi_model = scvi.model.SCANVI.load(res_dir + "scanvi_reference/", adata = adata_ref)
scvi.model.SCANVI.prepare_query_anndata(adata_query, scanvi_model)
scanvi_query = scvi.model.SCANVI.load_query_data(adata_query, scanvi_model)
scanvi_query.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0}, check_val_every_n_epoch=10,)
SCANVI_PREDICTIONS_KEY = "predictions_scanvi"

adata_query.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation()
adata_query.obs[SCANVI_PREDICTIONS_KEY] = scanvi_query.predict()

scanvi_query.save(res_dir + "scanvi_query", overwrite = True)
adata_query.write_h5ad(res_dir + "scanvi_query/adata_query.h5ad")

# In[]
# -----------------------------------------------------------------------------------------
#
# Downstream analysis: Load annotation result and marker information for downstream analysis
#
# -----------------------------------------------------------------------------------------
adata_query = sc.read_h5ad(res_dir + "scanvi_query/adata_query.h5ad")
SCANVI_LATENT_KEY = "X_scANVI"
SCANVI_PREDICTIONS_KEY = "predictions_scanvi"

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

# macrophage sub markers
markers["macrophage_m1"] = ["Ptprc", "Itgam", "Adgre1", "Cd68", "Nos2"]
markers["macrophage_m2"] = ["Ptprc", "Itgam", "Adgre1", "Mrc1", "Arg1"]

# useless, just for sanity check
# markers["sv"] = ["Pax2"]


# In[]
# -----------------------------------------------------------------------------------------
#
# NOTE: Downstream 1: Visualize the predicted labels
#
# -----------------------------------------------------------------------------------------
plt.ioff()
# NOTE: Sanity check, visualize using joint scANVI embedding, see if the query dataset distribution truly is aligned to reference
# scANVI can perform bad when there is mismatch
adata_full = anndata.concat([adata_ref, adata_query], join = "inner", label = "modality", keys = ["Ref", "Query"])
# use the ground truth annotation for the reference dataset and the predicted annotation for query dataset
adata_full.obs.loc[adata_query.obs.index,"annot"] = [x for x in adata_query.obs[SCANVI_PREDICTIONS_KEY].values]
# for comparison purpose, load the original annotation of query
# load scanvi model and obtain embedding
scanvi_query = scvi.model.SCANVI.load(res_dir + "scanvi_query/", adata = adata_full)
adata_full.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation(adata_full)
sc.pp.neighbors(adata_full, use_rep=SCANVI_LATENT_KEY)
sc.tl.umap(adata_full)


adata_full.obs["annot"] = adata_full.obs["annot"].astype(object)
adata_full.obs["annot_orig"] = adata_full.obs["annot_orig"].astype(object)

# Visualization using predicted annotation
# NOTE: X_umap is calculated based on scANVI embedding
fig = utils.plot_by_batch(adata_full.obsm["X_umap"], annos = np.array([x for x in adata_full.obs["annot"].values]),
                          batches = np.array([x for x in adata_full.obs["modality"].values]), figsize = (12, 6),
                          ncols = 2, s = 2, alpha = 0.7, markerscale = 10)
fig.savefig(res_dir + "scanvi_embed.png", bbox_inches = "tight", dpi = 150)

# for validation of annotation
fig = utils.plot_by_batch(adata_full.obsm["X_umap"], annos = np.array([x for x in adata_full.obs["annot_orig"].values]),
                          batches = np.array([x for x in adata_full.obs["modality"].values]), figsize = (15, 6),
                          ncols = 2, s = 2, alpha = 0.7, markerscale = 10)
fig.savefig(res_dir + "scanvi_embed_orig.png", bbox_inches = "tight", dpi = 150)

# TODO: confusion matrix?
from sklearn.metrics import confusion_matrix
import seaborn as sns
from matplotlib.colors import LogNorm
orig_annos = np.array([x for x in adata_full.obs.loc[adata_query.obs.index, "annot_orig"].values])
lt_annos = np.array([x for x in adata_full.obs.loc[adata_query.obs.index, "annot"].values])
lt_annos[(lt_annos == "Imm_B") | (lt_annos == "Imm_NK") | (lt_annos == "Imm_Tcell")] = "Lymphoid"
lt_annos[(lt_annos == "Imm_Macrophage")] = "Macrophage"
lt_annos[(lt_annos == "Imm_DC")] = "DC"
orig_annos[(orig_annos == "Luminal") | (orig_annos == "Luminal (Spink1+)")] = "Luminal"

# adata_ref.obs.loc[adata_ref.obs["annot"].isin(["Imm_Macrophage", "Immune-Imm_DC"])] = "Macrophage"
conf_mtx = confusion_matrix(y_true = orig_annos, y_pred = lt_annos, labels = np.unique(np.concatenate([orig_annos, lt_annos])) )
conf_mtx = pd.DataFrame(data = conf_mtx.astype(int), index = np.unique(np.concatenate([orig_annos, lt_annos])), columns = np.unique(np.concatenate([orig_annos, lt_annos])))

fig = plt.figure(figsize = (18, 16))
ax = fig.add_subplot()
sns.heatmap(conf_mtx + 1e-4, norm=LogNorm(), annot = True, ax = ax)
ax.set_xlabel("LT annot (scANVI)", fontsize = 20)
ax.set_ylabel("Orig annot", fontsize = 20)
ax.set_title("Confusion matrix", fontsize = 25)
fig.savefig(res_dir + "confusion_matrix.png", bbox_inches = "tight", dpi = 150)
# Three issue, according to the confusion matrix
# Original annotation miss assign DC into macrophage
# Label transfer miss assign monocytes as macrophage
# Original annotation Luminal is assigned as partially Basal and SV in LT
# Original Club epithelial is assigned as Luminal in LT
# Original annotation Basal is assigned as partially Luminal in LT
# The two annotations are highly consistent according to the remaining cells


# DC marker gene: CD11c (Itgax)
# for validation purpose
markers_dc = ["Itgax", "Ly75", "Cd14"]
marker_expr = pd.DataFrame(adata_full[adata_full.obs["modality"] == "Query", markers_dc].X.toarray(), columns = markers_dc, index = adata_full[adata_full.obs["modality"] == "Query"].obs.index.values)
fig = utils.plot_embeds_continuous(adata_full[adata_full.obs["modality"] == "Query"].obsm["X_umap"], annos = np.log1p(marker_expr), figsize = (10, 6), alpha = 0.7, markerscale = 10)
fig.savefig(res_dir + "scanvi_embed_dc_marker.png", bbox_inches = "tight", dpi = 150)
# TODO: boxplot

# update the original adata 
adata_query_annot = sc.read_h5ad("../dataset/data_scrnaseq/seurat_integration/adata_intact_seurat.h5ad")
# dc index
dc_idx = adata_full[(adata_full.obs["annot_orig"] == "Macrophage") & (adata_full.obs["annot"] == "Imm_DC") & (adata_full.obs["modality"] == "Query")].obs.index.values
adata_query_annot.obs["annot (lt_correct)"] = adata_query_annot.obs["annot"].astype(object)
adata_query_annot.obs.loc[dc_idx, "annot (lt_correct)"] = "DC"
adata_query_annot.obs["annot (lt_correct)"] = adata_query_annot.obs["annot (lt_correct)"].astype("category")

adata_full.obs.loc[dc_idx,"annot_orig"] = "DC"
fig = utils.plot_by_batch(adata_full.obsm["X_umap"], annos = np.array([x for x in adata_full.obs["annot_orig"].values]),
                          batches = np.array([x for x in adata_full.obs["modality"].values]), figsize = (15, 6),
                          ncols = 2, s = 2, alpha = 0.7, markerscale = 10)
fig.savefig(res_dir + "scanvi_embed_orig_correct.png", bbox_inches = "tight", dpi = 150)
adata_query_annot.write_h5ad("../dataset/data_scrnaseq/seurat_integration/adata_intact_seurat_lt_correct.h5ad")



# In[]
# -----------------------------------------------------------------------------------------
#
# NOTE: Downtream 3: Analysis on cell type composition
#
# -----------------------------------------------------------------------------------------
adata_intact = sc.read_h5ad("../dataset/data_scrnaseq/seurat_integration/adata_intact_seurat_lt_correct.h5ad")
adata_intact.obs["annot (lt_correct)"] = adata_intact.obs["annot (lt_correct)"].astype("str")
# immune cells count under 12wk
adata_ctr_12wk = adata_intact[(adata_intact.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata_intact.obs["age"] == "12wk"),:]
adata_ko_12wk = adata_intact[(adata_intact.obs["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata_intact.obs["age"] == "12wk"),:]
# immune cells count under 18wk
adata_ctr_18wk = adata_intact[(adata_intact.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata_intact.obs["age"] == "18wk"),:]
adata_ko_18wk = adata_intact[(adata_intact.obs["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata_intact.obs["age"] == "18wk"),:]


#. all cell type percentage
unique_ct_ctrl_12wk, counts_ct_ctrl_12wk = np.unique(adata_ctr_12wk.obs["annot (lt_correct)"].values, return_counts = True)
unique_ct_ctrl_18wk, counts_ct_ctrl_18wk = np.unique(adata_ctr_18wk.obs["annot (lt_correct)"].values, return_counts = True)
unique_ct_ko_12wk, counts_ct_ko_12wk = np.unique(adata_ko_12wk.obs["annot (lt_correct)"].values, return_counts = True)
unique_ct_ko_18wk, counts_ct_ko_18wk = np.unique(adata_ko_18wk.obs["annot (lt_correct)"].values, return_counts = True)

unique_ct = np.unique(np.concatenate([unique_ct_ctrl_12wk, unique_ct_ctrl_18wk, unique_ct_ko_12wk, unique_ct_ko_18wk]))
counts_df = pd.DataFrame(data = 0, columns = ["CTRL (12wk)", "CXCR7 KO (12wk)", "CTRL (18wk)", "CXCR7 KO (18wk)"], index = unique_ct)
counts_df.loc[unique_ct_ctrl_12wk, "CTRL (12wk)"] = counts_ct_ctrl_12wk
counts_df.loc[unique_ct_ko_12wk, "CXCR7 KO (12wk)"] = counts_ct_ko_12wk
counts_df.loc[unique_ct_ctrl_18wk, "CTRL (18wk)"] = counts_ct_ctrl_18wk
counts_df.loc[unique_ct_ko_18wk, "CXCR7 KO (18wk)"] = counts_ct_ko_18wk
pct_df = counts_df.copy()
pct_df.loc[:,:] = counts_df.values/np.sum(counts_df.values, axis = 0, keepdims = True) * 100

summary_df = pd.DataFrame(columns = ["ct", "condition", "count", "pct"])
summary_df["ct"] = np.concatenate([unique_ct] * 4)
summary_df["condition"] = np.array(["CTRL (12wk)"] * unique_ct.shape[0] + ["CXCR7 KO (12wk)"] * unique_ct.shape[0]
                          + ["CTRL (18wk)"] * unique_ct.shape[0] + ["CXCR7 KO (18wk)"] * unique_ct.shape[0])
summary_df["count"] = np.concatenate([counts_df["CTRL (12wk)"].values, counts_df["CXCR7 KO (12wk)"].values, 
                            counts_df["CTRL (18wk)"].values, counts_df["CXCR7 KO (18wk)"].values])
summary_df["pct"] = np.concatenate([pct_df["CTRL (12wk)"].values, pct_df["CXCR7 KO (12wk)"].values, 
                            pct_df["CTRL (18wk)"].values, pct_df["CXCR7 KO (18wk)"].values])

sns.set_theme()
fig = plt.figure(figsize = (17, 7))
ax = fig.add_subplot()
sns.barplot(summary_df, x = "ct", hue = "condition", y = "pct", ax = ax)
fig.savefig(res_dir + "cell_type_proportion.png")
# %%
