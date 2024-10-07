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
import matplotlib

from sklearn.metrics import confusion_matrix
import seaborn as sns
from matplotlib.colors import LogNorm


from sklearn.neighbors import NearestNeighbors

def class_infiltration(x_loc, annots, immune_annot, tumor_annot, n_neighbors = 40, thresh = 30):
    """\
        Description:
        ----------------
            Classify if the immune cell infiltrate the tumor or not
        
        Parameters:
        ----------------
            x_loc: the spatial locations of cells
            annots: the annotation of cells (including tumor and immune cell)
            immune_annot: the annotation of immune cell name
            tumor_annot: the annotation of tumor cell name
            n_neighbors: number of neighbors for knn graph construction
            thresh: knn classification threshold
    
    """
    nbrs_full = NearestNeighbors(n_neighbors = n_neighbors).fit(x_loc)
    distances, indices = nbrs_full.kneighbors(x_loc[annots == immune_annot, :])
    infiltration_annot = np.array([f"{immune_annot} (not infil)"] * indices.shape[0], dtype = object)
    for cell in range(indices.shape[0]):
        nbr_cts, nbr_ct_counts = np.unique(annots[indices[cell,:]], return_counts = True)
        if tumor_annot in nbr_cts:
            if nbr_ct_counts[nbr_cts == tumor_annot] > 30:
                infiltration_annot[cell] = f"{immune_annot} (infil)"
    annots[annots == immune_annot] = infiltration_annot

    return annots


# In[]
# -----------------------------------------------------------------------------------------
#
# load reference and query dataset
#
# -----------------------------------------------------------------------------------------
# NOTE: reference 2: public healthy data, with refined lymphoid annotations
# NOTE: the public dataset is already qc-ed
adata_ref1 =  sc.read_h5ad("../dataset/data_scrnaseq/data/qc_data/adata_qc_intact.h5ad")
adata_ref1.obs["annot"] = sc.read_h5ad("../dataset/data_scrnaseq/seurat_integration/adata_intact_seurat_lt_correct.h5ad").obs.loc[adata_ref1.obs.index,"annot (lt_correct)"].values
adata_ref1.obs["annot.level2"] = sc.read_h5ad("../dataset/data_scrnaseq/seurat_integration/adata_intact_seurat_lt_correct.h5ad").obs.loc[adata_ref1.obs.index,"annot.level2"].values

# select only lymphoid cells
adata_ref1 = adata_ref1[adata_ref1.obs["annot"] == "Lymphoid", :]
adata_ref1.layers["counts"] = adata_ref1.X.copy()

adata_ref2 = sc.read_h5ad("../dataset/reference_scrnaseq/reference_correct.h5ad")
# update the batch and annotation labels for LT
adata_ref2.obs["sample"] = adata_ref2.obs["batchID"].values
adata_ref2.obs["modality"] = "Ref"
# select only lymphoid cells
adata_ref2 = adata_ref2[adata_ref2.obs["annot"] == "Lymphoid", :]
adata_ref2.layers["counts"] = adata_ref2.X.copy()
# combine reference
adata_ref = anndata.concat([adata_ref1, adata_ref2], axis = 0, join = "inner", label = "modality", keys = ["Ref1", "Ref2"])

# result directory
# number of samples, maybe try different values for ref2, the population of ref2 is large, could run 500
dataset = "merfish0626"
nsamples_scanvi = 300
res_dir = f"../results_merfish/results_{dataset}_proseg/annot_{nsamples_scanvi}_combined/"
if not os.path.exists(res_dir):
    os.makedirs(res_dir)
    os.makedirs(res_dir + "plots_marker/")

# read in the combined annotation
adata_merfish = sc.read_h5ad(res_dir + "scanvi/adata_merfish_combfilter.h5ad")
# select only the annotation of lymphoid
adata_merfish_lymph = adata_merfish[adata_merfish.obs["annot_comb"] == "Lymphoid"].copy()

marker = ["Cd3e", "Cd4", "Cd8a", "Klrb1c", "Klrd1", "Nkg7", "Ms4a1"]

# In[]
# ------------------------------------------------------------------------------------------------------
#
# [Optional] Transfer lymphoid subtype annotation from reference 2 to merfish
#
# ------------------------------------------------------------------------------------------------------
# set the seed before running the model, for reproductivity
scvi.settings.seed = 0

# Train scVI on the reference dataset 
scvi.model.SCVI.setup_anndata(adata_ref, batch_key = "sample", layer = "counts", labels_key = "annot.level2")
scvi_ref = scvi.model.SCVI(adata_ref, use_layer_norm = "both", use_batch_norm = "none", encode_covariates = True, dropout_rate = 0.2, n_layers = 2)
scvi_ref.train()
# use the same label, layer, batch_label as in scvi
scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_ref, labels_key = "annot.level2", unlabeled_category="Unknown")
# default for n_samples_per_label is None, here we use the example value 100
scanvi_model.train(max_epochs = 20, n_samples_per_label = nsamples_scanvi)

# Label transfer
SCANVI_LATENT_KEY = "X_scANVI_lymph"
SCANVI_PREDICTIONS_KEY = "annot.level2"
scvi.model.SCANVI.prepare_query_anndata(adata_merfish_lymph, scanvi_model)
scanvi_query = scvi.model.SCANVI.load_query_data(adata_merfish_lymph, scanvi_model)
scanvi_query.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0}, check_val_every_n_epoch=10,)
adata_merfish_lymph.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation()
adata_merfish_lymph.obs[SCANVI_PREDICTIONS_KEY] = scanvi_query.predict()
scanvi_query.save(res_dir + "scanvi_lymph", overwrite = True)

# NOTE: visualize the integrated space
adata_merfish_lymph_pp = adata_merfish_lymph[adata_merfish_lymph.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"]
adata_merfish_lymph_ppc = adata_merfish_lymph[adata_merfish_lymph.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"]
adata_merfish_lymph_pp.obs["modality"] = "Merfish-pp"
adata_merfish_lymph_ppc.obs["modality"] = "Merfish-ppc"
adata_ref_merfish = anndata.concat([adata_ref, adata_merfish_lymph_pp, adata_merfish_lymph_ppc], join = "inner")
adata_ref_merfish.obs["modality"] = adata_ref_merfish.obs["modality"].astype("category") 
adata_ref_merfish.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation(adata_ref_merfish)
sc.pp.neighbors(adata_ref_merfish, use_rep=SCANVI_LATENT_KEY)
sc.tl.umap(adata_ref_merfish)
fig = utils.plot_by_batch(adata_ref_merfish.obsm["X_umap"], annos = np.array([x for x in adata_ref_merfish.obs[SCANVI_PREDICTIONS_KEY].values]),
                          batches = np.array([x for x in adata_ref_merfish.obs["modality"].values]), figsize = (10, 6),
                          ncols = 2, s = 5, alpha = 0.3, markerscale = 5)
fig.savefig(res_dir + "lymphoid_subtype_lt.png", bbox_inches = "tight", dpi = 150)
fig = utils.plot_embeds(adata_ref_merfish.obsm["X_umap"], annos = adata_ref_merfish.obs[["annot.level2", "modality"]].astype("category"), figsize = (9,6), markerscale = 5)
fig.savefig(res_dir + "lymphoid_subtype_lt2.png", bbox_inches = "tight", dpi = 150)
# # Leiden classification. 
# # NOTE: the annotation of immune subtypes are not accurate in reference scrna-seq dataset. 
# # need to redo the annotation using marker expression, separate CD4 T first from others using CD4 T markers 
# # (below is naive T cells without CD4-T marker). Then CD8 T and NK cells are mixed, CD3e can be partially expressed in NK, and NKG7 can be partially expressed in CD8
# sc.tl.leiden(adata_ref_merfish, resolution = 0.2)
# fig = utils.plot_by_batch(adata_ref_merfish.obsm["X_umap"], annos = np.array([x for x in adata_ref_merfish.obs.leiden.values]),
#                           batches = np.array([x for x in adata_ref_merfish.obs.modality.values]), figsize = (9,7), ncols = 3,
#                           markerscale = 5)
# fig.savefig(res_dir + "lymph_subtype_leiden.png", bbox_inches = "tight", dpi = 150)
# fig = utils.plot_embeds(adata_ref_merfish.obsm["X_umap"], annos = adata_ref_merfish.obs[["leiden", "modality"]],
#                         figsize = (9,7), ncols = 2, markerscale = 5)
# fig.savefig(res_dir + "lymph_subtype_leide2.png", bbox_inches = "tight", dpi = 150)

sc.pp.normalize_total(adata_ref_merfish, 1e4)
sc.pp.log1p(adata_ref_merfish)
# check marker gene expression for T cell, NK cell and B cell
marker = ["Cd3e", "Cd4", "Cd8a", "Klrb1c", "Klrd1", "Nkg7", "Ms4a1"]
gene_expr = pd.DataFrame(data = adata_ref_merfish[:, marker].X.toarray(), columns = marker, index = adata_ref_merfish.obs.index.values.squeeze())
vmax = [np.max(gene_expr[x].squeeze()) for x in marker]
for sample in ["Ref1", "Ref2", "Merfish-pp", "Merfish-ppc"]:
    fig = utils.plot_embeds_continuous(adata_ref_merfish[adata_ref_merfish.obs["modality"] == sample, :].obsm["X_umap"], 
                                                annos = gene_expr.loc[adata_ref_merfish.obs["modality"] == sample, :],
                                                colormap = utils.SUPER_MAGMA, vmax = vmax, figsize = (7,5), ncols = 4, alpha = 0.4, s = 2)

# fig = utils.plot_embeds_continuous(adata_ref_merfish.obsm["X_umap"], annos = gene_expr, colormap = utils.SUPER_MAGMA, vmax = vmax, figsize = (7,5), ncols = 4, alpha = 0.4, s = 2)
fig = utils.plot_embeds_continuous(adata_ref_merfish[adata_ref_merfish.obs["modality"].isin(["Merfish-pp", "Merfish-ppc"])].obsm["X_umap"], 
                                   annos = gene_expr.loc[adata_ref_merfish.obs["modality"].isin(["Merfish-pp", "Merfish-ppc"]),:], colormap = utils.SUPER_MAGMA, figsize = (7,5), ncols = 4, alpha = 0.4, s = 2)
fig.savefig(res_dir + f"plots_marker/lymphoid_markers_merfish.png", bbox_inches = "tight", dpi = 150)

adata_ref_merfish.write_h5ad(res_dir + "scanvi_lymph/adata_merfish_lymph.h5ad")

# In[]
# avoid re-running scanvi
adata_ref_merfish = sc.read_h5ad(res_dir + "scanvi_lymph/adata_merfish_lymph.h5ad")
# Incorporate the transferred annotation back to adata_merfish
adata_merfish.obs[SCANVI_PREDICTIONS_KEY] = adata_merfish.obs["annot_comb"].astype(object)
# NOTE: cancel reason, the lt result might not be accurate due to the inaccuracy of reference annotation
# adata_merfish.obs.loc[adata_merfish_lymph.obs.index.values, SCANVI_PREDICTIONS_KEY] = adata_merfish_lymph.obs[SCANVI_PREDICTIONS_KEY].values
# leiden annotation result
cell_idx = adata_ref_merfish.obs[~adata_ref_merfish.obs["modality"].isin(["Ref1", "Ref2"])].index.values.squeeze()
adata_merfish.obs.loc[cell_idx, SCANVI_PREDICTIONS_KEY] = adata_ref_merfish.obs.loc[cell_idx, "annot.level2"].values

# In[]
# ------------------------------------------------------------------------------------------------------
#
# Immune infiltration
#
# ------------------------------------------------------------------------------------------------------

# Epithelial cell can be treated as tumor: Luminal, Basal, Club epithelial
# Immune cell, first check macrophage
adata_merfish_pp = adata_merfish[adata_merfish.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"]
adata_merfish_ppc = adata_merfish[adata_merfish.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"]

immune_subtype = "Lymphoid"

annot_immune_infil = adata_merfish.obs[["annot_comb"]].astype(object)
annot_immune_infil.loc[~annot_immune_infil["annot_comb"].isin(["Club epithelia", "Basal", "Luminal"] + [immune_subtype]), "annot_comb"] = "1. Other"
annot_immune_infil.loc[annot_immune_infil["annot_comb"].isin(["Club epithelia", "Basal", "Luminal"]), "annot_comb"] = "2. Epithelial"
annot_immune_infil.loc[annot_immune_infil["annot_comb"].isin([immune_subtype]), "annot_comb"] = "3. " + immune_subtype
# cmap3colors = ListedColormap(colors = ["#e0e0e0", "#aacbeb", "#000000"])
# fig = utils.plot_embeds(adata_merfish_pp.obsm["spatial"], annos = annot_immune_infil[adata_merfish.obs["sample"] == "merfish0617-PP"].astype("category"), figsize = (20, 13), s = 3, alpha = 0.4, markerscale = 10, colormap = cmap3colors)
# fig = utils.plot_embeds(adata_merfish_ppc.obsm["spatial"], annos = annot_immune_infil[adata_merfish.obs["sample"] == "merfish0617-PPC"].astype("category"), figsize = (20, 13), s = 3, alpha = 0.4, markerscale = 10, colormap = cmap3colors)

annots_infil_pp = class_infiltration(x_loc = adata_merfish_pp.obsm["spatial"], annots = np.array([x for x in annot_immune_infil.loc[adata_merfish_pp.obs.index,"annot_comb"]], dtype = object),\
                                     immune_annot = "3. " + immune_subtype, tumor_annot = "2. Epithelial", n_neighbors = 40, thresh = 30)
        
annots_infil_ppc = class_infiltration(x_loc = adata_merfish_ppc.obsm["spatial"], annots = np.array([x for x in annot_immune_infil.loc[adata_merfish_ppc.obs.index,"annot_comb"]], dtype = object),\
                                     immune_annot = "3. " + immune_subtype, tumor_annot = "2. Epithelial", n_neighbors = 40, thresh = 30)

annot_immune_infil.loc[adata_merfish_pp.obs.index, "annot_comb"] = annots_infil_pp
annot_immune_infil.loc[adata_merfish_ppc.obs.index, "annot_comb"] = annots_infil_ppc
annot_immune_infil = annot_immune_infil.astype("category")
cmap4colors = ListedColormap(colors = ["#e0e0e0", "#aacbeb", "#000000", "#e888ab"])
fig = utils.plot_embeds(adata_merfish_pp.obsm["spatial"], annos = annot_immune_infil.loc[adata_merfish_pp.obs.index,:], figsize = (17, 10), s = 3, alpha = 0.4, markerscale = 10, colormap = cmap4colors)
fig.savefig(res_dir + f"{immune_subtype}_infiltration_pp.png", dpi = 150, bbox_inches = "tight")
fig = utils.plot_embeds(adata_merfish_ppc.obsm["spatial"], annos = annot_immune_infil.loc[adata_merfish_ppc.obs.index,:], figsize = (17, 10), s = 3, alpha = 0.4, markerscale = 10, colormap = cmap4colors)
fig.savefig(res_dir + f"{immune_subtype}_infiltration_ppc.png", dpi = 150, bbox_inches = "tight")


# In[]
# Compare infiltration across conditions
# lymphoid infiltration percentage
# NOTE: first step is remove SV cells, should not be included in the percentage calculation
lymph_status_pp, lymph_status_counts_pp = np.unique(annots_infil_pp[adata_merfish_pp.obs["annot_comb"] != "SV"], return_counts = True)
lymph_status_ppc, lymph_status_counts_ppc = np.unique(annots_infil_ppc[adata_merfish_ppc.obs["annot_comb"] != "SV"], return_counts = True)

lymph_status_pp = pd.DataFrame(index = lymph_status_pp, data = lymph_status_counts_pp[:, None], columns = ["counts"])
lymph_status_ppc = pd.DataFrame(index = lymph_status_ppc, data = lymph_status_counts_ppc[:, None], columns = ["counts"])
lymph_status_pp["pct"] = lymph_status_pp["counts"]/np.sum(lymph_status_pp["counts"])
lymph_status_ppc["pct"] = lymph_status_ppc["counts"]/np.sum(lymph_status_ppc["counts"])

lymph_counts_pp = lymph_status_pp.loc["3. Lymphoid (infil)", "counts"] + lymph_status_pp.loc["3. Lymphoid (not infil)", "counts"]
lymph_counts_ppc = lymph_status_ppc.loc["3. Lymphoid (infil)", "counts"] + lymph_status_ppc.loc["3. Lymphoid (not infil)", "counts"]
lymph_pct_pp = lymph_status_pp.loc["3. Lymphoid (infil)", "pct"] + lymph_status_pp.loc["3. Lymphoid (not infil)", "pct"]
lymph_pct_ppc = lymph_status_ppc.loc["3. Lymphoid (infil)", "pct"] + lymph_status_ppc.loc["3. Lymphoid (not infil)", "pct"]

lymph = pd.DataFrame(columns = ["counts", "pct", "cond"], index = ["PP", "PPC"], data = 0)
lymph["cond"] = np.array(["PP", "PPC"])
lymph["pct"] = np.array([lymph_pct_pp, lymph_pct_ppc])
lymph["counts"] = np.array([lymph_counts_pp, lymph_counts_ppc])

sns.set(font_scale=2)
fig = plt.figure(figsize = (7 * 2, 5))
ax = fig.subplots(nrows = 1, ncols = 2)
sns.barplot(data = lymph, x = "cond", y = "pct", ax = ax[0])
sns.barplot(data = lymph, x = "cond", y = "counts", ax = ax[1])
for i in ax[0].containers:
    ax[0].bar_label(i, fmt = "%.2f", fontsize = 15)
for i in ax[1].containers:
    ax[1].bar_label(i, fmt = "%d", fontsize = 15)
ax[0].set_xlabel(None)
ax[1].set_xlabel(None)
fig.suptitle("Lymphoid Percentage (no SV)", fontsize = 20)
plt.tight_layout()
fig.savefig(res_dir + "lymph_counts.png", dpi = 150, bbox_inches = "tight")

lymph_status_pct_pp = [lymph_status_pp.loc["3. Lymphoid (not infil)", "counts"]/np.sum(lymph_status_pp.loc["3. Lymphoid (not infil)", "counts"] + lymph_status_pp.loc["3. Lymphoid (infil)", "counts"]),
                       lymph_status_pp.loc["3. Lymphoid (infil)", "counts"]/np.sum(lymph_status_pp.loc["3. Lymphoid (not infil)", "counts"] + lymph_status_pp.loc["3. Lymphoid (infil)", "counts"])]
lymph_status_pct_ppc = [lymph_status_ppc.loc["3. Lymphoid (not infil)", "counts"]/np.sum(lymph_status_ppc.loc["3. Lymphoid (not infil)", "counts"] + lymph_status_ppc.loc["3. Lymphoid (infil)", "counts"]),
                       lymph_status_ppc.loc["3. Lymphoid (infil)", "counts"]/np.sum(lymph_status_ppc.loc["3. Lymphoid (not infil)", "counts"] + lymph_status_ppc.loc["3. Lymphoid (infil)", "counts"])]
lymph = pd.DataFrame(columns = ["counts", "pct", "cond", "infil_status"], index = ["PP", "PP", "PPC", "PPC"], data = 0)
lymph["cond"] = np.array(["PP", "PP", "PPC", "PPC"])
lymph["infil_status"] = np.array(["not infil", "infil", "not infil", "infil"])
lymph["pct"] = np.array(lymph_status_pct_pp + lymph_status_pct_ppc)
lymph["counts"] = np.array([lymph_status_pp.loc["3. Lymphoid (not infil)", "counts"],
                            lymph_status_pp.loc["3. Lymphoid (infil)", "counts"],
                            lymph_status_ppc.loc["3. Lymphoid (not infil)", "counts"],
                            lymph_status_ppc.loc["3. Lymphoid (infil)", "counts"]])

sns.set(font_scale=2)
fig = plt.figure(figsize = (7 * 2, 5))
ax = fig.subplots(nrows = 1, ncols = 2)
sns.barplot(data = lymph, x = "cond", hue = "infil_status", y = "pct", ax = ax[0])
sns.barplot(data = lymph, x = "cond", hue = "infil_status", y = "counts", ax = ax[1])
for i in ax[0].containers:
    ax[0].bar_label(i, fmt = "%.2f", fontsize = 15)
for i in ax[1].containers:
    ax[1].bar_label(i, fmt = "%d", fontsize = 15)
ax[0].set_xlabel(None)
ax[1].set_xlabel(None)
leg = ax[0].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
leg = ax[1].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))

fig.suptitle("Lymphoid infiltration percentage & count", fontsize = 20)
ax[0].set_title("Pct of Lymph (infil)", fontsize = 15)
ax[1].set_title("Count of Lymph (infil)", fontsize = 15)
plt.tight_layout()
fig.savefig(res_dir + "lymph_infil_counts.png", dpi = 150, bbox_inches = "tight")

# In[]
# calculate the pie chart of infiltrated lymphoid subtypes
annot_immune_infil["annot.sub.infil"] = adata_merfish.obs[SCANVI_PREDICTIONS_KEY].astype(object)
for barcode in annot_immune_infil.index.values:
    if annot_immune_infil.loc[barcode, "annot_comb"] == "1. Other":
        annot_immune_infil.loc[barcode, "annot.sub.infil"] = "1. Other"
    elif annot_immune_infil.loc[barcode, "annot_comb"] == "2. Epithelial":
        annot_immune_infil.loc[barcode, "annot.sub.infil"] = "2. Epithelial"
    elif annot_immune_infil.loc[barcode, "annot_comb"] == "3. Lymphoid (not infil)":
        annot_immune_infil.loc[barcode, "annot.sub.infil"] = "3. " + annot_immune_infil.loc[barcode, "annot.sub.infil"] + " not infil"
    elif annot_immune_infil.loc[barcode, "annot_comb"] == "3. Lymphoid (infil)":
        annot_immune_infil.loc[barcode, "annot.sub.infil"] = "3. " + annot_immune_infil.loc[barcode, "annot.sub.infil"] + " infil"

# NOTE: remove SV again for percentage calculation
lymph_status_pp = annot_immune_infil.loc[(adata_merfish.obs["annot_comb"] != "SV") & (adata_merfish.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"), "annot.sub.infil"].values
lymph_status_ppc = annot_immune_infil.loc[(adata_merfish.obs["annot_comb"] != "SV") & (adata_merfish.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"), "annot.sub.infil"].values


# proportion of infiltration for each subtypes
prop_df = pd.DataFrame(columns = ["prop", "ct", "gt"])
# select immune subtypes
prop_df["ct"] = ["CD4 T", "CD8 T", "DN T", "NK", "B"] * 2
prop_df["gt"] = ["PP"] * 5 + ["PPC"] * 5
prop_df["prop"] = 0
# NOTE: B population too small    
for ct in ["CD4 T", "CD8 T", "DN T", "NK"]:
    prop_df.loc[(prop_df["ct"] == ct) & (prop_df["gt"] == "PP"), "prop"] = np.sum(lymph_status_pp == f"3. {ct} infil")/np.sum((lymph_status_pp == f"3. {ct} infil") | (lymph_status_pp == f"3. {ct} not infil"))
    prop_df.loc[(prop_df["ct"] == ct) & (prop_df["gt"] == "PPC"), "prop"] = np.sum(lymph_status_ppc == f"3. {ct} infil")/np.sum((lymph_status_ppc == f"3. {ct} infil") | (lymph_status_ppc == f"3. {ct} not infil"))


fig = plt.figure(figsize = (7, 5))
ax = fig.subplots(nrows = 1, ncols = 1)
sns.barplot(prop_df, x = "ct", hue = "gt", y = "prop", ax = ax)
for i in ax.containers:
    ax.bar_label(i, fmt = "%.2f", fontsize = 15)
leg = ax.legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
ax.set_xlabel(None)
ax.set_ylabel("Proportion")
ax.set_title("Proportion of infiltration")

fig.savefig(res_dir + "lymph_sub_infil_counts.png", dpi = 150, bbox_inches = "tight")

# In[]
import matplotlib
matplotlib.rc_file_defaults()

cmap_infil = ListedColormap(colors = ["#aacbeb", "#e0e0e0", "#ffffff", "#d8789b", "#f898fb", "#cf4f0e", "#ff7f0e", "#bcbd22", "#2ca02c", "#17becf", "#3f97d4"])
annot_immune_infil["annot.sub.infil"] = annot_immune_infil["annot.sub.infil"].astype("category")
fig = utils.plot_embeds(adata_merfish_pp.obsm["spatial"], annos = annot_immune_infil.loc[adata_merfish_pp.obs.index, ["annot.sub.infil"]], figsize = (17, 10), s = 3, alpha = 1, markerscale = 10, colormap = cmap_infil)
fig = utils.plot_embeds(adata_merfish_ppc.obsm["spatial"], annos = annot_immune_infil.loc[adata_merfish_ppc.obs.index, ["annot.sub.infil"]], figsize = (17, 10), s = 3, alpha = 1, markerscale = 10, colormap = cmap_infil)



cts = ["3. B", "3. CD4 T", "3. CD8 T", "3. DN T", "3. NK"]
ct_df = pd.DataFrame(columns = cts, index = annot_immune_infil.index, data = "1. Other")
for ct in cts:
    ct_df.loc[annot_immune_infil["annot.sub.infil"] == ct + " infil", ct] = ct + " infil"
    ct_df.loc[annot_immune_infil["annot.sub.infil"] == ct + " not infil", ct] = ct + " not infil"
    ct_df.loc[annot_immune_infil["annot.sub.infil"] == "1. Other", ct] = "1. Other"
    ct_df.loc[annot_immune_infil["annot.sub.infil"] == "2. Epithelial", ct] = "2. Epithelial"
    ct_df[ct] = ct_df[ct].astype("category")
cmap4colors = ListedColormap(colors = ["#e0e0e0", "#aacbeb", "#000000", "#e888ab"])
fig = utils.plot_embeds(adata_merfish_pp.obsm["spatial"], annos = ct_df.loc[adata_merfish_pp.obs.index,:], figsize = (17, 10), s = 3, alpha = 0.4, markerscale = 10, colormap = cmap4colors)
fig.savefig(res_dir + f"{immune_subtype}_sub_infiltration_pp.png", dpi = 150, bbox_inches = "tight")
fig = utils.plot_embeds(adata_merfish_ppc.obsm["spatial"], annos = ct_df.loc[adata_merfish_ppc.obs.index,:], figsize = (17, 10), s = 3, alpha = 0.4, markerscale = 10, colormap = cmap4colors)
fig.savefig(res_dir + f"{immune_subtype}_sub_infiltration_ppc.png", dpi = 150, bbox_inches = "tight")



# %%
