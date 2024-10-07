# In[]
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
import sys, os
import matplotlib.pyplot as plt
sys.path.append("../")
import utils
import seaborn as sns
import matplotlib
from matplotlib.colors import ListedColormap
import scvi


# In[]
# -----------------------------------------------------------------------------------------
#
# load reference and query dataset
#
# -----------------------------------------------------------------------------------------
# NOTE: no refined annotations in the reference scRNA-seq dataset
# load reference scRNA-seq dataset, NOTE: use raw data, include 4 samples/batches
adata_ref1 =  sc.read_h5ad("../dataset/data_scrnaseq/data/qc_data/adata_qc_intact.h5ad")
adata_ref1.obs["annot"] = sc.read_h5ad("../dataset/data_scrnaseq/seurat_integration/adata_intact_seurat_lt_correct.h5ad").obs.loc[adata_ref1.obs.index,"annot (lt_correct)"].values
adata_ref1.obs["annot.level2"] = sc.read_h5ad("../dataset/data_scrnaseq/seurat_integration/adata_intact_seurat_lt_correct.h5ad").obs.loc[adata_ref1.obs.index,"annot.level2"].values

# select only the macrophage
adata_ref1 = adata_ref1[adata_ref1.obs["annot"] == "Macrophage", :]
adata_ref1.layers["counts"] = adata_ref1.X.astype(int)

adata_ref2 = sc.read_h5ad("../dataset/reference_scrnaseq/reference_correct.h5ad")
# update the batch and annotation labels for LT
adata_ref2.obs["sample"] = adata_ref2.obs["batchID"].values
# select only macrophage
adata_ref2 = adata_ref2[adata_ref2.obs["annot"] == "Macrophage", :]
adata_ref2.layers["counts"] = adata_ref2.X.astype(int)

# combine the reference annotation
adata_ref = anndata.concat([adata_ref1, adata_ref2], axis = 0, join = "inner", label = "modality", keys = ["Ref1", "Ref2"])

# NOTE: use the filtered anndata with higher quality immune cell annotation
dataset = "merfish0617"
nsamples_scanvi = 300
# NOTE: use proseg segmentation
res_dir = f"../results_merfish/results_{dataset}_proseg/annot_{nsamples_scanvi}_combined/"
adata_merfish = sc.read_h5ad(res_dir + "scanvi/adata_merfish_combfilter.h5ad")
# the prediction key used in the adata
SCANVI_PREDICTIONS_KEY = "annot_comb"
adata_merfish.X = adata_merfish.layers["counts"].copy()
# select only the macrophage cells
adata_merfish_macrophage = adata_merfish[adata_merfish.obs["annot_comb"] == "Macrophage"]

marker_m2 = ["Cd163", "Mrc1", "Ccl24", "Arg1", "Msr1", "Il10"]
marker_m1 = ["Il1b", "Il6", "Tnf", "Nos2", "Cxcl9"]


# In[]
# -----------------------------------------------------------------------------------------
#
# scanvi for label transfer, from annotated (M1 & M2) reference scRNA-seq to MerFISH
#
# -----------------------------------------------------------------------------------------
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

# Label transfer to merfish
SCANVI_LATENT_KEY = "X_scANVI_macrophage"
SCANVI_PREDICTIONS_KEY = "annot.level2"
scvi.model.SCANVI.prepare_query_anndata(adata_merfish_macrophage, scanvi_model)
scanvi_query = scvi.model.SCANVI.load_query_data(adata_merfish_macrophage, scanvi_model)
scanvi_query.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0}, check_val_every_n_epoch=10,)
adata_merfish_macrophage.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation()
adata_merfish_macrophage.obs[SCANVI_PREDICTIONS_KEY] = scanvi_query.predict()
scanvi_query.save(res_dir + "scanvi_macro", overwrite = True)

# visualize the integrated space
adata_merfish_macro_pp = adata_merfish_macrophage[adata_merfish_macrophage.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"]
adata_merfish_macro_ppc = adata_merfish_macrophage[adata_merfish_macrophage.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"]
adata_merfish_macro_pp.obs["modality"] = "Merfish-pp"
adata_merfish_macro_ppc.obs["modality"] = "Merfish-ppc"
adata_ref_merfish = anndata.concat([adata_ref, adata_merfish_macro_pp, adata_merfish_macro_ppc], join = "inner")
adata_ref_merfish.obs["modality"] = adata_ref_merfish.obs["modality"].astype("category") 
adata_ref_merfish.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation(adata_ref_merfish)
sc.pp.neighbors(adata_ref_merfish, use_rep=SCANVI_LATENT_KEY)
sc.tl.umap(adata_ref_merfish)
adata_ref_merfish.obs[SCANVI_PREDICTIONS_KEY] = adata_ref_merfish.obs[SCANVI_PREDICTIONS_KEY].astype("category")
fig = utils.plot_by_batch(adata_ref_merfish.obsm["X_umap"], annos = np.array([x for x in adata_ref_merfish.obs[SCANVI_PREDICTIONS_KEY].values]),
                          batches = np.array([x for x in adata_ref_merfish.obs["modality"].values]), figsize = (10, 6),
                          ncols = 4, s = 5, alpha = 0.3, markerscale = 5)
fig.savefig(res_dir + "macrophage_subtype_lt.png", bbox_inches = "tight", dpi = 150)
fig = utils.plot_embeds(adata_ref_merfish.obsm["X_umap"], annos = adata_ref_merfish.obs[[SCANVI_PREDICTIONS_KEY, "modality"]], figsize = (9,5), markerscale = 5)
fig.savefig(res_dir + "macrophage_subtype_lt2.png", bbox_inches = "tight", dpi = 150)

# check the marker gene expression, validate the annotation result 
adata_ref_merfish.X = adata_ref_merfish.layers["counts"].copy()
sc.pp.normalize_total(adata_ref_merfish, 10e4)
sc.pp.log1p(adata_ref_merfish)
gene_expr = pd.DataFrame(data = adata_ref_merfish[:, marker_m1 + marker_m2].X.toarray(), columns = marker_m1 + marker_m2, index = adata_ref_merfish.obs.index.values.squeeze())
# fig = utils.plot_embeds_continuous(adata_ref_merfish.obsm["X_umap"], annos = gene_expr, colormap = utils.SUPER_MAGMA, figsize = (7,5), ncols = 4, alpha = 0.4, s = 2)
# check only on merfish population
fig = utils.plot_embeds_continuous(adata_ref_merfish[adata_ref_merfish.obs["modality"].isin(["Merfish-pp", "Merfish-ppc"])].obsm["X_umap"], 
                                   annos = gene_expr.loc[adata_ref_merfish.obs["modality"].isin(["Merfish-pp", "Merfish-ppc"]),:], colormap = utils.SUPER_MAGMA, figsize = (7,5), ncols = 4, alpha = 0.4, s = 2)
fig.savefig(res_dir + "plots_marker/macrophage_markers_merfish.png", bbox_inches = "tight", dpi = 150)
adata_ref_merfish.write_h5ad(res_dir + "scanvi_macro/adata_merfish_macro.h5ad")


# In[]
# -----------------------------------------------------------------------------------------
#
# Macrophage infiltration analysis
#
# -----------------------------------------------------------------------------------------
# Epithelial cell can be treated as tumor: Luminal, Basal, Club epithelial
# Immune cell, first check macrophage
adata_merfish_pp = adata_merfish[adata_merfish.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"]
adata_merfish_ppc = adata_merfish[adata_merfish.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"]

immune_subtype = "Macrophage"
annot_immune_infil = adata_merfish.obs[["annot_comb"]].astype(object)
annot_immune_infil.loc[~annot_immune_infil["annot_comb"].isin(["Club epithelia", "Basal", "Luminal"] + [immune_subtype]), "annot_comb"] = "1. Other"
annot_immune_infil.loc[annot_immune_infil["annot_comb"].isin(["Club epithelia", "Basal", "Luminal"]), "annot_comb"] = "2. Epithelial"
annot_immune_infil.loc[annot_immune_infil["annot_comb"].isin([immune_subtype]), "annot_comb"] = "3. " + immune_subtype
# annot_immune_infil = annot_immune_infil.astype("category")
# cmap3colors = ListedColormap(colors = ["#e0e0e0", "#aacbeb", "#000000"])
# fig = utils.plot_embeds(adata_merfish_pp.obsm["spatial"], annos = annot_immune_infil[adata_merfish.obs["sample"] == "merfish0617-PP"], figsize = (20, 13), s = 3, alpha = 0.4, markerscale = 10, colormap = cmap3colors)
# fig = utils.plot_embeds(adata_merfish_ppc.obsm["spatial"], annos = annot_immune_infil[adata_merfish.obs["sample"] == "merfish0617-PPC"], figsize = (20, 13), s = 3, alpha = 0.4, markerscale = 10, colormap = cmap3colors)

annots_infil_pp = utils.class_infiltration(x_loc = adata_merfish_pp.obsm["spatial"], annots = np.array([x for x in annot_immune_infil.loc[adata_merfish_pp.obs.index,"annot_comb"]], dtype = object),\
                                     immune_annot = "3. " + immune_subtype, tumor_annot = "2. Epithelial", n_neighbors = 40, thresh = 30)
        
annots_infil_ppc = utils.class_infiltration(x_loc = adata_merfish_ppc.obsm["spatial"], annots = np.array([x for x in annot_immune_infil.loc[adata_merfish_ppc.obs.index,"annot_comb"]], dtype = object),\
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
# Macrophage infiltration percentage
# NOTE: For percentage calculation, always remove the SV cells (not in the tme)
macro_status_pp, macro_status_counts_pp = np.unique(annots_infil_pp[adata_merfish_pp.obs["annot_comb"] != "SV"], return_counts = True)
macro_status_ppc, macro_status_counts_ppc = np.unique(annots_infil_ppc[adata_merfish_ppc.obs["annot_comb"] != "SV"], return_counts = True)

macro_status_pp = pd.DataFrame(index = macro_status_pp, data = macro_status_counts_pp[:, None], columns = ["counts"])
macro_status_ppc = pd.DataFrame(index = macro_status_ppc, data = macro_status_counts_ppc[:, None], columns = ["counts"])
macro_status_pp["pct"] = macro_status_pp["counts"]/np.sum(macro_status_pp["counts"])
macro_status_ppc["pct"] = macro_status_ppc["counts"]/np.sum(macro_status_ppc["counts"])

macro_counts_pp = macro_status_pp.loc["3. Macrophage (infil)", "counts"] + macro_status_pp.loc["3. Macrophage (not infil)", "counts"]
macro_counts_ppc = macro_status_ppc.loc["3. Macrophage (infil)", "counts"] + macro_status_ppc.loc["3. Macrophage (not infil)", "counts"]
macro_pct_pp = macro_status_pp.loc["3. Macrophage (infil)", "pct"] + macro_status_pp.loc["3. Macrophage (not infil)", "pct"]
macro_pct_ppc = macro_status_ppc.loc["3. Macrophage (infil)", "pct"] + macro_status_ppc.loc["3. Macrophage (not infil)", "pct"]

macro = pd.DataFrame(columns = ["counts", "pct", "cond"], index = ["PP", "PPC"], data = 0)
macro["cond"] = np.array(["PP", "PPC"])
macro["pct"] = np.array([macro_pct_pp, macro_pct_ppc])
macro["counts"] = np.array([macro_counts_pp, macro_counts_ppc])

sns.set(font_scale=2)
fig = plt.figure(figsize = (7 * 2, 5))
ax = fig.subplots(nrows = 1, ncols = 2)
sns.barplot(data = macro, x = "cond", y = "pct", ax = ax[0])
sns.barplot(data = macro, x = "cond", y = "counts", ax = ax[1])
for i in ax[0].containers:
    ax[0].bar_label(i, fmt = "%.2f", fontsize = 15)
for i in ax[1].containers:
    ax[1].bar_label(i, fmt = "%d", fontsize = 15)
ax[0].set_xlabel(None)
ax[1].set_xlabel(None)
fig.suptitle("Macrophage Percentage (no SV)", fontsize = 20)
plt.tight_layout()
fig.savefig(res_dir + "macrophage_counts.png", dpi = 150, bbox_inches = "tight")

macro_status_pct_pp = [macro_status_pp.loc["3. Macrophage (not infil)", "counts"]/np.sum(macro_status_pp.loc["3. Macrophage (not infil)", "counts"] + macro_status_pp.loc["3. Macrophage (infil)", "counts"]),
                       macro_status_pp.loc["3. Macrophage (infil)", "counts"]/np.sum(macro_status_pp.loc["3. Macrophage (not infil)", "counts"] + macro_status_pp.loc["3. Macrophage (infil)", "counts"])]
macro_status_pct_ppc = [macro_status_ppc.loc["3. Macrophage (not infil)", "counts"]/np.sum(macro_status_ppc.loc["3. Macrophage (not infil)", "counts"] + macro_status_ppc.loc["3. Macrophage (infil)", "counts"]),
                       macro_status_ppc.loc["3. Macrophage (infil)", "counts"]/np.sum(macro_status_ppc.loc["3. Macrophage (not infil)", "counts"] + macro_status_ppc.loc["3. Macrophage (infil)", "counts"])]
macro = pd.DataFrame(columns = ["counts", "pct", "cond", "infil_status"], index = ["PP", "PP", "PPC", "PPC"], data = 0)
macro["cond"] = np.array(["PP", "PP", "PPC", "PPC"])
macro["infil_status"] = np.array(["not infil", "infil", "not infil", "infil"])
macro["pct"] = np.array(macro_status_pct_pp + macro_status_pct_ppc)
macro["counts"] = np.array([macro_status_pp.loc["3. Macrophage (not infil)", "counts"],
                            macro_status_pp.loc["3. Macrophage (infil)", "counts"],
                            macro_status_ppc.loc["3. Macrophage (not infil)", "counts"],
                            macro_status_ppc.loc["3. Macrophage (infil)", "counts"]])

sns.set(font_scale=2)
fig = plt.figure(figsize = (7 * 2, 5))
ax = fig.subplots(nrows = 1, ncols = 2)
sns.barplot(data = macro, x = "cond", hue = "infil_status", y = "pct", ax = ax[0])
sns.barplot(data = macro, x = "cond", hue = "infil_status", y = "counts", ax = ax[1])
for i in ax[0].containers:
    ax[0].bar_label(i, fmt = "%.2f", fontsize = 15)
for i in ax[1].containers:
    ax[1].bar_label(i, fmt = "%d", fontsize = 15)
ax[0].set_xlabel(None)
ax[1].set_xlabel(None)
leg = ax[0].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
leg = ax[1].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))

fig.suptitle("Macrophage infiltration percentage & count", fontsize = 20)
ax[0].set_title("Pct of Macro (infil)", fontsize = 15)
ax[1].set_title("Count of Macro (infil)", fontsize = 15)
plt.tight_layout()
fig.savefig(res_dir + "macro_infil_counts.png", dpi = 150, bbox_inches = "tight")

# In[]
# -----------------------------------------------------------------------------------------
#
# Macrophage subtypes infiltration analysis, 
# NOTE: start to use the subtype annotation generated from scanvi label transfer before
#
# -----------------------------------------------------------------------------------------
# calculate the pie chart of infiltrated macrophage subtypes
adata_merfish.obs["annot_comb"] = adata_merfish.obs["annot_comb"].astype(object)
adata_merfish.obs.loc[adata_merfish_macrophage.obs.index, "annot_comb"] = adata_ref_merfish.obs.loc[adata_merfish_macrophage.obs.index, "annot.level2"].values
annot_immune_infil["annot.sub.infil"] = adata_merfish.obs["annot_comb"].astype(object)
for barcode in annot_immune_infil.index.values:
    if annot_immune_infil.loc[barcode, "annot_comb"] == "1. Other":
        annot_immune_infil.loc[barcode, "annot.sub.infil"] = "1. Other"
    elif annot_immune_infil.loc[barcode, "annot_comb"] == "2. Epithelial":
        annot_immune_infil.loc[barcode, "annot.sub.infil"] = "2. Epithelial"
    elif annot_immune_infil.loc[barcode, "annot_comb"] == "3. Macrophage (not infil)":
        annot_immune_infil.loc[barcode, "annot.sub.infil"] = "3. " + annot_immune_infil.loc[barcode, "annot.sub.infil"] + " not infil"
    elif annot_immune_infil.loc[barcode, "annot_comb"] == "3. Macrophage (infil)":
        annot_immune_infil.loc[barcode, "annot.sub.infil"] = "3. " + annot_immune_infil.loc[barcode, "annot.sub.infil"] + " infil"

# NOTE: remove SV again for percentage calculation
macro_status_pp, macro_status_counts_pp = np.unique(annot_immune_infil.loc[(adata_merfish.obs["annot_comb"] != "SV") & (adata_merfish.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"), "annot.sub.infil"].values, return_counts = True)
macro_status_ppc, macro_status_counts_ppc = np.unique(annot_immune_infil.loc[(adata_merfish.obs["annot_comb"] != "SV") & (adata_merfish.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"), "annot.sub.infil"].values, return_counts = True)

# composition of the infiltrated cells
# fig = plt.figure(figsize = (20, 7))
# axs = fig.subplots(nrows = 1, ncols = 2)
# # axs[0].pie(lymph_status_counts_pp[1:4], labels = lymph_status_pp[1:4], autopct='%1.1f%%')
# axs[0].pie(lymph_status_counts_pp[4:-1], labels = lymph_status_pp[4:-1], autopct='%1.1f%%')
# axs[1].pie(lymph_status_counts_ppc[4:-1], labels = lymph_status_ppc[4:-1], autopct='%1.1f%%')

# proportion of infiltration for each subtypes
prop_df = pd.DataFrame(columns = ["prop", "ct", "gt"])
# select immune subtypes
prop_df["ct"] = ["M1 macro", "M2 macro"] * 2
prop_df["gt"] = ["PP"] * 2 + ["PPC"] * 2
prop_df["prop"] = 0
# NOTE: B population too small    
for ct in ["M1 macro", "M2 macro"]:
    prop_df.loc[(prop_df["ct"] == ct) & (prop_df["gt"] == "PP"), "prop"] = macro_status_counts_pp[macro_status_pp == f"3. {ct} infil"]/(macro_status_counts_pp[macro_status_pp == f"3. {ct} infil"] + macro_status_counts_pp[macro_status_pp == f"3. {ct} not infil"])
    prop_df.loc[(prop_df["ct"] == ct) & (prop_df["gt"] == "PPC"), "prop"] = macro_status_counts_ppc[macro_status_ppc == f"3. {ct} infil"]/(macro_status_counts_ppc[macro_status_ppc == f"3. {ct} infil"] + macro_status_counts_ppc[macro_status_ppc == f"3. {ct} not infil"])


fig = plt.figure(figsize = (7, 5))
ax = fig.subplots(nrows = 1, ncols = 1)
sns.barplot(prop_df, x = "ct", hue = "gt", y = "prop", ax = ax)
for i in ax.containers:
    ax.bar_label(i, fmt = "%.2f", fontsize = 15)
leg = ax.legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
ax.set_xlabel(None)
ax.set_ylabel("Proportion")
ax.set_title("Proportion of infiltration")

fig.savefig(res_dir + "macro_sub_infil_counts.png", dpi = 150, bbox_inches = "tight")


# In[]
matplotlib.rc_file_defaults()

# cmap_infil = ListedColormap(colors = ["#aacbeb", "#e0e0e0", "#ffffff", "#d8789b", "#f898fb", "#cf4f0e", "#ff7f0e", "#bcbd22", "#2ca02c", "#17becf", "#3f97d4"])
# annot_immune_infil["annot.sub.infil"] = annot_immune_infil["annot.sub.infil"].astype("category")
# fig = utils.plot_embeds(adata_merfish_pp.obsm["spatial"], annos = annot_immune_infil.loc[adata_merfish_pp.obs.index, ["annot.sub.infil"]], figsize = (17, 10), s = 3, alpha = 1, markerscale = 10, colormap = cmap_infil)
# fig = utils.plot_embeds(adata_merfish_ppc.obsm["spatial"], annos = annot_immune_infil.loc[adata_merfish_ppc.obs.index, ["annot.sub.infil"]], figsize = (17, 10), s = 3, alpha = 1, markerscale = 10, colormap = cmap_infil)

cts = ["3. M1 macro", "3. M2 macro"]
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

# In[]
# TODO: DE analysis for the infiltrated and not-infiltrated macrophage subtypes

# %%
