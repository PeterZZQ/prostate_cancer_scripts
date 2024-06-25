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
# 1. Reference: scRNA-seq dataset, include 4 samples/batches.
# 2. Query: merfish dataset, include 2 samples/batches: Merfish-PP15, Merfish-PPC15
#
# -----------------------------------------------------------------------------------------
#
# load query dataset, region 0: PPC15, region 1: PP15
# 1. load region 0
adata_merfish_ppc = sc.read_h5ad("../dataset/data_merfish/region_0/adata_qc.h5ad")
adata_merfish_ppc.obs["annot"] = pd.Categorical(["Unknown"] * adata_merfish_ppc.shape[0], categories = ["Unknown"])
adata_merfish_ppc.obs["genotype"] = pd.Categorical(["PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"] * adata_merfish_ppc.shape[0], categories = ["PbCre(+/-),Pten(-/-),P53(-/-)", "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"])
adata_merfish_ppc.obs["sample"] = pd.Categorical(["Merfish-PPC15"] * adata_merfish_ppc.shape[0], categories = ["Merfish-PP15", "Merfish-PPC15"])

# 2. load region 1
adata_merfish_pp = sc.read_h5ad("../dataset/data_merfish/region_1/adata_qc.h5ad")
adata_merfish_pp.obs["annot"] = pd.Categorical(["Unknown"] * adata_merfish_pp.shape[0], categories = ["Unknown"])
adata_merfish_pp.obs["genotype"] = pd.Categorical(["PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"] * adata_merfish_pp.shape[0], categories = ["PbCre(+/-),Pten(-/-),P53(-/-)", "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"])
adata_merfish_pp.obs["sample"] = pd.Categorical(["Merfish-PP15"] * adata_merfish_pp.shape[0], categories = ["Merfish-PP15", "Merfish-PPC15"])

# load reference scRNA-seq dataset, NOTE: use raw data, include 4 samples/batches
adata_ref =  sc.read_h5ad("../dataset/data_scrnaseq/data/qc_data/adata_qc_intact.h5ad")
adata_ref.obs["annot"] = sc.read_h5ad("../dataset/data_scrnaseq/seurat_integration/adata_intact_seurat.h5ad").obs.loc[adata_ref.obs.index,"annot"].values

# find overlapping features
gene_overlap = np.intersect1d(adata_merfish_pp.var.index.values, adata_ref.var.index.values)
adata_ref = adata_ref[:, gene_overlap]
adata_merfish_pp = adata_merfish_pp[:, gene_overlap]
adata_merfish_ppc = adata_merfish_ppc[:, gene_overlap]
adata_ref.layers["counts"] = adata_ref.X.copy()
adata_merfish_pp.layers["counts"] = adata_merfish_pp.X.copy()
adata_merfish_ppc.layers["counts"] = adata_merfish_ppc.X.copy()

# result directory
res_dir = "../results_merfish/annot/"
if not os.path.exists(res_dir):
    os.makedirs(res_dir)


# In[] 
# -----------------------------------------------------------------------------------------
#
# Label transfer -- scANVI
# obtain label annotation
#
# -----------------------------------------------------------------------------------------
# NOTE: Step 1, Train scVI on reference scRNA-seq dataset
# Train scVI on the reference dataset 
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

# In[]
# NOTE: Step 2, Train scANVI based on scVI model
# Reference
# use scvi initialized with the reference dataset for better result
scvi_model = scvi.model.SCVI.load(res_dir + "scvi_reference", adata = adata_ref)
# use the same label, layer, batch_label as in scvi
scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_model, labels_key = "annot", unlabeled_category="Unknown")
# default for n_samples_per_label is None, here we use the example value 100
scanvi_model.train(max_epochs = 20, n_samples_per_label = 500)
scanvi_model.save(res_dir + "scanvi_reference", overwrite=True)

SCANVI_LATENT_KEY = "X_scANVI"

# # Sanity check, visualize the trained reference result (scANVI)
# adata_ref.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation()
# sc.pp.neighbors(adata_ref, use_rep=SCANVI_LATENT_KEY)
# sc.tl.umap(adata_ref)
# sc.pl.umap(adata_ref, color=["sample", "annot"], frameon=False, ncols=1)
# utils.plot_latent(adata_ref.obsm["X_umap"], annos = np.array([x for x in adata_ref.obs["annot"].values]), batches = np.array([x for x in adata_ref.obs["sample"].values]), mode = "separate", figsize = (10, 25))

# 1. Query classification, obtain the PP annotation
scvi.model.SCANVI.prepare_query_anndata(adata_merfish_pp, scanvi_model)
scanvi_query_pp = scvi.model.SCANVI.load_query_data(adata_merfish_pp, scanvi_model)
scanvi_query_pp.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0}, check_val_every_n_epoch=10,)
SCANVI_PREDICTIONS_KEY = "predictions_scanvi"

adata_merfish_pp.obsm[SCANVI_LATENT_KEY] = scanvi_query_pp.get_latent_representation()
adata_merfish_pp.obs[SCANVI_PREDICTIONS_KEY] = scanvi_query_pp.predict()

scanvi_query_pp.save(res_dir + "scanvi_pp", overwrite = True)

adata_merfish_pp.write_h5ad(res_dir + "scanvi_pp/adata_merfish_pp.h5ad")
adata_merfish_pp.obs.to_csv(res_dir + "scanvi_pp/scanvi_annot_pp.csv")

# 2. Query classification, obtain the PPC annotation
scvi.model.SCANVI.prepare_query_anndata(adata_merfish_ppc, scanvi_model)
scanvi_query_ppc = scvi.model.SCANVI.load_query_data(adata_merfish_ppc, scanvi_model)
scanvi_query_ppc.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0}, check_val_every_n_epoch=10,)
SCANVI_PREDICTIONS_KEY = "predictions_scanvi"

adata_merfish_ppc.obsm[SCANVI_LATENT_KEY] = scanvi_query_ppc.get_latent_representation()
adata_merfish_ppc.obs[SCANVI_PREDICTIONS_KEY] = scanvi_query_ppc.predict()

scanvi_query_ppc.save(res_dir + "scanvi_ppc", overwrite = True)
adata_merfish_ppc.write_h5ad(res_dir + "scanvi_ppc/adata_merfish_ppc.h5ad")
adata_merfish_ppc.obs.to_csv(res_dir + "scanvi_ppc/scanvi_annot_ppc.csv")


# In[]
# -----------------------------------------------------------------------------------------
#
# Load annotation result and marker information for downstream analysis
#
# -----------------------------------------------------------------------------------------
adata_merfish_pp = sc.read_h5ad(res_dir + "scanvi_pp/adata_merfish_pp.h5ad")
adata_merfish_ppc = sc.read_h5ad(res_dir + "scanvi_ppc/adata_merfish_ppc.h5ad")
SCANVI_LATENT_KEY = "X_scANVI"
SCANVI_PREDICTIONS_KEY = "predictions_scanvi"
adata_merfish_pp.obsm["spatial"] = np.vstack([0.5 * (adata_merfish_pp.obs["min_x"].values + adata_merfish_pp.obs["max_x"].values), 0.5 * (adata_merfish_pp.obs["min_y"].values + adata_merfish_pp.obs["max_y"].values)]).T 
adata_merfish_ppc.obsm["spatial"] = np.vstack([0.5 * (adata_merfish_ppc.obs["min_x"].values + adata_merfish_ppc.obs["max_x"].values), 0.5 * (adata_merfish_ppc.obs["min_y"].values + adata_merfish_ppc.obs["max_y"].values)]).T 

# Read in the gene panel information
gene_panel_info = pd.read_csv("../dataset/data_merfish/MERSCOPE_gene_panel_info.csv", index_col = 0)
markers = {}
# picked only the one that are included in spatial data
markers["luminal"] = ["Ar", "Krt8"]
markers["basal"] = ["Trp63", "Krt5"]
markers["club_epithelia"] = ["Agr2", "Krt7"]
markers["endothelial"] = ["Ackr1", "Cldn5"]
markers["lymphoid"] = ["Cd3e", "Ms4a1", "Klrb1c"]
markers["myeloid"] = ["Ptprc", "Itgam"]
markers["monocytes"] = ["Ptprc", "Itgam", "Cd14", "S100a8", "S100a9"] 
markers["macrophage"] = ["Ptprc", "Itgam", "Adgre1"]
markers["macrophage_m1"] = ["Ptprc", "Itgam", "Adgre1", "Cd68", "Nos2"]
markers["macrophage_m2"] = ["Ptprc", "Itgam", "Adgre1", "Mrc1", "Arg1"]
markers["mesenchymal"] = ["Fgf10", "Wnt10a", "Wnt2"]
# useless, just for sanity check
# markers["sv"] = ["Pax2"]

# log-normalize the cells
adata_merfish = anndata.concat([adata_merfish_pp, adata_merfish_ppc])
sc.pp.normalize_total(adata_merfish, 1e4)
sc.pp.log1p(adata_merfish)

# In[]
# -----------------------------------------------------------------------------------------
#
# NOTE: Downstream 1: Visualize the predicted labels
#
# -----------------------------------------------------------------------------------------
# Visualize using joint scANVI embedding
# 1. PP
adata_ref_pp = anndata.concat([adata_ref, adata_merfish_pp])
adata_ref_pp.obs.loc[adata_merfish_pp.obs.index,"annot"] = [x for x in adata_merfish_pp.obs[SCANVI_PREDICTIONS_KEY].values]
# load scanvi model and obtain embedding
scanvi_query_pp = scvi.model.SCANVI.load(res_dir + "scanvi_pp/", adata = adata_ref_pp)
adata_ref_pp.obsm[SCANVI_LATENT_KEY] = scanvi_query_pp.get_latent_representation(adata_ref_pp)
sc.pp.neighbors(adata_ref_pp, use_rep=SCANVI_LATENT_KEY)
sc.tl.umap(adata_ref_pp)

# 2. PPC
adata_ref_ppc = anndata.concat([adata_ref, adata_merfish_ppc])
adata_ref_ppc.obs.loc[adata_merfish_ppc.obs.index,"annot"] = [x for x in adata_merfish_ppc.obs[SCANVI_PREDICTIONS_KEY].values]
scanvi_query_ppc = scvi.model.SCANVI.load(res_dir + "scanvi_ppc/", adata = adata_ref_ppc)
adata_ref_ppc.obsm[SCANVI_LATENT_KEY] = scanvi_query_ppc.get_latent_representation(adata_ref_ppc)
sc.pp.neighbors(adata_ref_ppc, use_rep=SCANVI_LATENT_KEY)
sc.tl.umap(adata_ref_ppc)

# Visualization using predicted annotation
sc.pl.umap(adata_ref_pp, color=["sample", "annot"], frameon=False, ncols=1)
utils.plot_latent(adata_ref_pp.obsm["X_umap"], annos = np.array([x for x in adata_ref_pp.obs["annot"].values]), batches = np.array([x for x in adata_ref_pp.obs["sample"].values]), mode = "separate", figsize = (15, 35), save = res_dir + "scanvi_pp/scanvi_umap.pdf")

sc.pl.umap(adata_ref_ppc, color=["sample", "annot"], frameon=False, ncols=1)
utils.plot_latent(adata_ref_ppc.obsm["X_umap"], annos = np.array([x for x in adata_ref_ppc.obs["annot"].values]), batches = np.array([x for x in adata_ref_ppc.obs["sample"].values]), mode = "separate", figsize = (15, 35), save = res_dir + "scanvi_ppc/scanvi_umap.pdf")

# In[]
# Visualize using UMAP embedding
# calculate the umap embedding
sc.tl.pca(adata_merfish, n_comps = 100)
sc.pp.neighbors(adata_merfish, use_rep = "X_pca")
sc.tl.umap(adata_merfish, min_dist = 0.1)

# Visualization using predicted annotation
sc.pl.umap(adata_merfish, color=["sample", SCANVI_PREDICTIONS_KEY], frameon=False, ncols=1)
# plot cell type on UMAP space
# plot cell type 
fig = utils.plot_by_batch(adata_merfish.obsm["X_umap"], annos = np.array([x for x in adata_merfish.obs[SCANVI_PREDICTIONS_KEY].values]), 
                          batches = np.array([x for x in adata_merfish.obs["sample"].values]), figsize = (10, 4), ncols = 2, s = 5, alpha = 0.4, markerscale = 5)
fig.savefig(res_dir + "umap_scanvi_merfish.pdf", bbox_inches = "tight")

fig = utils.plot_embeds(adata_merfish_pp.obsm["spatial"], annos = adata_merfish_pp.obs[[SCANVI_PREDICTIONS_KEY]], figsize = (25, 13), s = 1, alpha = 0.4, markerscale = 15)
fig.savefig(res_dir + "spatial_scanvi_merfish_pp.pdf", bbox_inches = "tight")
fig = utils.plot_embeds(adata_merfish_ppc.obsm["spatial"], annos = adata_merfish_ppc.obs[[SCANVI_PREDICTIONS_KEY]], figsize = (25, 13), s = 1, alpha = 0.4, markerscale = 15)
fig.savefig(res_dir + "spatial_scanvi_merfish_ppc.pdf", bbox_inches = "tight")


# In[]
# -----------------------------------------------------------------------------------------
#
# NOTE: Downstream 2: Plot marker gene expression heatmap on UMAP and spatial
#
# -----------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
for ct, marker in markers.items():
    nrows = len(marker)
    fig = plt.figure(figsize = (10 * 2, 7 * nrows))
    ax = fig.subplots(nrows = nrows, ncols = 2)
    # adata_intact is already log-normalized
    for idx, gene in enumerate(marker):
        vmax = np.max(adata_merfish[:, gene].X)
        vmin = 0.0
        ax[idx,0] = sc.pl.umap(adata_merfish[adata_merfish.obs["sample"] == "Merfish-PP15", :], color = gene, color_map = utils.SUPER_MAGMA, ax = ax[idx,0], show = False, vmin = vmin, vmax = vmax)
        ax[idx,1] = sc.pl.umap(adata_merfish[adata_merfish.obs["sample"] == "Merfish-PPC15", :], color = gene, color_map = utils.SUPER_MAGMA, ax = ax[idx,1], show = False, vmin = vmin, vmax = vmax)
        ax[idx,0].set_title(f"{gene} PP15", fontsize = 35)
        ax[idx,1].set_title(f"{gene} PPC15", fontsize = 35)

    fig.tight_layout()
    plt.show()
    fig.savefig(res_dir + f"{ct}_marker_umap.pdf", bbox_inches = "tight")

X_umap = adata_merfish.obsm["X_umap"].copy()
adata_merfish.obsm["X_umap"] = adata_merfish.obsm["spatial"]
import matplotlib.pyplot as plt
for ct, marker in markers.items():
    nrows = len(marker)
    fig = plt.figure(figsize = (10 * 2, 7 * nrows))
    ax = fig.subplots(nrows = nrows, ncols = 2)
    # adata_intact is already log-normalized
    for idx, gene in enumerate(marker):
        vmax = np.max(adata_merfish[:, gene].X)
        vmin = 0.0
        ax[idx,0] = sc.pl.umap(adata_merfish[adata_merfish.obs["sample"] == "Merfish-PP15", :], color = gene, color_map = utils.SUPER_MAGMA, ax = ax[idx,0], show = False, vmin = vmin, vmax = vmax)
        ax[idx,1] = sc.pl.umap(adata_merfish[adata_merfish.obs["sample"] == "Merfish-PPC15", :], color = gene, color_map = utils.SUPER_MAGMA, ax = ax[idx,1], show = False, vmin = vmin, vmax = vmax)
        ax[idx,0].set_title(f"{gene} PP15", fontsize = 35)
        ax[idx,1].set_title(f"{gene} PPC15", fontsize = 35)

    fig.tight_layout()
    plt.show()
    fig.savefig(res_dir + f"{ct}_marker_spatial.pdf", bbox_inches = "tight")

adata_merfish.obsm["X_umap"] = X_umap


# In[]
# -----------------------------------------------------------------------------------------
#
# NOTE: Downtream 3: Analysis on cell type composition
#
# -----------------------------------------------------------------------------------------
# Count cell type composition
annos_immune_pp = np.array([x for x in adata_merfish_pp.obs[SCANVI_PREDICTIONS_KEY]])
annos_immune_pp = np.where((annos_immune_pp == "Endothelial") | (annos_immune_pp == "SV"), "Other", annos_immune_pp)
annos_immune_ppc = np.array([x for x in adata_merfish_ppc.obs[SCANVI_PREDICTIONS_KEY]])
annos_immune_ppc = np.where((annos_immune_ppc == "Endothelial") | (annos_immune_ppc == "SV"), "Other", annos_immune_ppc)

ct_pp, ct_count_pp = np.unique(annos_immune_pp, return_counts = True)
ct_pct_pp = ct_count_pp/np.sum(ct_count_pp)
ct_ppc, ct_count_ppc = np.unique(annos_immune_ppc, return_counts = True)
ct_pct_ppc = ct_count_ppc/np.sum(ct_count_ppc)

unique_ct = ["Luminal", "Basal", "Monocytes", "Lymphoid", "Macrophage", "Luminal (Spink1+)", "Mesenchymal", "Club epithelia", "Other"]
count_ct = pd.DataFrame(data = 0.0, index = unique_ct, columns = ["PP (15wk)", "PPC (15wk)"])
count_ct.loc[ct_pp, "PP (15wk)"] = ct_count_pp
count_ct.loc[ct_pp, "PPC (15wk)"] = ct_count_ppc

pct_ct = pd.DataFrame(data = 0.0, index = unique_ct, columns = ["PP (15wk)", "PPC (15wk)"])
pct_ct.loc[ct_pp, "PP (15wk)"] = ct_pct_pp
pct_ct.loc[ct_pp, "PPC (15wk)"] = ct_pct_ppc

count_ct.to_csv(res_dir + "count_ct.csv")
pct_ct.to_csv(res_dir + "pct_ct.csv")

# NOTE: Combine the immune cell, for ease of pie chart plot
count_ct.loc["Immune", :] = count_ct.loc["Lymphoid", :].values + count_ct.loc["Macrophage", :].values + count_ct.loc["Monocytes", :].values
count_ct = count_ct.drop(["Lymphoid", "Macrophage", "Monocytes"], axis = 0)
pct_ct.loc["Immune", :] = pct_ct.loc["Lymphoid", :].values + pct_ct.loc["Macrophage", :].values + pct_ct.loc["Monocytes", :].values
pct_ct = pct_ct.drop(["Lymphoid", "Macrophage", "Monocytes"], axis = 0)

plt.rcParams["font.size"] = 10
fig = plt.figure(figsize = (20,10))
axs = fig.subplots(nrows = 1, ncols = 2)
axs[0].pie(pct_ct["PP (15wk)"].values.squeeze(), labels = pct_ct.index.values.squeeze(), autopct='%.1f%%', pctdistance = 0.8, labeldistance = 1.1)
axs[1].pie(pct_ct["PPC (15wk)"].values.squeeze(), labels = pct_ct.index.values.squeeze(), autopct='%.1f%%', pctdistance = 0.8, labeldistance = 1.1)
axs[0].set_title("PP 15wk")
axs[1].set_title("PPC 15wk")
fig.suptitle("Cell percentage", fontsize = 30)
fig.savefig(res_dir + "percent_cell.png", bbox_inches = "tight")

# In[]
# Plot the locations of the immune cells, include lymphoid, macrophage, monocytes
immune_cmap = ListedColormap(["#FAFAFA", "#98DF8A", "#AEC7E8", "#FF9896"])

annos_immune_pp = np.array([x for x in adata_merfish_pp.obs[SCANVI_PREDICTIONS_KEY]])
annos_immune_pp = np.where((annos_immune_pp != "Lymphoid") & (annos_immune_pp != "Macrophage") & (annos_immune_pp != "Monocytes"), " ", annos_immune_pp)
annos_immune_pp = pd.DataFrame(columns = [SCANVI_PREDICTIONS_KEY], data = annos_immune_pp[:,None])
fig = utils.plot_embeds(adata_merfish_pp.obsm["spatial"], annos = annos_immune_pp, figsize = (25, 13), s = 1, alpha = 0.7, markerscale = 15, colormap = immune_cmap)
fig.savefig(res_dir + "spatial_scanvi_immune_pp.pdf", bbox_inches = "tight")

annos_immune_ppc = np.array([x for x in adata_merfish_ppc.obs[SCANVI_PREDICTIONS_KEY]])
annos_immune_ppc = np.where((annos_immune_ppc != "Lymphoid") & (annos_immune_ppc != "Macrophage") & (annos_immune_ppc != "Monocytes"), " ", annos_immune_ppc)
annos_immune_ppc = pd.DataFrame(columns = [SCANVI_PREDICTIONS_KEY], data = annos_immune_ppc[:,None])
fig = utils.plot_embeds(adata_merfish_ppc.obsm["spatial"], annos = annos_immune_ppc, figsize = (25, 13), s = 1, alpha = 0.7, markerscale = 15, colormap = immune_cmap)
fig.savefig(res_dir + "spatial_scanvi_immune_ppc.pdf", bbox_inches = "tight")

# In[]
# --------------------------------------------------------------------------------
#
# NOTE: Downstream 4: Immune sub-type detection in space, assuming the lymphoid macrophage annotations are correct
#
# --------------------------------------------------------------------------------
if not os.path.exists(res_dir + "immune_subtypes/"):
    os.makedirs(res_dir + "immune_subtypes/")

three_cat = ListedColormap(["#FAFAFA", "#FF9896", "#AEC7E8"])

# NOTE: the markers used for immune subtype detection should be the same as the ones used in scRNA-seq and visium
markers_lymphoid = {}
# markers_lymphoid["lymphoid"] = ["Cd3e"]#, "Ms4a1", "Klrb1c"]
markers_lymphoid["CD4_T"] = ["Cd4"]
markers_lymphoid["CD8_T"] = ["Cd8a"]
markers_lymphoid["cytotoxicity T"] = ["Gzmb", "Prf1"]
markers_lymphoid["exhausted T"] = ["Tigit", "Havcr2"]
markers_lymphoid["NK"] = ["Klrd1"]
markers_lymphoid["B"] = ["Ms4a1"]

markers_macrophage = {}
# markers_macrophage["macrophage"] = ["Itgam", "Adgre1"]
# markers_macrophage["macrophage_m1"] = ["Cd80", "Cd86", "Nos2", "Tnf", "Il1b", "Il6", "Cxcl10"]
# markers_macrophage["macrophage_m2"] = ["Cd163", "Arg1", "Mrc1", "Il10"]
markers_macrophage["macrophage_m1"] = ["Cd80", "Cd86"]
markers_macrophage["macrophage_m2"] = ["Cd163", "Arg1", "Mrc1"]

markers_immune = {"Lymphoid": markers_lymphoid, "Macrophage": markers_macrophage}

# NOTE: calculate the percentage of immune cells expressing subtype markers
for subtype, markers in markers_immune.items():
    # Lymphoid and Macrophage
    for ct in markers.keys():
        # PP15
        marker_exprs = pd.DataFrame(index = adata_merfish_pp.obs.index.values)
        for gene in markers[ct]:
            # NOTE: the expression is already log-normed, see loading
            marker_expr = adata_merfish_pp[:,gene].X.toarray().squeeze()
            # binarize the expression value
            marker_bin = np.where(marker_expr > 0, "2. expr", "3. not expr")
            # set the expression of remaining cell types to be 0, remove confounders
            marker_bin[(adata_merfish_pp.obs[SCANVI_PREDICTIONS_KEY] != subtype).values.squeeze()] = "1. other"
            print(np.unique(marker_bin, return_counts = True))
            marker_exprs[gene] = marker_bin
        # Plot the spatial location of markers
        fig = utils.plot_embeds(embed = adata_merfish_pp.obsm["spatial"], annos = marker_exprs, colormap = three_cat, alpha = 1, s = 20, markerscale = 3)
        fig.suptitle(ct, fontsize = 25)
        fig.savefig(res_dir + f"immune_subtypes/{ct}_markers_pp.pdf", bbox_inches = "tight")
        marker_exprs.to_csv(res_dir + f"immune_subtypes/bin_{ct}_markers_pp.csv")

        # PPC15
        marker_exprs = pd.DataFrame(index = adata_merfish_ppc.obs.index.values)
        for gene in markers[ct]:
            # NOTE: the expression is already log-normed
            marker_expr = adata_merfish_ppc[:,gene].X.toarray().squeeze()
            # binarize the expression value
            marker_bin = np.where(marker_expr > 0, "2. expr", "3. not expr")
            # set the expression of remaining cell types to be 0, remove confounders
            marker_bin[(adata_merfish_ppc.obs[SCANVI_PREDICTIONS_KEY] != subtype).values.squeeze()] = "1. other"
            print(np.unique(marker_bin, return_counts = True))
            marker_exprs[gene] = marker_bin
        # Plot the spatial location of markers
        fig = utils.plot_embeds(embed = adata_merfish_ppc.obsm["spatial"], annos = marker_exprs, colormap = three_cat, alpha = 1, s = 20, markerscale = 3)
        fig.suptitle(ct, fontsize = 25)
        fig.savefig(res_dir + f"immune_subtypes/{ct}_markers_ppc.pdf", bbox_inches = "tight")
        marker_exprs.to_csv(res_dir + f"immune_subtypes/bin_{ct}_markers_ppc.csv")

# In[]
# NOTE: plot the pie chart of subtype composition, Lymphoid and Macrophage
for subtype, markers in markers_immune.items():
    for ct in markers.keys():
        pct_pos_df = pd.DataFrame(columns = ["PP 15wk", "PPC 15wk"], index = ["Pos", "Neg"], data = 0)
        marker_exprs = pd.read_csv(res_dir + f"immune_subtypes/bin_{ct}_markers_pp.csv", index_col = 0)
        # number of lymphoid cells
        num_total = (adata_merfish_pp.obs[SCANVI_PREDICTIONS_KEY] == subtype).values.squeeze().sum()
        # NOTE: here we use "or" instead of "and" in scRNA-seq, because "and" will filter out almost all cells in multi-marker cases
        # Merfish data is extremely shallow
        idx_pos = np.zeros((marker_exprs.shape[0],))
        for gene in marker_exprs.columns:
            idx_pos += np.where(marker_exprs[gene] == "2. expr", 1, 0) 
        num_pos = np.sum(idx_pos > 0)
        pct_pos_df.loc["Pos", "PP 15wk"] = num_pos
        pct_pos_df.loc["Neg", "PP 15wk"] = num_total - num_pos

        marker_exprs = pd.read_csv(res_dir + f"immune_subtypes/bin_{ct}_markers_ppc.csv", index_col = 0)
        # number of lymphoid cells
        num_total = (adata_merfish_ppc.obs[SCANVI_PREDICTIONS_KEY] == subtype).values.squeeze().sum()
        idx_pos = np.zeros((marker_exprs.shape[0],))
        for gene in marker_exprs.columns:
            idx_pos += np.where(marker_exprs[gene] == "2. expr", 1, 0) 
        num_pos = np.sum(idx_pos > 0)
        pct_pos_df.loc["Pos", "PPC 15wk"] = num_pos
        pct_pos_df.loc["Neg", "PPC 15wk"] = num_total - num_pos
        pct_pos_df.to_csv(res_dir + f"immune_subtypes/{ct}_pct.csv")


        plt.rcParams["font.size"] = 10
        fig = plt.figure(figsize = (20,10))
        axs = fig.subplots(nrows = 1, ncols = 2)
        axs[0].pie(pct_pos_df["PP 15wk"].values.squeeze(), labels = pct_pos_df.index.values.squeeze(), autopct='%.1f%%', pctdistance = 0.8, labeldistance = 1.1)
        axs[1].pie(pct_pos_df["PPC 15wk"].values.squeeze(), labels = pct_pos_df.index.values.squeeze(), autopct='%.1f%%', pctdistance = 0.8, labeldistance = 1.1)
        axs[0].set_title("PP 15wk")
        axs[1].set_title("PPC 15wk")
        fig.suptitle(f"{ct} percentage", fontsize = 30)
        fig.savefig(res_dir + f"immune_subtypes/{ct}_pct.png", bbox_inches = "tight")

# %%
# # -----------------------------------------------------------------------------------------
# #
# # Label transfer -- scArches
# # align the query dataset with the reference dataset
# #
# # -----------------------------------------------------------------------------------------
# scvi_ref_path = res_dir + "scvi_reference"
# scvi.model.SCVI.prepare_query_anndata(adata_merfish_pp, scvi_ref_path)
# scvi_query_pp = scvi.model.SCVI.load_query_data(adata_merfish_pp, scvi_ref_path)
# # NOTE: weight_decay to 0 to prevent the model from updating the latent embedding of the reference dataset
# scvi_query_pp.train(max_epochs=200, plan_kwargs={"weight_decay": 0.0})


# scvi.model.SCVI.prepare_query_anndata(adata_merfish_ppc, scvi_query_pp)
# scvi_query_ppc = scvi.model.SCVI.load_query_data(adata_merfish_ppc, scvi_query_pp)
# # NOTE: weight_decay to 0 to prevent the model from updating the latent embedding of the reference dataset
# scvi_query_ppc.train(max_epochs=200, plan_kwargs={"weight_decay": 0.0})

# scvi_query_ppc.save(res_dir + "scvi_final", overwrite = True)


# # In[]
# SCVI_LATENT_KEY = "X_scVI"
# adata_full = anndata.concat([adata_ref, adata_merfish_pp, adata_merfish_ppc])
# adata_full.obsm[SCVI_LATENT_KEY] = scvi_query_ppc.get_latent_representation(adata_full)

# sc.pp.neighbors(adata_full, use_rep=SCVI_LATENT_KEY)
# sc.tl.umap(adata_full)

# sc.pl.umap(adata_full, color=["sample", "annot"], frameon=False, ncols=1)
# utils.plot_latent(adata_full.obsm["X_umap"], annos = np.array([x for x in adata_full.obs["annot"].values]), batches = np.array([x for x in adata_full.obs["sample"].values]), mode = "separate", figsize = (10, 25))
