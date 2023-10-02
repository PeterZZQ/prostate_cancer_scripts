# In[]
import sys, os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from cell2location.utils import select_slide

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

result_dir = "./results_cell2loc/"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

# In[]
# Load Visium data
adata_pp12 = sc.read_visium(path = "spatial_data/Visium_python/225_PP12/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_pp18 = sc.read_visium(path = "spatial_data/Visium_python/1687_PP18/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_ppc12 = sc.read_visium(path = "spatial_data/Visium_python/1161_PPC/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_ppc18 = sc.read_visium(path = "spatial_data/Visium_python/1660_PPC18/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")

adata_pp12.obsm["spatial"] = adata_pp12.obsm["spatial"].astype(np.float64)
adata_pp18.obsm["spatial"] = adata_pp18.obsm["spatial"].astype(np.float64)
adata_ppc12.obsm["spatial"] = adata_ppc12.obsm["spatial"].astype(np.float64)
adata_ppc18.obsm["spatial"] = adata_ppc18.obsm["spatial"].astype(np.float64)

# adata_pp12.var_names_make_unique()
# adata_pp18.var_names_make_unique()
# adata_ppc12.var_names_make_unique()
# adata_ppc18.var_names_make_unique()

# rename gene names from SYMBOL to ENSEMBLE_ID
adata_pp12.var["SYMBOL"] = adata_pp12.var_names
adata_pp12.var.index = adata_pp12.var["gene_ids"].values
adata_pp18.var["SYMBOL"] = adata_pp18.var_names
adata_pp18.var.index = adata_pp18.var["gene_ids"].values
adata_ppc12.var["SYMBOL"] = adata_ppc12.var_names
adata_ppc12.var.index = adata_ppc12.var["gene_ids"].values
adata_ppc18.var["SYMBOL"] = adata_ppc18.var_names
adata_ppc18.var.index = adata_ppc18.var["gene_ids"].values

# In[]
# Load reference scRNA-seq data
adata = sc.read_h5ad("Cell_Ranger_output/adata_seurat.h5ad")

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
adata.obs.loc[adata.obs["seurat_clusters"] == "10", "annotation"] = "SV" # originally SV, Spink1+
adata.obs.loc[adata.obs["seurat_clusters"] == "11", "annotation"] = "Endothelial"
adata.obs.loc[adata.obs["seurat_clusters"] == "12", "annotation"] = "Luminal 2"
adata.obs.loc[adata.obs["seurat_clusters"] == "13", "annotation"] = "Luminal 3"
# Remove Luminal 2 and Luminal 3
adata_ref = adata[(adata.obs["annotation"] != "Luminal 2") & (adata.obs["annotation"] != "Luminal 3"),:]
adata_ref.var["SYMBOL"] = adata_ref.var.index.values
adata_ref.var.index = adata_ref.var["gene_ids"].values
adata_ref.obs["annotation"] = adata_ref.obs["annotation"].astype('category')

# In[]
# preprocessing
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
# filter the object
adata_ref = adata_ref[:, selected].copy()

# In[]
# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='sample',
                        # cell type, covariate used for constructing signatures
                        labels_key='annotation'
                       )

# create the regression model
mod = RegressionModel(adata_ref)
# view anndata_setup as a sanity check
mod.view_anndata_setup()
# train the model
mod.train(max_epochs=250, use_gpu=True)
# ELBO loss curve
mod.plot_history(20)

# In[]
# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)
# Save model
# mod.save(result_dir + "reference_signitures", overwrite=True)
# Save anndata object with results
# adata_ref.write(result_dir + "reference_signitures/sc.h5ad")

# adata_ref = sc.read_h5ad(result_dir + "reference_signitures/sc.h5ad")
# mod = cell2location.models.RegressionModel.load(result_dir + "reference_signitures", adata_ref)
# In[]
# 1. Reconstruction accuracy to assess if there are any issues with inference. 
# This 2D histogram plot should have most observations along a noisy diagonal.
# 
# 2. The estimated expression signatures are distinct from mean expression in 
# each cluster because of batch effects. For scRNA-seq datasets which do not 
# suffer from batch effect (this dataset does), cluster average expression can 
# be used instead of estimating signatures with a model. When this plot is very 
# different from a diagonal plot (e.g. very low values on Y-axis, density everywhere) 
# it indicates problems with signature estimation.
mod.plot_QC()

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
# inf_aver.iloc[0:5, 0:5]

# In[]
# find shared genes and subset both anndata and reference signatures
sample = "ppc18"
if sample == "pp12":
    adata_vis = adata_pp12
elif sample == "pp18":
    adata_vis = adata_pp18
elif sample == "ppc12":
    adata_vis = adata_ppc12
elif sample == "ppc18":
    adata_vis = adata_ppc18

intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# NOTE: the key parameters: N_cells_per_location, detection_alpha
# prepare anndata for cell2location model
N_cells_per_location = 10
detection_alpha = 20
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key=None)
# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    # Typical 10X visium: 5, default: 30
    N_cells_per_location=N_cells_per_location,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=detection_alpha
)
mod.view_anndata_setup()

# In[]
mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True,
         )

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);

# In[]
# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}", overwrite=True)

# mod = cell2location.models.Cell2location.load(result_dir + "cell2loc_map", adata_vis)

# Save anndata object with results
adata_vis.write( result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/sp.h5ad")

# In[]
mod.plot_QC()

# In[]
# Load data
sample = "ppc18"
N_cells_per_location = 30
detection_alpha = 20 
adata_vis = sc.read_h5ad(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/sp.h5ad")

# Extract deconvolution weight Wsf
select = "q05"
if select == "q05":
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
elif select == "mean":
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['means_cell_abundance_w_sf']

Wsf = adata_vis.obs[adata_vis.uns['mod']['factor_names']]
# Use the maximum Wsf
adata_vis.obs["max_posterior"] = Wsf.columns.values[np.argmax(Wsf.values, axis = 1)]
adata_vis.obs["max_posterior"] = pd.Categorical(adata_vis.obs["max_posterior"], categories=adata_ref.obs["annotation"].cat.categories.values)
adata_vis.obs["max_posterior"] = adata_vis.obs["max_posterior"].cat.rename_categories({"Luminal 1": "Luminal"})

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "max_posterior", ax = ax)
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_max_posterior_prediction.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_max_posterior_prediction.png", bbox_inches = "tight", dpi = 300)

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Spink1+ (Luminal)", ax = ax)
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_spink1.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_spink1.png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Luminal 1", ax = ax)
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Luminal.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Luminal.png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Basal", ax = ax)
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Basal.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Basal.png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Club epithelia (Luminal)", ax = ax)
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Club_epithelia_(Luminal).pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Club_epithelia_(Luminal).png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Endothelial", ax = ax)
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Endothelial.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Endothelial.png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Hillock epithelia (Basal)", ax = ax)
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Hillock_epithelia_(Basal).pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Hillock_epithelia_(Basal).png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Lymphoid (Total immune)", ax = ax)
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Lymphoid_(Total immune).pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Lymphoid_(Total immune).png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Macrophage (Myeloid, Total immune)", ax = ax)
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Macrophage_(Myeloid, Total immune).pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Macrophage_(Myeloid, Total immune).png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Mesenchymal", ax = ax)
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Mesenchymal.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Mesenchymal.png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Monocytes (Myeloid, Total immune)", ax = ax)
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Monocytes (Myeloid, Total immune).pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_Monocytes (Myeloid, Total immune).png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "SV", ax = ax)
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_SV.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"cell2loc_map_{sample}_{N_cells_per_location}_{detection_alpha}/{select}_SV.png", bbox_inches = "tight")

# %%
