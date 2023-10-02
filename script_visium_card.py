# In[]
import sys, os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs


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

adata_visium = {"pp12": adata_pp12, "pp18": adata_pp18, "ppc12": adata_ppc12, "ppc18": adata_ppc18}

celltypes_categories = np.array(['Basal', 'Club epithelia (Luminal)', 'Endothelial',
                        'Hillock epithelia (Basal)', 'Luminal',
                        'Lymphoid (Total immune)', 'Macrophage (Myeloid, Total immune)',
                        'Mesenchymal', 'Monocytes (Myeloid, Total immune)', 'SV',
                        'Spink1+ (Luminal)'])
# In[]
deconv_score_pp12 = pd.read_csv("results_card/deconv_score_pp12.csv", sep = "\t").loc[adata_visium["pp12"].obs.index.values,:]
deconv_score_pp18 = pd.read_csv("results_card/deconv_score_pp18.csv", sep = "\t").loc[adata_visium["pp18"].obs.index.values,:]
deconv_score_ppc12 = pd.read_csv("results_card/deconv_score_ppc12.csv", sep = "\t").loc[adata_visium["ppc12"].obs.index.values,:]
deconv_score_ppc18 = pd.read_csv("results_card/deconv_score_ppc18.csv", sep = "\t").loc[adata_visium["ppc18"].obs.index.values,:]

adata_visium["pp12"].obs = pd.concat([adata_visium["pp12"].obs, deconv_score_pp12], axis = 1)
adata_visium["pp18"].obs = pd.concat([adata_visium["pp18"].obs, deconv_score_pp18], axis = 1)
adata_visium["ppc12"].obs = pd.concat([adata_visium["ppc12"].obs, deconv_score_ppc12], axis = 1)
adata_visium["ppc18"].obs = pd.concat([adata_visium["ppc18"].obs, deconv_score_ppc18], axis = 1)

adata_visium["pp12"].obs["label_card"] = deconv_score_pp12.columns.values[np.argmax(deconv_score_pp12.values, axis = 1)]
adata_visium["pp12"].obs["label_card"] = pd.Categorical(adata_visium["pp12"].obs["label_card"], categories=celltypes_categories)

adata_visium["pp18"].obs["label_card"] = deconv_score_pp18.columns.values[np.argmax(deconv_score_pp18.values, axis = 1)]
adata_visium["pp18"].obs["label_card"] = pd.Categorical(adata_visium["pp18"].obs["label_card"], categories=celltypes_categories)

adata_visium["ppc12"].obs["label_card"] = deconv_score_ppc12.columns.values[np.argmax(deconv_score_ppc12.values, axis = 1)]
adata_visium["ppc12"].obs["label_card"] = pd.Categorical(adata_visium["ppc12"].obs["label_card"], categories=celltypes_categories)

adata_visium["ppc18"].obs["label_card"] = deconv_score_ppc18.columns.values[np.argmax(deconv_score_ppc18.values, axis = 1)]
adata_visium["ppc18"].obs["label_card"] = pd.Categorical(adata_visium["ppc18"].obs["label_card"], categories=celltypes_categories)

# In[]
# Load data
sample = "ppc18"
adata_vis = adata_visium[sample]
result_dir = "./results_card/results_" + sample + "/"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "label_card", ax = ax)
fig.savefig(result_dir + f"label_{sample}.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"label_{sample}.png", bbox_inches = "tight", dpi = 300)

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Spink1+ (Luminal)", ax = ax)
fig.savefig(result_dir + f"label_{sample}_spink1.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"label_{sample}_spink1.png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Luminal", ax = ax)
fig.savefig(result_dir + f"label_{sample}_Luminal.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"label_{sample}_Luminal.png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Basal", ax = ax)
fig.savefig(result_dir + f"label_{sample}_Basal.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"label_{sample}_Basal.png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Club epithelia (Luminal)", ax = ax)
fig.savefig(result_dir + f"label_{sample}_Club_epithelia_(Luminal).pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"label_{sample}_Club_epithelia_(Luminal).png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Endothelial", ax = ax)
fig.savefig(result_dir + f"label_{sample}_Endothelial.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"label_{sample}_Endothelial.png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Hillock epithelia (Basal)", ax = ax)
fig.savefig(result_dir + f"label_{sample}_Hillock_epithelia_(Basal).pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"label_{sample}_Hillock_epithelia_(Basal).png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Lymphoid (Total immune)", ax = ax)
fig.savefig(result_dir + f"label_{sample}_Lymphoid_(Total immune).pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"label_{sample}_Lymphoid_(Total immune).png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Macrophage (Myeloid, Total immune)", ax = ax)
fig.savefig(result_dir + f"label_{sample}_Macrophage_(Myeloid, Total immune).pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"label_{sample}_Macrophage_(Myeloid, Total immune).png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Mesenchymal", ax = ax)
fig.savefig(result_dir + f"label_{sample}_Mesenchymal.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"label_{sample}_Mesenchymal.png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "Monocytes (Myeloid, Total immune)", ax = ax)
fig.savefig(result_dir + f"label_{sample}_Monocytes (Myeloid, Total immune).pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"label_{sample}_Monocytes (Myeloid, Total immune).png", bbox_inches = "tight")

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot()
sc.pl.spatial(adata_vis, img_key = "hires", color = "SV", ax = ax)
fig.savefig(result_dir + f"label_{sample}_SV.pdf", bbox_inches = "tight")
fig.savefig(result_dir + f"label_{sample}_SV.png", bbox_inches = "tight")

# %%
