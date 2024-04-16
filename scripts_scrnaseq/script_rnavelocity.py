# In[]
import numpy as np 
import scanpy as sc 
import scvelo as scv
import matplotlib.pyplot as plt
import phate
import pandas as pd

import sys
sys.path.append(".")
import utils


# In[]
# Read the data
data_dir = "../dataset/data_scrnaseq/Cell_Ranger_output/loom_files/"
seurat_dir = "../dataset/data_scrnaseq/seurat_integration/"

# intact samples
M1416_adata = scv.read(data_dir + "M1416.loom", cache=True)
M1417_adata = scv.read(data_dir + "M1417.loom", cache=True)
M1436_adata = scv.read(data_dir + "M1436.loom", cache=True)
M1437_adata = scv.read(data_dir + "M1437.loom", cache=True)
# castrated sample, not needed
# M1476_adata = scv.read(data_dir + "M1476.loom", cache=True)
# M1477_adata = scv.read(data_dir + "M1477.loom", cache=True)

# assign cell types
adata_intact_seurat = sc.read_h5ad(seurat_dir + "adata_intact_seurat.h5ad")
# adata_castrated_seurat = sc.read_h5ad(seurat_dir + "adata_castrated_seurat.h5ad")

meta_m1416 = adata_intact_seurat.obs.loc[adata_intact_seurat.obs["sample"] == "M1416-PP18",:]
barcodes = ["M1416:" + x.split("-")[0] + "x" for x in meta_m1416.index.values.squeeze()]
meta_m1416.index = barcodes
# filter loom file accordingly
M1416_adata = M1416_adata[barcodes,:]
# M1416_adata.obs["annot"] = "Noisy"
M1416_adata.obs.loc[barcodes, "annot"] = meta_m1416["annot"].values
M1416_adata.obs["sample"] = "M1416-PP18"
M1416_adata.obsm = adata_intact_seurat[adata_intact_seurat.obs["sample"] == "M1416-PP18",:].obsm.copy()

meta_m1417 = adata_intact_seurat.obs.loc[adata_intact_seurat.obs["sample"] == "M1417-PP12",:]
barcodes = ["M1417:" + x.split("-")[0] + "x" for x in meta_m1417.index.values.squeeze()]
meta_m1417.index = barcodes
# filter loom file accordingly
M1417_adata = M1417_adata[barcodes,:]
# M1417_adata.obs["annot"] = "Noisy"
M1417_adata.obs.loc[barcodes,"annot"] = meta_m1417["annot"].values
M1417_adata.obs["sample"] = "M1417-PP12"
M1417_adata.obsm = adata_intact_seurat[adata_intact_seurat.obs["sample"] == "M1417-PP12",:].obsm.copy()

meta_m1436 = adata_intact_seurat.obs.loc[adata_intact_seurat.obs["sample"] == "M1436-PPC12",:]
barcodes = ["M1436:" + x.split("-")[0] + "x" for x in meta_m1436.index.values.squeeze()]
meta_m1436.index = barcodes
# filter loom file accordingly
M1436_adata = M1436_adata[barcodes,:]
# M1436_adata.obs["annot"] = "Noisy"
M1436_adata.obs.loc[barcodes,"annot"] = meta_m1436["annot"].values
M1436_adata.obs["sample"] = "M1436-PPC12"
M1436_adata.obsm = adata_intact_seurat[adata_intact_seurat.obs["sample"] == "M1436-PPC12",:].obsm.copy()

meta_m1437 = adata_intact_seurat.obs.loc[adata_intact_seurat.obs["sample"] == "M1437-PPC18",:]
barcodes = ["M1437:" + x.split("-")[0] + "x" for x in meta_m1437.index.values.squeeze()]
meta_m1437.index = barcodes
# filter loom file accordingly
M1437_adata = M1437_adata[barcodes,:]
# M1437_adata.obs["annot"] = "Noisy"
M1437_adata.obs.loc[barcodes,"annot"] = meta_m1437["annot"].values
M1437_adata.obs["sample"] = "M1437-PPC18"
M1437_adata.obsm = adata_intact_seurat[adata_intact_seurat.obs["sample"] == "M1437-PPC18",:].obsm.copy()

# In[]
# Preprocessing
M1416_adata.var_names_make_unique()
M1417_adata.var_names_make_unique()
M1436_adata.var_names_make_unique()
M1437_adata.var_names_make_unique()

adata = sc.concat([M1416_adata, M1417_adata, M1436_adata, M1437_adata])

# # redo umap
# scv.pp.filter_and_normalize(adata, min_shared_counts = 20, n_top_genes = 2000)
# sc.pp.neighbors(adata)
# sc.tl.umap(adata)
# M1416_adata.obsm["X_umap"] = adata[M1416_adata.obs.index.values].obsm["X_umap"].toarray()
# M1417_adata.obsm["X_umap"] = adata[M1417_adata.obs.index.values].obsm["X_umap"].toarray()
# M1436_adata.obsm["X_umap"] = adata[M1436_adata.obs.index.values].obsm["X_umap"].toarray()
# M1437_adata.obsm["X_umap"] = adata[M1437_adata.obs.index.values].obsm["X_umap"].toarray()

# check proportion of spliced and unspliced counts
scv.pl.proportions(adata)


# In[]
# Plot luminal markers
# Mesenchymal markers: Lum, Vim
# Luminal markers: Ar, Krt8, Cd24a, Krt18
# Stem cell markers: Cd44, Epcam, Sox2, Pou5f1, Nanog, Letm1 
# non expression of Sox2, Pou5f1, Nanog
with plt.rc_context({"figure.figsize": (7, 5), "figure.dpi": (300), "font.size": 15}):   
    fig = sc.pl.umap(M1416_adata, color = ["Ar", "Krt8", "Cd24a", "Krt18", "Lum", "Vim", "Cd44", "Epcam", "Letm1"], return_fig = True, color_map = utils.SUPER_MAGMA, ncols = 2)
    fig.savefig("../results_scrnaseq/results_emt/markers/m1416_pp18.png", bbox_inches = "tight")

with plt.rc_context({"figure.figsize": (7, 5), "figure.dpi": (300), "font.size": 15}):
    fig = sc.pl.umap(M1417_adata, color = ["Ar", "Krt8", "Cd24a", "Krt18", "Lum", "Vim", "Cd44", "Epcam", "Letm1"], return_fig = True, color_map = utils.SUPER_MAGMA, ncols = 2)
    fig.savefig("../results_scrnaseq/results_emt/markers/m1417_pp12.png", bbox_inches = "tight")

with plt.rc_context({"figure.figsize": (7, 5), "figure.dpi": (300), "font.size": 15}):
    fig = sc.pl.umap(M1436_adata, color = ["Ar", "Krt8", "Cd24a", "Krt18", "Lum", "Vim", "Cd44", "Epcam", "Letm1"], return_fig = True, color_map = utils.SUPER_MAGMA, ncols = 2)
    fig.savefig("../results_scrnaseq/results_emt/markers/m1436_ppc12.png", bbox_inches = "tight")

with plt.rc_context({"figure.figsize": (7, 5), "figure.dpi": (300), "font.size": 15}):
    fig = sc.pl.umap(M1437_adata, color = ["Ar", "Krt8", "Cd24a", "Krt18", "Lum", "Vim", "Cd44", "Epcam", "Letm1"], return_fig = True, color_map = utils.SUPER_MAGMA, ncols = 2)
    fig.savefig("../results_scrnaseq/results_emt/markers/m1437_ppc18.png", bbox_inches = "tight")


# In[]

# re-think the calculation step of RNA velocity, the RNA velocity should only exist in the trajectory, we should narrow down our focus
# why is the selection of highly variable genes affect the result that much??
# Select only the EMT population
M1416_adata_emt = M1416_adata[(M1416_adata.obs["annot"] == "Basal") | (M1416_adata.obs["annot"] == "Luminal") | (M1416_adata.obs["annot"] == "Mesenchymal")]
M1417_adata_emt = M1417_adata[(M1417_adata.obs["annot"] == "Basal") | (M1417_adata.obs["annot"] == "Luminal") | (M1417_adata.obs["annot"] == "Mesenchymal")]
M1436_adata_emt = M1436_adata[(M1436_adata.obs["annot"] == "Basal") | (M1436_adata.obs["annot"] == "Luminal") | (M1436_adata.obs["annot"] == "Mesenchymal")]
M1437_adata_emt = M1437_adata[(M1437_adata.obs["annot"] == "Basal") | (M1437_adata.obs["annot"] == "Luminal") | (M1437_adata.obs["annot"] == "Mesenchymal")]

# NOTE: filtering and hvg based on the EMT population, hard to find overlap if there are too few genes
scv.pp.filter_and_normalize(M1416_adata_emt, min_shared_counts=20, n_top_genes = 7000)
scv.pp.filter_and_normalize(M1417_adata_emt, min_shared_counts=20, n_top_genes = 7000)
scv.pp.filter_and_normalize(M1436_adata_emt, min_shared_counts=20, n_top_genes = 7000)
scv.pp.filter_and_normalize(M1437_adata_emt, min_shared_counts=20, n_top_genes = 7000)

# NOTE: filtering and hvg based on the whole population, 
# the velocity should not be calculated on the whole population as only emt is our focus
# filter and normalize the genes, and select only the top variable genes
# scv.pp.filter_and_normalize(M1416_adata, min_shared_counts=20, n_top_genes=2000)
# scv.pp.filter_and_normalize(M1417_adata, min_shared_counts=20, n_top_genes=2000)
# scv.pp.filter_and_normalize(M1436_adata, min_shared_counts=20, n_top_genes=2000)
# scv.pp.filter_and_normalize(M1437_adata, min_shared_counts=20, n_top_genes=2000)

# scv.pp.filter_and_normalize(M1416_adata, min_shared_counts=20, n_top_genes=5000)
# scv.pp.filter_and_normalize(M1417_adata, min_shared_counts=20, n_top_genes=5000)
# scv.pp.filter_and_normalize(M1436_adata, min_shared_counts=20, n_top_genes=5000)
# scv.pp.filter_and_normalize(M1437_adata, min_shared_counts=20, n_top_genes=5000)


# In[]
# calculate phate embedding on emt trajectory, select only EMT as phate is a trajectory illustration algorithm
# with 7000 hvg, there are 2605 overlap hvg
adata_emt = sc.concat([M1416_adata_emt, M1417_adata_emt, M1436_adata_emt, M1437_adata_emt], join = "inner")
phate_operator = phate.PHATE(n_jobs=-2)
X_phate_intact = phate_operator.fit_transform(adata_emt.X)
X_phate_intact = pd.DataFrame(data = X_phate_intact, index = adata_emt.obs.index, columns = ["X_phate1", "X_phate2"])
adata_emt.obsm["X_phate"] = X_phate_intact.values

# In[]
fig = plt.figure(figsize = (30, 20))
ax = fig.subplots(nrows = 2, ncols = 2)
plt.rcParams["font.size"] = 25
ax[0,0] = sc.pl.scatter(adata_emt[adata_emt.obs["sample"] == "M1417-PP12", :], basis = "phate", color = ["annot"], ax = ax[0,0], show = False)
ax[0,1] = sc.pl.scatter(adata_emt[adata_emt.obs["sample"] == "M1416-PP18", :], basis = "phate", color = ["annot"], ax = ax[0,1], show = False)
ax[1,0] = sc.pl.scatter(adata_emt[adata_emt.obs["sample"] == "M1436-PPC12", :], basis = "phate", color = ["annot"], ax = ax[1,0], show = False)
ax[1,1] = sc.pl.scatter(adata_emt[adata_emt.obs["sample"] == "M1437-PPC18", :], basis = "phate", color = ["annot"], ax = ax[1,1], show = False)
ax[0,0].set_title("PP-12wk", fontsize = 35)
ax[0,1].set_title("PP-18wk", fontsize = 35)
ax[1,0].set_title("PPC-12wk", fontsize = 35)
ax[1,1].set_title("PPC-18wk", fontsize = 35)

fig.tight_layout()
plt.show()
fig.savefig("../results_scrnaseq/annot_intact/" + "annotation_phate.png", bbox_inches = "tight", dpi = 150)

# In[]
# ---------------------------------------------------------------------------------------------
#
# RNA velocity inference on emt population [should not be calculated on the trajectory]
#
# ---------------------------------------------------------------------------------------------

# NOTE: Option 1: calculate the RNA velocity for each condition independently, 
# the hvgs that is used for the calculated is also selected independently (7000 unshared genes)
# to validate if the independent running shows the consistent velocity result
# and the result indeed shows strong consistency
M1416_adata_emt.obsm["X_phate"] = adata_emt[M1416_adata_emt.obs.index.values,:].obsm["X_phate"]
M1417_adata_emt.obsm["X_phate"] = adata_emt[M1417_adata_emt.obs.index.values,:].obsm["X_phate"]
M1436_adata_emt.obsm["X_phate"] = adata_emt[M1436_adata_emt.obs.index.values,:].obsm["X_phate"]
M1437_adata_emt.obsm["X_phate"] = adata_emt[M1437_adata_emt.obs.index.values,:].obsm["X_phate"]

scv.pp.moments(M1417_adata_emt, n_pcs=30, n_neighbors=30)
scv.pp.moments(M1416_adata_emt, n_pcs=30, n_neighbors=30)
scv.pp.moments(M1436_adata_emt, n_pcs=30, n_neighbors=30)
scv.pp.moments(M1437_adata_emt, n_pcs=30, n_neighbors=30)

mode = "dynamic"
if mode == "static":
    scv.tl.velocity(M1417_adata_emt)
    scv.tl.velocity(M1416_adata_emt)
    scv.tl.velocity(M1436_adata_emt)
    scv.tl.velocity(M1437_adata_emt)
    
elif mode == "dynamic":
    scv.tl.recover_dynamics(M1417_adata_emt, n_jobs = 8)
    scv.tl.recover_dynamics(M1416_adata_emt, n_jobs = 8)
    scv.tl.recover_dynamics(M1436_adata_emt, n_jobs = 8)
    scv.tl.recover_dynamics(M1437_adata_emt, n_jobs = 8)
    scv.tl.velocity(M1417_adata_emt, mode='dynamical')
    scv.tl.velocity(M1416_adata_emt, mode='dynamical')
    scv.tl.velocity(M1436_adata_emt, mode='dynamical')
    scv.tl.velocity(M1437_adata_emt, mode='dynamical')

scv.tl.velocity_graph(M1417_adata_emt)
scv.tl.velocity_graph(M1416_adata_emt)
scv.tl.velocity_graph(M1436_adata_emt)
scv.tl.velocity_graph(M1437_adata_emt)

# In[] project the velocity into phate

velo_dir = "../results_scrnaseq/results_emt/embed_rnavelo/"
import matplotlib.pyplot as plt
fig = plt.figure(figsize = (30, 30))
axs = fig.subplots(nrows = 2, ncols = 2)
scv.pl.velocity_embedding_stream(M1416_adata_emt, basis='phate', color = "annot", palette = "tab10", size = 60, title = "M1416-PP18", figsize = (15,10), fontsize = 20, ax = axs[1,0])
scv.pl.velocity_embedding_stream(M1417_adata_emt, basis='phate', color = "annot", palette = "tab10", size = 60, title = "M1417-PP12", figsize = (15,10), fontsize = 20, ax = axs[0,0])
scv.pl.velocity_embedding_stream(M1436_adata_emt, basis='phate', color = "annot", palette = "tab10", size = 60, title = "M1436-PPC12", figsize = (15,10), fontsize = 20, ax = axs[0,1])
scv.pl.velocity_embedding_stream(M1437_adata_emt, basis='phate', color = "annot", palette = "tab10", size = 60, title = "M1437-PPC18", figsize = (15,10), fontsize = 20, ax = axs[1,1])
if mode == "static":
    fig.savefig(velo_dir + "phate_embedding_emt_rnavelo.png", bbox_inches = "tight", dpi = 300)
elif mode == "dynamic":
    fig.savefig(velo_dir + "phate_embedding_emt_rnavelo_dyn.png", bbox_inches = "tight", dpi = 300)


fig = plt.figure(figsize = (30, 30))
axs = fig.subplots(nrows = 2, ncols = 2)
scv.pl.velocity_embedding_stream(M1416_adata_emt, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1416-PP18", figsize = (15,10), fontsize = 20, ax = axs[1,0])
scv.pl.velocity_embedding_stream(M1417_adata_emt, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1417-PP12", figsize = (15,10), fontsize = 20, ax = axs[0,0])
scv.pl.velocity_embedding_stream(M1436_adata_emt, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1436-PPC12", figsize = (15,10), fontsize = 20, ax = axs[0,1])
scv.pl.velocity_embedding_stream(M1437_adata_emt, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1437-PPC18", figsize = (15,10), fontsize = 20, ax = axs[1,1])
if mode == "static":
    fig.savefig(velo_dir + "umap_embedding_emt_rnavelo.png", bbox_inches = "tight", dpi = 300)
elif mode == "dynamic":
    fig.savefig(velo_dir + "umap_embedding_emt_rnavelo_dyn.png", bbox_inches = "tight", dpi = 300)


# In[]
fig = plt.figure(figsize = (30, 30))
axs = fig.subplots(nrows = 2, ncols = 2)
# scv.pl.velocity_embedding_stream(M1416_adata_emt, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1416-PP18", figsize = (15,10), fontsize = 20, ax = axs[1,0])
# scv.pl.velocity_embedding_stream(M1417_adata_emt, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1417-PP12", figsize = (15,10), fontsize = 20, ax = axs[0,0])
# scv.pl.velocity_embedding_stream(M1436_adata_emt, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1436-PPC12", figsize = (15,10), fontsize = 20, ax = axs[0,1])
# scv.pl.velocity_embedding_stream(M1437_adata_emt, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1437-PPC18", figsize = (15,10), fontsize = 20, ax = axs[1,1])

# axs[1,0].scatter(M1416_adata.obsm["X_umap"][:,0], M1416_adata.obsm["X_umap"][:,1], alpha = 0.5, s = 3, color = "gray")
# axs[0,0].scatter(M1417_adata.obsm["X_umap"][:,0], M1417_adata.obsm["X_umap"][:,1], alpha = 0.5, s = 3, color = "gray")
# axs[0,1].scatter(M1436_adata.obsm["X_umap"][:,0], M1436_adata.obsm["X_umap"][:,1], alpha = 0.5, s = 3, color = "gray")
# axs[1,1].scatter(M1437_adata.obsm["X_umap"][:,0], M1437_adata.obsm["X_umap"][:,1], alpha = 0.5, s = 3, color = "gray")

scv.pl.velocity_embedding_stream(M1416_adata_emt, basis='umap', palette = "tab10", size = 60, title = "M1416-PP18", figsize = (15,10), fontsize = 20, ax = axs[1,0], alpha = 0)
scv.pl.velocity_embedding_stream(M1417_adata_emt, basis='umap', palette = "tab10", size = 60, title = "M1417-PP12", figsize = (15,10), fontsize = 20, ax = axs[0,0], alpha = 0)
scv.pl.velocity_embedding_stream(M1436_adata_emt, basis='umap', palette = "tab10", size = 60, title = "M1436-PPC12", figsize = (15,10), fontsize = 20, ax = axs[0,1], alpha = 0)
scv.pl.velocity_embedding_stream(M1437_adata_emt, basis='umap', palette = "tab10", size = 60, title = "M1437-PPC18", figsize = (15,10), fontsize = 20, ax = axs[1,1], alpha = 0)

sc.pl.umap(M1416_adata, color = "annot", palette = "tab10", ax = axs[1,0], legend_loc = "on data", alpha = 0.2)
sc.pl.umap(M1417_adata, color = "annot", palette = "tab10", ax = axs[0,0], legend_loc = "on data", alpha = 0.2)
sc.pl.umap(M1436_adata, color = "annot", palette = "tab10", ax = axs[0,1], legend_loc = "on data", alpha = 0.2)
sc.pl.umap(M1437_adata, color = "annot", palette = "tab10", ax = axs[1,1], legend_loc = "on data", alpha = 0.2)

plt.tight_layout()

if mode == "static":
    fig.savefig(velo_dir + "umap_embedding_rnavelo.png", bbox_inches = "tight", dpi = 300)
elif mode == "dynamic":
    fig.savefig(velo_dir + "umap_embedding_rnavelo_dyn.png", bbox_inches = "tight", dpi = 300)

# In[]
# NOTE: Option 2: calculate the RNA velocity on the common set of hvgs, 
# there are totally 2605 genes, these genes are selected to be the overlapping genes across samples
# NOTE: the result does not vary a lot compared to the previous one, which shows that the inferred RNA velocity is robust

# M1417_adata_emt = adata_emt[adata_emt.obs["sample"] == "M1417-PP12", :]
# M1416_adata_emt = adata_emt[adata_emt.obs["sample"] == "M1416-PP18", :]
# M1436_adata_emt = adata_emt[adata_emt.obs["sample"] == "M1436-PPC12", :]
# M1437_adata_emt = adata_emt[adata_emt.obs["sample"] == "M1437-PPC18", :]

# scv.pp.moments(M1417_adata_emt, n_pcs=30, n_neighbors=30)
# scv.pp.moments(M1416_adata_emt, n_pcs=30, n_neighbors=30)
# scv.pp.moments(M1436_adata_emt, n_pcs=30, n_neighbors=30)
# scv.pp.moments(M1437_adata_emt, n_pcs=30, n_neighbors=30)

# mode = "static"
# if mode == "static":
#     scv.tl.velocity(M1417_adata_emt)
#     scv.tl.velocity(M1416_adata_emt)
#     scv.tl.velocity(M1436_adata_emt)
#     scv.tl.velocity(M1437_adata_emt)
    
# elif mode == "dynamic":
#     scv.tl.recover_dynamics(M1417_adata_emt, n_jobs = 8)
#     scv.tl.recover_dynamics(M1416_adata_emt, n_jobs = 8)
#     scv.tl.recover_dynamics(M1436_adata_emt, n_jobs = 8)
#     scv.tl.recover_dynamics(M1437_adata_emt, n_jobs = 8)
#     scv.tl.velocity(M1417_adata_emt, mode='dynamical')
#     scv.tl.velocity(M1416_adata_emt, mode='dynamical')
#     scv.tl.velocity(M1436_adata_emt, mode='dynamical')
#     scv.tl.velocity(M1437_adata_emt, mode='dynamical')

# scv.tl.velocity_graph(M1417_adata_emt)
# scv.tl.velocity_graph(M1416_adata_emt)
# scv.tl.velocity_graph(M1436_adata_emt)
# scv.tl.velocity_graph(M1437_adata_emt)

# # In[] project the velocity into phate

# velo_dir = "../results_scrnaseq/results_emt/embed_rnavelo/"
# import matplotlib.pyplot as plt
# fig = plt.figure(figsize = (30, 30))
# axs = fig.subplots(nrows = 2, ncols = 2)
# scv.pl.velocity_embedding_stream(M1416_adata_emt, basis='phate', color = "annot", palette = "tab10", size = 60, title = "M1416-PP18", figsize = (15,10), fontsize = 20, ax = axs[1,0])
# scv.pl.velocity_embedding_stream(M1417_adata_emt, basis='phate', color = "annot", palette = "tab10", size = 60, title = "M1417-PP12", figsize = (15,10), fontsize = 20, ax = axs[0,0])
# scv.pl.velocity_embedding_stream(M1436_adata_emt, basis='phate', color = "annot", palette = "tab10", size = 60, title = "M1436-PPC12", figsize = (15,10), fontsize = 20, ax = axs[0,1])
# scv.pl.velocity_embedding_stream(M1437_adata_emt, basis='phate', color = "annot", palette = "tab10", size = 60, title = "M1437-PPC18", figsize = (15,10), fontsize = 20, ax = axs[1,1])
# if mode == "static":
#     fig.savefig(velo_dir + "phate_embedding_emt_rnavelo2.png", bbox_inches = "tight", dpi = 300)
# elif mode == "dynamic":
#     fig.savefig(velo_dir + "phate_embedding_emt_rnavelo_dyn.png", bbox_inches = "tight", dpi = 300)


# fig = plt.figure(figsize = (30, 30))
# axs = fig.subplots(nrows = 2, ncols = 2)
# scv.pl.velocity_embedding_stream(M1416_adata_emt, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1416-PP18", figsize = (15,10), fontsize = 20, ax = axs[1,0])
# scv.pl.velocity_embedding_stream(M1417_adata_emt, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1417-PP12", figsize = (15,10), fontsize = 20, ax = axs[0,0])
# scv.pl.velocity_embedding_stream(M1436_adata_emt, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1436-PPC12", figsize = (15,10), fontsize = 20, ax = axs[0,1])
# scv.pl.velocity_embedding_stream(M1437_adata_emt, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1437-PPC18", figsize = (15,10), fontsize = 20, ax = axs[1,1])
# if mode == "static":
#     fig.savefig(velo_dir + "umap_embedding_emt_rnavelo2.png", bbox_inches = "tight", dpi = 300)
# elif mode == "dynamic":
#     fig.savefig(velo_dir + "umap_embedding_emt_rnavelo_dyn.png", bbox_inches = "tight", dpi = 300)




# %%
