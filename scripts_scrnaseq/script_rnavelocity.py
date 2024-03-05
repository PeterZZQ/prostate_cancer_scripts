# In[]
import numpy as np 
import scanpy as sc 
import scvelo as scv
import matplotlib.pyplot as plt

# In[]
# Read the data
data_dir = "../dataset/data_scrnaseq/Cell_Ranger_output/loom_files/"
seurat_dir = "../dataset/data_scrnaseq/seurat_integration/"

M1416_adata = scv.read(data_dir + "M1416.loom", cache=True)
M1417_adata = scv.read(data_dir + "M1417.loom", cache=True)
M1436_adata = scv.read(data_dir + "M1436.loom", cache=True)
M1437_adata = scv.read(data_dir + "M1437.loom", cache=True)
M1476_adata = scv.read(data_dir + "M1476.loom", cache=True)
M1477_adata = scv.read(data_dir + "M1477.loom", cache=True)

# assign cell types
adata_intact_seurat = sc.read_h5ad(seurat_dir + "adata_intact_seurat.h5ad")
adata_castrated_seurat = sc.read_h5ad(seurat_dir + "adata_castrated_seurat.h5ad")

meta_m1416 = adata_intact_seurat.obs.loc[adata_intact_seurat.obs["sample"] == "M1416-PP18",:]
barcodes = ["M1416:" + x.split("-")[0] + "x" for x in meta_m1416.index.values.squeeze()]
meta_m1416.index = barcodes
# filter loom file accordingly
# M1416_adata = M1416_adata[barcodes,:]
M1416_adata.obs["annot"] = "Noisy"
M1416_adata.obs.loc[barcodes, "annot"] = meta_m1416["annot"].values
M1416_adata.obs["sample"] = "M1416-PP18"
# M1416_adata.obsm = adata_intact_seurat[adata_intact_seurat.obs["sample"] == "M1416-PP18",:].obsm.copy()

meta_m1417 = adata_intact_seurat.obs.loc[adata_intact_seurat.obs["sample"] == "M1417-PP12",:]
barcodes = ["M1417:" + x.split("-")[0] + "x" for x in meta_m1417.index.values.squeeze()]
meta_m1417.index = barcodes
# filter loom file accordingly
# M1417_adata = M1417_adata[barcodes,:]
M1417_adata.obs["annot"] = "Noisy"
M1417_adata.obs.loc[barcodes,"annot"] = meta_m1417["annot"].values
M1417_adata.obs["sample"] = "M1417-PP12"
# M1417_adata.obsm = adata_intact_seurat[adata_intact_seurat.obs["sample"] == "M1417-PP12",:].obsm.copy()

meta_m1436 = adata_intact_seurat.obs.loc[adata_intact_seurat.obs["sample"] == "M1436-PPC12",:]
barcodes = ["M1436:" + x.split("-")[0] + "x" for x in meta_m1436.index.values.squeeze()]
meta_m1436.index = barcodes
# filter loom file accordingly
# M1436_adata = M1436_adata[barcodes,:]
M1436_adata.obs["annot"] = "Noisy"
M1436_adata.obs.loc[barcodes,"annot"] = meta_m1436["annot"].values
M1436_adata.obs["sample"] = "M1436-PPC12"
# M1436_adata.obsm = adata_intact_seurat[adata_intact_seurat.obs["sample"] == "M1436-PPC12",:].obsm.copy()

meta_m1437 = adata_intact_seurat.obs.loc[adata_intact_seurat.obs["sample"] == "M1437-PPC18",:]
barcodes = ["M1437:" + x.split("-")[0] + "x" for x in meta_m1437.index.values.squeeze()]
meta_m1437.index = barcodes
# filter loom file accordingly
# M1437_adata = M1437_adata[barcodes,:]
M1437_adata.obs["annot"] = "Noisy"
M1437_adata.obs.loc[barcodes,"annot"] = meta_m1437["annot"].values
M1437_adata.obs["sample"] = "M1437-PPC18"
# M1437_adata.obsm = adata_intact_seurat[adata_intact_seurat.obs["sample"] == "M1437-PPC18",:].obsm.copy()

meta_m1476 = adata_castrated_seurat.obs.loc[adata_castrated_seurat.obs["sample"] == "M1476-PP",:]
barcodes = ["M1476:" + x.split("-")[0] + "x" for x in meta_m1476.index.values.squeeze()]
# meta_m1476.index = barcodes
# filter loom file accordingly
# M1476_adata = M1476_adata[barcodes,:]
M1476_adata.obs["annot"] = "Noisy"
M1476_adata.obs.loc[barcodes,"annot"] = meta_m1476["annot"].values
M1476_adata.obs["sample"] = "M1476-PP"
# M1476_adata.obsm = adata_castrated_seurat[adata_castrated_seurat.obs["sample"] == "M1476-PP",:].obsm.copy()

meta_m1477 = adata_castrated_seurat.obs.loc[adata_castrated_seurat.obs["sample"] == "M1477-PPC",:]
barcodes = ["M1477:" + x.split("-")[0] + "x" for x in meta_m1477.index.values.squeeze()]
# meta_m1477.index = barcodes
# filter loom file accordingly
# M1477_adata = M1477_adata[barcodes,:]
M1477_adata.obs["annot"] = "Noisy"
M1477_adata.obs.loc[barcodes,"annot"] = meta_m1477["annot"].values
M1477_adata.obs["sample"] = "M1477-PPC"
# M1477_adata.obsm = adata_castrated_seurat[adata_castrated_seurat.obs["sample"] == "M1477-PPC",:].obsm.copy()

# In[]
# Preprocessing
M1416_adata.var_names_make_unique()
M1417_adata.var_names_make_unique()
M1436_adata.var_names_make_unique()
M1437_adata.var_names_make_unique()
M1476_adata.var_names_make_unique()
M1477_adata.var_names_make_unique()

adata = sc.concat([M1416_adata, M1417_adata, M1436_adata, M1437_adata, M1476_adata, M1477_adata])
# joint, for umap consistency
scv.pp.filter_and_normalize(adata, min_shared_counts = 20, n_top_genes = 2000)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

M1416_adata.obsm["X_umap"] = adata[M1416_adata.obs.index.values].obsm["X_umap"].toarray()
M1417_adata.obsm["X_umap"] = adata[M1417_adata.obs.index.values].obsm["X_umap"].toarray()
M1436_adata.obsm["X_umap"] = adata[M1436_adata.obs.index.values].obsm["X_umap"].toarray()
M1437_adata.obsm["X_umap"] = adata[M1437_adata.obs.index.values].obsm["X_umap"].toarray()
M1476_adata.obsm["X_umap"] = adata[M1476_adata.obs.index.values].obsm["X_umap"].toarray()
M1477_adata.obsm["X_umap"] = adata[M1477_adata.obs.index.values].obsm["X_umap"].toarray()

# check proportion of spliced and unspliced counts
scv.pl.proportions(adata)

# In[]
M1476_adata.var['mt'] = M1476_adata.var_names.str.startswith('mt-')
M1476_adata.obs["small count<200"] = np.where((M1476_adata.obs["n_genes_by_counts"].values <= 200), 1, 0)
sc.pp.calculate_qc_metrics(M1476_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


sc.pl.umap(M1476_adata, color = ["pct_counts_mt", "n_genes_by_counts", "total_counts", "small<200"])

# In[]
# Select only the EMT population
# M1416_adata = M1416_adata[(M1416_adata.obs["annot"] == "Basal") | (M1416_adata.obs["annot"] == "Luminal") | (M1416_adata.obs["annot"] == "Mesenchymal")]
# M1417_adata = M1417_adata[(M1417_adata.obs["annot"] == "Basal") | (M1417_adata.obs["annot"] == "Luminal") | (M1417_adata.obs["annot"] == "Mesenchymal")]
# M1436_adata = M1436_adata[(M1436_adata.obs["annot"] == "Basal") | (M1436_adata.obs["annot"] == "Luminal") | (M1436_adata.obs["annot"] == "Mesenchymal")]
# M1437_adata = M1437_adata[(M1437_adata.obs["annot"] == "Basal") | (M1437_adata.obs["annot"] == "Luminal") | (M1437_adata.obs["annot"] == "Mesenchymal")]
# M1476_adata = M1476_adata[(M1476_adata.obs["annot"] == "Basal") | (M1476_adata.obs["annot"] == "Luminal") | (M1476_adata.obs["annot"] == "Mesenchymal")]
# M1477_adata = M1477_adata[(M1477_adata.obs["annot"] == "Basal") | (M1477_adata.obs["annot"] == "Luminal") | (M1477_adata.obs["annot"] == "Mesenchymal")]

# In[]
# Plot luminal markers
# Mesenchymal markers: Lum, Vim
# Luminal markers: Ar, Krt8, Cd24a, Krt18
# Stem cell markers: Cd44, Epcam, Sox2, Pou5f1, Nanog, Letm1 
# non expression of Sox2, Pou5f1, Nanog
with plt.rc_context({"figure.figsize": (7, 5), "figure.dpi": (300), "font.size": 15}):   
    fig = sc.pl.umap(M1416_adata, color = ["Ar", "Krt8", "Cd24a", "Krt18", "Lum", "Vim", "Cd44", "Epcam", "Letm1"], return_fig = True, color_map = "Reds", ncols = 2)
    fig.savefig("../results_ti/castrated_basal_mesenchymal/markers/m1416_pp18.png", bbox_inches = "tight")

with plt.rc_context({"figure.figsize": (7, 5), "figure.dpi": (300), "font.size": 15}):
    fig = sc.pl.umap(M1417_adata, color = ["Ar", "Krt8", "Cd24a", "Krt18", "Lum", "Vim", "Cd44", "Epcam", "Letm1"], return_fig = True, color_map = "Reds", ncols = 2)
    fig.savefig("../results_ti/castrated_basal_mesenchymal/markers/m1417_pp12.png", bbox_inches = "tight")

with plt.rc_context({"figure.figsize": (7, 5), "figure.dpi": (300), "font.size": 15}):
    fig = sc.pl.umap(M1436_adata, color = ["Ar", "Krt8", "Cd24a", "Krt18", "Lum", "Vim", "Cd44", "Epcam", "Letm1"], return_fig = True, color_map = "Reds", ncols = 2)
    fig.savefig("../results_ti/castrated_basal_mesenchymal/markers/m1436_ppc12.png", bbox_inches = "tight")

with plt.rc_context({"figure.figsize": (7, 5), "figure.dpi": (300), "font.size": 15}):
    fig = sc.pl.umap(M1437_adata, color = ["Ar", "Krt8", "Cd24a", "Krt18", "Lum", "Vim", "Cd44", "Epcam", "Letm1"], return_fig = True, color_map = "Reds", ncols = 2)
    fig.savefig("../results_ti/castrated_basal_mesenchymal/markers/m1437_ppc18.png", bbox_inches = "tight")

with plt.rc_context({"figure.figsize": (7, 5), "figure.dpi": (300), "font.size": 15}):
    fig = sc.pl.umap(M1476_adata, color = ["Ar", "Krt8", "Cd24a", "Krt18", "Lum", "Vim", "Cd44", "Epcam", "Letm1"], return_fig = True, color_map = "Reds", ncols = 2)
    fig.savefig("../results_ti/castrated_basal_mesenchymal/markers/m1476_pp.png", bbox_inches = "tight")

with plt.rc_context({"figure.figsize": (7, 5), "figure.dpi": (300), "font.size": 15}):
    fig = sc.pl.umap(M1477_adata, color = ["Ar", "Krt8", "Cd24a", "Krt18", "Lum", "Vim", "Cd44", "Epcam", "Letm1"], return_fig = True, color_map = "Reds", ncols = 2)
    fig.savefig("../results_ti/castrated_basal_mesenchymal/markers/m1477_ppc.png", bbox_inches = "tight")

# In[]
scv.pp.filter_and_normalize(M1416_adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.filter_and_normalize(M1417_adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.filter_and_normalize(M1436_adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.filter_and_normalize(M1437_adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.filter_and_normalize(M1476_adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.filter_and_normalize(M1477_adata, min_shared_counts=20, n_top_genes=2000)


# In[]
scv.pp.moments(M1416_adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(M1417_adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(M1436_adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(M1437_adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(M1476_adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(M1477_adata, n_pcs=30, n_neighbors=30)

mode = "dynamic"
if mode == "static":
    scv.tl.velocity(M1416_adata)
    scv.tl.velocity(M1417_adata)
    scv.tl.velocity(M1436_adata)
    scv.tl.velocity(M1437_adata)
    scv.tl.velocity(M1476_adata)
    scv.tl.velocity(M1477_adata)
elif mode == "dynamic":
    scv.tl.recover_dynamics(M1416_adata, n_jobs = 8)
    scv.tl.recover_dynamics(M1417_adata, n_jobs = 8)
    scv.tl.recover_dynamics(M1436_adata, n_jobs = 8)
    scv.tl.recover_dynamics(M1437_adata, n_jobs = 8)
    scv.tl.recover_dynamics(M1476_adata, n_jobs = 8)
    scv.tl.recover_dynamics(M1477_adata, n_jobs = 8)
    scv.tl.velocity(M1416_adata, mode='dynamical')
    scv.tl.velocity(M1417_adata, mode='dynamical')
    scv.tl.velocity(M1436_adata, mode='dynamical')
    scv.tl.velocity(M1437_adata, mode='dynamical')
    scv.tl.velocity(M1476_adata, mode='dynamical')
    scv.tl.velocity(M1477_adata, mode='dynamical')

scv.tl.velocity_graph(M1416_adata)
scv.tl.velocity_graph(M1417_adata)
scv.tl.velocity_graph(M1436_adata)
scv.tl.velocity_graph(M1437_adata)
scv.tl.velocity_graph(M1476_adata)
scv.tl.velocity_graph(M1477_adata)

# In[]
velo_dir = "../results_ti/castrated_basal_mesenchymal/embed_rnavelo_filtered/"
import matplotlib.pyplot as plt
fig = plt.figure(figsize = (30, 30))
axs = fig.subplots(nrows = 3, ncols = 2)
scv.pl.velocity_embedding_stream(M1416_adata, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1416-PP18", figsize = (15,10), fontsize = 20, ax = axs[1,0])
scv.pl.velocity_embedding_stream(M1417_adata, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1417-PP12", figsize = (15,10), fontsize = 20, ax = axs[0,0])
scv.pl.velocity_embedding_stream(M1436_adata, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1436-PPC12", figsize = (15,10), fontsize = 20, ax = axs[0,1])
scv.pl.velocity_embedding_stream(M1437_adata, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1437-PPC18", figsize = (15,10), fontsize = 20, ax = axs[1,1])
scv.pl.velocity_embedding_stream(M1476_adata, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1476-PP", figsize = (15,10), fontsize = 20, ax = axs[2,0])
scv.pl.velocity_embedding_stream(M1477_adata, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1477-PPC", figsize = (15,10), fontsize = 20, ax = axs[2,1])
if mode == "static":
    fig.savefig(velo_dir + "umap_embedding_rnavelo.png", bbox_inches = "tight", dpi = 300)
elif mode == "dynamic":
    fig.savefig(velo_dir + "umap_embedding_rnavelo_dyn.png", bbox_inches = "tight", dpi = 300)

# fig = plt.figure(figsize = (30, 30))
# axs = fig.subplots(nrows = 3, ncols = 2)
# scv.pl.velocity_embedding_stream(M1416_adata, basis='pca', color = "annot", palette = "tab10", size = 60, title = "M1416-PP18", figsize = (15,10), fontsize = 20, ax = axs[1,0])
# scv.pl.velocity_embedding_stream(M1417_adata, basis='pca', color = "annot", palette = "tab10", size = 60, title = "M1417-PP12", figsize = (15,10), fontsize = 20, ax = axs[0,0])
# scv.pl.velocity_embedding_stream(M1436_adata, basis='pca', color = "annot", palette = "tab10", size = 60, title = "M1436-PPC12", figsize = (15,10), fontsize = 20, ax = axs[0,1])
# scv.pl.velocity_embedding_stream(M1437_adata, basis='pca', color = "annot", palette = "tab10", size = 60, title = "M1437-PPC18", figsize = (15,10), fontsize = 20, ax = axs[1,1])
# scv.pl.velocity_embedding_stream(M1476_adata, basis='pca', color = "annot", palette = "tab10", size = 60, title = "M1476-PP", figsize = (15,10), fontsize = 20, ax = axs[2,0])
# scv.pl.velocity_embedding_stream(M1477_adata, basis='pca', color = "annot", palette = "tab10", size = 60, title = "M1477-PPC", figsize = (15,10), fontsize = 20, ax = axs[2,1])
# if mode == "static":
#     fig.savefig(velo_dir + "pca_embedding_rnavelo.png", bbox_inches = "tight", dpi = 300)
# elif mode == "dynamic":
#     fig.savefig(velo_dir + "pca_embedding_rnavelo_dyn.png", bbox_inches = "tight", dpi = 300)


# In[]
import matplotlib.pyplot as plt
fig = plt.figure(figsize = (30, 30))
axs = fig.subplots(nrows = 3, ncols = 2)
scv.pl.velocity_embedding(M1416_adata, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1416-PP18", figsize = (15,10), fontsize = 20, ax = axs[1,0], arrow_length=4, arrow_size=1)
scv.pl.velocity_embedding(M1417_adata, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1417-PP12", figsize = (15,10), fontsize = 20, ax = axs[0,0], arrow_length=4, arrow_size=1)
scv.pl.velocity_embedding(M1436_adata, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1436-PPC12", figsize = (15,10), fontsize = 20, ax = axs[0,1], arrow_length=4, arrow_size=1)
scv.pl.velocity_embedding(M1437_adata, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1437-PPC18", figsize = (15,10), fontsize = 20, ax = axs[1,1], arrow_length=4, arrow_size=1)
scv.pl.velocity_embedding(M1476_adata, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1476-PP", figsize = (15,10), fontsize = 20, ax = axs[2,0], arrow_length=4, arrow_size=1)
scv.pl.velocity_embedding(M1477_adata, basis='umap', color = "annot", palette = "tab10", size = 60, title = "M1477-PPC", figsize = (15,10), fontsize = 20, ax = axs[2,1], arrow_length=4, arrow_size=1)
if mode == "static":
    fig.savefig(velo_dir + "umap_embedding_rnavelo_arrow.png", bbox_inches = "tight", dpi = 300)
elif mode == "dynamic":
    fig.savefig(velo_dir + "umap_embedding_rnavelo_dyn_arrow.png", bbox_inches = "tight", dpi = 300)

# fig = plt.figure(figsize = (30, 30))
# axs = fig.subplots(nrows = 3, ncols = 2)
# scv.pl.velocity_embedding(M1416_adata, basis='pca', color = "annot", palette = "tab10", size = 60, title = "M1416-PP18", figsize = (15,10), fontsize = 20, ax = axs[1,0], arrow_length=2, arrow_size=1)
# scv.pl.velocity_embedding(M1417_adata, basis='pca', color = "annot", palette = "tab10", size = 60, title = "M1417-PP12", figsize = (15,10), fontsize = 20, ax = axs[0,0], arrow_length=2, arrow_size=1)
# scv.pl.velocity_embedding(M1436_adata, basis='pca', color = "annot", palette = "tab10", size = 60, title = "M1436-PPC12", figsize = (15,10), fontsize = 20, ax = axs[0,1], arrow_length=2, arrow_size=1)
# scv.pl.velocity_embedding(M1437_adata, basis='pca', color = "annot", palette = "tab10", size = 60, title = "M1437-PPC18", figsize = (15,10), fontsize = 20, ax = axs[1,1], arrow_length=2, arrow_size=1)
# scv.pl.velocity_embedding(M1476_adata, basis='pca', color = "annot", palette = "tab10", size = 60, title = "M1476-PP", figsize = (15,10), fontsize = 20, ax = axs[2,0], arrow_length=2, arrow_size=1)
# scv.pl.velocity_embedding(M1477_adata, basis='pca', color = "annot", palette = "tab10", size = 60, title = "M1477-PPC", figsize = (15,10), fontsize = 20, ax = axs[2,1], arrow_length=2, arrow_size=1)
# if mode == "static":
#     fig.savefig(velo_dir + "pca_embedding_rnavelo_arrow.png", bbox_inches = "tight", dpi = 300)
# elif mode == "dynamic":
#     fig.savefig(velo_dir + "pca_embedding_rnavelo_dyn_arrow.png", bbox_inches = "tight", dpi = 300)


# %%
