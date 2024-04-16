# ------------------------------------------------------------ # 
# 
# Running RNA velocity on Luly's project data
#
# ------------------------------------------------------------ # 

# In[]
import numpy as np 
import scanpy as sc 
import scvelo as scv
import matplotlib.pyplot as plt
import phate
import pandas as pd


import sys, os
sys.path.append(".")
import utils


# In[]
# Read the data
data_dir = "../dataset/data_scrnaseq/FOXA1/"

M1385_adata = scv.read(data_dir + "M1385.loom", cache=True)
M1386_adata = scv.read(data_dir + "M1386.loom", cache=True)
M1387_adata = scv.read(data_dir + "M1387.loom", cache=True)
M1415_adata = scv.read(data_dir + "M1415.loom", cache=True)
M1418_adata = scv.read(data_dir + "M1418.loom", cache=True)
M1475_adata = scv.read(data_dir + "M1475.loom", cache=True)

# In[]
M1385_adata.var_names_make_unique()
M1386_adata.var_names_make_unique()
M1387_adata.var_names_make_unique()
M1415_adata.var_names_make_unique()
M1418_adata.var_names_make_unique()
M1475_adata.var_names_make_unique()

# included cells 
meta_cells = pd.read_csv(data_dir + "included_cell_barcode.csv", index_col = 0)

# In[]
# rename the index of meta cells
new_idx = []
samples = []
for barcode in meta_cells.index.values:
    if barcode.split("-")[1] == "1_1":
        new_idx.append("M1385:" + barcode.split("-")[0] + "x")
        samples.append("M1385")
    
    elif barcode.split("-")[1] == "1_2":
        new_idx.append("M1386:" + barcode.split("-")[0] + "x")
        samples.append("M1386")
    
    elif barcode.split("-")[1] == "1_3":
        new_idx.append("M1387:" + barcode.split("-")[0] + "x")
        samples.append("M1387")
    
    elif barcode.split("-")[1] == "1_4":
        new_idx.append("M1418:" + barcode.split("-")[0] + "x")
        samples.append("M1418")
    
    elif barcode.split("-")[1] == "1_5":
        new_idx.append("M1415:" + barcode.split("-")[0] + "x")
        samples.append("M1415")
    
    else:
        new_idx.append("M1475:" + barcode.split("-")[0] + "x")
        samples.append("M1475")


meta_cells.index = np.array(new_idx)
meta_cells["sample"] = np.array(samples)

# m1385_idx = np.array(["M1385:" + x.split("-")[0] + "x" for x in meta_cells.index.values if x.split("-")[1] == "1_1"])
# m1386_idx = np.array(["M1386:" + x.split("-")[0] + "x" for x in meta_cells.index.values if x.split("-")[1] == "1_2"])
# m1387_idx = np.array(["M1387:" + x.split("-")[0] + "x" for x in meta_cells.index.values if x.split("-")[1] == "1_3"])
# m1418_idx = np.array(["M1418:" + x.split("-")[0] + "x" for x in meta_cells.index.values if x.split("-")[1] == "1_4"])
# m1415_idx = np.array(["M1415:" + x.split("-")[0] + "x" for x in meta_cells.index.values if x.split("-")[1] == "1_5"])
# m1475_idx = np.array(["M1475:" + x.split("-")[0] + "x" for x in meta_cells.index.values if x.split("-")[1] == "1_6"])

# M1385_adata = M1385_adata[m1385_idx,:]
# M1386_adata = M1386_adata[m1386_idx,:]
# M1387_adata = M1387_adata[m1387_idx,:]
# M1418_adata = M1418_adata[m1418_idx,:]
# M1415_adata = M1415_adata[m1415_idx,:]
# M1475_adata = M1475_adata[m1475_idx,:]

adata_merge = sc.concat([M1385_adata, M1386_adata, M1387_adata, M1415_adata, M1418_adata, M1475_adata], join = "inner")
adata_merge = adata_merge[meta_cells.index.values,:]
adata_merge.obs[meta_cells.columns] = meta_cells

# use existing visualization
# adata_merge.obsm["X_umap"] = adata_merge.obs[["_X", "_Y"]].values

# In[]
scv.pl.proportions(adata_merge)

M1385_adata = adata_merge[adata_merge.obs["sample"] == "M1385",:]
M1386_adata = adata_merge[adata_merge.obs["sample"] == "M1386",:]
M1387_adata = adata_merge[adata_merge.obs["sample"] == "M1387",:]
M1415_adata = adata_merge[adata_merge.obs["sample"] == "M1415",:]
M1418_adata = adata_merge[adata_merge.obs["sample"] == "M1418",:]
M1475_adata = adata_merge[adata_merge.obs["sample"] == "M1475",:]

# separate according to samples
scv.pp.filter_and_normalize(M1385_adata, min_shared_counts=20, n_top_genes = 2000)
scv.pp.filter_and_normalize(M1386_adata, min_shared_counts=20, n_top_genes = 2000)
scv.pp.filter_and_normalize(M1387_adata, min_shared_counts=20, n_top_genes = 2000)
scv.pp.filter_and_normalize(M1415_adata, min_shared_counts=20, n_top_genes = 2000)
scv.pp.filter_and_normalize(M1418_adata, min_shared_counts=20, n_top_genes = 2000)
scv.pp.filter_and_normalize(M1475_adata, min_shared_counts=20, n_top_genes = 2000)

# In[]
# calculate umap visualization
sc.pp.normalize_total(adata_merge, target_sum=1e4)
sc.pp.log1p(adata_merge)

# automatically used for umap calculation
sc.pp.highly_variable_genes(adata_merge, n_top_genes = 2000)
adata_merge = adata_merge[:, adata_merge.var.highly_variable]

# calculate UMAP embedding
sc.pp.neighbors(adata_merge, n_neighbors = 50, n_pcs = 100)
sc.tl.umap(adata_merge, min_dist = 0.5)

M1385_adata.obsm["X_umap"] = adata_merge[M1385_adata.obs.index,:].obsm["X_umap"]
M1386_adata.obsm["X_umap"] = adata_merge[M1386_adata.obs.index,:].obsm["X_umap"]
M1387_adata.obsm["X_umap"] = adata_merge[M1387_adata.obs.index,:].obsm["X_umap"]
M1415_adata.obsm["X_umap"] = adata_merge[M1415_adata.obs.index,:].obsm["X_umap"]
M1418_adata.obsm["X_umap"] = adata_merge[M1418_adata.obs.index,:].obsm["X_umap"]
M1475_adata.obsm["X_umap"] = adata_merge[M1475_adata.obs.index,:].obsm["X_umap"]

# In[]
scv.pp.moments(M1385_adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(M1386_adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(M1387_adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(M1415_adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(M1418_adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(M1475_adata, n_pcs=30, n_neighbors=30)

mode = "dynamic"
if mode == "static":
    scv.tl.velocity(M1385_adata)
    scv.tl.velocity(M1386_adata)
    scv.tl.velocity(M1387_adata)
    scv.tl.velocity(M1415_adata)
    scv.tl.velocity(M1418_adata)
    scv.tl.velocity(M1475_adata)
    
elif mode == "dynamic":
    scv.tl.recover_dynamics(M1385_adata, n_jobs = 8)
    scv.tl.recover_dynamics(M1386_adata, n_jobs = 8)
    scv.tl.recover_dynamics(M1387_adata, n_jobs = 8)
    scv.tl.recover_dynamics(M1415_adata, n_jobs = 8)
    scv.tl.recover_dynamics(M1418_adata, n_jobs = 8)
    scv.tl.recover_dynamics(M1475_adata, n_jobs = 8)
    scv.tl.velocity(M1385_adata, mode='dynamical')
    scv.tl.velocity(M1386_adata, mode='dynamical')
    scv.tl.velocity(M1387_adata, mode='dynamical')
    scv.tl.velocity(M1415_adata, mode='dynamical')
    scv.tl.velocity(M1418_adata, mode='dynamical')
    scv.tl.velocity(M1475_adata, mode='dynamical')

scv.tl.velocity_graph(M1385_adata)
scv.tl.velocity_graph(M1386_adata)
scv.tl.velocity_graph(M1387_adata)
scv.tl.velocity_graph(M1415_adata)
scv.tl.velocity_graph(M1418_adata)
scv.tl.velocity_graph(M1475_adata)

# In[] project the velocity into phate

velo_dir = "../results_scrnaseq/results_emt/FOXA1/embed_rnavelo/"
if not os.path.exists(velo_dir):
    os.makedirs(velo_dir)

import matplotlib.pyplot as plt
fig = plt.figure(figsize = (30, 30))
axs = fig.subplots(nrows = 3, ncols = 2)
scv.pl.velocity_embedding_stream(M1385_adata, basis='umap', color = "celltypes_coarse", palette = "tab10", size = 60, title = "M1385", figsize = (15,10), fontsize = 20, ax = axs[0,0])
scv.pl.velocity_embedding_stream(M1386_adata, basis='umap', color = "celltypes_coarse", palette = "tab10", size = 60, title = "M1386", figsize = (15,10), fontsize = 20, ax = axs[0,1])
scv.pl.velocity_embedding_stream(M1387_adata, basis='umap', color = "celltypes_coarse", palette = "tab10", size = 60, title = "M1387", figsize = (15,10), fontsize = 20, ax = axs[1,0])
scv.pl.velocity_embedding_stream(M1415_adata, basis='umap', color = "celltypes_coarse", palette = "tab10", size = 60, title = "M1415", figsize = (15,10), fontsize = 20, ax = axs[1,1])
scv.pl.velocity_embedding_stream(M1418_adata, basis='umap', color = "celltypes_coarse", palette = "tab10", size = 60, title = "M1418", figsize = (15,10), fontsize = 20, ax = axs[2,0])
scv.pl.velocity_embedding_stream(M1475_adata, basis='umap', color = "celltypes_coarse", palette = "tab10", size = 60, title = "M1475", figsize = (15,10), fontsize = 20, ax = axs[2,1])

if mode == "static":
    fig.savefig(velo_dir + "umap_embedding_emt_rnavelo.png", bbox_inches = "tight", dpi = 300)
elif mode == "dynamic":
    fig.savefig(velo_dir + "umap_embedding_emt_rnavelo_dyn.png", bbox_inches = "tight", dpi = 300)


# %%
