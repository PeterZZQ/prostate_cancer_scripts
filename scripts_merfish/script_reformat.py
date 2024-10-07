# In[]
# --------------------------------------------------------
#
# Script adopted from the Vizgen example dataset:
# https://info.vizgen.com/mouse-liver-data?submissionGuid=5ca082a2-01f7-46c4-ad91-d005a6549c35 
#
# --------------------------------------------------------
import pandas as pd
import numpy as np 
from anndata import AnnData
import scanpy as sc
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import sys
sys.path.append("../")
import utils

# In[]
# --------------------------------------------------------
#
# Load data
#
# --------------------------------------------------------
# old merfish dataset
# base_path = "../dataset/data_merfish/"
# dataset_name = "region_0/" # PP15
# dataset_name = "region_1/" # PPC15

# merfish0617
# base_path = "../dataset/merfish0617/"
# dataset_name = "region_0/" # PPC15
# dataset_name = "region_1/" # PP15

# merfish0626
base_path = "../dataset/merfish0626/"
# dataset_name = "region_0/" # PPC19
dataset_name = "region_1/" # PP19

cell_by_gene = pd.read_csv(base_path + dataset_name + "cell_by_gene.csv", index_col=0)
meta_cell = pd.read_csv(base_path + dataset_name + "cell_metadata.csv", index_col=0)
meta_cell["barcodeCount"] = cell_by_gene.sum(axis=1)
# initialize meta_gene
meta_gene = pd.DataFrame(index=cell_by_gene.columns.tolist())
# drop blanks for single cell analysis
keep_genes = [x for x in cell_by_gene.columns.tolist() if 'Blank' not in x]
cell_by_gene = cell_by_gene[keep_genes]
meta_gene = meta_gene.loc[keep_genes]
meta_gene['expression'] = cell_by_gene.sum(axis=0)
vizgen_genes = meta_gene.index.tolist()

# create anndata
adata = AnnData(X = cell_by_gene, obs = meta_cell, var = meta_gene)
# assign spatial location of cells
adata.obsm["spatial"] = adata.obs[["center_x", "center_y"]].values
adata.obs.drop(columns=["center_x", "center_y"], inplace=True)

# # NOTE: only for pp in 0626
# selected_cells = np.array([str(x) for x in pd.read_csv("../dataset/merfish0626/selected_cellid.csv", header = None).values.squeeze()], dtype = object)
# adata = adata[selected_cells,:]
adata.write_h5ad(base_path + dataset_name + "adata.h5ad")

# In[]
# Basic preprocessing, adopted from:
# https://squidpy.readthedocs.io/en/latest/notebooks/tutorials/tutorial_vizgen_mouse_liver.html
# base_path = "../dataset/data_merfish/"
# dataset_name = "region_1/"
# adata = sc.read_h5ad(base_path + dataset_name + "adata.h5ad")
print(adata)

# check library size
libsize = adata.obs["barcodeCount"].values.squeeze()
plt.hist(libsize, bins = 30)

# no mitochondrial gene
adata.var_names_make_unique()
# fitler cells
sc.pp.filter_cells(adata, min_genes = 10)
print(adata)
# save the filtered anndata
adata.write_h5ad(base_path + dataset_name + "adata_qc.h5ad")

# totally 300 genes, no futher filter of genes
# sc.pp.filter_genes(adata, min_cells=10)
print("normalize total")
sc.pp.normalize_total(adata)
print("log transform")
sc.pp.log1p(adata)

# In[]
# Basic visualization
resolution = 1.5
print("PCA")
sc.tl.pca(adata, svd_solver="arpack")
print("neighbors")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
print("UMAP")
sc.tl.umap(adata)
# print("Leiden")
# sc.tl.leiden(adata, resolution=resolution)
sc.set_figure_params(figsize=(10, 10))
sc.pl.umap(adata)
# sc.pl.spatial(adata, shape=None, color="leiden", size=0.5, library_id="spatial", figsize=(10, 10))

# In[]
# --------------------------------------------------------
#
# Load Proseg data
#
# --------------------------------------------------------
merfish = "merfish0626"
if merfish == "merfish1":
    base_path = "../dataset/data_merfish/"
elif merfish == "merfish0617":
    base_path = "../dataset/merfish0617/"
elif merfish == "merfish0626":
    base_path = "../dataset/merfish0626/"

block = "region1"
if block == "region0":
    dataset_name = "region_0/proseg_results/"
else:
    dataset_name = "region_1/proseg_results/"

# In[]
cell_by_gene = pd.read_csv(base_path + dataset_name + "expected-counts.csv")
meta_cell = pd.read_csv(base_path + dataset_name + "cell-metadata.csv", index_col=0)
meta_cell["barcodeCount"] = cell_by_gene.sum(axis=1)
# initialize meta_gene
meta_gene = pd.DataFrame(index=cell_by_gene.columns.tolist())
# drop blanks for single cell analysis
keep_genes = [x for x in cell_by_gene.columns.tolist() if 'Blank' not in x]
cell_by_gene = cell_by_gene[keep_genes]
meta_gene = meta_gene.loc[keep_genes]
meta_gene['expression'] = cell_by_gene.sum(axis=0)
vizgen_genes = meta_gene.index.tolist()

# rename the cell barcode to make it unique
cell_by_gene.index = ["cell" + str(x) + "-" + block + "-merfish1" for x in cell_by_gene.index.values]
meta_cell.index = [x for x in cell_by_gene.index.values]
# create anndata
adata = AnnData(X = cell_by_gene, obs = meta_cell, var = meta_gene)
# the original count is fractional
adata.layers["counts.proseg.orig"] = csr_matrix(adata.X)
adata.X = csr_matrix(np.rint(adata.X))

# assign spatial location of cells
adata.obsm["X_spatial"] = adata.obs[["centroid_x", "centroid_y"]].values
adata.write_h5ad(base_path + dataset_name + "adata.h5ad")

# In[]
# QC on the proseg data
# Plot library size per cell
adata = sc.read_h5ad(base_path + dataset_name + "adata.h5ad")
print(adata)
# basic filtering
sc.pp.filter_cells(adata, min_genes = 10)
print(adata)
# calculate qc metrics automatically
sc.pp.calculate_qc_metrics(adata, percent_top = None, log1p=False, inplace = True)
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")
sc.pl.scatter(adata, basis = "spatial", color = "total_counts")
sc.pl.scatter(adata, basis = "spatial", color = "n_genes_by_counts")
# filter cells according to qc, mainly luminal
# merfish1 region 0
# adata = adata[adata.obs["n_genes_by_counts"] < 90, :]
# merfish1 region 1
# adata = adata[adata.obs["n_genes_by_counts"] < 80, :]
# merfish0617, merfish0626 region 0, 1
adata = adata[adata.obs["n_genes_by_counts"] < 140, :]


# In[]
adata.layers["counts.proseg.round"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.neighbors(adata, n_neighbors = 50)
sc.tl.umap(adata, min_dist = 0.1)
# move back the count
adata.X = adata.layers["counts.proseg.round"].copy()
adata.write_h5ad(base_path + dataset_name + "adata_qc.h5ad")

# cells are clustered in the space, but the clusters are fewer, try label transfer first
fig = plt.figure(figsize = (10, 7))
ax = fig.add_subplot()
sc.pl.umap(adata, color = "n_genes_by_counts", ax = ax)
fig.savefig(base_path + dataset_name + "umap_orig.png", dpi = 150, bbox_inches = "tight")


# In[]
assert merfish == "merfish0626"
assert block == "region1"
adata = sc.read_h5ad(base_path + dataset_name + "adata_qc.h5ad")
sc.pp.neighbors(adata, n_neighbors = 30, use_rep = "X_spatial")
sc.tl.leiden(adata, resolution = 0.01)
fig = utils.plot_embeds(embed = adata.obsm["X_spatial"], annos = adata.obs[["leiden"]], markerscale = 5)
# remove clusters 9, 3, 1
adata_sub = adata[~adata.obs["leiden"].isin(["9", "3", "1", "12"])]
adata_sub.write_h5ad(base_path + dataset_name + "adata_qc_f.h5ad")


# %%
