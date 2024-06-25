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


# In[]
# --------------------------------------------------------
#
# Load data
#
# --------------------------------------------------------

base_path = "../dataset/data_merfish/"
# PPC15
dataset_name = "region_0/"
# PP15
# dataset_name = "region_1/"

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

# In[]
# --------------------------------------------------------
#
# Create AnnData
#
# --------------------------------------------------------
# import squidpy as sq
# adata_sq = sq.read.vizgen(
#     base_path + dataset_name,
#     counts_file="cell_by_gene.csv",
#     meta_file="cell_metadata.csv",
# )

# create anndata
adata = AnnData(X = cell_by_gene, obs = meta_cell, var = meta_gene)
# assign spatial location of cells
adata.obsm["spatial"] = adata.obs[["center_x", "center_y"]].values
adata.obs.drop(columns=["center_x", "center_y"], inplace=True)

adata.write_h5ad(base_path + dataset_name + "adata.h5ad")



# In[]
# --------------------------------------------------------
#
# Basic preprocessing, adopted from:
# https://squidpy.readthedocs.io/en/latest/notebooks/tutorials/tutorial_vizgen_mouse_liver.html
#
# --------------------------------------------------------
base_path = "../dataset/data_merfish/"
dataset_name = "region_1/"
adata = sc.read_h5ad(base_path + dataset_name + "adata.h5ad")
print(adata)

# check library size
libsize = adata.obs["barcodeCount"].values.squeeze()
import matplotlib.pyplot as plt
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


# 1. Annotate cell type:
# 1.1. label transfer
# 1.2. annotate marker genes
# 1.3. adjust umap parameters

# 2. Region 0: a layer or multiple layers

# %%
