# In[]
import numpy as np
import pandas as pd
from scipy.io import mmread, mmwrite
from anndata import AnnData
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.sparse import csc_matrix, csr_matrix
import warnings
warnings.filterwarnings("ignore")

cellranger_dir = "../dataset/data_scrnaseq/Cell_Ranger_output/"
raw_dir = "../dataset/data_scrnaseq/data/"
datadir_m1416_pp18 = cellranger_dir + "M1416-PP18/"
datadir_m1417_pp12 = cellranger_dir + "M1417-PP12/"
datadir_m1436_ppc12 = cellranger_dir + "M1436-PPC12/"
datadir_m1437_ppc18 = cellranger_dir + "M1437-PPC18/"
datadir_m1476_pp = cellranger_dir + "M1476-XPP/"
datadir_m1477_ppc = cellranger_dir + "M1477-XPPC/"


# In[]
# ------------------------------------------------------------------------------------------
#
# 1. Processing raw data and generate raw anndata file
#
# ------------------------------------------------------------------------------------------
if False:
    # already saved
    for data_dir in [datadir_m1417_pp12, datadir_m1416_pp18, datadir_m1436_ppc12, datadir_m1437_ppc18, datadir_m1476_pp, datadir_m1477_ppc]:
        barcodes = pd.read_csv(data_dir + "filtered_feature_bc_matrix/barcodes.tsv", header = None, sep = "\t").values.squeeze()
        features = pd.read_csv(data_dir + "filtered_feature_bc_matrix/features.tsv", header = None, sep = "\t")
        matrix = mmread(data_dir + "filtered_feature_bc_matrix/matrix.mtx")

        adata = AnnData(X = csr_matrix(matrix.T))
        adata.obs.index = barcodes
        adata.var.index = features[1].values.squeeze()
        adata.var["gene_ids"] = features[0].values.squeeze()
        adata.var_names_make_unique()
        adata.write_h5ad(data_dir + "adata_raw.h5ad")

# In[]
# ------------------------------------------------------------------------------------------
#
# 2. Merge and redo the annotation files
# 
# ------------------------------------------------------------------------------------------
adata_m1416_pp18 = sc.read_h5ad(datadir_m1416_pp18 + "adata_raw.h5ad")
adata_m1417_pp12 = sc.read_h5ad(datadir_m1417_pp12 + "adata_raw.h5ad")
adata_m1436_ppc12 = sc.read_h5ad(datadir_m1436_ppc12 + "adata_raw.h5ad")
adata_m1437_ppc18 = sc.read_h5ad(datadir_m1437_ppc18 + "adata_raw.h5ad")
adata_m1476 = sc.read_h5ad(datadir_m1476_pp + "adata_raw.h5ad")
adata_m1477 = sc.read_h5ad(datadir_m1477_ppc + "adata_raw.h5ad")

adata_m1416_pp18.obs["sample"] = "M1416-PP18"
adata_m1416_pp18.obs["age"] = "18wk"
adata_m1416_pp18.obs["condition"] = "intact"
adata_m1416_pp18.obs["genotype"] = "PbCre(+/-),Pten(-/-),P53(-/-)"
# Batch/Submission mix name
adata_m1416_pp18.obs["batch"] = "JYMM93(PE150, FSU) & JYMM97(PE150, FSU)"

adata_m1417_pp12.obs["sample"] = "M1417-PP12"
adata_m1417_pp12.obs["age"] = "12wk"
adata_m1417_pp12.obs["condition"] = "intact"
adata_m1417_pp12.obs["genotype"] = "PbCre(+/-),Pten(-/-),P53(-/-)"
# Batch/Submission mix name
adata_m1417_pp12.obs["batch"] = "JYMM93(PE150, FSU) & JYMM97(PE150, FSU)"

adata_m1436_ppc12.obs["sample"] = "M1436-PPC12"
adata_m1436_ppc12.obs["age"] = "12wk"
adata_m1436_ppc12.obs["condition"] = "intact"
adata_m1436_ppc12.obs["genotype"] = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"
# Batch/Submission mix name
adata_m1436_ppc12.obs["batch"] = "JYMM98(PE150, FSU)"

adata_m1437_ppc18.obs["sample"] = "M1437-PPC18"
adata_m1437_ppc18.obs["age"] = "18wk"
adata_m1437_ppc18.obs["condition"] = "intact"
adata_m1437_ppc18.obs["genotype"] = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"
# Batch/Submission mix name
adata_m1437_ppc18.obs["batch"] = "JYMM98(PE150, FSU)"

adata_m1476.obs["sample"] = "M1476-PP"
adata_m1476.obs["condition"] = "castrated"
adata_m1476.obs["genotype"] = "PbCre(+/-),Pten(-/-),P53(-/-)"

adata_m1477.obs["sample"] = "M1477-PPC"
adata_m1477.obs["condition"] = "castrated"
adata_m1477.obs["genotype"] = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"

adata_intact = AnnData.concatenate(adata_m1416_pp18, adata_m1417_pp12, adata_m1436_ppc12, adata_m1437_ppc18, join = "inner")
adata_castrated = AnnData.concatenate(adata_m1476, adata_m1477, join = "inner")
adata_merge = AnnData.concatenate(adata_intact, adata_castrated, join = "inner")


# In[]
# write files, hard to read 
adata_intact.write_h5ad(raw_dir + "raw_data/adata_raw_intact.h5ad")
adata_castrated.write_h5ad(raw_dir + "raw_data/adata_raw_castrated.h5ad")
adata_merge.write_h5ad(raw_dir + "raw_data/adata_raw_merged.h5ad")

# easy to read by R
mmwrite(raw_dir + "raw_data/X_intact.mtx", adata_intact.X)
adata_intact.obs.to_csv(raw_dir + "raw_data/meta_intact.csv")
adata_intact.var.to_csv(raw_dir + "raw_data/gene_intact.csv")

mmwrite(raw_dir + "raw_data/X_castrated.mtx", adata_castrated.X)
adata_castrated.obs.to_csv(raw_dir + "raw_data/meta_castrated.csv")
adata_castrated.var.to_csv(raw_dir + "raw_data/gene_castrated.csv")


# In[]
# ------------------------------------------------------------------------------------------
#
# 3. Quality Control, Intact and Castrated samples are preprocessed separately. Intact:
#
# ------------------------------------------------------------------------------------------

print(adata_intact)
adata_intact.var_names_make_unique() # gene_ids are still unique
# # check the top expressed genes
# sc.pl.highest_expr_genes(adata_intact, n_top=20)

# basic filtering, remove genes and cells with low quality
sc.pp.filter_cells(adata_intact, min_genes=200)
sc.pp.filter_genes(adata_intact, min_cells=3)
print(adata_intact)

# check the mitochondrial genes
adata_intact.var['mt'] = adata_intact.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata_intact, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.scatter(adata_intact, x='total_counts', y='pct_counts_mt')
# n_genes_by_counts: the number of genes expressed in the count matrix
sc.pl.scatter(adata_intact, x='total_counts', y='n_genes_by_counts')
sc.pl.violin(adata_intact, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

# In[]
# remove cells with too many mitochondrial genes or too many total counts
adata_intact = adata_intact[adata_intact.obs.n_genes_by_counts < 10000, :]
adata_intact = adata_intact[adata_intact.obs.pct_counts_mt < 20, :]
adata_intact.write_h5ad(raw_dir + "qc_data/adata_qc_intact.h5ad")

mmwrite(raw_dir + "qc_data/X_intact.mtx", adata_intact.X)
adata_intact.obs.to_csv(raw_dir + "qc_data/meta_intact.csv")
adata_intact.var.to_csv(raw_dir + "qc_data/gene_intact.csv")

# In[]
# ------------------------------------------------------------------------------------------
#
# 4. Quality Control, Intact and Castrated samples are preprocessed separately. Castrated:
#
# ------------------------------------------------------------------------------------------

print(adata_castrated)
adata_castrated.var_names_make_unique() # gene_ids are still unique
# # check the top expressed genes
# sc.pl.highest_expr_genes(adata_castrated, n_top=20)

# basic filtering, remove genes and cells with low quality
sc.pp.filter_cells(adata_castrated, min_genes=200)
sc.pp.filter_genes(adata_castrated, min_cells=3)
print(adata_castrated)

# check the mitochondrial genes
adata_castrated.var['mt'] = adata_castrated.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata_castrated, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.scatter(adata_castrated, x='total_counts', y='pct_counts_mt')
# n_genes_by_counts: the number of genes expressed in the count matrix
sc.pl.scatter(adata_castrated, x='total_counts', y='n_genes_by_counts')
sc.pl.violin(adata_castrated, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

# In[]
# remove cells with too many mitochondrial genes or too many total counts
adata_castrated = adata_castrated[adata_castrated.obs.n_genes_by_counts < 10000, :]
adata_castrated = adata_castrated[adata_castrated.obs.pct_counts_mt < 20, :]
adata_castrated.write_h5ad(raw_dir + "qc_data/adata_qc_castrated.h5ad")

mmwrite(raw_dir + "qc_data/X_castrated.mtx", adata_castrated.X)
adata_castrated.obs.to_csv(raw_dir + "qc_data/meta_castrated.csv")
adata_castrated.var.to_csv(raw_dir + "qc_data/gene_castrated.csv")

# # In[]
# # ---------------------------------------------------------------
# #
# # Sanity check
# #
# # ---------------------------------------------------------------
# adata_merge = AnnData.concatenate(adata_intact, adata_castrated, join = "inner")

# # log-normalize the gene expression data
# sc.pp.normalize_total(adata_merge, target_sum=1e4)
# sc.pp.log1p(adata_merge)

# # keep only highly variable genes
# sc.pp.highly_variable_genes(adata_merge, n_top_genes = 2000)
# adata_merge.raw = adata_merge
# adata_merge = adata_merge[:, adata_merge.var.highly_variable]

# # In[]
# sc.pp.neighbors(adata_merge, n_neighbors=30, n_pcs=100)
# sc.tl.umap(adata_merge)
# sc.pl.umap(adata_merge, color = ["sample"])

# # In[]
# # load cell type annotations
# adata_intact_seurat = sc.read_h5ad(raw_dir + "adata_seurat.h5ad")
# annot_intact = adata_intact_seurat.obs[["annot"]]
# adata_intact.obs["annot"] = annot_intact.loc[adata_intact.obs.index,"annot"]

# sc.pl.umap(adata_intact, color = ["sample", "annot"])

# adata_castrated_seurat = sc.read_h5ad(raw_dir + "adata_castrated_seurat.h5ad")
# annot_merge = np.concatenate([adata_intact_seurat.obs["annot"].values, adata_castrated_seurat.obs["annot_transfer"].values], axis = 0)
# annot_merge = pd.DataFrame(index = np.concatenate([adata_intact_seurat.obs.index.values, adata_castrated_seurat.obs.index.values], axis = 0), data = annot_merge, columns = ["annot"])


# %%
