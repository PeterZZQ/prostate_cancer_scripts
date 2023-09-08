# In[]
import scanpy as sc
import scvi
import numpy as np 
import pandas as pd 
from scvi.model.utils import mde

# In[]
adata = sc.read_h5ad("Cell_Ranger_output/adata_merge.h5ad")
# sc.pp.filter_cells(adata, min_genes = 200)
sc.pp.filter_genes(adata, min_cells = 3)
# keep the original count before normalization and log transform
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
# store the normalized values in `.raw` before filtering
adata.raw = adata

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()

# In[]
'''
sc.pp.neighbors(adata, n_neighbors = 30)
sc.tl.umap(adata)

sc.pl.umap(adata, color = ["sample"])
sc.pl.umap(adata, color = ["age"])
sc.pl.umap(adata, color = ["genotype"])
'''

# In[]
adata.obsm["X_scVI"] = vae.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.leiden(adata)

# In[]
adata.obsm["X_mde"] = mde(adata.obsm["X_scVI"])
sc.pl.embedding(
    adata,
    basis="X_mde",
    color=["sample", "leiden"],
    frameon=False,
    ncols=1,
)
# %%
