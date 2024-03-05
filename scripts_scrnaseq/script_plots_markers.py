# In[]
import sys
sys.path.append(".")
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc
import anndata
import scipy as sci
import utils
import warnings
warnings.filterwarnings("ignore")

raw_dir = "../dataset/data_scrnaseq/data/qc_data/"
seurat_dir = "../dataset/data_scrnaseq/seurat_integration/"
vis_intact_dir = "../results_scrnaseq/annot_intact/"
vis_cas_dir = "../results_scrnaseq/annot_cas/"

# In[]
# ---------------------------------------------------------------- #
#
# Preprocessing, read in the intact and castrated dataset
#
# ---------------------------------------------------------------- #

# Load anndata for intact
adata_intact = sc.read_h5ad(raw_dir + "adata_qc_intact.h5ad")

# Load seurat embeddings
cca_umap_intact = pd.read_csv(seurat_dir + "cca_umap_intact.csv", sep = "\t")
cca_pca_intact = pd.read_csv(seurat_dir + "cca_pca_intact.csv", sep = "\t")
meta_intact_seurat = pd.read_csv(seurat_dir + "meta_intact.csv", sep = "\t")

adata_intact.obsm["X_seurat_pca"] = cca_pca_intact.loc[adata_intact.obs.index,:].values
adata_intact.obsm["X_seurat_umap"] = cca_umap_intact.loc[adata_intact.obs.index,:].values
adata_intact.obs["seurat_cluster"] = meta_intact_seurat.loc[adata_intact.obs.index,"seurat_clusters"].astype("category")

# castrated sample is not in the scope
# # Load anndata for castrated
# adata_castrated = sc.read_h5ad(raw_dir + "adata_qc_castrated.h5ad")
# # adata_castrated.obs["label_old"] = sc.read_h5ad("../dataset/data_scrnaseq/data/adata_castrated_seurat.h5ad").obs.loc[adata_castrated.obs.index.values,"annot_transfer"].values
# # adata_castrated.write_h5ad(raw_dir + "adata_qc_castrated.h5ad")
# meta_castrated_seurat = pd.read_csv(seurat_dir + "meta_castrated_transfer.csv", sep = "\t")
# adata_castrated.obs["seurat_cluster"] = meta_castrated_seurat.loc[:, "seurat_cluster"].astype("category")


# In[]
# ---------------------------------------------------------------- #
#
# Log-transformation and UMAP calculation
#
# ---------------------------------------------------------------- #

# not using castrated sample
# adata_merge = anndata.concat([adata_intact, adata_castrated], join = "inner")
adata_merge = adata_intact.copy()

sc.pp.normalize_total(adata_merge, target_sum=1e4)
sc.pp.log1p(adata_merge)
# automatically used for umap calculation
sc.pp.highly_variable_genes(adata_merge, n_top_genes = 2000)
sc.pp.neighbors(adata_merge, n_neighbors = 100, n_pcs = 100)
sc.tl.umap(adata_merge, min_dist = 0.5)

X_umap_intact = adata_merge[adata_merge.obs["condition"] == "intact",:].obsm["X_umap"]
X_umap_intact = pd.DataFrame(data = X_umap_intact, index = adata_merge.obs[adata_merge.obs["condition"] == "intact"].index, columns = ["X_umap1", "X_umap2"])
X_umap_intact.to_csv(seurat_dir + "umap_intact.csv", sep = "\t")

# not using castrated sample
# X_umap_cas = adata_merge[adata_merge.obs["condition"] == "castrated",:].obsm["X_umap"]
# X_umap_cas = pd.DataFrame(data = X_umap_cas, index = adata_merge.obs[adata_merge.obs["condition"] == "castrated"].index, columns = ["X_umap1", "X_umap2"])
# X_umap_cas.to_csv(seurat_dir + "umap_cas.csv", sep = "\t")
# # ---------------------------------------------------------------- #
# #
# # Log-normalization
# #
# # ---------------------------------------------------------------- #
# adata_intact_raw = adata_intact.copy()
# adata_castrated_raw = adata_castrated.copy()
# adata_intact_raw.obsm["X_umap"] = X_umap_intact.loc[adata_intact_raw.obs.index,:].values
# adata_castrated_raw.obsm["X_umap"] = X_umap_cas.loc[adata_castrated_raw.obs.index,:].values

# sc.pp.normalize_total(adata_intact, target_sum=1e4)
# sc.pp.log1p(adata_intact)

# sc.pp.normalize_total(adata_castrated, target_sum=1e4)
# sc.pp.log1p(adata_castrated)


# In[]
# ---------------------------------------------------------------- #
#
# Plot the marker gene expression [Log-normalized]
#
# ---------------------------------------------------------------- #

# read the cell markers
# markers = pd.read_csv("../dataset/markers_info/Cell annotation gene markers.csv", sep = ",")

# selected markers
markers = {}
markers["luminal"] = ["Ar", "Krt8", "Cd24a", "Krt18", "Spink1"]
markers["basal"] = ["Trp63", "Krt5", "Krt14"]
markers["club_epithelia"] = ["Agr2", "Krt7"]
markers["endothelial"] = ["Ackr1", "Cldn5"]
markers["lymphoid"] = ["Cd3e", "Ms4a1", "Klrb1c"]
markers["myeloid"] = ["Ptprc", "Itgam"]
markers["monocytes"] = ["Ptprc", "Itgam", "Cd14", "S100a8", "S100a9"] 
markers["macrophage"] = ["Ptprc", "Itgam", "Adgre1"]
markers["macrophage_m1"] = ["Ptprc", "Itgam", "Adgre1", "Cd68", "Nos2"]
markers["macrophage_m2"] = ["Ptprc", "Itgam", "Adgre1", "Mrc1", "Arg1"]
markers["mesenchymal"] = ["Fgf10", "Rorb", "Rspo1", "Sult1e1", "Wnt10a", "Wnt2", "Wnt6"]
markers["sv"] = ["Pate4", "Pax2", "Svs2"]

# Plot the heatmap of gene expression on seurat integrated space
# temporary: for ease of plot umap
adata_intact.obsm["X_umap"] = adata_intact.obsm["X_seurat_umap"]

for ct in markers.keys():
    fig = plt.figure(figsize = (len(markers[ct]) * 7, 5))
    axs = fig.subplots(nrows = 1, ncols = len(markers[ct]))
    for idx, marker in enumerate(markers[ct]):
        sc.pl.umap(adata_intact, color = marker, color_map = utils.SUPER_MAGMA, ax = axs[idx], show = False)
    # fig.suptitle(ct, fontsize = 30)
    fig.savefig(vis_intact_dir + f"markers_{ct}.png", bbox_inches = "tight", dpi = 150)


# # run UMAP for castrated
# sc.pp.neighbors(adata_castrated)
# sc.tl.umap(adata_castrated)
# for ct in markers.keys():
#     fig = plt.figure(figsize = (len(markers[ct]) * 7, 10))
#     axs = fig.subplots(nrows = 2, ncols = len(markers[ct]))
#     for idx, marker in enumerate(markers[ct]):
#         sc.pl.umap(adata_castrated[adata_castrated.obs["sample"] == "M1476-PP",:], color = marker, color_map = utils.SUPER_MAGMA, ax = axs[0,idx], show = False)
#         axs[0,idx].set_title(f"{marker} (M1476-PP)")
#         sc.pl.umap(adata_castrated[adata_castrated.obs["sample"] == "M1477-PPC",:], color = marker, color_map = utils.SUPER_MAGMA, ax = axs[1,idx], show = False)
#         axs[1,idx].set_title(f"{marker} (M1477-PPC)")
#     # fig.suptitle(ct, fontsize = 30)
#     fig.savefig(vis_cas_dir + f"markers_{ct}.png", bbox_inches = "tight", dpi = 150)


# In[]
# ---------------------------------------------------------------- #
#
# Annotate cell types
#
# ---------------------------------------------------------------- #

# Summary:
# Epithelial, Luminal (Spink1-): cluster 0, 2, 3, 5, 6, 8, 9, 11?, 17, 18?
# Epithelial, Luminal (Spink1+): cluster 9
# Epithelial, Basal: cluster 1, 11? 6? (cluster 11 is a mixture of basal and luminal)
# Epithelial, Club Epithelia: 18
# Endothelial: 15
# Lymphoid: 13
# Myeloid: 4, 7, 12
# Monocytes: 4
# Macrophage: 7, 12
# M1 macrophage:
# M2 macrophage:
# Mesenchymal: 10
# SV: 14


for adata in [adata_intact]:
    adata.obs["annot"] = "Other"
    adata.obs.loc[adata.obs["seurat_cluster"].isin([0, 2, 3, 6, 8, 9, 17, 18]), "annot"] = "Luminal"
    adata.obs.loc[adata.obs["seurat_cluster"] == 9, "annot"] = "Luminal (Spink1+)"
    adata.obs.loc[adata.obs["seurat_cluster"] == 1, "annot"] = "Basal"
    # adata.obs.loc[adata.obs["seurat_cluster"] == 11, "annot"] = "Luminal & Basal"
    adata.obs.loc[adata.obs["seurat_cluster"] == 11, "annot"] = "Luminal"
    adata.obs.loc[adata.obs["seurat_cluster"] == 5, "annot"] = "Club epithelia"
    adata.obs.loc[adata.obs["seurat_cluster"] == 15, "annot"] = "Endothelial"
    adata.obs.loc[adata.obs["seurat_cluster"] == 13, "annot"] = "Lymphoid"
    adata.obs.loc[adata.obs["seurat_cluster"] == 4, "annot"] = "Monocytes"
    adata.obs.loc[adata.obs["seurat_cluster"].isin([7,12]), "annot"] = "Macrophage"
    adata.obs.loc[adata.obs["seurat_cluster"] == 10, "annot"] = "Mesenchymal"
    adata.obs.loc[adata.obs["seurat_cluster"].isin([14, 16]), "annot"] = "SV"
    adata.obs["annot"].astype("category")


# In[]
adata_intact.obsm["X_umap"] = X_umap_intact.loc[adata_intact.obs.index.values,:].values
# adata_castrated.obsm["X_umap"] = X_umap_cas.loc[adata_castrated.obs.index.values,:].values
# adata_intact_raw.obs["annot"] = adata_intact.obs["annot"].values
# adata_castrated_raw.obs["annot"] = adata_castrated.obs["annot"].values

adata_intact_raw.write_h5ad(seurat_dir + "adata_intact_seurat.h5ad")
adata_castrated_raw.write_h5ad(seurat_dir + "adata_castrated_seurat.h5ad")


# In[]
# seurat_labels_old = pd.read_csv("../dataset/data_scrnaseq/data/seurat_meta.txt", sep = "\t")
# adata_intact.obs["label_old"] = seurat_labels_old.loc[adata_intact.obs.index,"annot"].values
plt.rcParams["font.size"] = 20
fig = plt.figure(figsize = (45, 10))
ax = fig.subplots(nrows = 1, ncols = 3)
ax[0] = sc.pl.scatter(adata_intact, basis = "seurat_umap", color = ["sample"], ax = ax[0], show = False)
ax[1] = sc.pl.scatter(adata_intact, basis = "seurat_umap", color = ["seurat_cluster"], ax = ax[1], show = False)
# ax[2] = sc.pl.scatter(adata_intact, basis = "seurat_umap", color = ["label_old"], ax = ax[2], show = False)
ax[2] = sc.pl.scatter(adata_intact, basis = "seurat_umap", color = ["annot"], ax = ax[2], show = False)

fig.tight_layout()
plt.show()
fig.savefig(vis_intact_dir + "annotation_seurat.png", bbox_inches = "tight", dpi = 150)

# In[]
fig = plt.figure(figsize = (40, 20))
ax = fig.subplots(nrows = 2, ncols = 2)
plt.rcParams["font.size"] = 25
ax[0,0] = sc.pl.scatter(adata_intact[adata_intact.obs["sample"] == "M1417-PP12", :], basis = "umap", color = ["annot"], ax = ax[0,0], show = False)
ax[0,1] = sc.pl.scatter(adata_intact[adata_intact.obs["sample"] == "M1416-PP18", :], basis = "umap", color = ["annot"], ax = ax[0,1], show = False)
ax[1,0] = sc.pl.scatter(adata_intact[adata_intact.obs["sample"] == "M1436-PPC12", :], basis = "umap", color = ["annot"], ax = ax[1,0], show = False)
ax[1,1] = sc.pl.scatter(adata_intact[adata_intact.obs["sample"] == "M1437-PPC18", :], basis = "umap", color = ["annot"], ax = ax[1,1], show = False)
ax[0,0].set_title("PP-12wk", fontsize = 35)
ax[0,1].set_title("PP-18wk", fontsize = 35)
ax[1,0].set_title("PPC-12wk", fontsize = 35)
ax[1,1].set_title("PPC-18wk", fontsize = 35)

fig.tight_layout()
plt.show()
fig.savefig(vis_intact_dir + "annotation_umap.png", bbox_inches = "tight", dpi = 150)


# In[]
adata_castrated.obs["annot"] = adata_castrated.obs["annot"].astype("category")
adata_castrated_pp = adata_castrated[adata_castrated.obs["sample"] == "M1476-PP",:]
adata_castrated_ppc = adata_castrated[adata_castrated.obs["sample"] == "M1477-PPC",:]

adata_castrated_pp.obs["label_old"] = adata_castrated_pp.obs["label_old"].cat.set_categories(adata_castrated.obs["label_old"].cat.categories)
adata_castrated_ppc.obs["label_old"] = adata_castrated_ppc.obs["label_old"].cat.set_categories(adata_castrated.obs["label_old"].cat.categories)
adata_castrated_pp.obs["annot"] = adata_castrated_pp.obs["annot"].cat.set_categories(adata_castrated.obs["annot"].cat.categories)
adata_castrated_ppc.obs["annot"] = adata_castrated_ppc.obs["annot"].cat.set_categories(adata_castrated.obs["annot"].cat.categories)

plt.rcParams["font.size"] = 20
fig = plt.figure(figsize = (35, 20))
ax = fig.subplots(nrows = 2, ncols = 2)
ax[0,0] = sc.pl.scatter(adata_castrated_pp, basis = "umap", color = ["seurat_cluster"], ax = ax[0,0], show = False)
# ax[0,1] = sc.pl.scatter(adata_castrated_pp, basis = "umap", color = ["label_old"], ax = ax[0,1], show = False)
ax[1,0] = sc.pl.scatter(adata_castrated_pp, basis = "umap", color = ["annot"], ax = ax[1,0], show = False)
ax[0,0].set_title("Seurat_cluster: M1476-PP", fontsize = 35)
# ax[0,1].set_title("label_old: M1476-PP")
ax[1,0].set_title("Annot: M1476-PP", fontsize = 35)


ax[0,1] = sc.pl.scatter(adata_castrated_ppc, basis = "umap", color = ["seurat_cluster"], ax = ax[0,1], show = False)
# ax[1,1] = sc.pl.scatter(adata_castrated_ppc, basis = "umap", color = ["label_old"], ax = ax[1,1], show = False)
ax[1,1] = sc.pl.scatter(adata_castrated_ppc, basis = "umap", color = ["annot"], ax = ax[1,1], show = False)
ax[0,1].set_title("Seurat_cluster: M1477-PPC", fontsize = 35)
# ax[1,1].set_title("label_old: M1477-PPC")
ax[1,1].set_title("Annot: M1477-PPC", fontsize = 35)

fig.tight_layout()
plt.show()
fig.savefig(vis_cas_dir + "annotation_umap.png", bbox_inches = "tight", dpi = 150)


# In[]
# ---------------------------------------------------------------- #
#
# Additional interest genes: ``interest genes.csv'', include both scrna-seq genes and visium genes
#
# ---------------------------------------------------------------- #
#
# Trp53, Pten, Ackr3: perturbed genes
# Cxcr4: interacts
# Cxcl12, Mif: potential ligand
# Spink1: affected
# Mki67: proliferation

# markers for other cell types
# Ptprc: immune cells
# Cd4: CD4+ T cells
# Cd8a: CD8+ T cells
# Foxp3: Treg T cells
# Ctla4: T cells
# Itgax: Myeloid cells (DCs)
# Il2, Il4, Il6: Interleukin
# Cxcl2, Cxcl12: Chemokine Ligand
# Tmprss4: Transmembrane Serine Protease
# Syp, Chga, Sox2, Mycn: NE 
# Dcn, FN1: Fibroblasts
# Nkx3-1, Pbsn: L1
# LGR5: Stem cells
# Ccl21: Endothelial cells
# Cxcr5, Ms4a1: B-cells 
# Klrd1, Klrb1c: NK cell
# Gzmb, Prf1, Tnfa, Ifng, Tgfb, Il6: Cytotoxic

interest_genes_intact_dir = "../results_scrnaseq/interest_genes_intact/"
interest_genes_cas_dir = "../results_scrnaseq/interest_genes_cas/"

interest_genes = pd.read_csv("../dataset/markers_info/interest_genes.csv", index_col = 0).index.values
functions = pd.read_csv("../dataset/markers_info/interest_genes.csv", index_col = 0)["Note"].values
interest_genes = [x for x in interest_genes if (x != "Cxcl11") & (x != "Tnfa") & (x != "Tgfb") & (x != "Il12") & (x != "Pd1") & (x != "Tim3") & (x != "Ccl21")]
full_names = [x + "_" + y for x,y in zip(interest_genes, functions) if (x != "Cxcl11") & (x != "Tnfa") & (x != "Tgfb") & (x != "Il12") & (x != "Pd1") & (x != "Tim3") & (x != "Ccl21")]

for gene, file in zip(interest_genes, full_names):
    fig = plt.figure(figsize = (7, 5))
    ax = fig.subplots(nrows = 1, ncols = 1)
    sc.pl.umap(adata_intact, color = gene, color_map = utils.SUPER_MAGMA, ax = ax, show = False)
    fig.savefig(interest_genes_intact_dir + f"{file}.png", bbox_inches = "tight", dpi = 150)


for gene, file in zip(interest_genes, full_names):
    if gene != "Il13":
        fig = plt.figure(figsize = (14, 5))
        ax = fig.subplots(nrows = 1, ncols = 2)
        sc.pl.umap(adata_castrated_pp, color = gene, color_map = utils.SUPER_MAGMA, ax = ax[0], show = False)
        sc.pl.umap(adata_castrated_ppc, color = gene, color_map = utils.SUPER_MAGMA, ax = ax[1], show = False)
        fig.savefig(interest_genes_cas_dir + f"{file}.png", bbox_inches = "tight", dpi = 150)

# %%
