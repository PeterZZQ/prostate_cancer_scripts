# In[]
import sys
sys.path.append("../")
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc
import anndata
import scipy as sci
import seaborn as sns

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
# 1. Preprocessing, read in the intact dataset [analysis of castrated dataset already removed]
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


# In[]
# ---------------------------------------------------------------- #
#
# 2. Log-transformation and UMAP calculation
#
# ---------------------------------------------------------------- #
sc.pp.normalize_total(adata_intact, target_sum=1e4)
sc.pp.log1p(adata_intact)

# automatically used for umap calculation
sc.pp.highly_variable_genes(adata_intact, n_top_genes = 2000)
adata_merge = adata_intact[:, adata_intact.var.highly_variable]

# calculate UMAP embedding
sc.pp.neighbors(adata_merge, n_neighbors = 50, n_pcs = 100)
sc.tl.umap(adata_merge, min_dist = 0.5)
X_umap_intact = pd.DataFrame(data = adata_merge.obsm["X_umap"], index = adata_merge.obs.index, columns = ["X_umap1", "X_umap2"])
# X_umap_intact.to_csv(seurat_dir + "umap_intact.csv", sep = "\t")

# In[]
# ---------------------------------------------------------------- #
#
# 3. Plot the marker gene expression [Log-normalized]
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

# Plot the heatmap of gene expression on seurat integrated space, temporary: for ease of plot umap
key = "seurat"
if key == "seurat":
    adata_intact.obsm["X_umap"] = adata_intact.obsm["X_seurat_umap"]
else:
    adata_intact.obsm["X_umap"] = X_umap_intact.values

for ct in markers.keys():
    fig = plt.figure(figsize = (len(markers[ct]) * 7, 5))
    axs = fig.subplots(nrows = 1, ncols = len(markers[ct]))
    for idx, marker in enumerate(markers[ct]):
        # adata_intact is already log-normalized
        sc.pl.umap(adata_intact, color = marker, color_map = utils.SUPER_MAGMA, ax = axs[idx], show = False)
    # fig.suptitle(ct, fontsize = 30)
    fig.savefig(vis_intact_dir + f"markers_{ct}.png", bbox_inches = "tight", dpi = 150)


# In[]
# ---------------------------------------------------------------- #
#
# 4. Annotate cell types
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


adata_intact.obsm["X_umap"] = X_umap_intact.values
adata_intact.write_h5ad(seurat_dir + "adata_intact_seurat.h5ad")

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
fig = plt.figure(figsize = (30, 20))
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
# ---------------------------------------------------------------- #
#
# 5. Additional interest genes: ``interest genes.csv'', include both scrna-seq genes and visium genes
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

# Heatmap

interest_genes = pd.read_csv("../dataset/markers_info/interest_genes.csv", index_col = 0).index.values
functions = pd.read_csv("../dataset/markers_info/interest_genes.csv", index_col = 0)["Note"].values
interest_genes = [x for x in interest_genes if (x != "Cxcl11") & (x != "Tnfa") & (x != "Tgfb") & (x != "Il12") & (x != "Pd1") & (x != "Tim3") & (x != "Ccl21")]
full_names = [x + "_" + y for x,y in zip(interest_genes, functions) if (x != "Cxcl11") & (x != "Tnfa") & (x != "Tgfb") & (x != "Il12") & (x != "Pd1") & (x != "Tim3") & (x != "Ccl21")]

# In[]
interest_genes_intact_dir = "../results_scrnaseq/interest_genes_intact/heatmap/"

for gene, file in zip(interest_genes, full_names):
    fig = plt.figure(figsize = (30, 20))
    ax = fig.subplots(nrows = 2, ncols = 2)
    # adata_intact is already log-normalized
    
    vmax = np.max(adata_intact[:, gene].X.todense())
    vmin = 0.0

    ax[0,0] = sc.pl.umap(adata_intact[adata_intact.obs["sample"] == "M1417-PP12", :], color = gene, color_map = utils.SUPER_MAGMA, ax = ax[0,0], show = False, vmin = vmin, vmax = vmax)
    ax[0,1] = sc.pl.umap(adata_intact[adata_intact.obs["sample"] == "M1416-PP18", :], color = gene, color_map = utils.SUPER_MAGMA, ax = ax[0,1], show = False, vmin = vmin, vmax = vmax)
    ax[1,0] = sc.pl.umap(adata_intact[adata_intact.obs["sample"] == "M1436-PPC12", :], color = gene, color_map = utils.SUPER_MAGMA, ax = ax[1,0], show = False, vmin = vmin, vmax = vmax)
    ax[1,1] = sc.pl.umap(adata_intact[adata_intact.obs["sample"] == "M1437-PPC18", :], color = gene, color_map = utils.SUPER_MAGMA, ax = ax[1,1], show = False, vmin = vmin, vmax = vmax)
    ax[0,0].set_title("PP-12wk", fontsize = 35)
    ax[0,1].set_title("PP-18wk", fontsize = 35)
    ax[1,0].set_title("PPC-12wk", fontsize = 35)
    ax[1,1].set_title("PPC-18wk", fontsize = 35)


    fig.tight_layout()
    plt.show()

    fig.savefig(interest_genes_intact_dir + f"{file}.png", bbox_inches = "tight", dpi = 150)


# In[]
# Boxplot of expression level
interest_genes_intact_dir = "../results_scrnaseq/interest_genes_intact/boxplot/"

for gene, file in zip(interest_genes, full_names):
    fig = plt.figure(figsize = (20, 14))
    ax = fig.subplots(nrows = 2, ncols = 1)

    # only on luminal cells?
    expr_df = pd.DataFrame(columns = ["Expr", "annot", "genotype", "age"])
    expr_df["Expr"] = adata_intact[:,gene].X.toarray().squeeze()
    expr_df[["annot", "genotype", "age"]] = adata_intact.obs[["annot", "genotype", "age"]].values
    sns.boxplot(data = expr_df[expr_df["age"] == "12wk"], x = "annot", y = "Expr", hue = "genotype", ax =ax[0])
    leg = ax[0].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
    ax[0].set_title("12wk", fontsize = 35)
    ax[0].tick_params(axis='x', labelrotation=45)
    ax[0].set_xlabel(None)

    sns.boxplot(data = expr_df[expr_df["age"] == "18wk"], x = "annot", y = "Expr", hue = "genotype", ax =ax[1])
    leg = ax[1].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
    ax[1].set_title("18wk", fontsize = 35)
    ax[1].tick_params(axis='x', labelrotation=45)
    ax[1].set_xlabel(None)

    fig.suptitle(gene, fontsize = 40)

    fig.tight_layout()
    plt.show()

    fig.savefig(interest_genes_intact_dir + f"{file}.png", bbox_inches = "tight", dpi = 150)
    
# In[]
# Percentage of cells with expressed genes 
interest_genes_intact_dir = "../results_scrnaseq/interest_genes_intact/percent_expr/"

plt.rcParams["font.size"] = 15

for gene, file in zip(interest_genes, full_names):

    # only on luminal cells?
    expr_df = pd.DataFrame(columns = ["Expr", "annot", "genotype", "age"])
    expr_df["Expr"] = adata_intact[:,gene].X.toarray().squeeze()
    expr_df[["annot", "genotype", "age"]] = adata_intact.obs[["annot", "genotype", "age"]].values

    # calculate the percentage of cells expressing the gene
    expr_df_pp12 = expr_df[(expr_df["age"] == "12wk") & (expr_df["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")]
    expr_df_ppc12 = expr_df[(expr_df["age"] == "12wk") & (expr_df["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)")]
    expr_df_pp18 = expr_df[(expr_df["age"] == "18wk") & (expr_df["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")]
    expr_df_ppc18 = expr_df[(expr_df["age"] == "18wk") & (expr_df["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)")]

    percent_expr_pp12 = np.sum(expr_df_pp12["Expr"].values != 0)/expr_df_pp12.shape[0]
    percent_expr_ppc12 = np.sum(expr_df_ppc12["Expr"].values != 0)/expr_df_ppc12.shape[0]
    percent_expr_pp18 = np.sum(expr_df_pp18["Expr"].values != 0)/expr_df_pp18.shape[0]
    percent_expr_ppc18 = np.sum(expr_df_ppc18["Expr"].values != 0)/expr_df_ppc18.shape[0]

    # calculate the percentage of immune cells expressing the gene
    expr_df = expr_df[expr_df["annot"].isin(["Lymphoid", "Monocytes", "Macrophage"])]
    expr_df_pp12 = expr_df[(expr_df["age"] == "12wk") & (expr_df["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")]
    expr_df_ppc12 = expr_df[(expr_df["age"] == "12wk") & (expr_df["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)")]
    expr_df_pp18 = expr_df[(expr_df["age"] == "18wk") & (expr_df["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")]
    expr_df_ppc18 = expr_df[(expr_df["age"] == "18wk") & (expr_df["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)")]

    percent_expr_immune_pp12 = np.sum(expr_df_pp12["Expr"].values != 0)/expr_df_pp12.shape[0]
    percent_expr_immune_ppc12 = np.sum(expr_df_ppc12["Expr"].values != 0)/expr_df_ppc12.shape[0]
    percent_expr_immune_pp18 = np.sum(expr_df_pp18["Expr"].values != 0)/expr_df_pp18.shape[0]
    percent_expr_immune_ppc18 = np.sum(expr_df_ppc18["Expr"].values != 0)/expr_df_ppc18.shape[0]

    fig = plt.figure(figsize = (12,12))
    axs = fig.subplots(nrows = 2, ncols = 2)
    axs[0,0].pie([percent_expr_pp12, 1-percent_expr_pp12], labels = ["Expressed", None], autopct='%.2f%%')
    axs[0,1].pie([percent_expr_ppc12, 1-percent_expr_ppc12], labels = ["Expressed", None], autopct='%.2f%%')
    axs[1,0].pie([percent_expr_pp18, 1-percent_expr_pp18], labels = ["Expressed", None], autopct='%.2f%%')
    axs[1,1].pie([percent_expr_ppc18, 1-percent_expr_ppc18], labels = ["Expressed", None], autopct='%.2f%%')

    axs[0,0].set_title("CTRL (12wk)")
    axs[0,1].set_title("CXCR7 KO (12wk)")
    axs[1,0].set_title("CTRL (18wk)")
    axs[1,1].set_title("CXCR7 KO (18wk)")

    fig.suptitle(f"{gene}: percentage of expression", fontsize = 30)
    fig.tight_layout()
    fig.savefig(interest_genes_intact_dir + f"all/{file}.png", bbox_inches = "tight", dpi = 150)

    fig = plt.figure(figsize = (12,12))
    axs = fig.subplots(nrows = 2, ncols = 2)
    axs[0,0].pie([percent_expr_immune_pp12, 1-percent_expr_immune_pp12], labels = ["Expressed", None], autopct='%.2f%%')
    axs[0,1].pie([percent_expr_immune_ppc12, 1-percent_expr_immune_ppc12], labels = ["Expressed", None], autopct='%.2f%%')
    axs[1,0].pie([percent_expr_immune_pp18, 1-percent_expr_immune_pp18], labels = ["Expressed", None], autopct='%.2f%%')
    axs[1,1].pie([percent_expr_immune_ppc18, 1-percent_expr_immune_ppc18], labels = ["Expressed", None], autopct='%.2f%%')

    axs[0,0].set_title("CTRL (12wk)")
    axs[0,1].set_title("CXCR7 KO (12wk)")
    axs[1,0].set_title("CTRL (18wk)")
    axs[1,1].set_title("CXCR7 KO (18wk)")

    fig.suptitle(f"{gene}: percentage of immune expression", fontsize = 30)
    fig.tight_layout()
    fig.savefig(interest_genes_intact_dir + f"immune/immune_{file}.png", bbox_inches = "tight", dpi = 150)


# In[]
# ---------------------------------------------------------------- #
#
# 6. Plot cell percentage
#
# ---------------------------------------------------------------- #
adata_intact = sc.read_h5ad(seurat_dir + "adata_intact_seurat.h5ad")
adata_intact.obs["annot"] = adata_intact.obs["annot"].astype("str")
adata_intact.obs.loc[adata_intact.obs["annot"].isin(["Lymphoid", "Macrophage", "Monocytes"]),"annot"] = "Immune"
adata_intact.obs.loc[adata_intact.obs["annot"].isin(["Endothelial", "SV"]),"annot"] = "Other"
# immune cells count under 12wk
adata_ctr_12wk = adata_intact[(adata_intact.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata_intact.obs["age"] == "12wk"),:]
adata_ko_12wk = adata_intact[(adata_intact.obs["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata_intact.obs["age"] == "12wk"),:]
# immune cells count under 18wk
adata_ctr_18wk = adata_intact[(adata_intact.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata_intact.obs["age"] == "18wk"),:]
adata_ko_18wk = adata_intact[(adata_intact.obs["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata_intact.obs["age"] == "18wk"),:]



# In[]
# 1. Immune cell percentage
n_immune_ctr_12wk = np.sum(adata_ctr_12wk.obs["annot"] == "Immune")
n_other_ctr_12wk = adata_ctr_12wk.shape[0] - n_immune_ctr_12wk
n_immune_ko_12wk = np.sum(adata_ko_12wk.obs["annot"] == "Immune")
n_other_ko_12wk = adata_ko_12wk.shape[0] - n_immune_ko_12wk
n_immune_ctr_18wk = np.sum(adata_ctr_18wk.obs["annot"] == "Immune")
n_other_ctr_18wk = adata_ctr_18wk.shape[0] - n_immune_ctr_18wk
n_immune_ko_18wk = np.sum(adata_ko_18wk.obs["annot"] == "Immune")
n_other_ko_18wk = adata_ko_18wk.shape[0] - n_immune_ko_18wk


plt.rcParams["font.size"] = 15
fig = plt.figure(figsize = (12,12))
axs = fig.subplots(nrows = 2, ncols = 2)
axs[0,0].pie([n_immune_ctr_12wk, n_other_ctr_12wk], labels = ["Immune", "Other"], autopct='%.0f%%')
axs[0,1].pie([n_immune_ko_12wk, n_other_ko_12wk], labels = ["Immune", "Other"], autopct='%.0f%%')
axs[1,0].pie([n_immune_ctr_18wk, n_other_ctr_18wk], labels = ["Immune", "Other"], autopct='%.0f%%')
axs[1,1].pie([n_immune_ko_18wk, n_other_ko_18wk], labels = ["Immune", "Other"], autopct='%.0f%%')
axs[0,0].set_title("CTRL (12wk)")
axs[0,1].set_title("CXCR7 KO (12wk)")
axs[1,0].set_title("CTRL (18wk)")
axs[1,1].set_title("CXCR7 KO (18wk)")
fig.suptitle("Immune cell percentage", fontsize = 30)

fig.savefig("../results_scrnaseq/percent_immune_cell.png", bbox_inches = "tight")

# In[]
# 2. all cell type percentage
unique_ct_ctrl_12wk, counts_ct_ctrl_12wk = np.unique(adata_ctr_12wk.obs["annot"].values, return_counts = True)
unique_ct_ctrl_18wk, counts_ct_ctrl_18wk = np.unique(adata_ctr_18wk.obs["annot"].values, return_counts = True)
unique_ct_ko_12wk, counts_ct_ko_12wk = np.unique(adata_ko_12wk.obs["annot"].values, return_counts = True)
unique_ct_ko_18wk, counts_ct_ko_18wk = np.unique(adata_ko_18wk.obs["annot"].values, return_counts = True)
unique_ct = np.unique(np.concatenate([unique_ct_ctrl_12wk, unique_ct_ctrl_18wk, unique_ct_ko_12wk, unique_ct_ko_18wk], axis = 0))
unique_ct = ["Luminal", "Basal", "Immune", "Luminal (Spink1+)", "Mesenchymal", "Club epithelia", "Other"]

plt.rcParams["font.size"] = 10
fig = plt.figure(figsize = (12,12))
axs = fig.subplots(nrows = 2, ncols = 2)
axs[0,0].pie(counts_ct_ctrl_12wk, labels = unique_ct_ctrl_12wk, autopct='%.0f%%')
axs[0,1].pie(counts_ct_ko_12wk, labels = unique_ct_ko_12wk, autopct='%.0f%%')
axs[1,0].pie(counts_ct_ctrl_18wk, labels = unique_ct_ctrl_18wk, autopct='%.0f%%')
axs[1,1].pie(counts_ct_ko_18wk, labels = unique_ct_ko_18wk, autopct='%.0f%%')
axs[0,0].set_title("CTRL (12wk)", fontsize = 20)
axs[0,1].set_title("CXCR7 KO (12wk)", fontsize = 20)
axs[1,0].set_title("CTRL (18wk)", fontsize = 20)
axs[1,1].set_title("CXCR7 KO (18wk)", fontsize = 20)
fig.suptitle("Cell percentage", fontsize = 30)
fig.tight_layout()

fig.savefig("../results_scrnaseq/precent_celltype.png", bbox_inches = "tight")

# In[]
# stack bar plot
import seaborn as sns
from matplotlib.pylab import cm
import matplotlib.patches as mpatches


cmap = cm.get_cmap("tab10")
plt.rcParams["font.size"] = 15
fig = plt.figure(figsize = (10,7))
ax = fig.subplots(nrows = 1, ncols = 1)

cum_ctrl_12wk = 100
cum_ctrl_18wk = 100
cum_ko_12wk = 100
cum_ko_18wk = 100
legend_ct = []

for idx, ct in enumerate(unique_ct[::-1]):
    # count
    if ct in unique_ct_ctrl_12wk:
        count_ctrl_12wk = counts_ct_ctrl_12wk[np.where(unique_ct_ctrl_12wk == ct)[0][0]]
    else:
        count_ctrl_12wk = 0

    if ct in unique_ct_ctrl_18wk:
        count_ctrl_18wk = counts_ct_ctrl_18wk[np.where(unique_ct_ctrl_18wk == ct)[0][0]]
    else:
        count_ctrl_18wk = 0
    
    if ct in unique_ct_ko_12wk:
        count_ko_12wk = counts_ct_ko_12wk[np.where(unique_ct_ko_12wk == ct)[0][0]]
    else:
        count_ko_12wk = 0
    
    if ct in unique_ct_ko_18wk:
        count_ko_18wk = counts_ct_ko_18wk[np.where(unique_ct_ko_18wk == ct)[0][0]]
    else:
        count_ko_18wk = 0

    # percentage
    pct_ctrl_12wk = count_ctrl_12wk/np.sum(counts_ct_ctrl_12wk) * 100
    pct_ctrl_18wk = count_ctrl_18wk/np.sum(counts_ct_ctrl_18wk) * 100
    pct_ko_12wk = count_ko_12wk/np.sum(counts_ct_ko_12wk) * 100
    pct_ko_18wk = count_ko_18wk/np.sum(counts_ct_ko_18wk) * 100


    data = pd.DataFrame(columns = ["condition", "percentage"])
    data["Condition"] = ["CTRL (12wk)", "CXCR7 KO (12wk)", "CTRL (18wk)", "CXCR7 (18wk)"]
    data["Percentage"] = [cum_ctrl_12wk, cum_ko_12wk, cum_ctrl_18wk, cum_ko_18wk]
    # plot
    sns.barplot(data, x = "Condition", y = "Percentage", ax = ax, color = cmap(idx))

    cum_ctrl_12wk = cum_ctrl_12wk - pct_ctrl_12wk
    cum_ko_12wk = cum_ko_12wk - pct_ko_12wk
    cum_ctrl_18wk = cum_ctrl_18wk - pct_ctrl_18wk
    cum_ko_18wk = cum_ko_18wk - pct_ko_18wk

    legend_ct.append(mpatches.Patch(color = cmap(idx), label = ct))


leg = ax.legend(loc='upper left', prop={'size': 15}, frameon = False, handles = legend_ct, bbox_to_anchor=(1.04, 1))
fig.savefig("../results_scrnaseq/barplot_celltype.png", bbox_inches = "tight")


# In[]
# Save data frame
adata_intact = sc.read_h5ad(seurat_dir + "adata_intact_seurat.h5ad")
adata_intact.obs["annot"] = adata_intact.obs["annot"].astype("str")
adata_intact.obs.loc[adata_intact.obs["annot"].isin(["Endothelial", "SV"]),"annot"] = "Other"
# immune cells count under 12wk
adata_ctr_12wk = adata_intact[(adata_intact.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata_intact.obs["age"] == "12wk"),:]
adata_ko_12wk = adata_intact[(adata_intact.obs["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata_intact.obs["age"] == "12wk"),:]
# immune cells count under 18wk
adata_ctr_18wk = adata_intact[(adata_intact.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata_intact.obs["age"] == "18wk"),:]
adata_ko_18wk = adata_intact[(adata_intact.obs["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata_intact.obs["age"] == "18wk"),:]


counts_df = pd.DataFrame(data = 0, columns = ["CTRL (12wk)", "CXCR7 KO (12wk)", "CTRL (18wk)", "CXCR7 KO (18wk)"], index = ["Luminal", "Basal", "Lymphoid", "Macrophage", "Monocytes", "Luminal (Spink1+)", "Mesenchymal", "Club epithelia", "Other"])
ct_ctr_12wk, count_ctr_12wk = np.unique(adata_ctr_12wk.obs["annot"].values, return_counts = True)
ct_ko_12wk, count_ko_12wk = np.unique(adata_ko_12wk.obs["annot"].values, return_counts = True)
ct_ctr_18wk, count_ctr_18wk = np.unique(adata_ctr_18wk.obs["annot"].values, return_counts = True)
ct_ko_18wk, count_ko_18wk = np.unique(adata_ko_18wk.obs["annot"].values, return_counts = True)

counts_df.loc[ct_ctr_12wk, "CTRL (12wk)"] = count_ctr_12wk
counts_df.loc[ct_ko_12wk, "CXCR7 KO (12wk)"] = count_ko_12wk
counts_df.loc[ct_ctr_18wk, "CTRL (18wk)"] = count_ctr_18wk
counts_df.loc[ct_ko_18wk, "CXCR7 KO (18wk)"] = count_ko_18wk

counts_df.to_csv("../results_scrnaseq/counts_celltype.csv")

pct_df = counts_df.copy()
pct_df.loc[:,:] = counts_df.values/np.sum(counts_df.values, axis = 0, keepdims = True) * 100
pct_df.to_csv("../results_scrnaseq/pct_celltype.csv")

# %%
