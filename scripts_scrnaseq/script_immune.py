# In[]
import sys
sys.path.append("..")
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc
import anndata
import seaborn as sns
import scvi

import utils
import warnings
warnings.filterwarnings("ignore")

import matplotlib

# # NOTE: DE analysis on M1 and M2 macrophages in reference scRNA-seq, only used the generated scrna-seq dataset (ref1), batch effect in ref2 is severe
# adata_ref1 = adata_ref[adata_ref.obs["modality"] == "Ref1"]
# adata_ref2 = adata_ref[adata_ref.obs["modality"] == "Ref2"]
# sc.tl.rank_genes_groups(adata_ref1, groupby = "annot.leiden", method = "wilcoxon", reference = "rest", tie_correct = True, pts = True, rankby_abs = False)
# DE_ct_results = adata_ref1.uns["rank_genes_groups"]
# for ct in adata_ref1.obs["annot.leiden"].cat.categories:
#     if ct != "M1 & M2 macro":
#         # the genes are ordered according to the z-score, then used to calculated the pvalues
#         DE_df = pd.DataFrame(data = 0, columns = ["zscore", "pval", "pval (adj)", "log_fold", "pts", "pts_rest"], index = DE_ct_results["names"][ct])
#         DE_df["zscore"] = DE_ct_results["scores"][ct]
#         DE_df["pval"] = DE_ct_results["pvals"][ct]
#         DE_df["pval (adj)"] = DE_ct_results["pvals_adj"][ct]
#         DE_df["log_fold"] = DE_ct_results["logfoldchanges"][ct]
#         DE_df["pts"] = DE_ct_results["pts"][ct]
#         DE_df["pts_rest"] = DE_ct_results["pts_rest"][ct]

#         # need to separate according to the enriched and depleted genes
#         # the enriched and depleted genes are separated according to the mean log-fold change
#         DE_df_enrich = DE_df[(DE_df["log_fold"] >= 0) & (DE_df["pts"] >= 0.3)]
#         DE_df_deplete = DE_df[(DE_df["log_fold"] < 0) & (DE_df["pts_rest"] >= 0.3)]
#         DE_df_enrich.to_csv(res_dir + f"plots_marker/DE_{ct}_enrich.csv", sep = ",")
#         # DE_df_deplete.to_csv(res_dir + f"plots_marker/DE_{ct}_deplete.csv", sep = ",")

#         # filtering based on the percentage of expression
#         # DE enriched should have more than 20% of cells expressed in the cluster
#         DE_df = pd.concat([DE_df_enrich, DE_df_deplete], axis = 0)
#         fig, ax = utils.plot_volcano(DE_df, x = "log_fold", y = "pval (adj)", x_cutoffs = [-np.inf, 1], y_cutoff = 0.01, gene_name = True, xlim = 3, ylim = 1e-5)
#         ax.set_title(f"DE: {ct}", fontsize = 25)
#         fig.savefig(res_dir + f"plots_marker/DE_volcano_{ct}.png", bbox_inches = "tight", dpi = 150)
    

# In[]
# -----------------------------------------------------------------------------------------
#
# sub type annotation of macrophage
#
# -----------------------------------------------------------------------------------------
# NOTE: no refined annotations in the reference scRNA-seq dataset
# load reference scRNA-seq dataset, NOTE: use raw data, include 4 samples/batches
adata_ref1 =  sc.read_h5ad("../dataset/data_scrnaseq/data/qc_data/adata_qc_intact.h5ad")
# The smallest cell cluster in reference: Lymphoid, include 349 cells
# lt_correct is the corrected version of the label 
adata_ref1.obs["annot"] = sc.read_h5ad("../dataset/data_scrnaseq/seurat_integration/adata_intact_seurat_lt_correct.h5ad").obs.loc[adata_ref1.obs.index,"annot (lt_correct)"].values
# select only the macrophage
adata_ref1_macropphage = adata_ref1[adata_ref1.obs["annot"] == "Macrophage", :]
adata_ref1_macropphage.layers["counts"] = adata_ref1_macropphage.X.astype(int)

adata_ref2 = sc.read_h5ad("../dataset/reference_scrnaseq/reference_raw.h5ad")
# update the batch and annotation labels for LT
adata_ref2.obs["sample"] = adata_ref2.obs["batchID"].values
# select only macrophage
adata_ref2_macrophage = adata_ref2[adata_ref2.obs["FullTypeMerged"] == "Imm_Macrophage", :]
adata_ref2_macrophage.layers["counts"] = adata_ref2_macrophage.X.astype(int)

# combine the reference annotation
adata_ref_macrophage = anndata.concat([adata_ref1_macropphage, adata_ref2_macrophage], axis = 0, join = "inner", label = "modality", keys = ["Ref1", "Ref2"])

marker_m2 = ["Cd163", "Mrc1", "Ccl24", "Arg1", "Msr1", "Il10"]
marker_m1 = ["Il1b", "Il6", "Tnf", "Nos2", "Cxcl9"]

scvi.settings.seed = 0
scvi.model.SCVI.setup_anndata(adata_ref_macrophage, layer = "counts", batch_key = "sample")
scvi_model = scvi.model.SCVI(adata_ref_macrophage, use_layer_norm = "both", use_batch_norm = "none", encode_covariates = True, dropout_rate = 0.2, n_layers = 2, gene_likelihood = "nb")
scvi_model.train()
SCVI_LATENT_KEY = "X_scVI"
adata_ref_macrophage.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation()
sc.pp.neighbors(adata_ref_macrophage, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata_ref_macrophage, resolution = 0.2)
sc.tl.umap(adata_ref_macrophage)

sc.pp.normalize_total(adata_ref_macrophage, 10e4)
sc.pp.log1p(adata_ref_macrophage)
gene_expr = pd.DataFrame(data = adata_ref_macrophage[:, marker_m1 + marker_m2].X.toarray(), columns = marker_m1 + marker_m2, index = adata_ref_macrophage.obs.index.values.squeeze())
fig = utils.plot_embeds_continuous(adata_ref_macrophage.obsm["X_umap"], annos = gene_expr, colormap = utils.SUPER_MAGMA, figsize = (7,5), ncols = 4, alpha = 0.4, s = 2)
fig.savefig("../results_scrnaseq/annot_immune_subtype/macrophage_markers_Ref.png", bbox_inches = "tight")

# annotate macrophage by m1 & m2
adata_ref_macrophage.obs["annot.leiden"] = "M1 macro"
adata_ref_macrophage.obs.loc[adata_ref_macrophage.obs["leiden"].isin(["1"]), "annot.leiden"] = "M2 macro"
adata_ref_macrophage.obs["annot.leiden"] = adata_ref_macrophage.obs["annot.leiden"].astype("category")
# fig = utils.plot_embeds(adata_ref_macrophage[adata_ref_macrophage.obs["modality"] == "Ref1"].obsm["X_umap"], annos = adata_ref_macrophage[adata_ref_macrophage.obs["modality"] == "Ref1"].obs[["annot.leiden", "modality"]], figsize = (9,5), markerscale = 5)
fig = utils.plot_embeds(adata_ref_macrophage.obsm["X_umap"], annos = adata_ref_macrophage.obs[["leiden", "annot.leiden", "modality"]], figsize = (9,5), markerscale = 5)
fig.savefig("../results_scrnaseq/annot_immune_subtype/macrophage_leiden.png", bbox_inches = "tight")

# In[]
# -----------------------------------------------------------------------------------------
#
# sub type annotation of lymphoid
#
# -----------------------------------------------------------------------------------------
adata_ref1 =  sc.read_h5ad("../dataset/data_scrnaseq/data/qc_data/adata_qc_intact.h5ad")
adata_ref1.obs["annot"] = sc.read_h5ad("../dataset/data_scrnaseq/seurat_integration/adata_intact_seurat_lt_correct.h5ad").obs.loc[adata_ref1.obs.index,"annot (lt_correct)"].values
# select only lymphoid cells
adata_ref1_lymphoid = adata_ref1[adata_ref1.obs["annot"] == "Lymphoid", :]
adata_ref1_lymphoid.layers["counts"] = adata_ref1_lymphoid.X.copy()

adata_ref2 = sc.read_h5ad("../dataset/reference_scrnaseq/reference_raw.h5ad")
# update the batch and annotation labels for LT
adata_ref2.obs["sample"] = adata_ref2.obs["batchID"].values
adata_ref2.obs["modality"] = "Ref"
# select only lymphoid cells
adata_ref2_lymphoid = adata_ref2[adata_ref2.obs["FullTypeMerged"].isin(["Imm_B", "Imm_NK", "Imm_Tcell"])]
adata_ref2_lymphoid.layers["counts"] = adata_ref2_lymphoid.X.copy()
# combine reference
adata_ref_lymphoid = anndata.concat([adata_ref1_lymphoid, adata_ref2_lymphoid], axis = 0, join = "inner", label = "modality", keys = ["Ref1", "Ref2"])


marker = ["Cd3e", "Cd4", "Cd8a", "Klrb1c", "Klrd1", "Nkg7", "Ms4a1"]

scvi.settings.seed = 0
scvi.model.SCVI.setup_anndata(adata_ref_lymphoid, layer = "counts", batch_key = "sample")
scvi_model = scvi.model.SCVI(adata_ref_lymphoid, use_layer_norm = "both", use_batch_norm = "none", encode_covariates = True, dropout_rate = 0.2, n_layers = 2, gene_likelihood = "nb")
scvi_model.train()
SCVI_LATENT_KEY = "X_scVI"
adata_ref_lymphoid.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation()
sc.pp.neighbors(adata_ref_lymphoid, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata_ref_lymphoid, resolution = 1)
sc.tl.umap(adata_ref_lymphoid)

sc.pp.normalize_total(adata_ref_lymphoid, 10e4)
sc.pp.log1p(adata_ref_lymphoid)
gene_expr = pd.DataFrame(data = adata_ref_lymphoid[:, marker].X.toarray(), columns = marker, index = adata_ref_lymphoid.obs.index.values.squeeze())
fig = utils.plot_embeds_continuous(adata_ref_lymphoid.obsm["X_umap"], annos = gene_expr, colormap = utils.SUPER_MAGMA, figsize = (7,5), ncols = 4, alpha = 0.4, s = 2)
fig.savefig("../results_scrnaseq/annot_immune_subtype/lymphoid_markers_Ref.png", bbox_inches = "tight")

# annotate macrophage by m1 & m2
adata_ref_lymphoid.obs["annot.leiden"] = "Lymphoid"
adata_ref_lymphoid.obs.loc[adata_ref_lymphoid.obs["leiden"].isin(["0"]), "annot.leiden"] = "NK"
# 1 is a subset of CD8 T population that also express NK markers, including NTG7 AND KLRD1, differentiate from the remaining CD8+ T cell population
adata_ref_lymphoid.obs.loc[adata_ref_lymphoid.obs["leiden"].isin(["1"]), "annot.leiden"] = "CD8 T"
adata_ref_lymphoid.obs.loc[adata_ref_lymphoid.obs["leiden"].isin(["5", "9", "13"]), "annot.leiden"] = "CD8 T"
adata_ref_lymphoid.obs.loc[adata_ref_lymphoid.obs["leiden"].isin(["3", "7", "11", "10"]), "annot.leiden"] = "CD4 T"
adata_ref_lymphoid.obs.loc[adata_ref_lymphoid.obs["leiden"].isin(["4", "6", "2", "12"]), "annot.leiden"] = "DN T"
adata_ref_lymphoid.obs.loc[adata_ref_lymphoid.obs["leiden"].isin(["8"]), "annot.leiden"] = "B"

adata_ref_lymphoid.obs["annot.leiden"] = adata_ref_lymphoid.obs["annot.leiden"].astype("category")
fig = utils.plot_embeds(adata_ref_lymphoid.obsm["X_umap"], annos = adata_ref_lymphoid.obs[["annot.leiden", "modality"]], figsize = (9,5), markerscale = 5)
fig.savefig("../results_scrnaseq/annot_immune_subtype/lymphoid_leiden.png", bbox_inches = "tight")



# In[]
# refine the annotation of ref2 (scrnaseq)
adata_ref = sc.read_h5ad("../dataset/data_scrnaseq/seurat_integration/adata_intact_seurat_lt_correct.h5ad")
adata_ref.obs["annot.level2"] = adata_ref.obs["annot (lt_correct)"].values.astype(object)
adata_ref.obs.loc[adata_ref_lymphoid[adata_ref_lymphoid.obs["modality"] == "Ref1"].obs.index, "annot.level2"] = adata_ref_lymphoid[adata_ref_lymphoid.obs["modality"] == "Ref1"].obs["annot.leiden"].values
adata_ref.obs.loc[adata_ref_macrophage[adata_ref_macrophage.obs["modality"] == "Ref1"].obs.index, "annot.level2"] = adata_ref_macrophage[adata_ref_macrophage.obs["modality"] == "Ref1"].obs["annot.leiden"].values
adata_ref.write_h5ad("../dataset/data_scrnaseq/seurat_integration/adata_intact_seurat_lt_correct.h5ad")


adata_ref = sc.read_h5ad("../dataset/reference_scrnaseq/reference_raw.h5ad")
adata_ref.obs["annot"] = adata_ref.obs["FullTypeMerged"].values.astype(object)
# remove the predicted doublets
adata_ref = adata_ref[~adata_ref.obs["annot"].isin(["PredDoublet_Epi_Imm", "PredDoublet_Str_Epi", "PredDoublet_Str_Imm"]),:]
# combine the SV cells, not of interest, better to remove SV cells when calculating the cell type proportion
adata_ref.obs.loc[adata_ref.obs["annot"].isin(["Epi_SV_Basal", "Epi_SV_Ionocyte", "Epi_SV_Luminal"]), "annot"] = "SV"
# combine endothelium
adata_ref.obs.loc[adata_ref.obs["annot"].isin(["Str_Endothelium_Lymphatic", "Str_Endothelium_Vascular"]), "annot"] = "Endothelial"
# combine luminal
adata_ref.obs.loc[adata_ref.obs["annot"].isin(["Epi_Luminal", "Epi_Luminal_2Psca", "Epi_Luminal_3Foxi1"]), "annot"] = "Luminal"
# combine lymphoid, NOTE: the lymphoid subtype annotation is not very accurate according to marker gene checking
adata_ref.obs.loc[adata_ref.obs["annot"].isin(["Imm_B", "Imm_NK", "Imm_Tcell"]), "annot"] = "Lymphoid"
# change name
adata_ref.obs.loc[adata_ref.obs["annot"].isin(["Str_Mesenchymal"]), "annot"] = "Mesenchymal"
adata_ref.obs.loc[adata_ref.obs["annot"].isin(["Epi_Basal"]), "annot"] = "Basal"
adata_ref.obs.loc[adata_ref.obs["annot"].isin(["Imm_Macrophage"]), "annot"] = "Macrophage"
adata_ref.obs.loc[adata_ref.obs["annot"].isin(["Imm_DC"]), "annot"] = "DC"
adata_ref.obs["annot"] = adata_ref.obs["annot"].astype("category")
adata_ref.obs["annot.level2"] = adata_ref.obs["annot"].astype(object)
adata_ref.obs.loc[adata_ref_lymphoid[adata_ref_lymphoid.obs["modality"] == "Ref2"].obs.index, "annot.level2"] = adata_ref_lymphoid[adata_ref_lymphoid.obs["modality"] == "Ref2"].obs["annot.leiden"].values
adata_ref.obs.loc[adata_ref_macrophage[adata_ref_macrophage.obs["modality"] == "Ref2"].obs.index, "annot.level2"] = adata_ref_macrophage[adata_ref_macrophage.obs["modality"] == "Ref2"].obs["annot.leiden"].values
adata_ref.write_h5ad("../dataset/reference_scrnaseq/reference_correct.h5ad")

# In[]
annot_sub = sc.read_h5ad("../dataset/data_scrnaseq/seurat_integration/adata_intact_seurat_lt_correct.h5ad").obs[["annot (lt_correct)","annot.level2", "genotype", "age", "sample"]]
annot_lymphoid = annot_sub[annot_sub["annot (lt_correct)"] == "Lymphoid"]
annot_macrophage = annot_sub[annot_sub["annot (lt_correct)"] == "Macrophage"]


count_lymphoid_sub = pd.DataFrame(columns = ["CellType", "Genotype", "Counts", "Sample", "Age", "Pct"])
lymphoid_subtypes = ["CD4 T", "CD8 T", "DN T", "NK", "B"]
ncell_pp12 = np.sum(annot_sub["sample"] == "M1417-PP12")
ncell_ppc12 = np.sum(annot_sub["sample"] == "M1436-PPC12")
ncell_pp18 = np.sum(annot_sub["sample"] == "M1416-PP18")
ncell_ppc18 = np.sum(annot_sub["sample"] == "M1437-PPC18")

for lymphoid_subtype in lymphoid_subtypes:
    counts_pp12 = np.sum(annot_lymphoid.loc[annot_lymphoid["sample"] == "M1417-PP12", "annot.level2"].values == lymphoid_subtype)
    counts_ppc12 = np.sum(annot_lymphoid.loc[annot_lymphoid["sample"] == "M1436-PPC12", "annot.level2"].values == lymphoid_subtype)
    counts_pp18 = np.sum(annot_lymphoid.loc[annot_lymphoid["sample"] == "M1416-PP18", "annot.level2"].values == lymphoid_subtype)
    counts_ppc18 = np.sum(annot_lymphoid.loc[annot_lymphoid["sample"] == "M1437-PPC18", "annot.level2"].values == lymphoid_subtype)
    count_lymphoid_sub.loc[len(count_lymphoid_sub), ["CellType", "Genotype", "Counts", "Sample", "Age", "Pct"]] = lymphoid_subtype, "PbCre(+/-),Pten(-/-),P53(-/-)", counts_pp12, "PP12", "12wk", counts_pp12/ncell_pp12
    count_lymphoid_sub.loc[len(count_lymphoid_sub), ["CellType", "Genotype", "Counts", "Sample", "Age", "Pct"]] = lymphoid_subtype, "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", counts_ppc12, "PPC12", "12wk", counts_ppc12/ncell_ppc12
    count_lymphoid_sub.loc[len(count_lymphoid_sub), ["CellType", "Genotype", "Counts", "Sample", "Age", "Pct"]] = lymphoid_subtype, "PbCre(+/-),Pten(-/-),P53(-/-)", counts_pp18, "PP18", "18wk", counts_pp18/ncell_pp18
    count_lymphoid_sub.loc[len(count_lymphoid_sub), ["CellType", "Genotype", "Counts", "Sample", "Age", "Pct"]] = lymphoid_subtype, "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", counts_ppc18, "PPC18", "18wk", counts_ppc18/ncell_ppc18


count_macrophage_sub = pd.DataFrame(columns = ["CellType", "Genotype", "Counts", "Sample", "Age", "Pct"])
macrophage_subtypes = ["M1 macro", "M2 macro"]
for macrophage_subtype in macrophage_subtypes:
    counts_pp12 = np.sum(annot_macrophage.loc[annot_macrophage["sample"] == "M1417-PP12", "annot.level2"].values == macrophage_subtype)
    counts_ppc12 = np.sum(annot_macrophage.loc[annot_macrophage["sample"] == "M1436-PPC12", "annot.level2"].values == macrophage_subtype)
    counts_pp18 = np.sum(annot_macrophage.loc[annot_macrophage["sample"] == "M1416-PP18", "annot.level2"].values == macrophage_subtype)
    counts_ppc18 = np.sum(annot_macrophage.loc[annot_macrophage["sample"] == "M1437-PPC18", "annot.level2"].values == macrophage_subtype)
    count_macrophage_sub.loc[len(count_macrophage_sub), ["CellType", "Genotype", "Counts", "Sample", "Age", "Pct"]] = macrophage_subtype, "PbCre(+/-),Pten(-/-),P53(-/-)", counts_pp12, "PP12", "12wk", counts_pp12/ncell_pp12
    count_macrophage_sub.loc[len(count_macrophage_sub), ["CellType", "Genotype", "Counts", "Sample", "Age", "Pct"]] = macrophage_subtype, "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", counts_ppc12, "PPC12", "12wk", counts_ppc12/ncell_ppc12
    count_macrophage_sub.loc[len(count_macrophage_sub), ["CellType", "Genotype", "Counts", "Sample", "Age", "Pct"]] = macrophage_subtype, "PbCre(+/-),Pten(-/-),P53(-/-)", counts_pp18, "PP18", "18wk", counts_pp18/ncell_pp18
    count_macrophage_sub.loc[len(count_macrophage_sub), ["CellType", "Genotype", "Counts", "Sample", "Age", "Pct"]] = macrophage_subtype, "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", counts_ppc18, "PPC18", "18wk", counts_ppc18/ncell_ppc18

sns.set_theme()
fig = plt.figure(figsize = (25, 14))
ax = fig.subplots(nrows = 2, ncols = 2)
ax[0,0] = sns.barplot(data = count_lymphoid_sub, x = "CellType", y = "Counts", hue = "Genotype", ax = ax[0,0], width=0.7, estimator=sum, ci=None)
sns.move_legend(ax[0,0], "upper left", bbox_to_anchor=(1.04, 1), frameon = False, fontsize = 15, title = None)
ax[0,1] = sns.barplot(data = count_macrophage_sub, x = "CellType", y = "Counts", hue = "Genotype", ax = ax[0,1], width=0.7, estimator=sum, ci=None)
sns.move_legend(ax[0,1], "upper left", bbox_to_anchor=(1.04, 1), frameon = False, fontsize = 15, title = None)
ax[0,0].set_title("Lymphoid Subtype", fontsize = 20)
ax[0,0].set_ylabel("#cells", fontsize = 20)
ax[0,1].set_title("Macrophage Subtype", fontsize = 20)
ax[0,1].set_ylabel("#cells", fontsize = 20)
for i in ax[0,0].containers:
    ax[0,0].bar_label(i, fmt = "%d", fontsize = 15)
for i in ax[0,1].containers:
    ax[0,1].bar_label(i, fmt = "%d", fontsize = 15)

ax[1,0] = sns.barplot(data = count_lymphoid_sub, x = "CellType", y = "Counts", hue = "Sample", ax = ax[1,0], width=0.7)
sns.move_legend(ax[1,0], "upper left", bbox_to_anchor=(1.04, 1), frameon = False, fontsize = 15, title = None)
ax[1,1] = sns.barplot(data = count_macrophage_sub, x = "CellType", y = "Counts", hue = "Sample", ax = ax[1,1], width=0.7)
sns.move_legend(ax[1,1], "upper left", bbox_to_anchor=(1.04, 1), frameon = False, fontsize = 15, title = None)
ax[1,0].set_title("Lymphoid Subtype", fontsize = 20)
ax[1,0].set_ylabel("#cells", fontsize = 20)
ax[1,1].set_title("Macrophage Subtype", fontsize = 20)
ax[1,1].set_ylabel("#cells", fontsize = 20)
for i in ax[1,0].containers:
    ax[1,0].bar_label(i, fmt = "%d", fontsize = 15)
for i in ax[1,1].containers:
    ax[1,1].bar_label(i, fmt = "%d", fontsize = 15)
plt.tight_layout()
fig.savefig("../results_scrnaseq/annot_immune_subtype/immune_subtype_counts.png", bbox_inches = "tight", dpi = 150)

fig = plt.figure(figsize = (30, 7))
ax = fig.subplots(nrows = 1, ncols = 2)
ax[0] = sns.barplot(data = count_lymphoid_sub, x = "CellType", y = "Pct", hue = "Sample", ax = ax[0], width=0.7)
sns.move_legend(ax[0], "upper left", bbox_to_anchor=(1.04, 1), frameon = False, fontsize = 15, title = None)
ax[1] = sns.barplot(data = count_macrophage_sub, x = "CellType", y = "Pct", hue = "Sample", ax = ax[1], width=0.7)
sns.move_legend(ax[1], "upper left", bbox_to_anchor=(1.04, 1), frameon = False, fontsize = 15, title = None)
ax[0].set_title("Lymphoid Subtype", fontsize = 20)
ax[0].set_ylabel("Pct/all cells", fontsize = 20)
ax[1].set_title("Macrophage Subtype", fontsize = 20)
ax[1].set_ylabel("Pct/all cells", fontsize = 20)
for i in ax[0].containers:
    ax[0].bar_label(i, fmt = "%.3f", fontsize = 15)
for i in ax[1].containers:
    ax[1].bar_label(i, fmt = "%.3f", fontsize = 15)
plt.tight_layout()
fig.savefig("../results_scrnaseq/annot_immune_subtype/immune_subtype_pct.png", bbox_inches = "tight", dpi = 150)

matplotlib.rc_file_defaults()

# In[]
# NOTE: Lymphoid
fig = plt.figure(figsize = (14,22))
axs = fig.subplots(nrows = 3, ncols = 2)
axs[0,0].pie([x for x in count_lymphoid_sub.loc[count_lymphoid_sub["Sample"] == "PP12", "Counts"]], labels = [x for x in count_lymphoid_sub.loc[count_lymphoid_sub["Sample"] == "PP12", "CellType"]], autopct='%.2f%%')
axs[0,1].pie([x for x in count_lymphoid_sub.loc[count_lymphoid_sub["Sample"] == "PPC12", "Counts"]], labels = [x for x in count_lymphoid_sub.loc[count_lymphoid_sub["Sample"] == "PPC12", "CellType"]], autopct='%.2f%%')
axs[1,0].pie([x for x in count_lymphoid_sub.loc[count_lymphoid_sub["Sample"] == "PP18", "Counts"]], labels = [x for x in count_lymphoid_sub.loc[count_lymphoid_sub["Sample"] == "PP18", "CellType"]], autopct='%.2f%%')
axs[1,1].pie([x for x in count_lymphoid_sub.loc[count_lymphoid_sub["Sample"] == "PPC18", "Counts"]], labels = [x for x in count_lymphoid_sub.loc[count_lymphoid_sub["Sample"] == "PPC18", "CellType"]], autopct='%.2f%%')
axs[2,0].pie([x for x in count_lymphoid_sub.loc[count_lymphoid_sub["Sample"] == "PP12", "Counts"].values + count_lymphoid_sub.loc[count_lymphoid_sub["Sample"] == "PPC12", "Counts"].values],
             labels = [x for x in count_lymphoid_sub.loc[count_lymphoid_sub["Sample"] == "PP12", "CellType"]], autopct='%.2f%%')
axs[2,1].pie([x for x in count_lymphoid_sub.loc[count_lymphoid_sub["Sample"] == "PP18", "Counts"].values + count_lymphoid_sub.loc[count_lymphoid_sub["Sample"] == "PPC18", "Counts"].values],
             labels = [x for x in count_lymphoid_sub.loc[count_lymphoid_sub["Sample"] == "PPC12", "CellType"]], autopct='%.2f%%')

axs[0,0].set_title("PP (12wk)", fontsize = 20)
axs[0,1].set_title("PPC (12wk)", fontsize = 20)
axs[1,0].set_title("PP (18wk)", fontsize = 20)
axs[1,1].set_title("PPC (18wk)", fontsize = 20)
axs[2,0].set_title("PP (Total)", fontsize = 20)
axs[2,1].set_title("PPC (Total)", fontsize = 20)
fig.suptitle(f"Lymphoid subtype / all lymphoid", fontsize = 25)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
fig.savefig("../results_scrnaseq/annot_immune_subtype/lymphoid_subtype_pcts.png", bbox_inches = "tight", dpi = 150)

# NOTE: Macrophage
fig = plt.figure(figsize = (14,22))
axs = fig.subplots(nrows = 3, ncols = 2)
axs[0,0].pie([x for x in count_macrophage_sub.loc[count_macrophage_sub["Sample"] == "PP12", "Counts"]], labels = [x for x in count_macrophage_sub.loc[count_macrophage_sub["Sample"] == "PP12", "CellType"]], autopct='%.2f%%')
axs[0,1].pie([x for x in count_macrophage_sub.loc[count_macrophage_sub["Sample"] == "PPC12", "Counts"]], labels = [x for x in count_macrophage_sub.loc[count_macrophage_sub["Sample"] == "PPC12", "CellType"]], autopct='%.2f%%')
axs[1,0].pie([x for x in count_macrophage_sub.loc[count_macrophage_sub["Sample"] == "PP18", "Counts"]], labels = [x for x in count_macrophage_sub.loc[count_macrophage_sub["Sample"] == "PP18", "CellType"]], autopct='%.2f%%')
axs[1,1].pie([x for x in count_macrophage_sub.loc[count_macrophage_sub["Sample"] == "PPC18", "Counts"]], labels = [x for x in count_macrophage_sub.loc[count_macrophage_sub["Sample"] == "PPC18", "CellType"]], autopct='%.2f%%')
axs[2,0].pie([x for x in count_macrophage_sub.loc[count_macrophage_sub["Sample"] == "PP12", "Counts"].values + count_macrophage_sub.loc[count_macrophage_sub["Sample"] == "PPC12", "Counts"].values],
             labels = [x for x in count_macrophage_sub.loc[count_macrophage_sub["Sample"] == "PP12", "CellType"]], autopct='%.2f%%')
axs[2,1].pie([x for x in count_macrophage_sub.loc[count_macrophage_sub["Sample"] == "PP18", "Counts"].values + count_macrophage_sub.loc[count_macrophage_sub["Sample"] == "PPC18", "Counts"].values],
             labels = [x for x in count_macrophage_sub.loc[count_macrophage_sub["Sample"] == "PPC12", "CellType"]], autopct='%.2f%%')
axs[0,0].set_title("PP (12wk)", fontsize = 20)
axs[0,1].set_title("PPC (12wk)", fontsize = 20)
axs[1,0].set_title("PP (18wk)", fontsize = 20)
axs[1,1].set_title("PPC (18wk)", fontsize = 20)
axs[2,0].set_title("PP (Total)", fontsize = 20)
axs[2,1].set_title("PPC (Total)", fontsize = 20)

fig.suptitle(f"Macrophage subtype / all macrophage", fontsize = 25)
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
fig.savefig("../results_scrnaseq/annot_immune_subtype/macrophage_subtype_pcts.png", bbox_inches = "tight", dpi = 150)


# In[]
# # ---------------------------------------------------------------- #
# #
# # 4. NOTE: subtype annotation (method 2)
# #  * count the number of cells expression the T cell subtype markers (percentage) under two different genotypes
# #  * boxplot showing the nonzero marker gene expression (remove zeros) level under two different genotypes
# #
# # ---------------------------------------------------------------- #

# # 1. Lymphoid cell subtype comparison
# # extract expression
# adata_lymphoid_total = adata_lymphoid.copy()
# lymphoid_sub_df = pd.DataFrame(columns = ["Subtype", "GT", "Percentage", "Mean expr (nonzero)"])
# for age in ["12wk", "18wk"]:
#     adata_lymphoid = adata_lymphoid_total[adata_lymphoid_total.obs["age"] == age]
#     x_pp_cd4T = adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)", "Cd4"].X.toarray().squeeze()
#     x_ppc_cd4T = adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", "Cd4"].X.toarray().squeeze()
#     x_pp_cd8T = adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)", "Cd8a"].X.toarray().squeeze()
#     x_ppc_cd8T = adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", "Cd8a"].X.toarray().squeeze()
#     x_pp_cytoT = adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)", ["Gzmb", "Prf1"]].X.toarray().squeeze()
#     x_ppc_cytoT = adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ["Gzmb", "Prf1"]].X.toarray().squeeze()
#     x_pp_exhaustT = adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)", ["Tigit", "Havcr2"]].X.toarray().squeeze()
#     x_ppc_exhaustT = adata_lymphoid[adata_lymphoid.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ["Tigit", "Havcr2"]].X.toarray().squeeze()

#     # calculate percentage
#     pct_pp_cd4T = np.sum((x_pp_cd4T > 0))/x_pp_cd4T.shape[0]
#     pct_ppc_cd4T = np.sum((x_ppc_cd4T > 0))/x_ppc_cd4T.shape[0]
#     pct_pp_cd8T = np.sum((x_pp_cd8T > 0))/x_pp_cd8T.shape[0]
#     pct_ppc_cd8T = np.sum((x_ppc_cd8T > 0))/x_ppc_cd8T.shape[0]
#     pct_pp_cytoT = np.sum((x_pp_cytoT[:,0] > 0) & (x_pp_cytoT[:,1] > 0))/x_pp_cytoT.shape[0]
#     pct_ppc_cytoT = np.sum((x_ppc_cytoT[:,0] > 0) & (x_ppc_cytoT[:,1] > 0))/x_ppc_cytoT.shape[0]
#     pct_pp_exhaustT = np.sum((x_pp_exhaustT[:,0] > 0) & (x_pp_exhaustT[:,1] > 0))/x_pp_exhaustT.shape[0]
#     pct_ppc_exhaustT = np.sum((x_ppc_exhaustT[:,0] > 0) & (x_ppc_exhaustT[:,1] > 0))/x_ppc_exhaustT.shape[0]

#     expr_pp_cd4T = np.mean(x_pp_cd4T[x_pp_cd4T > 0])
#     expr_ppc_cd4T = np.mean(x_ppc_cd4T[x_ppc_cd4T > 0])
#     expr_pp_cd8T = np.mean(x_pp_cd8T[x_pp_cd8T > 0])
#     expr_ppc_cd8T = np.mean(x_ppc_cd8T[x_ppc_cd8T > 0])
#     expr_pp_cytoT = np.mean(x_pp_cytoT[(x_pp_cytoT[:,0] > 0) & (x_pp_cytoT[:,1] > 0)])
#     expr_ppc_cytoT = np.mean(x_ppc_cytoT[(x_ppc_cytoT[:,0] > 0) & (x_ppc_cytoT[:,1] > 0)])
#     expr_pp_exhaustT = np.mean(x_pp_exhaustT[(x_pp_exhaustT[:,0] > 0) & (x_pp_exhaustT[:,1] > 0)])
#     expr_ppc_exhaustT = np.mean(x_ppc_exhaustT[(x_ppc_exhaustT[:,0] > 0) & (x_ppc_exhaustT[:,1] > 0)])

#     lymphoid_sub_df = lymphoid_sub_df._append({"Subtype": "CD4-T", "GT": "PP " + age, "Percentage": pct_pp_cd4T, "Mean expr (nonzero)": expr_pp_cd4T}, ignore_index = True)
#     lymphoid_sub_df = lymphoid_sub_df._append({"Subtype": "CD4-T", "GT": "PPC " + age, "Percentage": pct_ppc_cd4T, "Mean expr (nonzero)": expr_ppc_cd4T}, ignore_index = True)
#     lymphoid_sub_df = lymphoid_sub_df._append({"Subtype": "CD8-T", "GT": "PP " + age, "Percentage": pct_pp_cd8T, "Mean expr (nonzero)": expr_pp_cd8T}, ignore_index = True)
#     lymphoid_sub_df = lymphoid_sub_df._append({"Subtype": "CD8-T", "GT": "PPC " + age, "Percentage": pct_ppc_cd8T, "Mean expr (nonzero)": expr_ppc_cd8T}, ignore_index = True)
#     lymphoid_sub_df = lymphoid_sub_df._append({"Subtype": "CytoT", "GT": "PP " + age, "Percentage": pct_pp_cytoT, "Mean expr (nonzero)": expr_pp_cytoT}, ignore_index = True)
#     lymphoid_sub_df = lymphoid_sub_df._append({"Subtype": "CytoT", "GT": "PPC " + age, "Percentage": pct_ppc_cytoT, "Mean expr (nonzero)": expr_ppc_cytoT}, ignore_index = True)
#     lymphoid_sub_df = lymphoid_sub_df._append({"Subtype": "ExhaustT", "GT": "PP " + age, "Percentage": pct_pp_exhaustT, "Mean expr (nonzero)": expr_pp_exhaustT}, ignore_index = True)
#     lymphoid_sub_df = lymphoid_sub_df._append({"Subtype": "ExhaustT", "GT": "PPC " + age, "Percentage": pct_ppc_exhaustT, "Mean expr (nonzero)": expr_ppc_exhaustT}, ignore_index = True)

# # 2. Macrophage subtypes
# adata_macrophage_total = adata_macrophage.copy()
# macrophage_sub_df = pd.DataFrame(columns = ["Subtype", "GT", "Percentage", "Mean expr (nonzero)"])
# for age in ["12wk", "18wk"]:
#     adata_macrophage = adata_macrophage_total[adata_macrophage_total.obs["age"] == age]
#     x_pp_m1 = adata_macrophage[adata_macrophage.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)", ["Cd80", "Cd86"]].X.toarray().squeeze() # NOTE: expression of Nos2 is too few
#     x_ppc_m1 = adata_macrophage[adata_macrophage.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ["Cd80", "Cd86"]].X.toarray().squeeze()
#     x_pp_m2 = adata_macrophage[adata_macrophage.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)", ["Arg1", "Mrc1", "Cd163"]].X.toarray().squeeze()
#     x_ppc_m2 = adata_macrophage[adata_macrophage.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ["Arg1", "Mrc1", "Cd163"]].X.toarray().squeeze()

#     pct_pp_m1 = np.sum((x_pp_m1[:,0] > 0) & (x_pp_m1[:,1] > 0))/x_pp_m1.shape[0]
#     pct_ppc_m1 = np.sum((x_ppc_m1[:,0] > 0) & (x_ppc_m1[:,1] > 0))/x_ppc_m1.shape[0]
#     pct_pp_m2 = np.sum((x_pp_m2[:,0] > 0) & (x_pp_m2[:,1] > 0) & (x_pp_m2[:,2] > 0))/x_pp_m2.shape[0]
#     pct_ppc_m2 = np.sum((x_ppc_m2[:,0] > 0) & (x_ppc_m2[:,1] > 0) & (x_ppc_m2[:,2] > 0))/x_ppc_m2.shape[0]

#     expr_pp_m1 = np.mean(x_pp_m1[(x_pp_m1[:,0] > 0) & (x_pp_m1[:,1] > 0)])
#     expr_ppc_m1 = np.mean(x_ppc_m1[(x_ppc_m1[:,0] > 0) & (x_ppc_m1[:,1] > 0)])
#     expr_pp_m2 = np.mean(x_pp_m2[(x_pp_m2[:,0] > 0) & (x_pp_m2[:,1] > 0) & (x_pp_m2[:,2] > 0)])
#     expr_ppc_m2 = np.mean(x_ppc_m2[(x_ppc_m2[:,0] > 0) & (x_ppc_m2[:,1] > 0) & (x_ppc_m2[:,2] > 0)])

#     macrophage_sub_df = macrophage_sub_df._append({"Subtype": "M1 Macro", "GT": "PP " + age, "Percentage": pct_pp_m1, "Mean expr (nonzero)": expr_pp_m1}, ignore_index = True)
#     macrophage_sub_df = macrophage_sub_df._append({"Subtype": "M1 Macro", "GT": "PPC " + age, "Percentage": pct_ppc_m1, "Mean expr (nonzero)": expr_ppc_m1}, ignore_index = True)
#     macrophage_sub_df = macrophage_sub_df._append({"Subtype": "M2 Macro", "GT": "PP " + age, "Percentage": pct_pp_m2, "Mean expr (nonzero)": expr_pp_m2}, ignore_index = True)
#     macrophage_sub_df = macrophage_sub_df._append({"Subtype": "M2 Macro", "GT": "PPC " + age, "Percentage": pct_ppc_m2, "Mean expr (nonzero)": expr_ppc_m2}, ignore_index = True)

# # In[]
# sns.set_theme()
# fig = plt.figure(figsize = (25, 14))
# ax = fig.subplots(nrows = 2, ncols = 2)
# ax[0,0] = sns.barplot(data = lymphoid_sub_df, x = "Subtype", y = "Percentage", hue = "GT", ax = ax[0,0], width=0.8)
# sns.move_legend(ax[0,0], "upper left", bbox_to_anchor=(1.04, 1), frameon = False, fontsize = 15, title = None)
# ax[1,0] = sns.barplot(data = lymphoid_sub_df, x = "Subtype", y = "Mean expr (nonzero)", hue = "GT", ax = ax[1,0], width=0.8)
# sns.move_legend(ax[1,0], "upper left", bbox_to_anchor=(1.04, 1), frameon = False, fontsize = 15, title = None)
# ax[0,0].set_title("Lymphoid Subtype", fontsize = 15)
# ax[0,0].set_ylabel("Pct of expr", fontsize = 15)
# ax[1,0].set_title("Lymphoid Subtype", fontsize = 15)
# ax[1,0].set_ylabel("Mean expr", fontsize = 15)
# for i in ax[0,0].containers:
#     ax[0,0].bar_label(i, fmt = "%.2f", fontsize = 15)
# for i in ax[1,0].containers:
#     ax[1,0].bar_label(i, fmt = "%.2f", fontsize = 15)
# plt.tight_layout()

# ax[0,1] = sns.barplot(data = macrophage_sub_df, x = "Subtype", y = "Percentage", hue = "GT", ax = ax[0,1], width=0.4)
# sns.move_legend(ax[0,1], "upper left", bbox_to_anchor=(1.04, 1), frameon = False, fontsize = 15, title = None)
# ax[1,1] = sns.barplot(data = macrophage_sub_df, x = "Subtype", y = "Mean expr (nonzero)", hue = "GT", ax = ax[1,1], width=0.4)
# sns.move_legend(ax[1,1], "upper left", bbox_to_anchor=(1.04, 1), frameon = False, fontsize = 15, title = None)
# ax[0,1].set_title("Macrophage Subtype", fontsize = 15)
# ax[0,1].set_ylabel("Pct of exp", fontsize = 15)
# ax[1,1].set_title("Macrophage Subtype", fontsize = 15)
# ax[1,1].set_ylabel("Mean expr", fontsize = 15)
# for i in ax[0,1].containers:
#     ax[0,1].bar_label(i, fmt = "%.3f", fontsize = 15)
# for i in ax[1,1].containers:
#     ax[1,1].bar_label(i, fmt = "%.3f", fontsize = 15)
# plt.tight_layout()
# # sns.reset_orig()
# fig.savefig(res_dir + "subtypes_expr.png", bbox_inches = "tight", dpi = 150)

# %%
