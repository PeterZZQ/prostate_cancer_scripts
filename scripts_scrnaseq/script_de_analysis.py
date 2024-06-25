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
import seaborn as sns

import utils
import warnings
warnings.filterwarnings("ignore")


# Volcano plots
import matplotlib.pyplot as plt
import seaborn as sns
# optional
from adjustText import adjust_text

def plot_volcano(df, x = "logFC", y = "adj.P.Val", x_cutoff = 2, y_cutoff = 0.01, gene_name = True, ylim = None, xlim = None):
    fig = plt.figure(figsize = (10, 7))
    ax = fig.add_subplot()
    if ylim is None:
        ylim = np.inf
    if xlim is None:
        xlim = np.inf
    else:
        df.loc[df[x] > xlim, x] = xlim
        df.loc[df[x] < -xlim, x] = -xlim

    ax.scatter(x = df[x], y = df[y].apply(lambda x:-np.log10(max(x, ylim))), s = 1, color = "gray")#, label = "Not significant")
    
    # highlight down- or up- regulated genes
    down = df[(df[x] <= -x_cutoff) & (df[y] <= y_cutoff)]
    up = df[(df[x] >= x_cutoff) & (df[y] <= y_cutoff)]
    ax.scatter(x = down[x], y = down[y].apply(lambda x:-np.log10(max(x, ylim))), s = 3, label = "Down-regulated", color = "blue")
    ax.scatter(x = up[x], y = up[y].apply(lambda x:-np.log10(max(x, ylim))), s = 3, label = "Up-regulated", color = "red")

    # add legends
    ax.set_xlabel("logFC", fontsize = 15)
    ax.set_ylabel("-logPVal", fontsize = 15)
    ax.set_xlim([-np.max(np.abs(df[x].values)) - 0.5, np.max(np.abs(df[x].values)) + 0.5])
    # ax.set_ylim(-3, 50)
    ax.axvline(-x_cutoff, color = "grey", linestyle = "--")
    ax.axvline(x_cutoff, color = "grey", linestyle = "--")
    ax.axhline(-np.log10(y_cutoff), color = "grey", linestyle = "--")
    leg = ax.legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1), markerscale = 3)
    for lh in leg.legendHandles: 
        lh.set_alpha(1)

    # add gene names
    if gene_name:
        texts = []
        for i,r in down.iterrows():
            texts.append(plt.text(x = r[x], y = -np.log10(max(r[y], ylim)), s = i, fontsize = 7))

        for i,r in up.iterrows():
            texts.append(plt.text(x = r[x], y = -np.log10(max(r[y], ylim)), s = i, fontsize = 7))
        # # optional, adjust text
        adjust_text(texts)#,arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
    
    return fig, ax



# In[]
# ---------------------------------------------------------------- #
#
# 1. Data loading
#
# ---------------------------------------------------------------- #
raw_dir = "../dataset/data_scrnaseq/data/qc_data/"
seurat_dir = "../dataset/data_scrnaseq/seurat_integration/"
vis_intact_dir = "../results_scrnaseq/annot_intact/"

adata_intact_raw = sc.read_h5ad(raw_dir + "adata_qc_intact.h5ad")
# log-normalized data
adata_intact_seurat = sc.read_h5ad(seurat_dir + "adata_intact_seurat.h5ad")
# save raw count in the seurat adata
adata_intact_raw.obs = adata_intact_seurat.obs[["annot", "sample", "age", "condition", "genotype", "batch"]]
# load the interest genes
interest_genes = pd.read_csv("../dataset/markers_info/interest_genes.csv", index_col = 0).index.values
interest_genes = [x for x in interest_genes if (x != "Cxcl11") & (x != "Tnfa") & (x != "Tgfb") & (x != "Il12") & (x != "Pd1") & (x != "Tim3") & (x != "Ccl21")]
adata_intact_seurat.var["interest_gene"] = False
adata_intact_raw.var["interest_gene"] = False
adata_intact_seurat.var.loc[interest_genes, "interest_gene"] = True
adata_intact_raw.var.loc[interest_genes, "interest_gene"] = True

# In[]
# ---------------------------------------------------------------- #
#
# 2. Differential expression analysis across different cell types
# aim to find cell-type-specific DE genes, which includes all genotypes
# Wilcoxon Rank Sum test is used
#
# ---------------------------------------------------------------- #

# NOTE: extract only highly-variable genes for de analysis, this is only valid for between cell type de (whole cell population is used)
# the highly-variable genes where calculated in `script_plots_markers.py` on the whole intact cell population
# incorporate the important marker genes too
adata_intact_seurat.var.loc[interest_genes,"highly_variable"] = True
adata_intact_seurat = adata_intact_seurat[:, adata_intact_seurat.var["highly_variable"]]

# For each cell type, compare the gene expression between the cell type and the remaining cells
# NOTE: tie_correct, correct the tied scores (wilcoxon), pts: bool, compute the fraction of cells expressing the gene
# rankby_abs: rank by the absolute value or signed value of scores
sc.tl.rank_genes_groups(adata_intact_seurat, groupby="annot", method="wilcoxon", reference = "rest", tie_correct = True, pts = True, rankby_abs = False)
# returns
DE_ct_dir = "../results_scrnaseq/DE_CT/"
DE_ct_results = adata_intact_seurat.uns["rank_genes_groups"]
for ct in adata_intact_seurat.obs["annot"].cat.categories:
    # the genes are ordered according to the z-score, then used to calculated the pvalues
    DE_df = pd.DataFrame(data = 0, columns = ["zscore", "pval", "pval (adj)", "log_fold", "pts", "pts_rest"], index = DE_ct_results["names"][ct])
    DE_df["zscore"] = DE_ct_results["scores"][ct]
    DE_df["pval"] = DE_ct_results["pvals"][ct]
    DE_df["pval (adj)"] = DE_ct_results["pvals_adj"][ct]
    DE_df["log_fold"] = DE_ct_results["logfoldchanges"][ct]
    DE_df["pts"] = DE_ct_results["pts"][ct]
    DE_df["pts_rest"] = DE_ct_results["pts_rest"][ct]
    # need to separate according to the enriched and depleted genes
    # the enriched and depleted genes are separated according to the mean log-fold change
    DE_df_enrich = DE_df[DE_df["log_fold"] >= 0]
    DE_df_deplete = DE_df[DE_df["log_fold"] < 0]

    DE_df_enrich.to_csv(DE_ct_dir + f"DE_{ct}_enrich.csv", sep = ",")
    DE_df_deplete.to_csv(DE_ct_dir + f"DE_{ct}_deplete.csv", sep = ",")

# In[]
# 2.1. select the marker genes according to the logfold change and pvalues
# draw the volcano plot
DE_ct_dir = "../results_scrnaseq/DE_CT/"
for ct in adata_intact_seurat.obs["annot"].cat.categories:
    DE_df_enrich = pd.read_csv(DE_ct_dir + f"DE_{ct}_enrich.csv", sep = ",", index_col = 0)
    DE_df_deplete = pd.read_csv(DE_ct_dir + f"DE_{ct}_deplete.csv", sep = ",", index_col = 0)
    DE_df = pd.concat([DE_df_enrich, DE_df_deplete], axis = 0, ignore_index = False)
    
    # # [removed] drop the genes that are not expressed at all in ct/background, possibly technical noise
    # DE_df = DE_df[(DE_df["pts"] > 0) & (DE_df["pts_rest"] > 0)]

    fig, ax = plot_volcano(df = DE_df, x = "log_fold", y = "pval (adj)", x_cutoff = 2, y_cutoff = 0.01, gene_name = False, ylim = 10e-15, xlim = 15)
    ax.set_title(ct, fontsize = 20)
    fig.savefig(DE_ct_dir + f"volcano_{ct}.png", bbox_inches = "tight", dpi = 150)

# In[]
# 2.2. Heatmap of the enriched genes for each cluster 
DEs_enriched = []
org_cts = ['Luminal', 'Basal', 'Luminal (Spink1+)', 'Club epithelia', 'Endothelial', 'Lymphoid', 'Macrophage', 'Mesenchymal', 'Monocytes', 'SV']
for ct in org_cts:
    DE_df_enrich = pd.read_csv(DE_ct_dir + f"DE_{ct}_enrich.csv", sep = ",", index_col = 0)
    # # [removed] filter all-zero genes
    # DE_df_enrich = DE_df_enrich[(DE_df_enrich["pts"] > 0) & (DE_df_enrich["pts_rest"] > 0)]
    # filter according to the log-fold change and pvalues
    DE_df_enrich = DE_df_enrich[(DE_df_enrich["log_fold"].values.squeeze() > 2) & (DE_df_enrich["pval (adj)"].values.squeeze() < 0.01)]
    # thresholding extremely small pvalues
    DE_df_enrich.loc[DE_df_enrich["pval (adj)"] < 10e-15, "pval (adj)"] = 10e-15
    # sort according to the pvalues & log-fold change, extract the top-10 genes
    orders = np.lexsort((-DE_df_enrich["log_fold"].values.squeeze(), DE_df_enrich["pval (adj)"].values.squeeze()))[:10]
    DEs_enriched.append(DE_df_enrich.iloc[orders,:].index.values.squeeze())


# re-organize the count by seurat clusters and DE_enriched genes
adata_clusters = []
for ct in org_cts:
    adata_clusters.append(adata_intact_seurat[adata_intact_seurat.obs["annot"] == ct,:])
adata_merge = anndata.AnnData.concatenate(*adata_clusters, join = "inner")
counts = adata_merge[:,np.concatenate(DEs_enriched, axis = 0)].X.toarray()
# z-score
# counts = (counts - np.mean(counts, axis = 1, keepdims = True))/np.sqrt(np.mean((counts - np.mean(counts, axis = 1, keepdims = True))**2, axis = 1, keepdims = True) + 1e-9)

# plot heatmap
plt.rcParams["font.size"] = 10
lut = plt.cm.get_cmap("tab20b", 14)
yticks = []
yticklabels = []
count = 0
for ct in org_cts:

    yticks.append(count + int(np.sum(adata_merge.obs["annot"] == ct)/2))
    count += int(np.sum(adata_intact_seurat.obs["annot"] == ct))
    yticklabels.append(ct)
yticks = np.array(yticks)
yticklabels = np.array(yticklabels)

fig = plt.figure(figsize = (70, 60))
ax = fig.add_subplot()
ax = sns.heatmap(counts, ax = ax)
ax.set_yticks(yticks)
ax.set_xticks(np.arange(np.concatenate(DEs_enriched, axis = 0).shape[0]))
ax.tick_params(axis='y', which='major', pad=30, length=0) # extra padding to leave room for the row colors
ax.set_yticklabels(yticklabels, rotation=45)
ax.set_xticklabels(np.concatenate(DEs_enriched, axis = 0), rotation = 45)

count = 0
for i, ct in enumerate(org_cts):
    ax.add_patch(plt.Rectangle(xy=(-0.007, count), width=0.007, height=int(np.sum(adata_intact_seurat.obs["annot"] == ct)), 
                         color=lut(i), lw=0, transform=ax.get_yaxis_transform(), clip_on=False))
    count += int(np.sum(adata_intact_seurat.obs["annot"] == ct))

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=40)
fig.savefig(DE_ct_dir + "heatmap_annotation.png", bbox_inches = "tight", dpi = 100)

# In[]
# ---------------------------------------------------------------- #
#
# 3. Differential expression analysis across different genotype
# aim to find the genotype-associated DE genes for each cell type separately
# optional, remove the cell type with not enough number of cells
#
# ---------------------------------------------------------------- #
use_interest_genes = True
if use_interest_genes:
    DE_gt_dir = "../results_scrnaseq/DE_GT/with_interest_genes/"
else:
    DE_gt_dir = "../results_scrnaseq/DE_GT/original/"

for ct in adata_intact_seurat.obs["annot"].cat.categories:
    print(ct)
    adata_ct = adata_intact_raw[adata_intact_raw.obs["annot"] == ct]
    print(f"number of cells: {adata_ct.shape[0]}")
    ncells_pp = adata_ct[adata_ct.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"].shape[0]
    ncells_ppc = adata_ct[adata_ct.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"].shape[0]
    # smallest population, endothelial, only 41/75 cells
    print(f"number of cells PP: {ncells_pp}")
    print(f"number of cells PPC: {ncells_ppc}")
    
    # NOTE: redo the filtering and hvg selection on the selected cell population, before de analysis across genotypes
    # [removed]
    # gene_subset1, _ = sc.pp.filter_genes(adata_ct[adata_ct.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"], min_cells = int(ncells_pp/10), inplace = False)
    # gene_subset2, _ = sc.pp.filter_genes(adata_ct[adata_ct.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"], min_cells = int(ncells_ppc/10), inplace = False)
    # gene_subset = (gene_subset1 & gene_subset2) 

    # should be filtering the gene on both condition jointly
    gene_subset, _ = sc.pp.filter_genes(adata_ct, min_cells = int((ncells_pp + ncells_ppc)/10), inplace = False)
    if use_interest_genes:
        gene_subset |= adata_ct.var["interest_gene"].values.squeeze()
    adata_ct = adata_ct[:, gene_subset]
    print("filter genes, number of genes kept: {:d}".format(adata_ct.shape[1]))
    sc.pp.normalize_per_cell(adata_ct)
    sc.pp.log1p(adata_ct)
    sc.pp.highly_variable_genes(adata_ct, n_top_genes = 2000)
    if use_interest_genes:
        adata_ct.var.loc[interest_genes,"highly_variable"] = True
    adata_ct = adata_ct[:, adata_ct.var["highly_variable"]]

    # conduct cross-genotype comparison for both 12wk and 18wk
    # compare PPC with PP, the background is PP, foreground is PPC
    sc.tl.rank_genes_groups(adata_ct, groupby = "genotype", groups = ["PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"], key_added = ct + " (GT)", method = "wilcoxon", 
                            reference = "PbCre(+/-),Pten(-/-),P53(-/-)", tie_correct = True, pts = True, rankby_abs = False)

    DE_results = adata_ct.uns[ct + " (GT)"]
    DE_df = pd.DataFrame(data = 0, columns = ["zscore", "pval", "pval (adj)", "log_fold", "pts", "pts_rest"], index = DE_results["names"]["PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"])
    DE_df["zscore"] = DE_results["scores"]["PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"]
    DE_df["pval"] = DE_results["pvals"]["PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"]
    DE_df["pval (adj)"] = DE_results["pvals_adj"]["PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"]
    DE_df["log_fold"] = DE_results["logfoldchanges"]["PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"]
    DE_df["pts"] = DE_results["pts"]["PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"]

    # # SANITY CHECK: percentage of expression
    # # separation at around 15 according to the volcano plot
    # log_fold_cutoff = 10
    # gene_largelogfold = DE_df[(DE_df["log_fold"] >= log_fold_cutoff) | (DE_df["log_fold"] <= -log_fold_cutoff)].index.values.squeeze()
    # gene_rest = DE_df[(DE_df["log_fold"] <= log_fold_cutoff) & (DE_df["log_fold"] >= -log_fold_cutoff)].index.values.squeeze()
    # pct_df = pd.DataFrame(columns = ["category", "percentage of expr 1", "percentage of expr 2", "smaller pct"], index = np.concatenate([gene_largelogfold, gene_rest], axis = 0), data = 0)
    # pct_df.loc[gene_largelogfold, "category"] = "large_logfold"
    # pct_df.loc[gene_largelogfold, "percentage of expr 1"] = np.sum(adata_ct[adata_ct.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)", gene_largelogfold].X.toarray().squeeze() != 0, axis = 0)/ncells_pp
    # pct_df.loc[gene_largelogfold, "percentage of expr 2"] = np.sum(adata_ct[adata_ct.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", gene_largelogfold].X.toarray().squeeze() != 0, axis = 0)/ncells_ppc
    # pct_df.loc[gene_rest, "category"] = "control"
    # pct_df.loc[gene_rest, "percentage of expr 1"] = np.sum(adata_ct[adata_ct.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)", gene_rest].X.toarray().squeeze() != 0, axis = 0)/ncells_pp
    # pct_df.loc[gene_rest, "percentage of expr 2"] = np.sum(adata_ct[adata_ct.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", gene_rest].X.toarray().squeeze() != 0, axis = 0)/ncells_ppc    
    # group1_large = (pct_df["percentage of expr 1"].values.squeeze() - pct_df["percentage of expr 2"].values.squeeze() >= 0)
    # pct_df.loc[group1_large, "smaller pct"] = pct_df.loc[group1_large, "percentage of expr 2"].values
    # pct_df.loc[~group1_large, "smaller pct"] = pct_df.loc[~group1_large, "percentage of expr 1"].values 
    # fig = plt.figure(figsize = (14,5))
    # ax = fig.subplots(nrows = 1, ncols = 3)
    # sns.boxplot(data = pct_df, x = "category", y = "percentage of expr 1", ax = ax[0])
    # sns.boxplot(data = pct_df, x = "category", y = "percentage of expr 2", ax = ax[1])
    # sns.boxplot(data = pct_df, x = "category", y = "smaller pct", ax = ax[2])
    # ax[0].set_ylim([0,0.4])
    # ax[0].set_title(ct)
    # ax[1].set_ylim([0,0.4])
    # ax[1].set_title(ct)
    # ax[2].set_ylim([0,0.4])
    # ax[2].set_title(ct)
    # print("max percentage of expression (large log-fold): {:.4f}".format(np.max(pct_df.loc[pct_df["category"] == "large_logfold", "smaller pct"].values)))

    # only if reference is "all"
    # DE_df["pts_rest"] = DE_results["pts_rest"]["PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"]
    # need to separate according to the enriched and depleted genes
    # the enriched and depleted genes are separated according to the mean log-fold change
    DE_df_enrich = DE_df[DE_df["log_fold"] >= 0]
    DE_df_deplete = DE_df[DE_df["log_fold"] < 0]

    if use_interest_genes:
        DE_df_enrich.to_csv(DE_gt_dir + f"DE_{ct}_PPC_enrich_interest_genes.csv", sep = ",")
        DE_df_deplete.to_csv(DE_gt_dir + f"DE_{ct}_PPC_deplete_interest_genes.csv", sep = ",")   
    else:
        DE_df_enrich.to_csv(DE_gt_dir + f"DE_{ct}_PPC_enrich.csv", sep = ",")
        DE_df_deplete.to_csv(DE_gt_dir + f"DE_{ct}_PPC_deplete.csv", sep = ",")

# In[]
# draw volcano plot, considering both pval and logfold-change
for ct in adata_intact_seurat.obs["annot"].cat.categories:
    if use_interest_genes:
        DE_df_enrich = pd.read_csv(DE_gt_dir + f"DE_{ct}_PPC_enrich_interest_genes.csv", sep = ",", index_col = 0)
        DE_df_deplete = pd.read_csv(DE_gt_dir + f"DE_{ct}_PPC_deplete_interest_genes.csv", sep = ",", index_col = 0)
    else:
        DE_df_enrich = pd.read_csv(DE_gt_dir + f"DE_{ct}_PPC_enrich.csv", sep = ",", index_col = 0)
        DE_df_deplete = pd.read_csv(DE_gt_dir + f"DE_{ct}_PPC_deplete.csv", sep = ",", index_col = 0)

    DE_df = pd.concat([DE_df_enrich, DE_df_deplete], axis = 0, ignore_index = False)
    fig, ax = plot_volcano(df = DE_df, x = "log_fold", y = "pval (adj)", x_cutoff = 1.5, y_cutoff = 0.01, gene_name = True, ylim = 10e-15, xlim = 15)
    ax.set_title(ct, fontsize = 20)
    if use_interest_genes:
        fig.savefig(DE_gt_dir + f"volcano_{ct}_interest_genes.png", bbox_inches = "tight", dpi = 150)
    else:
        fig.savefig(DE_gt_dir + f"volcano_{ct}.png", bbox_inches = "tight", dpi = 150)


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      # %%
