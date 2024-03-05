# In[]
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc
from anndata import AnnData
import scipy as sci
import warnings
from matplotlib.colors import LinearSegmentedColormap
warnings.filterwarnings("ignore")
from scipy.stats import ranksums



# In[]
# ---------------------------------------------------------------- #
#
# Preprocessing, read in the intact and castrated dataset
#
# ---------------------------------------------------------------- #

# Intact
adata = sc.read_h5ad("Cell_Ranger_output/adata_seurat.h5ad")
# adata.var['MT_gene'] = [gene.startswith('mt-') for gene in adata.var.index]
# adata = adata[:, ~adata.var['MT_gene'].values]
# sc.pp.filter_genes(adata, min_cells = 1)
# adata.write_h5ad("Cell_Ranger_output/adata_seurat.h5ad")
# sc.pp.filter_cells(adata, min_genes = 200)
sc.pp.filter_genes(adata, min_cells = 3)
# keep the original count before normalization and log transform
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.neighbors(adata, n_neighbors = 30)
sc.tl.umap(adata)

# Castrated
# adata = sc.read_h5ad("Cell_Ranger_output/adata_castrated.h5ad")
# meta_castrated = pd.read_csv("Cell_Ranger_output/meta_castrated.csv", index_col = 0, sep = "\t")
# adata.obs["annot_transfer"] = meta_castrated.loc[adata.obs.index.values, "annot.transfer"]
# adata.X = adata.X.astype(np.float64)
# adata.write_h5ad("Cell_Ranger_output/adata_castrated_seurat.h5ad")

adata_castrated = sc.read_h5ad("Cell_Ranger_output/adata_castrated_seurat.h5ad")
seurat_pca = pd.read_csv("Cell_Ranger_output/xpca_castrated_integrated.csv", sep = "\t", index_col = 0)
seurat_umap = pd.read_csv("Cell_Ranger_output/xumap_castrated_integrated.csv", sep = "\t", index_col = 0)

sc.pp.filter_genes(adata_castrated, min_cells = 3)
# keep the original count before normalization and log transform
adata_castrated.layers["counts"] = adata_castrated.X.copy()

sc.pp.normalize_per_cell(adata_castrated)
sc.pp.log1p(adata_castrated)
sc.pp.neighbors(adata_castrated, n_neighbors = 30)
sc.tl.umap(adata_castrated)

# In[]
# ---------------------------------------------------------------- #
#
# Plot the UMAP visualization
#
# ---------------------------------------------------------------- #
plot_latent(adata.obsm["X_umap"], annos = adata.obs["transfer_label"].squeeze(), batches = adata.obs["sample"], mode = "separate", figsize = (10,20), s = 1, markerscale = 6, axis_label = "UMAP", save = "results_seurat_scrnaseq/xumap_transferlabel_separate.png")
plot_latent(adata.obsm["X_umap"], annos = adata.obs["transfer_label"].squeeze(), batches = adata.obs["sample"], mode = "annos", figsize = (10,5), s = 1, markerscale = 6, axis_label = "UMAP", save = "results_seurat_scrnaseq/xumap_transferlabel.png")
plot_latent(adata.obsm["X_umap"], annos = np.array([eval(x) for x in adata.obs["seurat_clusters"].squeeze()]), batches = adata.obs["sample"], mode = "separate", figsize = (10,30), s = 1, markerscale = 6, axis_label = "UMAP", save = "results_seurat_scrnaseq/xumap_clusters_separate.png", label_inplace = True)
plot_latent(adata.obsm["X_umap"], annos = np.array([eval(x) for x in adata.obs["seurat_clusters"].squeeze()]), batches = adata.obs["sample"], mode = "annos", figsize = (10,7), s = 1, markerscale = 6, axis_label = "UMAP", save = "results_seurat_scrnaseq/xumap_clusters.png", label_inplace = True)
plot_latent(adata.obsm["X_umap"], annos = adata.obs["annot"].squeeze(), batches = adata.obs["sample"], mode = "separate", figsize = (10,20), s = 1, markerscale = 6, axis_label = "UMAP", save = "results_seurat_scrnaseq/xumap_annotation.png")

plot_latent(adata_castrated.obsm["X_umap"], annos = adata_castrated.obs["annot_transfer"].squeeze(), batches = adata_castrated.obs["sample"], mode = "separate", figsize = (10,10), s = 1, markerscale = 6, axis_label = "UMAP", save = "results_seurat_scrnaseq/xumap_castrated_separate.png")
plot_latent(adata_castrated.obsm["X_umap"], annos = adata_castrated.obs["annot_transfer"].squeeze(), batches = adata_castrated.obs["sample"], mode = "annos", figsize = (10,5), s = 1, markerscale = 6, axis_label = "UMAP", save = "results_seurat_scrnaseq/xumap_castrated.png")

plot_latent(seurat_umap.values, annos = adata_castrated.obs["annot_transfer"].squeeze(), batches = adata_castrated.obs["sample"], mode = "separate", figsize = (10,10), s = 1, markerscale = 6, axis_label = "UMAP", save = "results_seurat_scrnaseq/intumap_castrated_separate.png")
plot_latent(seurat_umap.values, annos = adata_castrated.obs["annot_transfer"].squeeze(), batches = adata_castrated.obs["sample"], mode = "annos", figsize = (10,5), s = 1, markerscale = 6, axis_label = "UMAP", save = "results_seurat_scrnaseq/intumap_castrated.png")

# Macrophage stages: M1: Mrc1, Arg1, M2: Cd86, Nos2
adata_castrated_pp = adata_castrated[adata_castrated.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"]
adata_castrated_ppc = adata_castrated[adata_castrated.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"]
fig = sc.pl.umap(adata_castrated_pp, color = ["Mrc1", "Arg1", "Cd86", "Nos2"], vmin = 0, vmax = 6, return_fig = True)
fig.savefig("results_seurat_scrnaseq/castrated_pp_macrophage_markers.png", bbox_inches = "tight", dpi = 200)
fig2 = sc.pl.umap(adata_castrated_ppc, color = ["Mrc1", "Arg1", "Cd86", "Nos2"], vmin = 0, vmax = 6, return_fig = True)
fig2.savefig("results_seurat_scrnaseq/castrated_ppc_macrophage_markers.png", bbox_inches = "tight", dpi = 200)

# In[]
# ---------------------------------------------------------------- # 
#
# Trajectory inference from Basal cell to Luminal cell to Mesenchymal cells
#
# ---------------------------------------------------------------- #
adata_cas_traj = adata_castrated[(adata_castrated.obs["annot_transfer"] == "Basal") | (adata_castrated.obs["annot_transfer"] == "Luminal") | (adata_castrated.obs["annot_transfer"] == "Mesenchymal"),:]
# adata_cas_traj.uns['iroot'] = 3765 # np.flatnonzero(adata_cas_traj.obs['annot_transfer'] == 'Basal')[0]
# sc.tl.dpt(adata_cas_traj, n_dcs=10, n_branchings=0, min_group_size=0.01, allow_kendall_tau_shift=True, neighbors_key=None, copy=False)

pt_slingshot = pd.read_csv("results_ti/castrated_basal_mesenchymal/epithelial_mesenchymal_pt.csv", sep = "\t")
adata_cas_traj.obs.loc[pt_slingshot.index.values.squeeze(), "pt_slingshot"] = pt_slingshot["pt_slingshot"].values
# adata_cas_traj.obsm["X_umap"] = seurat_umap.loc[adata_cas_traj.obs.index.values,:].values
with plt.rc_context({"figure.figsize": (12, 8), "figure.dpi": (300)}):
    fig = sc.pl.umap(adata_cas_traj, color = ["pt_slingshot", "annot_transfer", "genotype"], return_fig = True)
fig.savefig("results_ti/castrated_basal_mesenchymal/vis_traj.png", bbox_inches = "tight")
fig = sc.pl.umap(adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)",:], color = ["pt_slingshot", "annot_transfer"], return_fig = True)
fig.savefig("results_ti/castrated_basal_mesenchymal/vis_traj_pp.png", bbox_inches = "tight")
fig = sc.pl.umap(adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)",:], color = ["pt_slingshot","annot_transfer"], return_fig = True)
fig.savefig("results_ti/castrated_basal_mesenchymal/vis_traj_ppc.png", bbox_inches = "tight")

# In[]
# EMT marker list: 
# 1. Mesenchymal:
emt_mesen = ["Vim", "Twist", "Sparc", "Zeb1", "Cd44", "Smad4", "Fn1", "Zeb2", "Mmp2", "Cdh2", "Mmp14", "Mmp16",\
             "Tgfb1", "Fgf2", "Cdh11", "Smad2", "Snai2", "Fgfr1", "Tnc", "Col3a1", "Foxc2", "Snai1", "Gsc", "Itbg6",\
                "Mmp3", "Mmp9", "Sox10", "Dpysl3", "Alcam", "Sdc1", "Nt5e"]

emt_mesen = ["Vim", "Sparc", "Zeb1", "Cd44", "Smad4", "Fn1", "Zeb2", "Mmp2", "Cdh2", "Mmp14", "Mmp16",\
             "Tgfb1", "Fgf2", "Cdh11", "Smad2", "Snai2", "Fgfr1", "Tnc", "Col3a1", "Foxc2", "Snai1", "Gsc",\
                "Mmp3", "Mmp9", "Sox10", "Dpysl3", "Alcam", "Sdc1", "Nt5e"]
# Epithelial
emt_epith = ["Cldn7", "Cldn4", "Krt8", "Jup", "Cdh1", "Ocln", "Krt18", "Ctnnd1", "Krt14", "Dsp", "Tjp1", "Ovol1",\
             "Esrp1", "Lsr", "S100a14", "C1orf116", "Cldn1", "Dsg3", "Dsg2", "Krt19"]

emt_epith = ["Cldn7", "Cldn4", "Krt8", "Jup", "Cdh1", "Ocln", "Krt18", "Ctnnd1", "Krt14", "Dsp", "Tjp1", "Ovol1",\
             "Esrp1", "Lsr", "S100a14", "Cldn1", "Dsg3", "Dsg2", "Krt19"]

adata_cas_traj = adata_castrated[(adata_castrated.obs["annot_transfer"] == "Basal") | (adata_castrated.obs["annot_transfer"] == "Luminal") | (adata_castrated.obs["annot_transfer"] == "Mesenchymal"),:]
adata_intact_traj = adata[(adata.obs["annot"] == "Basal") | (adata.obs["annot"] == "Luminal") | (adata.obs["annot"] == "Mesenchymal"),:]



with plt.rc_context({"figure.figsize": (12, 10), "figure.dpi": (300)}):
    fig = sc.pl.umap(adata_cas_traj, color = ["annot_transfer"] + [x for x in emt_mesen], return_fig = True)
fig.savefig("results_ti/castrated_basal_mesenchymal/vis_emt_mesen.pdf", bbox_inches = "tight")

with plt.rc_context({"figure.figsize": (12, 10), "figure.dpi": (300)}):
    fig = sc.pl.umap(adata_cas_traj, color = ["annot_transfer"] + [x for x in emt_epith], return_fig = True)
fig.savefig("results_ti/castrated_basal_mesenchymal/vis_emt_epith.pdf", bbox_inches = "tight")

# separate, mesenchymal
with plt.rc_context({"figure.figsize": (12, 10), "figure.dpi": (300)}):
    fig = sc.pl.umap(adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)",:], color = ["annot_transfer"] + [x for x in emt_mesen], return_fig = True)
fig.savefig("results_ti/castrated_basal_mesenchymal/vis_emt_mesen_pp.pdf", bbox_inches = "tight")

with plt.rc_context({"figure.figsize": (12, 10), "figure.dpi": (300)}):
    fig = sc.pl.umap(adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)",:], color = ["annot_transfer"] + [x for x in emt_mesen], return_fig = True)
fig.savefig("results_ti/castrated_basal_mesenchymal/vis_emt_mesen_ppc.pdf", bbox_inches = "tight")

# separate, epithelial
with plt.rc_context({"figure.figsize": (12, 10), "figure.dpi": (300)}):
    fig = sc.pl.umap(adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)",:], color = ["annot_transfer"] + [x for x in emt_epith], return_fig = True)
fig.savefig("results_ti/castrated_basal_mesenchymal/vis_emt_epith_pp.pdf", bbox_inches = "tight")

with plt.rc_context({"figure.figsize": (12, 10), "figure.dpi": (300)}):
    fig = sc.pl.umap(adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)",:], color = ["annot_transfer"] + [x for x in emt_epith], return_fig = True)
fig.savefig("results_ti/castrated_basal_mesenchymal/vis_emt_epith_ppc.pdf", bbox_inches = "tight")

# In[]
for gene in emt_mesen:
    X_cas_pp = adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)",gene].X.toarray().squeeze()
    X_cas_ppc = adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)",gene].X.toarray().squeeze()
    meta_cas_pp = adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)",:].obs["annot_transfer"].values.astype(str)
    meta_cas_ppc = adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)",:].obs["annot_transfer"].values.astype(str)
    meta_cas_pp[(meta_cas_pp == "Basal") | (meta_cas_pp == "Luminal")] = "Epith"
    meta_cas_ppc[(meta_cas_ppc == "Basal") | (meta_cas_ppc == "Luminal")] = "Epith"
    meta_cas_pp[meta_cas_pp == "Mesenchymal"] = "Mesen"
    meta_cas_ppc[meta_cas_ppc == "Mesenchymal"] = "Mesen"


    X_intact_pp12 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1417-PP12", gene].X.toarray().squeeze()
    X_intact_pp18 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1416-PP18", gene].X.toarray().squeeze()
    meta_intact_pp12 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1417-PP12",:].obs["annot"].values.astype(str)
    meta_intact_pp18 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1416-PP18",:].obs["annot"].values.astype(str)
    meta_intact_pp12[(meta_intact_pp12 == "Basal") | (meta_intact_pp12 == "Luminal")] = "Epith"
    meta_intact_pp18[(meta_intact_pp18 == "Basal") | (meta_intact_pp18 == "Luminal")] = "Epith"
    meta_intact_pp12[meta_intact_pp12 == "Mesenchymal"] = "Mesen"
    meta_intact_pp18[meta_intact_pp18 == "Mesenchymal"] = "Mesen"

    X_intact_ppc12 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1436-PPC12", gene].X.toarray().squeeze()
    X_intact_ppc18 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1437-PPC18", gene].X.toarray().squeeze()
    meta_intact_ppc12 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1436-PPC12",:].obs["annot"].values.astype(str)
    meta_intact_ppc18 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1437-PPC18",:].obs["annot"].values.astype(str)
    meta_intact_ppc12[(meta_intact_ppc12 == "Basal") | (meta_intact_ppc12 == "Luminal")] = "Epith"
    meta_intact_ppc18[(meta_intact_ppc18 == "Basal") | (meta_intact_ppc18 == "Luminal")] = "Epith"
    meta_intact_ppc12[meta_intact_ppc12 == "Mesenchymal"] = "Mesen"
    meta_intact_ppc18[meta_intact_ppc18 == "Mesenchymal"] = "Mesen"

    df_pp = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_pp["expr"] = X_cas_pp
    df_pp["genotype"] = "PP"
    df_pp["annot"] = meta_cas_pp

    df_ppc = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_ppc["expr"] = X_cas_ppc
    df_ppc["genotype"] = "PPC"
    df_ppc["annot"] = meta_cas_ppc

    df_pp12 = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_pp12["expr"] = X_intact_pp12
    df_pp12["genotype"] = "PP"
    df_pp12["annot"] = meta_intact_pp12

    df_pp18 = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_pp18["expr"] = X_intact_pp18
    df_pp18["genotype"] = "PP"
    df_pp18["annot"] = meta_intact_pp18

    df_ppc12 = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_ppc12["expr"] = X_intact_ppc12
    df_ppc12["genotype"] = "PPC"
    df_ppc12["annot"] = meta_intact_ppc12

    df_ppc18 = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_ppc18["expr"] = X_intact_ppc18
    df_ppc18["genotype"] = "PPC"
    df_ppc18["annot"] = meta_intact_ppc18

    df_cas = pd.concat([df_pp, df_ppc], axis = 0)
    df_12 = pd.concat([df_pp12, df_ppc12], axis = 0)
    df_18 = pd.concat([df_pp18, df_ppc18], axis = 0)

    res_cas = ranksums(x = df_pp.loc[df_pp["annot"] == "Mesen", "expr"].values.squeeze(), y = df_ppc.loc[df_ppc["annot"] == "Mesen", "expr"].values.squeeze(), alternative = "two-sided")
    res_12 = ranksums(x = df_pp12.loc[df_pp12["annot"] == "Mesen", "expr"].values.squeeze(), y = df_ppc12.loc[df_ppc12["annot"] == "Mesen", "expr"].values.squeeze(), alternative = "two-sided")
    res_18 = ranksums(x = df_pp18.loc[df_pp18["annot"] == "Mesen", "expr"].values.squeeze(), y = df_ppc18.loc[df_ppc18["annot"] == "Mesen", "expr"].values.squeeze(), alternative = "two-sided")

    plt.rcParams["font.size"] = 20
    fig = plt.figure(figsize = (12,5))
    ax = fig.subplots(nrows = 1, ncols = 3)
    sns.boxplot(data = df_cas, x = "annot", y = "expr", hue = "genotype", ax = ax[0], order = ["Epith", "Mesen"], fliersize = 1)
    sns.boxplot(data = df_12, x = "annot", y = "expr", hue = "genotype", ax = ax[1], order = ["Epith", "Mesen"], fliersize = 1)
    sns.boxplot(data = df_18, x = "annot", y = "expr", hue = "genotype", ax = ax[2], order = ["Epith", "Mesen"], fliersize = 1)
    fig.suptitle(gene)
    leg = ax[0].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
    leg = ax[1].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
    leg = ax[2].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
    ax[0].set_xlabel(None)
    ax[1].set_xlabel(None)
    ax[2].set_xlabel(None)
    ax[0].set_title("cast\npvalue: {:.2e}".format(res_cas.pvalue))
    ax[1].set_title("12wk\npvalue: {:.2e}".format(res_12.pvalue))
    ax[2].set_title("18wk\npvalue: {:.2e}".format(res_18.pvalue))
    
    fig.tight_layout()
    fig.savefig(f"results_ti/castrated_basal_mesenchymal/emt_markers/emt_mesen_{gene}.png", bbox_inches = "tight")

for gene in emt_epith:
    X_cas_pp = adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)",gene].X.toarray().squeeze()
    X_cas_ppc = adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)",gene].X.toarray().squeeze()
    meta_cas_pp = adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)",:].obs["annot_transfer"].values.astype(str)
    meta_cas_ppc = adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)",:].obs["annot_transfer"].values.astype(str)
    meta_cas_pp[(meta_cas_pp == "Basal") | (meta_cas_pp == "Luminal")] = "Epith"
    meta_cas_ppc[(meta_cas_ppc == "Basal") | (meta_cas_ppc == "Luminal")] = "Epith"
    meta_cas_pp[meta_cas_pp == "Mesenchymal"] = "Mesen"
    meta_cas_ppc[meta_cas_ppc == "Mesenchymal"] = "Mesen"

    X_intact_pp12 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1417-PP12", gene].X.toarray().squeeze()
    X_intact_pp18 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1416-PP18", gene].X.toarray().squeeze()
    meta_intact_pp12 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1417-PP12",:].obs["annot"].values.astype(str)
    meta_intact_pp18 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1416-PP18",:].obs["annot"].values.astype(str)
    meta_intact_pp12[(meta_intact_pp12 == "Basal") | (meta_intact_pp12 == "Luminal")] = "Epith"
    meta_intact_pp18[(meta_intact_pp18 == "Basal") | (meta_intact_pp18 == "Luminal")] = "Epith"
    meta_intact_pp12[meta_intact_pp12 == "Mesenchymal"] = "Mesen"
    meta_intact_pp18[meta_intact_pp18 == "Mesenchymal"] = "Mesen"

    X_intact_ppc12 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1436-PPC12", gene].X.toarray().squeeze()
    X_intact_ppc18 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1437-PPC18", gene].X.toarray().squeeze()
    meta_intact_ppc12 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1436-PPC12",:].obs["annot"].values.astype(str)
    meta_intact_ppc18 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1437-PPC18",:].obs["annot"].values.astype(str)
    meta_intact_ppc12[(meta_intact_ppc12 == "Basal") | (meta_intact_ppc12 == "Luminal")] = "Epith"
    meta_intact_ppc18[(meta_intact_ppc18 == "Basal") | (meta_intact_ppc18 == "Luminal")] = "Epith"
    meta_intact_ppc12[meta_intact_ppc12 == "Mesenchymal"] = "Mesen"
    meta_intact_ppc18[meta_intact_ppc18 == "Mesenchymal"] = "Mesen"

    df_pp = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_pp["expr"] = X_cas_pp
    df_pp["genotype"] = "PP"
    df_pp["annot"] = meta_cas_pp

    df_ppc = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_ppc["expr"] = X_cas_ppc
    df_ppc["genotype"] = "PPC"
    df_ppc["annot"] = meta_cas_ppc

    df_pp12 = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_pp12["expr"] = X_intact_pp12
    df_pp12["genotype"] = "PP"
    df_pp12["annot"] = meta_intact_pp12

    df_pp18 = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_pp18["expr"] = X_intact_pp18
    df_pp18["genotype"] = "PP"
    df_pp18["annot"] = meta_intact_pp18

    df_ppc12 = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_ppc12["expr"] = X_intact_ppc12
    df_ppc12["genotype"] = "PPC"
    df_ppc12["annot"] = meta_intact_ppc12

    df_ppc18 = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_ppc18["expr"] = X_intact_ppc18
    df_ppc18["genotype"] = "PPC"
    df_ppc18["annot"] = meta_intact_ppc18   

    df_cas = pd.concat([df_pp, df_ppc], axis = 0)
    df_12 = pd.concat([df_pp12, df_ppc12], axis = 0)
    df_18 = pd.concat([df_pp18, df_ppc18], axis = 0)

    res_cas = ranksums(x = df_pp.loc[df_pp["annot"] == "Epith", "expr"].values.squeeze(), y = df_ppc.loc[df_ppc["annot"] == "Epith", "expr"].values.squeeze(), alternative = "two-sided")
    res_12 = ranksums(x = df_pp12.loc[df_pp12["annot"] == "Epith", "expr"].values.squeeze(), y = df_ppc12.loc[df_ppc12["annot"] == "Epith", "expr"].values.squeeze(), alternative = "two-sided")
    res_18 = ranksums(x = df_pp18.loc[df_pp18["annot"] == "Epith", "expr"].values.squeeze(), y = df_ppc18.loc[df_ppc18["annot"] == "Epith", "expr"].values.squeeze(), alternative = "two-sided")

    plt.rcParams["font.size"] = 20
    fig = plt.figure(figsize = (12,5))
    ax = fig.subplots(nrows = 1, ncols = 3)
    sns.boxplot(data = df_cas, x = "annot", y = "expr", hue = "genotype", ax = ax[0], order = ["Epith", "Mesen"], fliersize = 1)
    sns.boxplot(data = df_12, x = "annot", y = "expr", hue = "genotype", ax = ax[1], order = ["Epith", "Mesen"], fliersize = 1)
    sns.boxplot(data = df_18, x = "annot", y = "expr", hue = "genotype", ax = ax[2], order = ["Epith", "Mesen"], fliersize = 1)
    fig.suptitle(gene)
    leg = ax[0].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
    leg = ax[1].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
    leg = ax[2].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
    ax[0].set_xlabel(None)
    ax[1].set_xlabel(None)
    ax[2].set_xlabel(None)
    ax[0].set_title("cast\npvalue: {:.2e}".format(res_cas.pvalue))
    ax[1].set_title("12wk\npvalue: {:.2e}".format(res_12.pvalue))
    ax[2].set_title("18wk\npvalue: {:.2e}".format(res_18.pvalue))

    fig.tight_layout()
    fig.savefig(f"results_ti/castrated_basal_mesenchymal/emt_markers/emt_epith_{gene}.png", bbox_inches = "tight")


# In[]
# Compare Ackr3 expression for conditions 
genes = ["Ackr3"]
for gene in genes:
    X_cas_pp = adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)",gene].X.toarray().squeeze()
    X_cas_ppc = adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)",gene].X.toarray().squeeze()
    meta_cas_pp = adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)",:].obs["annot_transfer"].values.astype(str)
    meta_cas_ppc = adata_cas_traj[adata_cas_traj.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)",:].obs["annot_transfer"].values.astype(str)
    meta_cas_pp[(meta_cas_pp == "Basal") | (meta_cas_pp == "Luminal")] = "Epith"
    meta_cas_ppc[(meta_cas_ppc == "Basal") | (meta_cas_ppc == "Luminal")] = "Epith"
    meta_cas_pp[meta_cas_pp == "Mesenchymal"] = "Mesen"
    meta_cas_ppc[meta_cas_ppc == "Mesenchymal"] = "Mesen"

    X_intact_pp12 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1417-PP12", gene].X.toarray().squeeze()
    X_intact_pp18 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1416-PP18", gene].X.toarray().squeeze()
    meta_intact_pp12 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1417-PP12",:].obs["annot"].values.astype(str)
    meta_intact_pp18 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1416-PP18",:].obs["annot"].values.astype(str)
    meta_intact_pp12[(meta_intact_pp12 == "Basal") | (meta_intact_pp12 == "Luminal")] = "Epith"
    meta_intact_pp18[(meta_intact_pp18 == "Basal") | (meta_intact_pp18 == "Luminal")] = "Epith"
    meta_intact_pp12[meta_intact_pp12 == "Mesenchymal"] = "Mesen"
    meta_intact_pp18[meta_intact_pp18 == "Mesenchymal"] = "Mesen"

    X_intact_ppc12 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1436-PPC12", gene].X.toarray().squeeze()
    X_intact_ppc18 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1437-PPC18", gene].X.toarray().squeeze()
    meta_intact_ppc12 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1436-PPC12",:].obs["annot"].values.astype(str)
    meta_intact_ppc18 = adata_intact_traj[adata_intact_traj.obs["sample"] == "M1437-PPC18",:].obs["annot"].values.astype(str)
    meta_intact_ppc12[(meta_intact_ppc12 == "Basal") | (meta_intact_ppc12 == "Luminal")] = "Epith"
    meta_intact_ppc18[(meta_intact_ppc18 == "Basal") | (meta_intact_ppc18 == "Luminal")] = "Epith"
    meta_intact_ppc12[meta_intact_ppc12 == "Mesenchymal"] = "Mesen"
    meta_intact_ppc18[meta_intact_ppc18 == "Mesenchymal"] = "Mesen"

    df_pp = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_pp["expr"] = X_cas_pp
    df_pp["genotype"] = "PP"
    df_pp["annot"] = meta_cas_pp

    df_ppc = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_ppc["expr"] = X_cas_ppc
    df_ppc["genotype"] = "PPC"
    df_ppc["annot"] = meta_cas_ppc

    df_pp12 = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_pp12["expr"] = X_intact_pp12
    df_pp12["genotype"] = "PP"
    df_pp12["annot"] = meta_intact_pp12

    df_pp18 = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_pp18["expr"] = X_intact_pp18
    df_pp18["genotype"] = "PP"
    df_pp18["annot"] = meta_intact_pp18

    df_ppc12 = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_ppc12["expr"] = X_intact_ppc12
    df_ppc12["genotype"] = "PPC"
    df_ppc12["annot"] = meta_intact_ppc12

    df_ppc18 = pd.DataFrame(columns = ["expr", "genotype", "annot"])
    df_ppc18["expr"] = X_intact_ppc18
    df_ppc18["genotype"] = "PPC"
    df_ppc18["annot"] = meta_intact_ppc18

    df_cas = pd.concat([df_pp, df_ppc], axis = 0)
    df_12 = pd.concat([df_pp12, df_ppc12], axis = 0)
    df_18 = pd.concat([df_pp18, df_ppc18], axis = 0)

    res = ranksums(x = df_cas.loc[df_cas["annot"] == "Mesen", "expr"].values.squeeze(), y = np.concatenate([df_12.loc[df_12["annot"] == "Mesen", "expr"].values.squeeze(), 
                                                                                                                df_18.loc[df_18["annot"] == "Mesen", "expr"].values.squeeze()], axis = 0), alternative = "two-sided")

    plt.rcParams["font.size"] = 20
    fig = plt.figure(figsize = (12,5))
    ax = fig.subplots(nrows = 1, ncols = 3)
    sns.boxplot(data = df_cas, x = "annot", y = "expr", hue = "genotype", ax = ax[0], order = ["Epith", "Mesen"], fliersize = 1)
    sns.boxplot(data = df_12, x = "annot", y = "expr", hue = "genotype", ax = ax[1], order = ["Epith", "Mesen"], fliersize = 1)
    sns.boxplot(data = df_18, x = "annot", y = "expr", hue = "genotype", ax = ax[2], order = ["Epith", "Mesen"], fliersize = 1)
    fig.suptitle(gene)
    leg = ax[0].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
    leg = ax[1].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
    leg = ax[2].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
    ax[0].set_title("Cas")
    ax[1].set_title("Int-12")
    ax[2].set_title("Int-18")
    ax[0].set_xlabel(None)
    ax[1].set_xlabel(None)
    ax[2].set_xlabel(None)
    ax[0].set_ylim(-0.3, 4.5)
    ax[1].set_ylim(-0.3, 4.5)
    ax[2].set_ylim(-0.3, 4.5)
    fig.suptitle("p value: {:.3e}".format(res.pvalue))

    fig.tight_layout()
    fig.savefig(f"results_ti/castrated_basal_mesenchymal/emt_markers/{gene}.png", bbox_inches = "tight")

# In[]
# ---------------------------------------------------------------- #
#
# Heatmap: differentially expressed gene 
#
# ---------------------------------------------------------------- #
# read in the enriched and depleted genes of each cluster
DEs_enriched = []
DEs_depleted = []
DEs = []
for idx in adata.obs["annot"].cat.categories:
    DE_enriched = pd.read_csv(f"results_seurat_scrnaseq/DE/DE_{idx}_enriched.csv").index.values[:10]
    DE_depleted = pd.read_csv(f"results_seurat_scrnaseq/DE/DE_{idx}_depleted.csv").index.values[:10]
    DEs_enriched.append(DE_enriched)
    DEs_depleted.append(DE_depleted)
    DEs.append(np.concatenate([DE_enriched, DE_depleted], axis = 0))

# group by seurat clusters
adata_clusters = []
for idx in adata.obs["annot"].cat.categories:
    adata_clusters.append(adata[adata.obs["annot"] == idx,:])
adata_merge = AnnData.concatenate(*adata_clusters, join = "inner")
counts = adata_merge[:,np.concatenate(DEs_enriched, axis = 0)].X.toarray()

# z-score
counts = (counts - np.mean(counts, axis = 1, keepdims = True))/np.sqrt(np.mean((counts - np.mean(counts, axis = 1, keepdims = True))**2, axis = 1, keepdims = True) + 1e-9)

# plot heatmap
plt.rcParams["font.size"] = 10
lut = plt.cm.get_cmap("tab20b", 14)
yticks = []
yticklabels = []
count = 0
for idx in adata.obs["annot"].cat.categories:

    yticks.append(count + int(np.sum(adata_merge.obs["annot"] == idx)/2))
    count += int(np.sum(adata.obs["annot"] == idx))
    yticklabels.append(idx)
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
for i, idx in enumerate(adata.obs["annot"].cat.categories):
    ax.add_patch(plt.Rectangle(xy=(-0.007, count), width=0.007, height=int(np.sum(adata.obs["annot"] == idx)), 
                         color=lut(i), lw=0, transform=ax.get_yaxis_transform(), clip_on=False))
    count += int(np.sum(adata.obs["annot"] == idx))

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=40)
fig.savefig("results_seurat_scrnaseq/DE/heatmap_annotation.png", bbox_inches = "tight", dpi = 100)


# In[]
# ---------------------------------------------------------------- #
#
# Plot spatial interest genes with UMAPs
#
# ---------------------------------------------------------------- #
plt.rcParams["font.size"] = 15

interest_genes = pd.read_csv("spatial_data/interest_genes.csv", index_col = 0)
gene_names =[x for x in interest_genes["Genes"] if (x != "Cxcl11") & (x != "Tnfa") & (x != "Tgfb") & (x != "Il12") & (x != "Pd1") & (x != "Tim3")]

SUPER_MAGMA = LinearSegmentedColormap.from_list('super_magma', colors=['#e0e0e0', '#dedede', '#fff68f', '#ffec8b', '#ffc125', '#ee7600', '#ee5c42', '#cd3278', '#c71585', '#68228b'], N=500)

fig = plt.figure(figsize = (20, 10))
ax = fig.subplots(nrows = 2, ncols = 2)
adata_pp12 = adata[adata.obs["sample"] == "M1417-PP12",:]
adata_pp18 = adata[adata.obs["sample"] == "M1416-PP18",:]
adata_ppc12 = adata[adata.obs["sample"] == "M1436-PPC12",:]
adata_ppc18 = adata[adata.obs["sample"] == "M1437-PPC18",:]

sc.pl.umap(adata_pp12, color = "annot",save = None, ax = ax[0,0], show=False)
sc.pl.umap(adata_ppc12, color = "annot", save = None, ax = ax[0,1], show=False)
sc.pl.umap(adata_pp18, color = "annot",save = None, ax = ax[1,0], show=False)
sc.pl.umap(adata_ppc18, color = "annot", save = None, ax = ax[1,1], show=False)
# ax[0,0].tick_params(axis='both', which='major', labelsize=10)
# ax[0,1].tick_params(axis='both', which='major', labelsize=10)
# ax[1,0].tick_params(axis='both', which='major', labelsize=10)
# ax[1,1].tick_params(axis='both', which='major', labelsize=10)

ax[0,0].set_title("PP12")
ax[0,1].set_title("PP18")
ax[1,0].set_title("PPC12")
ax[1,1].set_title("PPC18")
fig.tight_layout()

fig.savefig(f"results_seurat_scrnaseq/figure_markers/important_genes/annot.png", dpi = 200, bbox_inches = "tight")

for gene in gene_names:
    fig = plt.figure(figsize = (15, 10))
    ax = fig.subplots(nrows = 2, ncols = 2)

    vmin = min(np.min(adata_pp12[:, gene].X), np.min(adata_ppc12[:, gene].X), np.min(adata_pp18[:, gene].X), np.min(adata_ppc18[:, gene].X))
    vmax = max(np.max(adata_pp12[:, gene].X), np.max(adata_ppc12[:, gene].X), np.max(adata_pp18[:, gene].X), np.max(adata_ppc18[:, gene].X))

    sc.pl.umap(adata_pp12, color = gene, save = None, ax = ax[0,0], show=False, vmin = vmin, vmax = vmax, color_map = SUPER_MAGMA)
    sc.pl.umap(adata_ppc12, color = gene, save = None, ax = ax[0,1], show=False, vmin = vmin, vmax = vmax, color_map = SUPER_MAGMA)
    sc.pl.umap(adata_pp18, color = gene, save = None, ax = ax[1,0], show=False, vmin = vmin, vmax = vmax, color_map = SUPER_MAGMA)
    sc.pl.umap(adata_ppc18, color = gene, save = None, ax = ax[1,1], show=False, vmin = vmin, vmax = vmax, color_map = SUPER_MAGMA)

    fig.suptitle(gene, fontsize = 20)
    ax[0,0].set_title("PP12")
    ax[0,1].set_title("PP18")
    ax[1,0].set_title("PPC12")
    ax[1,1].set_title("PPC18")
    fig.tight_layout()

    fig.savefig(f"results_seurat_scrnaseq/figure_markers/important_genes/{gene}.png", dpi = 200, bbox_inches = "tight")

with plt.rc_context({"figure.figsize": (7, 5), "figure.dpi": (300), "font.size": 15}):


    fig = sc.pl.umap(adata[(adata.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)") & (adata.obs["age"] == "12wk"),:], color = [x for x in gene_names if x != "Il13"], save = None, return_fig = True, ncols = 6, color_map = SUPER_MAGMA)
    fig.tight_layout()
    fig.savefig("results_seurat_scrnaseq/scrnaseq_pp12_visium_interest_genes.pdf", bbox_inches = "tight")

with plt.rc_context({"figure.figsize": (7, 5), "figure.dpi": (300), "font.size": 15}):
    fig = sc.pl.umap(adata[(adata.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)") & (adata.obs["age"] == "18wk"),:], color = gene_names, save =  None, return_fig = True, ncols = 6, color_map = SUPER_MAGMA)
    fig.tight_layout()
    fig.savefig("results_seurat_scrnaseq/scrnaseq_pp18_visium_interest_genes.pdf", bbox_inches = "tight")

with plt.rc_context({"figure.figsize": (7, 5), "figure.dpi": (300), "font.size": 15}):
    fig = sc.pl.umap(adata[(adata.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") & (adata.obs["age"] == "12wk"),:], color = gene_names, save =  None, return_fig = True, ncols = 6, color_map = SUPER_MAGMA)
    fig.tight_layout()
    fig.savefig("results_seurat_scrnaseq/scrnaseq_ppc12_visium_interest_genes.pdf", bbox_inches = "tight")

with plt.rc_context({"figure.figsize": (7, 5), "figure.dpi": (300), "font.size": 15}):
    fig = sc.pl.umap(adata[(adata.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") & (adata.obs["age"] == "18wk"),:], color = gene_names, save =  None, return_fig = True, ncols = 6, color_map = SUPER_MAGMA)
    fig.tight_layout()
    fig.savefig("results_seurat_scrnaseq/scrnaseq_ppc18_visium_interest_genes.pdf", bbox_inches = "tight")

# In[]
# ---------------------------------------------------------------- #
#
# Boxplots showing the expression change of genes across genotypes
#
# ---------------------------------------------------------------- #

interest_genes = pd.read_csv("spatial_data/interest_genes.csv", index_col = 0)
interest_genes2 = pd.read_excel("markers_info/Gene list for cell count in genotype.xlsx")
gene_names2 =[x for x in interest_genes2["Genes"] if (x!="Tnfa")]
gene_names = list(set(gene_names + gene_names2))
gene_names =[x for x in interest_genes["Genes"] if (x != "Cxcl11") & (x != "Tnfa") & (x != "Tgfb") & (x != "Il12") & (x != "Pd1") & (x != "Tim3") & (x != "Il13")]

for gene in gene_names:
    # intact tissue
    gene_expr = pd.DataFrame(columns = ["expr", "genotype", "annot", "age"])
    gene_expr["Expr"] = adata[:, gene].X.toarray().squeeze()
    gene_expr["genotype"] = adata.obs["genotype"].values
    gene_expr["age"] = adata.obs["age"].values
    gene_expr["annot"] = adata.obs["annot"].values

    # castrated tissue
    gene_expr_cas = pd.DataFrame(columns = ["expr", "genotype", "annot", "age"])
    gene_expr_cas["Expr"] = adata_castrated[:, gene].X.toarray().squeeze()
    gene_expr_cas["genotype"] = adata_castrated.obs["genotype"].values
    gene_expr_cas["annot"] = adata_castrated.obs["annot_transfer"].values

    fig = plt.figure(figsize = (30,10))
    ax = fig.subplots(nrows = 3, ncols = 1)
    sns.boxplot(data = gene_expr[gene_expr["age"] == "12wk"], x = "annot", hue = "genotype", y = "Expr", ax = ax[0])
    leg = ax[0].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
    ax[0].set_xlabel(None)
    ax[0].set_title(gene + " (12wk)", fontsize = 20)
    sns.boxplot(data = gene_expr[gene_expr["age"] == "18wk"], x = "annot", hue = "genotype", y = "Expr", ax = ax[1])
    leg = ax[1].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
    ax[1].set_xlabel(None)
    ax[1].set_title(gene + " (18wk)", fontsize = 20)
    sns.boxplot(data = gene_expr_cas, x = "annot", hue = "genotype", y = "Expr", ax = ax[2])
    leg = ax[2].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
    ax[2].set_xlabel(None)
    ax[2].set_title(gene + " (castrated)", fontsize = 20)
    fig.tight_layout()
    fig.savefig(f"results_seurat_scrnaseq/DE_genotype/interest_genes_genotype_compare/castrated/{gene}_genotype_compare.png", bbox_inches = "tight")

# interest_genes2 = pd.read_excel("markers_info/Gene list for cell count in genotype.xlsx")
# gene_names2 =[x for x in interest_genes2["Genes"] if (x!="Tnfa")]
# for gene in gene_names2:
#     gene_expr = pd.DataFrame(columns = ["expr", "genotype", "annot", "age"])
#     x_genes = adata[:, gene].X.toarray().squeeze()
#     gene_expr["Expr"] = x_genes
#     gene_expr["genotype"] = adata.obs["genotype"].values
#     gene_expr["age"] = adata.obs["age"].values
#     gene_expr["annot"] = adata.obs["annot"].values

#     fig = plt.figure(figsize = (20,10))
#     ax = fig.subplots(nrows = 2, ncols = 1)
#     sns.boxplot(data = gene_expr[gene_expr["age"] == "12wk"], x = "annot", hue = "genotype", y = "Expr", ax = ax[0])
#     leg = ax[0].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
#     ax[0].set_xlabel(None)
#     ax[0].set_title(gene + " (12wk)", fontsize = 20)
#     sns.boxplot(data = gene_expr[gene_expr["age"] == "18wk"], x = "annot", hue = "genotype", y = "Expr", ax = ax[1])
#     leg = ax[1].legend(loc='upper left', prop={'size': 15}, frameon = False, bbox_to_anchor=(1.04, 1))
#     ax[1].set_xlabel(None)
#     ax[1].set_title(gene + " (18wk)", fontsize = 20)
#     fig.savefig(f"results_seurat_scrnaseq/DE_genotype/interest_genes_genotype_compare/{gene}_genotype_compare.png", bbox_inches = "tight")


# In[]
# ---------------------------------------------------------------- #
#
# Plot cell composition pie charts
#
# ---------------------------------------------------------------- #
adata_castrated = sc.read_h5ad("Cell_Ranger_output/adata_castrated_seurat.h5ad")
adata_castrated.obs["annot_transfer"] = adata_castrated.obs["annot_transfer"].astype("str")
adata_castrated.obs.loc[adata_castrated.obs["annot_transfer"].isin(["Lymphoid", "Macrophage", "Monocytes"]),"annot_transfer"] = "Immune"
adata_ctr_castrated = adata_castrated[(adata_castrated.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"),:]
adata_ko_castrated = adata_castrated[(adata_castrated.obs["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)"),:]

adata = sc.read_h5ad("Cell_Ranger_output/adata_seurat.h5ad")
adata.obs["annot"] = adata.obs["annot"].astype("str")
adata.obs.loc[adata.obs["annot"].isin(["Lymphoid", "Macrophage", "Monocytes"]),"annot"] = "Immune"
# immune cells count under 12wk
adata_ctr_12wk = adata[(adata.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata.obs["age"] == "12wk"),:]
adata_ko_12wk = adata[(adata.obs["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata.obs["age"] == "12wk"),:]
# immune cells count under 18wk
adata_ctr_18wk = adata[(adata.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata.obs["age"] == "18wk"),:]
adata_ko_18wk = adata[(adata.obs["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata.obs["age"] == "18wk"),:]


n_immune_ctr_12wk = np.sum(adata_ctr_12wk.obs["annot"] == "Immune")
n_other_ctr_12wk = adata_ctr_12wk.shape[0] - n_immune_ctr_12wk
n_immune_ko_12wk = np.sum(adata_ko_12wk.obs["annot"] == "Immune")
n_other_ko_12wk = adata_ko_12wk.shape[0] - n_immune_ko_12wk
n_immune_ctr_18wk = np.sum(adata_ctr_18wk.obs["annot"] == "Immune")
n_other_ctr_18wk = adata_ctr_18wk.shape[0] - n_immune_ctr_18wk
n_immune_ko_18wk = np.sum(adata_ko_18wk.obs["annot"] == "Immune")
n_other_ko_18wk = adata_ko_18wk.shape[0] - n_immune_ko_18wk
n_immune_ctr_castrated = np.sum(adata_ctr_castrated.obs["annot_transfer"] == "Immune")
n_other_ctr_castrated = np.sum(adata_ctr_castrated.obs["annot_transfer"] != "Immune")
n_immune_ko_castrated = np.sum(adata_ko_castrated.obs["annot_transfer"] == "Immune")
n_other_ko_castrated = np.sum(adata_ko_castrated.obs["annot_transfer"] != "Immune")


plt.rcParams["font.size"] = 15
fig = plt.figure(figsize = (12,18))
axs = fig.subplots(nrows = 3, ncols = 2)
axs[0,0].pie([n_immune_ctr_12wk, n_other_ctr_12wk], labels = ["Immune", "Other"], autopct='%.0f%%')
axs[0,1].pie([n_immune_ko_12wk, n_other_ko_12wk], labels = ["Immune", "Other"], autopct='%.0f%%')
axs[1,0].pie([n_immune_ctr_18wk, n_other_ctr_18wk], labels = ["Immune", "Other"], autopct='%.0f%%')
axs[1,1].pie([n_immune_ko_18wk, n_other_ko_18wk], labels = ["Immune", "Other"], autopct='%.0f%%')
axs[2,0].pie([n_immune_ctr_castrated, n_other_ctr_castrated], labels = ["Immune", "Other"], autopct='%.0f%%')
axs[2,1].pie([n_immune_ko_castrated, n_other_ko_castrated], labels = ["Immune", "Other"], autopct='%.0f%%')
axs[0,0].set_title("CTRL (12wk)")
axs[0,1].set_title("CXCR7 KO (12wk)")
axs[1,0].set_title("CTRL (18wk)")
axs[1,1].set_title("CXCR7 KO (18wk)")
axs[2,0].set_title("CTRL (castrated)")
axs[2,1].set_title("CXCR7 KO (castrated)")

fig.savefig("results_seurat_scrnaseq/immune_percentage.png", bbox_inches = "tight")

# In[]
# ---------------------------------------------------------------- #
#
# Plot cell composition table
#
# ---------------------------------------------------------------- #
adata_castrated = sc.read_h5ad("Cell_Ranger_output/adata_castrated_seurat.h5ad")
adata_castrated.obs["annot_transfer"] = adata_castrated.obs["annot_transfer"].astype("str")
adata_castrated.obs.loc[adata_castrated.obs["annot_transfer"].isin(["Lymphoid", "Macrophage", "Monocytes"]),"annot_transfer"] = "Immune"
adata_ctr_castrated = adata_castrated[(adata_castrated.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"),:]
adata_ko_castrated = adata_castrated[(adata_castrated.obs["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)"),:]

adata = sc.read_h5ad("Cell_Ranger_output/adata_seurat.h5ad")
adata.obs["annot"] = adata.obs["annot"].astype("str")
adata.obs.loc[adata.obs["annot"].isin(["Lymphoid", "Macrophage", "Monocytes"]),"annot"] = "Immune"
# immune cells count under 12wk
adata_ctr_12wk = adata[(adata.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata.obs["age"] == "12wk"),:]
adata_ko_12wk = adata[(adata.obs["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata.obs["age"] == "12wk"),:]
# immune cells count under 18wk
adata_ctr_18wk = adata[(adata.obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata.obs["age"] == "18wk"),:]
adata_ko_18wk = adata[(adata.obs["genotype"] != "PbCre(+/-),Pten(-/-),P53(-/-)")&(adata.obs["age"] == "18wk"),:]

adata_ctr_castrated.obs["annot_transfer"].value_counts()
adata_ko_castrated.obs["annot_transfer"].value_counts()
adata_ctr_12wk.obs["annot_transfer"].value_counts()
adata_ko_12wk.obs["annot_transfer"].value_counts()
adata_ctr_18wk.obs["annot_transfer"].value_counts()
adata_ko_18wk.obs["annot_transfer"].value_counts()


# %%
