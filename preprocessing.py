# In[]
import numpy as np
import pandas as pd
from scipy.io import mmread
from anndata import AnnData
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import scanpy as sc

datadir_m1416_pp18 = "Cell_Ranger_output/M1416-PP18/filtered_feature_bc_matrix/"
datadir_m1417_pp12 = "Cell_Ranger_output/M1417-PP12/filtered_feature_bc_matrix/"
datadir_m1436_ppc12 = "Cell_Ranger_output/M1436-PPC12/filtered_feature_bc_matrix/"
datadir_m1437_ppc18 = "Cell_Ranger_output/M1437-PPC18/filtered_feature_bc_matrix/"



def plot_latent(zs, annos = None, batches = None, mode = "annos", save = None, figsize = (20,10), axis_label = "Latent", label_inplace = False, legend = True, **kwargs):
    """\
    Description
        Plot latent space
    Parameters
        z1
            the latent space of first data batch, of the shape (n_samples, n_dimensions)
        z2
            the latent space of the second data batch, of the shape (n_samples, n_dimensions)
        anno1
            the cluster annotation of the first data batch, of the  shape (n_samples,)
        anno2
            the cluster annotation of the second data batch, of the  shape (n_samples,)
        mode
            "joint": plot two latent spaces(from two batches) into one figure
            "separate" plot two latent spaces separately
        save
            file name for the figure
        figsize
            figure size
    """
    _kwargs = {
        "s": 10,
        "alpha": 0.7,
        "markerscale": 1,
        "text_size": "large",
        "colormap": None
    }
    _kwargs.update(kwargs)

    fig = plt.figure(figsize = figsize, dpi = 300, constrained_layout=True)
    if (mode == "annos") | (mode == "batches"):
        ax = fig.add_subplot()
        if mode == "annos":
            unique_cluster = np.unique(annos)
            if _kwargs["colormap"] is None:
                colormap = plt.cm.get_cmap("tab20", len(unique_cluster))
            else:
                colormap = _kwargs["colormap"]
        else:
            unique_cluster = np.unique(batches)
            if _kwargs["colormap"] is None:
                colormap = plt.cm.get_cmap("tab20", len(unique_cluster))
            else:
                colormap = _kwargs["colormap"]

        texts = []
        for i, cluster_type in enumerate(unique_cluster):
            if mode == "annos":
                index = np.where(annos == cluster_type)[0]
            else:
                index = np.where(batches == cluster_type)[0]
            z_clust = zs[index,:]
            ax.scatter(z_clust[:,0], z_clust[:,1], color = colormap(i), label = cluster_type, s = _kwargs["s"], alpha = _kwargs["alpha"])
            # text on plot
            if label_inplace:
                texts.append(ax.text(np.median(z_clust[:,0]), np.median(z_clust[:,1]), color = "black", s = unique_cluster[i], fontsize = _kwargs["text_size"], weight = 'semibold', in_layout = True))
        
        if legend:
            leg = ax.legend(loc='upper left', prop={'size': 15}, frameon = False, ncol = (len(unique_cluster) // 15) + 1, bbox_to_anchor=(1.04, 1), markerscale = _kwargs["markerscale"])
            for lh in leg.legendHandles: 
                lh.set_alpha(1)

        ax.tick_params(axis = "both", which = "major", labelsize = 15)

        ax.set_xlabel(axis_label + " 1", fontsize = 19)
        ax.set_ylabel(axis_label + " 2", fontsize = 19)
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
        ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)  
        # adjust position
        # if label_inplace:
        #     adjust_text(texts, only_move={'points':'xy', 'texts':'xy'})

    elif mode == "separate":
        unique_batch = np.unique(batches)
        unique_cluster = np.unique(annos) 
        axs = fig.subplots(len(unique_batch),1)
        if _kwargs["colormap"] is None:
            colormap = plt.cm.get_cmap("tab20", len(unique_cluster))
        else:
            colormap = _kwargs["colormap"]

        for i, batch in enumerate(unique_batch):
            zs_batch = zs[batches == batch]
            annos_batch = annos[batches == batch]
            z_clust = []
            texts = []
            for j, cluster_type in enumerate(unique_cluster):
                index = np.where(annos_batch == cluster_type)[0]
                if len(index) > 0:
                    axs[i].scatter(zs_batch[index,0], zs_batch[index,1], color = colormap(j), label = cluster_type, s = _kwargs["s"], alpha = _kwargs["alpha"])
                    # text on plot
                    if label_inplace:
                        # if exist cells
                        if zs_batch[index,0].shape[0] > 0:
                            texts.append(axs[i].text(np.median(zs_batch[index,0]), np.median(zs_batch[index,1]), color = "black", s = cluster_type, fontsize = _kwargs["text_size"], weight = 'semibold', in_layout = True))
            
            if legend:
                leg = axs[i].legend(loc='upper left', prop={'size': 15}, frameon = False, ncol = (len(unique_cluster) // 15) + 1, bbox_to_anchor=(1.04, 1), markerscale = _kwargs["markerscale"])
                for lh in leg.legendHandles: 
                    lh.set_alpha(1)

                
            axs[i].set_title(batch, fontsize = 25)

            axs[i].tick_params(axis = "both", which = "major", labelsize = 15)

            axs[i].set_xlabel(axis_label + " 1", fontsize = 19)
            axs[i].set_ylabel(axis_label + " 2", fontsize = 19)

            axs[i].spines['right'].set_visible(False)
            axs[i].spines['top'].set_visible(False)  

            axs[i].set_xlim(np.min(zs[:,0]), np.max(zs[:,0]))
            axs[i].set_ylim(np.min(zs[:,1]), np.max(zs[:,1]))
            axs[i].xaxis.set_major_locator(plt.MaxNLocator(4))
            axs[i].yaxis.set_major_locator(plt.MaxNLocator(4))
            axs[i].xaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
            axs[i].yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))

            # if label_inplace:
            #     adjust_text(texts, only_move={'points':'xy', 'texts':'xy'})        
            plt.tight_layout()
    if save:
        fig.savefig(save, bbox_inches = "tight")


# In[]
data_dir = datadir_m1437_ppc18
barcodes = pd.read_csv(data_dir + "barcodes.tsv", header = None, sep = "\t").values.squeeze()
features = pd.read_csv(data_dir + "features.tsv", header = None, sep = "\t")
matrix = mmread(data_dir + "matrix.mtx")

adata = AnnData(X = csr_matrix(matrix.T))
adata.obs.index = barcodes
adata.var.index = features[1].values.squeeze()
adata.var["gene_ids"] = features[0].values.squeeze()
adata.var_names_make_unique()

adata.write_h5ad(data_dir + "adata_raw.h5ad")

# In[]

adata_m1416_pp18 = sc.read_h5ad(datadir_m1416_pp18 + "adata_raw.h5ad")
adata_m1417_pp12 = sc.read_h5ad(datadir_m1417_pp12 + "adata_raw.h5ad")
adata_m1436_ppc12 = sc.read_h5ad(datadir_m1436_ppc12 + "adata_raw.h5ad")
adata_m1437_ppc18 = sc.read_h5ad(datadir_m1437_ppc18 + "adata_raw.h5ad")

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

adata_merge = AnnData.concatenate(adata_m1416_pp18, adata_m1417_pp12, adata_m1436_ppc12, adata_m1437_ppc18, join = "inner")

adata_merge.write_h5ad("Cell_Ranger_output/adata_merge.h5ad")


# In[]
# original adata file that includes the seurat cluster result and predicted label
adata = sc.read_h5ad("Cell_Ranger_output/adata_seurat.h5ad")
# sc.pp.filter_cells(adata, min_genes = 200)
sc.pp.filter_genes(adata, min_cells = 3)
# keep the original count before normalization and log transform
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
# store the normalized values in `.raw` before filtering
adata.raw = adata

sc.pp.neighbors(adata, n_neighbors = 30)
sc.tl.umap(adata)

# In[]
# sc.pl.umap(adata, color = ["sample"])
# sc.pl.umap(adata, color = ["age"])
# sc.pl.umap(adata, color = ["genotype"])
plot_latent(adata.obsm["X_umap"], annos = adata.obs["predicted_id"].squeeze(), batches = adata.obs["sample"], mode = "separate", figsize = (10,20), s = 1, markerscale = 6, axis_label = "UMAP", save = "transferred_label_separate.png")
plot_latent(adata.obsm["X_umap"], annos = adata.obs["predicted_id"].squeeze(), batches = adata.obs["sample"], mode = "annos", figsize = (10,5), s = 1, markerscale = 6, axis_label = "UMAP", save = "transferred_label.png")
plot_latent(adata.obsm["X_umap"], annos = np.array([eval(x) for x in adata.obs["seurat_clusters"].squeeze()]), batches = adata.obs["sample"], mode = "separate", figsize = (10,30), s = 1, markerscale = 6, axis_label = "UMAP", save = "seurat_clusters_separate.png", label_inplace = True)
plot_latent(adata.obsm["X_umap"], annos = np.array([eval(x) for x in adata.obs["seurat_clusters"].squeeze()]), batches = adata.obs["sample"], mode = "annos", figsize = (10,7), s = 1, markerscale = 6, axis_label = "UMAP", save = "seurat_clusters.png", label_inplace = True)

adata.obs["annotation"] = "Other"
adata.obs.loc[(adata.obs["seurat_clusters"] == "0") | (adata.obs["seurat_clusters"] == "1") | (adata.obs["seurat_clusters"] == "12") | (adata.obs["seurat_clusters"] == "13"), "annotation"] = "Luminal"
adata.obs.loc[adata.obs["seurat_clusters"] == "2", "annotation"] = "Basal"
adata.obs.loc[adata.obs["seurat_clusters"] == "3", "annotation"] = "Club epithelia (Luminal)"
adata.obs.loc[adata.obs["seurat_clusters"] == "4", "annotation"] = "Macrophages (Myeloid)"
adata.obs.loc[adata.obs["seurat_clusters"] == "5", "annotation"] = "Monocyte (Myeloid)"
adata.obs.loc[adata.obs["seurat_clusters"] == "6", "annotation"] = "Fibroblast (mesenchymal)"
adata.obs.loc[adata.obs["seurat_clusters"] == "7", "annotation"] = "Luminal (Spink1+)"
adata.obs.loc[adata.obs["seurat_clusters"] == "8", "annotation"] = "Basal"
adata.obs.loc[adata.obs["seurat_clusters"] == "9", "annotation"] = "Lymphoid"
adata.obs.loc[adata.obs["seurat_clusters"] == "10", "annotation"] = "SV (Spink1+)"
adata.obs.loc[adata.obs["seurat_clusters"] == "11", "annotation"] = "Endothelial"

plot_latent(adata.obsm["X_umap"], annos = adata.obs["annotation"].squeeze(), batches = adata.obs["sample"], mode = "separate", figsize = (10,20), s = 1, markerscale = 6, axis_label = "UMAP", save = "seurat_annotation.png")


# In[]
# original adata file that includes the seurat cluster result and predicted label
adata = sc.read_h5ad("Cell_Ranger_output/adata_seurat.h5ad")
# sc.pp.filter_cells(adata, min_genes = 200)
sc.pp.filter_genes(adata, min_cells = 3)
# keep the original count before normalization and log transform
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
# store the normalized values in `.raw` before filtering
adata.raw = adata

# separate by batch
adata_m1416_pp18 = adata[adata.obs["sample"] == "M1416-PP18",]
adata_m1417_pp12 = adata[adata.obs["sample"] == "M1417-PP12",]
adata_m1436_ppc12 = adata[adata.obs["sample"] == "M1436-PPC12",]
adata_m1437_ppc18 = adata[adata.obs["sample"] == "M1437-PPC18",]

sc.pp.neighbors(adata_m1416_pp18, n_neighbors = 30)
sc.tl.umap(adata_m1416_pp18)

sc.pp.neighbors(adata_m1417_pp12, n_neighbors = 30)
sc.tl.umap(adata_m1417_pp12)

sc.pp.neighbors(adata_m1436_ppc12, n_neighbors = 30)
sc.tl.umap(adata_m1436_ppc12)

sc.pp.neighbors(adata_m1437_ppc18, n_neighbors = 30)
sc.tl.umap(adata_m1437_ppc18)

# In[]

adata_merge = AnnData.concatenate(adata_m1416_pp18, adata_m1417_pp12, adata_m1436_ppc12, adata_m1437_ppc18, join = "inner")

plot_latent(adata_merge.obsm["X_umap"], annos = np.array([eval(x) for x in adata_merge.obs["seurat_clusters"].squeeze()]), batches = adata_merge.obs["sample"], mode = "separate", figsize = (10,30), s = 1, markerscale = 6, axis_label = "UMAP", save = "seurat_clusters_separate2.png", label_inplace = True)



# %%
