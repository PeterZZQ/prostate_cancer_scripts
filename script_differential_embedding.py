# In[]
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances
from anndata import AnnData
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from umap import UMAP
import matplotlib.pyplot as plt

import leidenalg

try:
    from leidenalg.VertexPartition import MutableVertexPartition
except ImportError:
    class MutableVertexPartition: pass
    MutableVertexPartition.__module__ = 'leidenalg.VertexPartition'

def get_igraph_from_adjacency(adjacency, directed=None):
    """Get igraph graph from adjacency matrix."""
    import igraph as ig
    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets]
    if isinstance(weights, np.matrix):
        weights = weights.A1
    g = ig.Graph(directed=directed)
    g.add_vertices(adjacency.shape[0])  # this adds adjacency.shape[0] vertices
    g.add_edges(list(zip(sources, targets)))
    try:
        g.es['weight'] = weights
    except:
        pass
    if g.vcount() != adjacency.shape[0]:
        print( 'Your adjacency matrix contained redundant nodes.' )
    return g


def cluster_cells_leiden(
    X, 
    resolution = 30.0,
    random_state = 0,
    n_iterations: int = -1,
    k_neighs = 30,
    sigma = 1,
    **partition_kwargs):

    from sklearn.neighbors import NearestNeighbors

    partition_kwargs = dict(partition_kwargs)

    neighbor = NearestNeighbors(n_neighbors=k_neighs)
    neighbor.fit(X)
    # get test connectivity result 0-1 adj_matrix, mode = 'connectivity' by default
    adj_matrix = neighbor.kneighbors_graph(X).toarray()
    dist_matrix = neighbor.kneighbors_graph(X, mode='distance').toarray()

    adj_matrix += adj_matrix.T
    # change 2 back to 1
    adj_matrix[adj_matrix.nonzero()[0],adj_matrix.nonzero()[1]] = 1

    affin = np.exp(- (dist_matrix - np.min(dist_matrix, axis = 1)[:,None]) ** 2/sigma)
    affin = adj_matrix * affin
    affin = (affin + affin.T)/2
        
    partition_type = leidenalg.RBConfigurationVertexPartition
    g = get_igraph_from_adjacency(affin, directed = True)

    partition_kwargs['n_iterations'] = n_iterations
    partition_kwargs['seed'] = random_state
    partition_kwargs['resolution_parameter'] = resolution

    partition_kwargs['weights'] = np.array(g.es['weight']).astype(np.float64)
    part = leidenalg.find_partition(g, partition_type, **partition_kwargs)
    groups = np.array(part.membership)

    return groups

def plot_latent(zs, annos = None, mode = "discrete", save = None, figsize = (20,10), axis_label = "Latent", label_inplace = False, legend = True, title = None, **kwargs):
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
        "alpha": 0.9,
        "markerscale": 1,
        "text_size": "large",
        "colormap": None
    }
    _kwargs.update(kwargs)

    fig = plt.figure(figsize = figsize, dpi = 300, constrained_layout=True)

    ax = fig.add_subplot()
    if mode == "discrete":
        unique_cluster = np.unique(annos)
        if _kwargs["colormap"] is None:
            colormap = plt.cm.get_cmap("tab20b", len(unique_cluster))
        else:
            colormap = _kwargs["colormap"]

        texts = []
        for i, cluster_type in enumerate(unique_cluster):
            index = np.where(annos == cluster_type)[0]
            z_clust = zs[index,:]
            ax.scatter(z_clust[:,0], z_clust[:,1], color = colormap(i), label = cluster_type, s = _kwargs["s"], alpha = _kwargs["alpha"])
            # text on plot
            if label_inplace:
                texts.append(ax.text(np.median(z_clust[:,0]), np.median(z_clust[:,1]), color = "black", s = unique_cluster[i], fontsize = _kwargs["text_size"], weight = 'semibold', in_layout = True))
        
        if legend:
            leg = ax.legend(loc='upper left', prop={'size': 15}, frameon = False, ncol = (len(unique_cluster) // 15) + 1, bbox_to_anchor=(1.04, 1), markerscale = _kwargs["markerscale"])
            for lh in leg.legendHandles: 
                lh.set_alpha(1)
    
    else:
        points = ax.scatter(zs[:,0], zs[:,1], c = annos, s = _kwargs["s"], alpha = _kwargs["alpha"], cmap = plt.get_cmap("PuBu"))
        fig.colorbar(points)
    ax.tick_params(axis = "both", which = "major", labelsize = 15)
    if title is not None:
        ax.set_title(title, fontsize = 20)

    ax.set_xlabel(axis_label + " 1", fontsize = 19)
    ax.set_ylabel(axis_label + " 2", fontsize = 19)
    ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    ax.yaxis.set_major_locator(plt.MaxNLocator(4))
    ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False) 
    if save is not None:
        fig.savefig(save)

# In[]
seurat_pca = pd.read_csv("Cell_Ranger_output/seurat_pca.txt", sep = "\t")
meta = pd.read_csv("Cell_Ranger_output/seurat_meta.txt", sep = "\t")
adata_seurat = sc.read_h5ad("Cell_Ranger_output/adata_seurat.h5ad")

celltypes = meta["seurat_clusters"].values
celltypes = np.where(celltypes == 0, "Luminal 1", celltypes)
celltypes = np.where(celltypes == "1", "Luminal 1", celltypes)
celltypes = np.where(celltypes == "2", "Hillock epithelia (Basal)", celltypes)
celltypes = np.where(celltypes == "3", "Club epithelia (Luminal)", celltypes)
celltypes = np.where(celltypes == "4", "Macrophage (Myeloid, Total immune)", celltypes)
celltypes = np.where(celltypes == "5", "Monocytes (Myeloid, Total immune)", celltypes)
celltypes = np.where(celltypes == "6", "Mesenchymal", celltypes)
celltypes = np.where(celltypes == "7", "Spink1+ (Luminal)", celltypes)
celltypes = np.where(celltypes == "8", "Basal", celltypes)
celltypes = np.where(celltypes == "9", "Lymphoid (Total immune)", celltypes)
celltypes = np.where(celltypes == "10", "SV, Spink1+", celltypes)
celltypes = np.where(celltypes == "11", "Endothelial", celltypes)
celltypes = np.where(celltypes == "12", "Luminal 2", celltypes)
celltypes = np.where(celltypes == "13", "Luminal 3", celltypes)
meta["seurat_celltypes"] = celltypes

# In[]
sc.pp.filter_genes(adata_seurat, min_cells = 3)
sc.pp.normalize_per_cell(adata_seurat)
sc.pp.log1p(adata_seurat)
adata_seurat = adata_seurat[meta.index.values,:]
counts = adata_seurat.X.toarray()
# In[]
x_pca_ctrl_12wk = seurat_pca.loc[meta.index.values[(meta["age"] == "12wk") & (meta["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")], :]
counts_ctrl_12wk = counts[((meta["age"] == "12wk") & (meta["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")).values, :]
meta_ctrl_12wk = meta.loc[(meta["age"] == "12wk") & (meta["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"),:]

x_pca_stim_12wk = seurat_pca.loc[meta.index.values[(meta["age"] == "12wk") & (meta["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)")], :]
counts_stim_12wk = counts[((meta["age"] == "12wk") & (meta["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)")).values, :]
meta_stim_12wk = meta.loc[(meta["age"] == "12wk") & (meta["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"),:]

x_pca_ctrl_18wk = seurat_pca.loc[meta.index.values[(meta["age"] == "18wk") & (meta["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")], :]
counts_ctrl_18wk = counts[((meta["age"] == "18wk") & (meta["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")).values, :]
meta_ctrl_18wk = meta.loc[(meta["age"] == "18wk") & (meta["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"),:]


x_pca_stim_18wk = seurat_pca.loc[meta.index.values[(meta["age"] == "18wk") & (meta["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)")], :]
counts_stim_18wk = counts[((meta["age"] == "18wk") & (meta["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)")).values, :]
meta_stim_18wk = meta.loc[(meta["age"] == "18wk") & (meta["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"),:]

# In[]
# find knn across conditions

# calculate pairwise distance
# 1. 12wk, ctrl by stim
pdist_12wk = pairwise_distances(X = x_pca_ctrl_12wk.values, Y = x_pca_stim_12wk.values)
# 2. 18wk, ctrl by stim
pdist_18wk = pairwise_distances(X = x_pca_ctrl_18wk.values, Y = x_pca_stim_18wk.values)

# find knn
K = 20
knn_index = np.argpartition(pdist_12wk, kth = 20-1, axis = 1)[:, (K-1)]
kth_dist = np.take_along_axis(pdist_12wk, knn_index[:,None], axis = 1)


x_diffs_12wk = []
counts_stim_12wk_mapped = []
for celli in range(pdist_12wk.shape[0]):
    idx = pdist_12wk[celli,:] <= kth_dist[celli,0]
    x_perturb = np.mean(counts_stim_12wk[idx,:], axis = 0)
    x_orig = counts_ctrl_12wk[[celli],:]
    x_diff = x_perturb - x_orig
    x_diffs_12wk.append(x_diff)
    counts_stim_12wk_mapped.append(x_perturb[None,:])

x_diffs_12wk = np.concatenate(x_diffs_12wk, axis = 0)
counts_stim_12wk_mapped = np.concatenate(counts_stim_12wk_mapped, axis = 0)

K = 20
knn_index = np.argpartition(pdist_18wk, kth = 20-1, axis = 1)[:, (K-1)]
kth_dist = np.take_along_axis(pdist_18wk, knn_index[:,None], axis = 1)

x_diffs_18wk = []
counts_stim_18wk_mapped = []
for celli in range(pdist_18wk.shape[0]):
    idx = pdist_18wk[celli,:] <= kth_dist[celli,0]
    x_perturb = np.mean(counts_stim_18wk[idx,:], axis = 0)
    x_orig = counts_ctrl_18wk[[celli],:]
    x_diff = x_perturb - x_orig
    x_diffs_18wk.append(x_diff)
    counts_stim_18wk_mapped.append(x_perturb[None,:])

x_diffs_18wk = np.concatenate(x_diffs_18wk, axis = 0)
counts_stim_18wk_mapped = np.concatenate(counts_stim_18wk_mapped, axis = 0)

# In[]
x_pca_diffs_12wk = PCA(n_components = 30).fit_transform(StandardScaler().fit_transform(x_diffs_12wk))
x_umap_diffs_12wk = UMAP().fit_transform(x_pca_diffs_12wk)

x_pca_diffs_18wk = PCA(n_components = 30).fit_transform(StandardScaler().fit_transform(x_diffs_18wk))
x_umap_diffs_18wk = UMAP().fit_transform(x_pca_diffs_18wk)



# In[]
# mean perturbed expression
mean_diff_12wk = np.mean(x_diffs_12wk, axis = 0)
mean_diff_18wk = np.mean(x_diffs_18wk, axis = 0)
# variation of perturbation
var_diff_12wk = np.mean((x_diffs_12wk - mean_diff_12wk) ** 2, axis = 0) 
var_diff_18wk = np.mean((x_diffs_18wk - mean_diff_18wk) ** 2, axis = 0) 

mean_diff_12wk = pd.DataFrame(data = np.abs(mean_diff_12wk), index = adata_seurat.var.index.values)
mean_diff_18wk = pd.DataFrame(data = np.abs(mean_diff_18wk), index = adata_seurat.var.index.values)
var_diff_12wk = pd.DataFrame(data = var_diff_12wk, index = adata_seurat.var.index.values)
var_diff_18wk = pd.DataFrame(data = var_diff_18wk, index = adata_seurat.var.index.values)

# mean_diff_12wk = mean_diff_12wk.loc[mean_diff_12wk.values > 1e-4,:]
# mean_diff_18wk = mean_diff_18wk.loc[mean_diff_18wk.values > 1e-4,:]
# var_diff_12wk = var_diff_12wk.loc[var_diff_12wk.values > 1e-4,:]
# var_diff_18wk = var_diff_18wk.loc[var_diff_18wk.values > 1e-4,:]
mean_diff_12wk = mean_diff_12wk.sort_values(by = [0])[:-2001:-1]
mean_diff_18wk = mean_diff_18wk.sort_values(by = [0])[:-2001:-1]
var_diff_12wk = var_diff_12wk.sort_values(by = [0])[:-2001:-1]
var_diff_18wk = var_diff_18wk.sort_values(by = [0])[:-2001:-1]

mean_diff = pd.concat([mean_diff_12wk, mean_diff_18wk], axis = 0)
mean_diff["time"] = ["12wk"] * mean_diff_12wk.shape[0] + ["18wk"] * mean_diff_18wk.shape[0]
mean_diff.columns = ["Mean", "time"]
var_diff = pd.concat([var_diff_12wk, var_diff_18wk], axis = 0)
var_diff["time"] = ["12wk"] * mean_diff_12wk.shape[0] + ["18wk"] * mean_diff_18wk.shape[0]
var_diff.columns = ["Variance", "time"]

# In[]
# # important genes
# important_genes = ["Trp53", "Pten", "Ackr3", "Spink1"]
# mean_diff_important = mean_diff.loc[important_genes, :]

import seaborn as sns
from adjustText import adjust_text
fig = plt.figure(figsize = (10, 7))
ax = fig.add_subplot()
sns.violinplot(mean_diff, x = "time", y = "Mean", ax = ax)
top_mean_diff = pd.concat([mean_diff.iloc[:10,:], mean_diff.iloc[2000:2010,:]])
g = sns.stripplot(x="time", y = "Mean", data = top_mean_diff, color="red", edgecolor="gray", size = 6)
texts_red = []
for i in range(top_mean_diff.shape[0]):
    if i <= 20:
        texts_red.append(g.text(x=top_mean_diff["time"].values[i], y=top_mean_diff["Mean"].values[i]+0.001, s=top_mean_diff.index.values[i], horizontalalignment='right', size='medium', color='black', fontsize = 10))
adjust_text(texts_red, only_move={'points':'x', 'texts':'x'})
fig.savefig("mean_perturbed.png", bbox_inches = "tight")


fig = plt.figure(figsize = (10, 7))
ax = fig.add_subplot()
sns.violinplot(var_diff, x = "time", y = "Variance", ax = ax)
top_var_diff = pd.concat([var_diff.iloc[:10,:], var_diff.iloc[2000:2010,:]])
g = sns.stripplot(x="time", y = "Variance", data = top_var_diff, color="red", edgecolor="gray", size = 6)
texts_red = []
for i in range(top_var_diff.shape[0]):
    if i <= 20:
        texts_red.append(g.text(x=top_var_diff["time"].values[i], y=top_var_diff["Variance"].values[i]+0.001, s=top_var_diff.index.values[i], horizontalalignment='right', size='medium', color='black', fontsize = 10))
adjust_text(texts_red, only_move={'points':'x', 'texts':'x'})
fig.savefig("var_perturbed.png", bbox_inches = "tight")




# In[]
plot_latent(x_umap_diffs_12wk, annos = meta_ctrl_12wk["seurat_celltypes"].values, s = 12, markerscale = 4, save = "expr_change_12wk.png",figsize = (15, 10))
plot_latent(x_umap_diffs_18wk, annos = meta_ctrl_18wk["seurat_celltypes"].values, s = 12, markerscale = 4, save = "expr_change_18wk.png",figsize = (15, 10))

# re-run clustering algorithm on the latent space, obtain cell clusters based on the perturbation effect
cluster_12wk = cluster_cells_leiden(x_pca_diffs_12wk, resolution = 0.15)
cluster_18wk = cluster_cells_leiden(x_pca_diffs_18wk, resolution = 0.15)

plot_latent(x_umap_diffs_12wk, annos = cluster_12wk, mode = "discrete", s = 12, markerscale = 4, alpha = 0.5)
plot_latent(x_umap_diffs_18wk, annos = cluster_18wk, mode = "discrete", s = 12, markerscale = 4, alpha = 0.5)
plot_latent(x_umap_diffs_12wk, annos = np.exp(x_diffs_12wk[:, adata_seurat.var.index.values == "Spink1"].squeeze())-1, mode = "continuous", s = 12, markerscale = 4, alpha = 0.5, title = "Spink1")
plot_latent(x_umap_diffs_18wk, annos = np.exp(x_diffs_18wk[:, adata_seurat.var.index.values == "Spink1"].squeeze())-1, mode = "continuous", s = 12, markerscale = 4, alpha = 0.5, title = "Spink1")
plot_latent(x_umap_diffs_12wk, annos = np.exp(x_diffs_12wk[:, adata_seurat.var.index.values == "Pten"].squeeze())-1, mode = "continuous", s = 12, markerscale = 4, alpha = 0.5, title = "Pten")
plot_latent(x_umap_diffs_18wk, annos = np.exp(x_diffs_18wk[:, adata_seurat.var.index.values == "Pten"].squeeze())-1, mode = "continuous", s = 12, markerscale = 4, alpha = 0.5, title = "Pten")

plot_latent(x_umap_diffs_12wk, annos = np.exp(x_diffs_12wk[:, adata_seurat.var.index.values == "Hspb1"].squeeze())-1, mode = "continuous", s = 12, markerscale = 4, alpha = 0.5, title = "Hspb1")
plot_latent(x_umap_diffs_18wk, annos = np.exp(x_diffs_18wk[:, adata_seurat.var.index.values == "Hspb1"].squeeze())-1, mode = "continuous", s = 12, markerscale = 4, alpha = 0.5, title = "Hspb1")
plot_latent(x_umap_diffs_12wk, annos = np.exp(x_diffs_12wk[:, adata_seurat.var.index.values == "Hspa1a"].squeeze())-1, mode = "continuous", s = 12, markerscale = 4, alpha = 0.5, title = "Hspa1a")
plot_latent(x_umap_diffs_18wk, annos = np.exp(x_diffs_18wk[:, adata_seurat.var.index.values == "Hspa1a"].squeeze())-1, mode = "continuous", s = 12, markerscale = 4, alpha = 0.5, title = "Hspa1a")

# In[]
# adata_diffs_12wk = AnnData(X = x_diffs_12wk, obs = meta_ctrl_12wk)
# adata_diffs_18wk = AnnData(X = x_diffs_18wk, obs = meta_ctrl_18wk)

# sc.pp.neighbors(adata_diffs_12wk)
# sc.tl.leiden(adata_diffs_12wk, resolution = 0.1)
# sc.tl.umap(adata_diffs_12wk)
# sc.pl.umap(adata_diffs_12wk, color = "seurat_celltypes")


# In[]
# Check the change of proliferation marker genes "Mki67" and "Pcna"
# x_diffs_12wk & meta_ctrl_12wk
plt.rcParams["font.size"] = 15
luminal_idx_12wk = meta_ctrl_12wk["seurat_celltypes"].isin(["Luminal 1"]).values.squeeze()
luminal_idx_18wk = meta_ctrl_18wk["seurat_celltypes"].isin(["Luminal 1"]).values.squeeze()
mki67_idx = np.where(adata_seurat.var.index == "Mki67")[0][0]
pcna_idx = np.where(adata_seurat.var.index == "Pcna")[0][0]

fig = plt.figure(figsize = (14,7))
ax = fig.subplots(nrows = 1, ncols = 2)
ax[0].scatter(counts_ctrl_12wk[luminal_idx_12wk, mki67_idx], counts_stim_12wk_mapped[luminal_idx_12wk, mki67_idx], c = "r", label = "Mki67 (12wk)", s = 2)
ax[0].scatter(counts_ctrl_18wk[luminal_idx_18wk, mki67_idx], counts_stim_18wk_mapped[luminal_idx_18wk, mki67_idx], c = "b", label = "Mki67 (18wk)", s = 2)
ax[0].plot(np.arange(0, 4, 0.1), np.arange(0, 4, 0.1), "g--")
ax[0].set_xlim([0, 4])
ax[0].set_ylim([0, 4])
ax[0].legend()
ax[0].set_xlabel("PbCre(+/-),Pten(-/-),P53(-/-)")
ax[0].set_ylabel("PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)")
ax[0].set_title("Mki67")

ax[1].scatter(counts_ctrl_12wk[luminal_idx_12wk, pcna_idx], counts_stim_12wk_mapped[luminal_idx_12wk, pcna_idx], c = "r", label = "Pcna (12wk)", s = 2)
ax[1].scatter(counts_ctrl_18wk[luminal_idx_18wk, pcna_idx], counts_stim_18wk_mapped[luminal_idx_18wk, pcna_idx], c = "b", label = "Pcna (18wk)", s = 2)
ax[1].plot(np.arange(0, 4, 0.1), np.arange(0, 4, 0.1), "g--")
ax[1].set_xlim([0, 4])
ax[1].set_ylim([0, 4])
ax[1].legend()
ax[1].set_xlabel("PbCre(+/-),Pten(-/-),P53(-/-)")
ax[1].set_ylabel("PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)")
ax[1].set_title("Pcna")
plt.suptitle("Log-normalized counts (proliferation marker)")
fig.savefig("proliferation_marker_luminal.png", bbox_inches = "tight")

plt.rcParams["font.size"] = 15
luminal_idx_12wk = meta_ctrl_12wk["seurat_celltypes"].isin(["Basal"]).values.squeeze()
luminal_idx_18wk = meta_ctrl_18wk["seurat_celltypes"].isin(["Basal"]).values.squeeze()
mki67_idx = np.where(adata_seurat.var.index == "Mki67")[0][0]
pcna_idx = np.where(adata_seurat.var.index == "Pcna")[0][0]

fig = plt.figure(figsize = (14,7))
ax = fig.subplots(nrows = 1, ncols = 2)
ax[0].scatter(counts_ctrl_12wk[luminal_idx_12wk, mki67_idx], counts_stim_12wk_mapped[luminal_idx_12wk, mki67_idx], c = "r", label = "Mki67 (12wk)", s = 2)
ax[0].scatter(counts_ctrl_18wk[luminal_idx_18wk, mki67_idx], counts_stim_18wk_mapped[luminal_idx_18wk, mki67_idx], c = "b", label = "Mki67 (18wk)", s = 2)
ax[0].plot(np.arange(0, 4, 0.1), np.arange(0, 4, 0.1), "g--")
ax[0].set_xlim([0, 4])
ax[0].set_ylim([0, 4])
ax[0].legend()
ax[0].set_xlabel("PbCre(+/-),Pten(-/-),P53(-/-)")
ax[0].set_ylabel("PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)")
ax[0].set_title("Mki67")

ax[1].scatter(counts_ctrl_12wk[luminal_idx_12wk, pcna_idx], counts_stim_12wk_mapped[luminal_idx_12wk, pcna_idx], c = "r", label = "Pcna (12wk)", s = 2)
ax[1].scatter(counts_ctrl_18wk[luminal_idx_18wk, pcna_idx], counts_stim_18wk_mapped[luminal_idx_18wk, pcna_idx], c = "b", label = "Pcna (18wk)", s = 2)
ax[1].plot(np.arange(0, 4, 0.1), np.arange(0, 4, 0.1), "g--")
ax[1].set_xlim([0, 4])
ax[1].set_ylim([0, 4])
ax[1].legend()
ax[1].set_xlabel("PbCre(+/-),Pten(-/-),P53(-/-)")
ax[1].set_ylabel("PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)")
ax[1].set_title("Pcna")
plt.suptitle("Log-normalized counts (proliferation marker)")
fig.savefig("proliferation_marker_basal.png", bbox_inches = "tight")
# %%
