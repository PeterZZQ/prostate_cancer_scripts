
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# ColorMap
# Define the colormap from grey to white to purple
colors = [(0.7, 0.7, 0.7), (1, 1, 1), (1, 0, 0)]  # Grey to White to Purple
cmap_name = 'grey_to_white_to_purple'
cmap_expr = LinearSegmentedColormap.from_list(cmap_name, colors)
green_red = LinearSegmentedColormap.from_list('green_to_red', [(0, 0.7, 0), (1, 1, 1), (1, 0, 0), (0.5, 0, 0.5)], N=500)
SUPER_MAGMA = LinearSegmentedColormap.from_list('super_magma', colors=['#e0e0e0', '#dedede', '#fff68f', '#ffec8b', '#ffc125', '#ee7600', '#ee5c42', '#cd3278', '#c71585', '#68228b'], N=500)




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

