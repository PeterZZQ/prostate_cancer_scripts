# In[]
import os,csv,re, time
import pickle
import random
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
from scipy import stats
from scipy.sparse import issparse
import scanpy as sc
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import cv2
import TESLA as tesla
from IPython.display import Image


def coord2pixel(coords, scale_factor):
    """
    coords: the spatial coordinates of cells in visium data
    scale_factor: the scale factors of cells in visium data
    """
    return coords * scale_factor


# In[]
#Read in gene expression and spatial location
adata_pp12 = sc.read_visium(path = "../dataset/data_visium/spatial_data/Visium_python/225_PP12/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_pp18 = sc.read_visium(path = "../dataset/data_visium/spatial_data/Visium_python/1687_PP18/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_ppc12 = sc.read_visium(path = "../dataset/data_visium/spatial_data/Visium_python/1161_PPC/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_ppc18 = sc.read_visium(path = "../dataset/data_visium/spatial_data/Visium_python/1660_PPC18/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")

#Read in hitology image, 
# Pixels in two ways are the same, rgb color scales are different
img_pp12 = cv2.imread("../dataset/data_visium/spatial_data/Visium_python/225_PP12/spatial/tissue_hires_image.png")
img_pp18 = cv2.imread("../dataset/data_visium/spatial_data/Visium_python/1687_PP18/spatial/tissue_hires_image.png")
img_ppc12 = cv2.imread("../dataset/data_visium/spatial_data/Visium_python/1161_PPC/spatial/tissue_hires_image.png")
img_ppc18 = cv2.imread("../dataset/data_visium/spatial_data/Visium_python/1660_PPC18/spatial/tissue_hires_image.png")
# values: rgb color, index: pixel
# img_pp12 = adata_pp12.uns["spatial"]["255"]["images"]["hires"]
# img_pp18 = adata_pp18.uns["spatial"]["1687PP"]["images"]["hires"]
# img_ppc12 = adata_ppc12.uns["spatial"]["1161"]["images"]["hires"]
# img_ppc18 = adata_ppc18.uns["spatial"]["1660PPC"]["images"]["hires"]

# scale the cell spatial coordinates into pixel
scale_factor = adata_pp12.uns["spatial"]["255"]["scalefactors"]["tissue_hires_scalef"]
coords = adata_pp12.obsm["spatial"].astype(np.double)
adata_pp12.obs[["pixel_x", "pixel_y"]] = coord2pixel(coords, scale_factor)

scale_factor = adata_pp18.uns["spatial"]["1687PP"]["scalefactors"]["tissue_hires_scalef"]
coords = adata_pp18.obsm["spatial"].astype(np.double)
adata_pp18.obs[["pixel_x", "pixel_y"]] = coord2pixel(coords, scale_factor)

scale_factor = adata_ppc12.uns["spatial"]["1161"]["scalefactors"]["tissue_hires_scalef"]
coords = adata_ppc12.obsm["spatial"].astype(np.double)
adata_ppc12.obs[["pixel_x", "pixel_y"]] = coord2pixel(coords, scale_factor)

scale_factor = adata_ppc18.uns["spatial"]["1660PPC"]["scalefactors"]["tissue_hires_scalef"]
coords = adata_ppc18.obsm["spatial"].astype(np.double)
adata_ppc18.obs[["pixel_x", "pixel_y"]] = coord2pixel(coords, scale_factor)

# # Test
# fig = plt.figure(figsize = (10, 10))
# ax = fig.add_subplot()
# ax.scatter(adata_pp12.obs["pixel_x"].values, adata_pp12.obs["pixel_y"].values)
# ax.imshow(img_pp12)

adata_dicts = {"pp12": adata_pp12, "pp18": adata_pp18, "ppc12": adata_ppc12, "ppc18": adata_ppc18}
img_dicts = {"pp12": img_pp12, "pp18": img_pp18, "ppc12": img_ppc12, "ppc18": img_ppc18}

# In[]
# Preprocessing
# 1. resize image
sample = "ppc18"
img = img_dicts[sample]
counts = adata_dicts[sample]

# make directory
result_dir = f"../results_visium/results_tesla/{sample}/"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

resize_factor=1000/np.min(img.shape[0:2])
resize_width=int(img.shape[1]*resize_factor)
resize_height=int(img.shape[0]*resize_factor)

# capitalize the gene name?
counts.var.index=[i.upper() for i in counts.var.index]
counts.var_names_make_unique()
counts.raw=counts
# impute on log scale
sc.pp.log1p(counts) 
if issparse(counts.X):counts.X=counts.X.A.copy()

# In[]
# Contour detection
cnt=tesla.cv2_detect_contour(img, apertureSize=5,L2gradient = True)

binary=np.zeros((img.shape[0:2]), dtype=np.uint8)
cv2.drawContours(binary, [cnt], -1, (1), thickness=-1)
#Enlarged filter
cnt_enlarged = tesla.scale_contour(cnt, 1.05)
binary_enlarged = np.zeros(img.shape[0:2])
cv2.drawContours(binary_enlarged, [cnt_enlarged], -1, (1), thickness=-1)
img_new = img.copy()
cv2.drawContours(img_new, [cnt], -1, (255), thickness=50)
img_new=cv2.resize(img_new, ((resize_width, resize_height)))

cv2.imwrite(result_dir + 'cnt.jpg', img_new)
Image(filename=result_dir + 'cnt.jpg')


# In[]
# Main function: Gene expression enhancement
#Set size of superpixel
res=10
# Note, if the numer of superpixels is too large and take too long, you can increase the res to 100
enhanced_exp_adata=tesla.imputation(img=img, raw=counts, cnt=cnt, genes=counts.var.index.tolist(), shape="None", res=res, s=1, k=2, num_nbs=10)
enhanced_exp_adata.write_h5ad(result_dir + "enhanced_exp.h5ad")


# Sanity check: plot the imputed gene expression data
cnt_color = clr.LinearSegmentedColormap.from_list('magma', ["#000003",  "#3b0f6f",  "#8c2980",   "#f66e5b", "#fd9f6c", "#fbfcbf"], N=256)
g="AR"
enhanced_exp_adata.obs[g]=enhanced_exp_adata.X[:,enhanced_exp_adata.var.index==g]
fig=sc.pl.scatter(enhanced_exp_adata,alpha=1,x="y",y="x",color=g,color_map=cnt_color,show=False,size=10)
fig.set_aspect('equal', 'box')
fig.invert_yaxis()
plt.gcf().set_dpi(600)
fig.figure.show()

# In[]
# TESLA cell type annotation
markers = {}
markers["luminal"] = ["Ar", "Krt8"] #, "Cd24a", "Krt18", "Spink1"]
markers["basal"] = ["Trp63", "Krt5"] #, "Krt14"]
# markers["club_epithelia"] = ["Agr2", "Krt7"]
# markers["endothelial"] = ["Ackr1", "Cldn5"]
# markers["lymphoid"] = ["Cd3e", "Ms4a1", "Klrb1c"]
# markers["myeloid"] = ["Ptprc", "Itgam"]
# markers["monocytes"] = ["Ptprc", "Itgam", "Cd14", "S100a8", "S100a9"] 
# markers["macrophage"] = ["Ptprc", "Itgam", "Adgre1"]
# markers["macrophage_m1"] = ["Ptprc", "Itgam", "Adgre1", "Cd68", "Nos2"]
# markers["macrophage_m2"] = ["Ptprc", "Itgam", "Adgre1", "Mrc1", "Arg1"]
# markers["mesenchymal"] = ["Fgf10", "Rorb", "Rspo1", "Sult1e1", "Wnt10a", "Wnt2", "Wnt6"]
# markers["sv"] = ["Pate4", "Pax2", "Svs2"]

# load the imputed data to prevent calculation [time consuming]
enhanced_exp_adata = sc.read_h5ad(result_dir + "enhanced_exp.h5ad")

# make capital, following before
for ct in markers.keys():
    markers[ct] = [i.upper() for i in markers[ct]]
    genes=list(set([i.upper() for i in markers[ct] if i.upper() in enhanced_exp_adata.var.index ]))


# In[]
for ct in markers.keys():
    #target_size can be set to "small" or "large".
    # num_required: the number of marker that should be used to identify the cell type
    # all the markers should be used, default was 1
    pred_refined, target_clusters, c_m=tesla.annotation(img=img, 
                                                        binary=binary,
                                                        sudo_adata=enhanced_exp_adata, 
                                                        genes=markers[ct], 
                                                        resize_factor=resize_factor,
                                                        # some bugs? why minus 1?
                                                        num_required=len(markers[ct]) - 1, 
                                                        target_size="small")
    #Plot
    ret_img=tesla.visualize_annotation(img=img, 
                                binary=binary, 
                                resize_factor=resize_factor,
                                pred_refined=pred_refined, 
                                target_clusters=target_clusters, 
                                c_m=c_m)

    cv2.imwrite(result_dir + ct + ".png", ret_img)
    Image(filename=result_dir + ct + ".png")

    #Save
    np.save(result_dir + ct + "_annotation.npy", pred_refined)
    print("Target_clusters: ", target_clusters, "\n")
    #Save the cluster density information
    c_d={i[0]:i[1] for i in c_m[0:len(target_clusters)]}
    print("Cluster_density : ", c_d)
    with open(result_dir + ct + '_annotation_c_d.pkl', 'wb') as f: pickle.dump(c_d, f)



# %%
