# ---------------------------------------------------------------------------------------------- # 
# 
# Spatial analysis provided by SquidPy:
# 1. Neighborhood analysis, find the tendency of neighboring cell types for each cell type
# 2. Ligand-receptor interaction: used cellPhoneDB
# 3. Spatially variable genes (used: Moran's I, other: SPARK, Spatial DE, trendsceek, HMRF)
# 
# ---------------------------------------------------------------------------------------------- # 


# In[]
import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc
import squidpy as sq

sc.logging.print_header()
print(f"squidpy=={sq.__version__}")


# In[]
# ---------------------------------------------------------------------
#
# Load Visium data
#
# ---------------------------------------------------------------------

adata_pp12 = sc.read_visium(path = "../dataset/data_visium/spatial_data/Visium_python/225_PP12/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_pp18 = sc.read_visium(path = "../dataset/data_visium/spatial_data/Visium_python/1687_PP18/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_ppc12 = sc.read_visium(path = "../dataset/data_visium/spatial_data/Visium_python/1161_PPC/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_ppc18 = sc.read_visium(path = "../dataset/data_visium/spatial_data/Visium_python/1660_PPC18/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")


adata_pp12.obsm["spatial"] = adata_pp12.obsm["spatial"].astype(np.float64)
adata_pp18.obsm["spatial"] = adata_pp18.obsm["spatial"].astype(np.float64)
adata_ppc12.obsm["spatial"] = adata_ppc12.obsm["spatial"].astype(np.float64)
adata_ppc18.obsm["spatial"] = adata_ppc18.obsm["spatial"].astype(np.float64)

# adata_pp12.var_names_make_unique()
# adata_pp18.var_names_make_unique()
# adata_ppc12.var_names_make_unique()
# adata_ppc18.var_names_make_unique()

# rename gene names from SYMBOL to ENSEMBLE_ID
adata_pp12.var["SYMBOL"] = adata_pp12.var_names
adata_pp12.var.index = adata_pp12.var["gene_ids"].values
adata_pp18.var["SYMBOL"] = adata_pp18.var_names
adata_pp18.var.index = adata_pp18.var["gene_ids"].values
adata_ppc12.var["SYMBOL"] = adata_ppc12.var_names
adata_ppc12.var.index = adata_ppc12.var["gene_ids"].values
adata_ppc18.var["SYMBOL"] = adata_ppc18.var_names
adata_ppc18.var.index = adata_ppc18.var["gene_ids"].values


# In[]
# ---------------------------------------------------------------------
#
# Cell type assignment
#
# ---------------------------------------------------------------------
# Use CARD decomposition result and Cell2Loc result, use consistent
dir_card = "../results_visium/results_card/"
dir_cell2loc = "../results_visium/results_cell2loc/"

adata_pp12_cell2loc = sc.read_h5ad(dir_cell2loc + "cell2loc_map_pp12_10_20/sp.h5ad")
deconv_pp12_cell2loc = adata_pp12_cell2loc.obsm["means_cell_abundance_w_sf"]
ct_pp12_cell2loc = np.array([x.split("_")[-1] for x in adata_pp12_cell2loc.obsm["means_cell_abundance_w_sf"].columns])
adata_pp12.obs["max annot (cell2loc)"] = ct_pp12_cell2loc[np.argmax(deconv_pp12_cell2loc.values, axis = 1)]

deconv_pp12_card = pd.read_csv(dir_card + "deconv_score_pp12.csv", index_col = 0, sep = "\t")
adata_pp12.obs["max annot (CARD)"] = deconv_pp12_card.columns[np.argmax(deconv_pp12_card.values, axis = 1)]

adata_pp12.obs["max annot"] = [x if x == y else "other" for x,y in zip(adata_pp12.obs["max annot (CARD)"].values.squeeze(), adata_pp12.obs["max annot (cell2loc)"].values.squeeze())]

adata_pp12.obs["max annot (CARD)"] = adata_pp12.obs["max annot (CARD)"].astype('category')
adata_pp12.obs["max annot (cell2loc)"] = adata_pp12.obs["max annot (cell2loc)"].astype('category')
adata_pp12.obs["max annot"] = adata_pp12.obs["max annot"].astype('category')


adata_pp18_cell2loc = sc.read_h5ad(dir_cell2loc + "cell2loc_map_pp18_10_20/sp.h5ad")
deconv_pp18_cell2loc = adata_pp18_cell2loc.obsm["means_cell_abundance_w_sf"]
ct_pp18_cell2loc = np.array([x.split("_")[-1] for x in adata_pp18_cell2loc.obsm["means_cell_abundance_w_sf"].columns])
adata_pp18.obs["max annot (cell2loc)"] = ct_pp18_cell2loc[np.argmax(deconv_pp18_cell2loc.values, axis = 1)]

deconv_pp18_card = pd.read_csv(dir_card + "deconv_score_pp18.csv", index_col = 0, sep = "\t")
adata_pp18.obs["max annot (CARD)"] = deconv_pp18_card.columns[np.argmax(deconv_pp18_card.values, axis = 1)]

adata_pp18.obs["max annot"] = [x if x == y else "other" for x,y in zip(adata_pp18.obs["max annot (CARD)"].values.squeeze(), adata_pp18.obs["max annot (cell2loc)"].values.squeeze())]

adata_pp18.obs["max annot (CARD)"] = adata_pp18.obs["max annot (CARD)"].astype('category')
adata_pp18.obs["max annot (cell2loc)"] = adata_pp18.obs["max annot (cell2loc)"].astype('category')
adata_pp18.obs["max annot"] = adata_pp18.obs["max annot"].astype('category')


adata_ppc12_cell2loc = sc.read_h5ad(dir_cell2loc + "cell2loc_map_ppc12_10_20/sp.h5ad")
deconv_ppc12_cell2loc = adata_ppc12_cell2loc.obsm["means_cell_abundance_w_sf"]
ct_ppc12_cell2loc = np.array([x.split("_")[-1] for x in adata_ppc12_cell2loc.obsm["means_cell_abundance_w_sf"].columns])
adata_ppc12.obs["max annot (cell2loc)"] = ct_ppc12_cell2loc[np.argmax(deconv_ppc12_cell2loc.values, axis = 1)]

deconv_ppc12 = pd.read_csv(dir_card + "deconv_score_ppc12.csv", index_col = 0, sep = "\t")
adata_ppc12.obs["max annot (CARD)"] = deconv_ppc12.columns[np.argmax(deconv_ppc12.values, axis = 1)]

adata_ppc12.obs["max annot"] = [x if x == y else "other" for x,y in zip(adata_ppc12.obs["max annot (CARD)"].values.squeeze(), adata_ppc12.obs["max annot (cell2loc)"].values.squeeze())]

adata_ppc12.obs["max annot (CARD)"] = adata_ppc12.obs["max annot (CARD)"].astype('category')
adata_ppc12.obs["max annot (cell2loc)"] = adata_ppc12.obs["max annot (cell2loc)"].astype('category')
adata_ppc12.obs["max annot"] = adata_ppc12.obs["max annot"].astype('category')


adata_ppc18_cell2loc = sc.read_h5ad(dir_cell2loc + "cell2loc_map_ppc18_10_20/sp.h5ad")
deconv_ppc18_cell2loc = adata_ppc18_cell2loc.obsm["means_cell_abundance_w_sf"]
ct_ppc18_cell2loc = np.array([x.split("_")[-1] for x in adata_ppc18_cell2loc.obsm["means_cell_abundance_w_sf"].columns])
adata_ppc18.obs["max annot (cell2loc)"] = ct_ppc18_cell2loc[np.argmax(deconv_ppc18_cell2loc.values, axis = 1)]

deconv_ppc18 = pd.read_csv(dir_card + "deconv_score_ppc18.csv", index_col = 0, sep = "\t")
adata_ppc18.obs["max annot (CARD)"] = deconv_ppc18.columns[np.argmax(deconv_ppc18.values, axis = 1)]

adata_ppc18.obs["max annot"] = [x if x == y else "other" for x,y in zip(adata_ppc18.obs["max annot (CARD)"].values.squeeze(), adata_ppc18.obs["max annot (cell2loc)"].values.squeeze())]

adata_ppc18.obs["max annot (CARD)"] = adata_ppc18.obs["max annot (CARD)"].astype('category')
adata_ppc18.obs["max annot (cell2loc)"] = adata_ppc18.obs["max annot (cell2loc)"].astype('category')
adata_ppc18.obs["max annot"] = adata_ppc18.obs["max annot"].astype('category')



# In[]
# ---------------------------------------------------------------------
#
# Neighborhood Analysis
#
# ---------------------------------------------------------------------

# [Optional] remove the non-relevant cell types
# adata_pp12 = adata_pp12[(adata_pp12.obs["max annot"] != "SV") & (adata_pp12.obs["max annot"] != "other"),:]
# adata_pp18 = adata_pp18[(adata_pp18.obs["max annot"] != "SV") & (adata_pp18.obs["max annot"] != "other"),:]
# adata_ppc12 = adata_ppc12[(adata_ppc12.obs["max annot"] != "SV") & (adata_ppc12.obs["max annot"] != "other"),:]
# adata_ppc18 = adata_ppc18[(adata_ppc18.obs["max annot"] != "SV") & (adata_ppc18.obs["max annot"] != "other"),:]

# calculate neighborhood
sq.gr.spatial_neighbors(adata_pp12)
sq.gr.nhood_enrichment(adata_pp12, cluster_key="max annot")
sq.gr.spatial_neighbors(adata_pp18)
sq.gr.nhood_enrichment(adata_pp18, cluster_key="max annot")
sq.gr.spatial_neighbors(adata_ppc12)
sq.gr.nhood_enrichment(adata_ppc12, cluster_key="max annot")
sq.gr.spatial_neighbors(adata_ppc18)
sq.gr.nhood_enrichment(adata_ppc18, cluster_key="max annot")

unified_ct = np.unique(np.concatenate([adata_pp12.obs["max annot"].cat.categories.values, adata_pp18.obs["max annot"].cat.categories, adata_ppc12.obs["max annot"].cat.categories, adata_ppc18.obs["max annot"].cat.categories], axis = 0))

nbr_zscore_pp12 = pd.DataFrame(data = 0, index = unified_ct, columns = unified_ct)
nbr_zscore_pp18 = pd.DataFrame(data = 0, index = unified_ct, columns = unified_ct)
nbr_zscore_ppc12 = pd.DataFrame(data = 0, index = unified_ct, columns = unified_ct)
nbr_zscore_ppc18 = pd.DataFrame(data = 0, index = unified_ct, columns = unified_ct)
nbr_zscore_pp12.loc[adata_pp12.obs["max annot"].cat.categories, adata_pp12.obs["max annot"].cat.categories] = adata_pp12.uns["max annot_nhood_enrichment"]["zscore"]
nbr_zscore_pp18.loc[adata_pp18.obs["max annot"].cat.categories, adata_pp18.obs["max annot"].cat.categories] = adata_pp18.uns["max annot_nhood_enrichment"]["zscore"]
nbr_zscore_ppc12.loc[adata_ppc12.obs["max annot"].cat.categories, adata_ppc12.obs["max annot"].cat.categories] = adata_ppc12.uns["max annot_nhood_enrichment"]["zscore"]
nbr_zscore_ppc18.loc[adata_ppc18.obs["max annot"].cat.categories, adata_ppc18.obs["max annot"].cat.categories] = adata_ppc18.uns["max annot_nhood_enrichment"]["zscore"]

nbr_counts_pp12 = pd.DataFrame(data = 0, index = unified_ct, columns = unified_ct)
nbr_counts_pp18 = pd.DataFrame(data = 0, index = unified_ct, columns = unified_ct)
nbr_counts_ppc12 = pd.DataFrame(data = 0, index = unified_ct, columns = unified_ct)
nbr_counts_ppc18 = pd.DataFrame(data = 0, index = unified_ct, columns = unified_ct)
nbr_counts_pp12.loc[adata_pp12.obs["max annot"].cat.categories, adata_pp12.obs["max annot"].cat.categories] = adata_pp12.uns["max annot_nhood_enrichment"]["count"]
nbr_counts_pp18.loc[adata_pp18.obs["max annot"].cat.categories, adata_pp18.obs["max annot"].cat.categories] = adata_pp18.uns["max annot_nhood_enrichment"]["count"]
nbr_counts_ppc12.loc[adata_ppc12.obs["max annot"].cat.categories, adata_ppc12.obs["max annot"].cat.categories] = adata_ppc12.uns["max annot_nhood_enrichment"]["count"]
nbr_counts_ppc18.loc[adata_ppc18.obs["max annot"].cat.categories, adata_ppc18.obs["max annot"].cat.categories] = adata_ppc18.uns["max annot_nhood_enrichment"]["count"]


adata_pp12.uns["max annot_nhood_enrichment"]["zscore"] = nbr_zscore_pp12.values
adata_pp18.uns["max annot_nhood_enrichment"]["zscore"] = nbr_zscore_pp18.values
adata_ppc12.uns["max annot_nhood_enrichment"]["zscore"] = nbr_zscore_ppc12.values
adata_ppc18.uns["max annot_nhood_enrichment"]["zscore"] = nbr_zscore_ppc18.values

adata_pp12.uns["max annot_nhood_enrichment"]["count"] = nbr_counts_pp12.values
adata_pp18.uns["max annot_nhood_enrichment"]["count"] = nbr_counts_pp18.values
adata_ppc12.uns["max annot_nhood_enrichment"]["count"] = nbr_counts_ppc12.values
adata_ppc18.uns["max annot_nhood_enrichment"]["count"] = nbr_counts_ppc18.values

adata_pp12.obs["max annot"] = adata_pp12.obs["max annot"].cat.set_categories(unified_ct)
adata_pp18.obs["max annot"] = adata_pp18.obs["max annot"].cat.set_categories(unified_ct)
adata_ppc12.obs["max annot"] = adata_ppc12.obs["max annot"].cat.set_categories(unified_ct)
adata_ppc18.obs["max annot"] = adata_ppc18.obs["max annot"].cat.set_categories(unified_ct)


import matplotlib.pyplot as plt
plt.rcParams["font.size"] = 20
fig = plt.figure(figsize = (40, 30))
ax = fig.subplots(nrows = 2, ncols = 2)
sq.pl.nhood_enrichment(adata_pp12, cluster_key="max annot", annotate = True, ax = ax[0,0], title = "PP12", mode = "zscore")
sq.pl.nhood_enrichment(adata_pp18, cluster_key="max annot", annotate = True, ax = ax[0,1], title = "PP18", mode = "zscore")
sq.pl.nhood_enrichment(adata_ppc12, cluster_key="max annot", annotate = True, ax = ax[1,0], title = "PPC12", mode = "zscore")
sq.pl.nhood_enrichment(adata_ppc18, cluster_key="max annot", annotate = True, ax = ax[1,1], title = "PPC18", mode = "zscore")

fig.savefig("../results_visium/results_squidpy/neighborhood_analysis.png", bbox_inches = "tight", dpi = 200)

fig = plt.figure(figsize = (40, 30))
ax = fig.subplots(nrows = 2, ncols = 2)
sq.pl.nhood_enrichment(adata_pp12, cluster_key="max annot", annotate = True, ax = ax[0,0], title = "PP12", mode = "count")
sq.pl.nhood_enrichment(adata_pp18, cluster_key="max annot", annotate = True, ax = ax[0,1], title = "PP18", mode = "count")
sq.pl.nhood_enrichment(adata_ppc12, cluster_key="max annot", annotate = True, ax = ax[1,0], title = "PPC12", mode = "count")
sq.pl.nhood_enrichment(adata_ppc18, cluster_key="max annot", annotate = True, ax = ax[1,1], title = "PPC18", mode = "count")

fig.savefig("../results_visium/results_squidpy/neighborhood_analysis_count.png", bbox_inches = "tight", dpi = 200)


# In[]
# ---------------------------------------------------------------------
#
# Visualize annotations
#
# ---------------------------------------------------------------------
fig = plt.figure(figsize = (40,30))
ax = fig.subplots(nrows = 2, ncols = 2)
sc.pl.spatial(adata_pp12, img_key = "hires", color = "max annot", ax = ax[0,0])
sc.pl.spatial(adata_pp18, img_key = "hires", color = "max annot", ax = ax[0,1])
sc.pl.spatial(adata_ppc12, img_key = "hires", color = "max annot", ax = ax[1,0])
sc.pl.spatial(adata_ppc18, img_key = "hires", color = "max annot", ax = ax[1,1])

fig.savefig("../results_visium/results_squidpy/annot.png", bbox_inches = "tight", dpi = 300)

# In[]
# ---------------------------------------------------------------------
#
# Spatial DE analysis
#
# ---------------------------------------------------------------------
import NaiveDE
import SpatialDE

counts_pp12 = pd.DataFrame(data = adata_pp12.X.toarray(), columns = adata_pp12.var.index.values, index = adata_pp12.obs.index.values)
counts_pp18 = pd.DataFrame(data = adata_pp18.X.toarray(), columns = adata_pp18.var.index.values, index = adata_pp18.obs.index.values)
counts_ppc12 = pd.DataFrame(data = adata_ppc12.X.toarray(), columns = adata_ppc12.var.index.values, index = adata_ppc12.obs.index.values)
counts_ppc18 = pd.DataFrame(data = adata_ppc18.X.toarray(), columns = adata_ppc18.var.index.values, index = adata_ppc18.obs.index.values)

sample_info_pp12 = pd.DataFrame(data = adata_pp12.obsm["spatial"], index = adata_pp12.obs.index.values, columns = ["x", "y"])
sample_info_pp18 = pd.DataFrame(data = adata_pp18.obsm["spatial"], index = adata_pp18.obs.index.values, columns = ["x", "y"])
sample_info_ppc12 = pd.DataFrame(data = adata_ppc12.obsm["spatial"], index = adata_ppc12.obs.index.values, columns = ["x", "y"])
sample_info_ppc18 = pd.DataFrame(data = adata_ppc18.obsm["spatial"], index = adata_ppc18.obs.index.values, columns = ["x", "y"])

sample_info_pp12["total_counts"] = np.sum(counts_pp12, axis = 1)
sample_info_pp18["total_counts"] = np.sum(counts_pp18, axis = 1)
sample_info_ppc12["total_counts"] = np.sum(counts_ppc12, axis = 1)
sample_info_ppc18["total_counts"] = np.sum(counts_ppc18, axis = 1)


# SpatialDE
norm_expr_pp12 = NaiveDE.stabilize(counts_pp12.T).T
resid_expr_pp12 = NaiveDE.regress_out(sample_info_pp12, norm_expr_pp12.T, 'np.log(total_counts)').T
X = sample_info_pp12[['x', 'y']]
results_pp12 = SpatialDE.run(X, resid_expr_pp12)
results_pp12 = results_pp12.sort_values('qval')[["g", "l", "pval", "qval"]]
results_pp12.index = adata_pp12.var.loc[results_pp12["g"], "SYMBOL"].values
results_pp12.to_csv("../results_visium/results_squidpy/spatialde/de_pp12.csv")


norm_expr_pp18 = NaiveDE.stabilize(counts_pp18.T).T
resid_expr_pp18 = NaiveDE.regress_out(sample_info_pp18, norm_expr_pp18.T, 'np.log(total_counts)').T
X = sample_info_pp18[['x', 'y']]
results_pp18 = SpatialDE.run(X, resid_expr_pp18)
results_pp18 = results_pp18.sort_values('qval')[["g", "l", "pval", "qval"]]
results_pp18.index = adata_pp18.var.loc[results_pp18["g"], "SYMBOL"].values
results_pp18.to_csv("../results_visium/results_squidpy/spatialde/de_pp18.csv")


norm_expr_ppc12 = NaiveDE.stabilize(counts_ppc12.T).T
resid_expr_ppc12 = NaiveDE.regress_out(sample_info_ppc12, norm_expr_ppc12.T, 'np.log(total_counts)').T
X = sample_info_ppc12[['x', 'y']]
results_ppc12 = SpatialDE.run(X, resid_expr_ppc12)
results_ppc12 = results_ppc12.sort_values('qval')[["g", "l", "pval", "qval"]]
results_ppc12.index = adata_ppc12.var.loc[results_ppc12["g"], "SYMBOL"].values
results_ppc12.to_csv("../results_visium/results_squidpy/spatialde/de_ppc12.csv")


norm_expr_ppc18 = NaiveDE.stabilize(counts_ppc18.T).T
resid_expr_ppc18 = NaiveDE.regress_out(sample_info_ppc18, norm_expr_ppc18.T, 'np.log(total_counts)').T
X = sample_info_ppc18[['x', 'y']]
results_ppc18 = SpatialDE.run(X, resid_expr_ppc18)
results_ppc18 = results_ppc18.sort_values('qval')[["g", "l", "pval", "qval"]]
results_ppc18.index = adata_ppc18.var.loc[results_ppc18["g"], "SYMBOL"].values
results_ppc18.to_csv("../results_visium/results_squidpy/spatialde/de_ppc18.csv")


assert False
# In[]
# ---------------------------------------------------------------------
#
# Image-morphology-based clustering, 
#
# ---------------------------------------------------------------------
# Extract image
img_hires_pp12 = sq.im.ImageContainer(adata_pp12.uns["spatial"]["255"]["images"]["hires"], scale = adata_pp12.uns["spatial"]["255"]["scalefactors"]["tissue_hires_scalef"])
img_hires_pp18 = sq.im.ImageContainer(adata_pp18.uns["spatial"]["1687PP"]["images"]["hires"], scale = adata_pp18.uns["spatial"]["1687PP"]["scalefactors"]["tissue_hires_scalef"])
img_hires_ppc12 = sq.im.ImageContainer(adata_ppc12.uns["spatial"]["1161"]["images"]["hires"], scale = adata_ppc12.uns["spatial"]["1161"]["scalefactors"]["tissue_hires_scalef"])
img_hires_ppc18 = sq.im.ImageContainer(adata_ppc18.uns["spatial"]["1660PPC"]["images"]["hires"], scale = adata_ppc18.uns["spatial"]["1660PPC"]["scalefactors"]["tissue_hires_scalef"])

# Calculate image features, from multiple scales
# calculate features for different scales (higher value means more context) 
# 1.0: only image below the spot
# 2.0: only image below 2 times the spot
adata_list = [adata_pp12, adata_pp18, adata_ppc12, adata_ppc18]
img_list = [img_hires_pp12, img_hires_pp18, img_hires_ppc12, img_hires_ppc18]
for adata, img in zip(adata_list, img_list):
    # use only scale = 1.0
    scale = 1.0
    feature_name = f"features"
    sq.im.calculate_image_features(
        adata,
        img.compute(),
        # feature calculation methods
        features="summary",
        key_added=feature_name,
        n_jobs=1,
        # size of the image used
        scale=scale,
    )

    # # combine features in one dataframe
    # adata.obsm["features"] = pd.concat(
    #     [adata.obsm[f] for f in adata.obsm.keys() if "features_summary" in f],
    #     axis="columns",
    # )

    # # make sure that we have no duplicated feature names in the combined table
    # adata.obsm["features"].columns = ad.utils.make_index_unique(
    #     adata.obsm["features"].columns
    # )
# %%
