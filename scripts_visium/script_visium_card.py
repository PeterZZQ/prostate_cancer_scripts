# In[]
import sys, os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs
import seaborn as sns
from matplotlib.pyplot import cm 
import matplotlib.patches as mpatches

sys.path.append("../")
import utils

# In[]
# --------------------------------------------------------------------------------
#
# Load Visium data
#
# --------------------------------------------------------------------------------
adata_pp12 = sc.read_visium(path = "../dataset/data_visium/spatial_data/Visium_python/225_PP12/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_pp18 = sc.read_visium(path = "../dataset/data_visium/spatial_data/Visium_python/1687_PP18/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_ppc12 = sc.read_visium(path = "../dataset/data_visium/spatial_data/Visium_python/1161_PPC/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")
adata_ppc18 = sc.read_visium(path = "../dataset/data_visium/spatial_data/Visium_python/1660_PPC18/", count_file = "filtered_feature_bc_matrix.h5", source_image_path = "spatial")

adata_pp12.var_names_make_unique()
adata_pp18.var_names_make_unique()
adata_ppc12.var_names_make_unique()
adata_ppc18.var_names_make_unique()

adata_pp12.obsm["spatial"] = adata_pp12.obsm["spatial"].astype(np.float64)
adata_pp18.obsm["spatial"] = adata_pp18.obsm["spatial"].astype(np.float64)
adata_ppc12.obsm["spatial"] = adata_ppc12.obsm["spatial"].astype(np.float64)
adata_ppc18.obsm["spatial"] = adata_ppc18.obsm["spatial"].astype(np.float64)

# # rename gene names from SYMBOL to ENSEMBLE_ID
# adata_pp12.var["SYMBOL"] = adata_pp12.var_names
# adata_pp12.var.index = adata_pp12.var["gene_ids"].values
# adata_pp18.var["SYMBOL"] = adata_pp18.var_names
# adata_pp18.var.index = adata_pp18.var["gene_ids"].values
# adata_ppc12.var["SYMBOL"] = adata_ppc12.var_names
# adata_ppc12.var.index = adata_ppc12.var["gene_ids"].values
# adata_ppc18.var["SYMBOL"] = adata_ppc18.var_names
# adata_ppc18.var.index = adata_ppc18.var["gene_ids"].values

adata_visium = {"pp12": adata_pp12, "pp18": adata_pp18, "ppc12": adata_ppc12, "ppc18": adata_ppc18}
celltypes_categories = np.array(['Luminal', 'Luminal (Spink1+)', 'Basal', 'Club epithelia', 'Endothelial', 
                                 'Lymphoid', 'Monocytes', 'Macrophage', 'Mesenchymal', 'SV'])

res_dir = "../results_visium/results_card/"

# In[]
# --------------------------------------------------------------------------------
#
# Load deconvolution scores of cards, and calculate the cell type composition
#
# --------------------------------------------------------------------------------
deconv_score_pp12 = pd.read_csv(res_dir + "deconv_score_pp12.csv", sep = "\t").loc[adata_visium["pp12"].obs.index.values,:]
deconv_score_pp18 = pd.read_csv(res_dir + "deconv_score_pp18.csv", sep = "\t").loc[adata_visium["pp18"].obs.index.values,:]
deconv_score_ppc12 = pd.read_csv(res_dir + "deconv_score_ppc12.csv", sep = "\t").loc[adata_visium["ppc12"].obs.index.values,:]
deconv_score_ppc18 = pd.read_csv(res_dir + "deconv_score_ppc18.csv", sep = "\t").loc[adata_visium["ppc18"].obs.index.values,:]

adata_visium["pp12"].obs = pd.concat([adata_visium["pp12"].obs, deconv_score_pp12], axis = 1)
adata_visium["pp18"].obs = pd.concat([adata_visium["pp18"].obs, deconv_score_pp18], axis = 1)
adata_visium["ppc12"].obs = pd.concat([adata_visium["ppc12"].obs, deconv_score_ppc12], axis = 1)
adata_visium["ppc18"].obs = pd.concat([adata_visium["ppc18"].obs, deconv_score_ppc18], axis = 1)

# calculate cell type counts
ct_pp12 = deconv_score_pp12.sum(axis = 0)
ct_pp18 = deconv_score_pp18.sum(axis = 0)
ct_ppc12 = deconv_score_ppc12.sum(axis = 0)
ct_ppc18 = deconv_score_ppc18.sum(axis = 0)

ct_composition = pd.concat([ct_pp12, ct_pp18, ct_ppc12, ct_ppc18], axis = 1)
ct_composition.columns = ["CTRL (12wk)", "CXCR7 KO (12wk)", "CTRL (18wk)", "CXCR7 KO (18wk)"]
ct_composition.to_csv(res_dir + "count_celltype.csv")

# normalize the "count" into percentage
ct_composition.loc[:,:] = ct_composition.values/np.sum(ct_composition.values, axis = 0, keepdims = True) * 100
ct_composition.to_csv(res_dir + "pct_celltype.csv")

# In[]
# --------------------------------------------------------------------------------
#
# Plot the cell type composition 
#
# --------------------------------------------------------------------------------
# 1. Plot the pie chart of cell type composition, merge the immune and other cells
ct_composition.loc["Other",:] = ct_composition.loc[["SV", "Endothelial"],:].sum(axis = 0)
ct_composition.loc["Immune",:] = ct_composition.loc[["Lymphoid", "Monocytes", "Macrophage"],:].sum(axis = 0)
# ct_composition = ct_composition.loc[["Luminal", "Basal", "Mesenchymal", "Club epithelia", "Luminal (Spink1+)", "Other", "Immune"],:]
ct_composition = ct_composition.loc[["Luminal", "Basal", "Immune", "Luminal (Spink1+)", "Mesenchymal", "Club epithelia", "Other"],:]
ct_pp12 = ct_composition["CTRL (12wk)"]
ct_ppc12 = ct_composition["CXCR7 KO (12wk)"]
ct_pp18 = ct_composition["CTRL (18wk)"]
ct_ppc18 = ct_composition["CXCR7 KO (18wk)"]

plt.rcParams["font.size"] = 10
fig = plt.figure(figsize = (12,12))
axs = fig.subplots(nrows = 2, ncols = 2)
axs[0,0].pie(ct_pp12.values, labels = ct_pp12.index.values, autopct='%.0f%%')
axs[0,1].pie(ct_ppc12.values, labels = ct_ppc12.index.values, autopct='%.0f%%')
axs[1,0].pie(ct_pp18.values, labels = ct_pp18.index.values, autopct='%.0f%%')
axs[1,1].pie(ct_ppc18.values, labels = ct_ppc18.index.values, autopct='%.0f%%')

axs[0,0].set_title("CTRL (12wk)", fontsize = 20)
axs[0,1].set_title("CXCR7 KO (12wk)", fontsize = 20)
axs[1,0].set_title("CTRL (18wk)", fontsize = 20)
axs[1,1].set_title("CXCR7 KO (18wk)", fontsize = 20)
fig.suptitle("Cell percentage", fontsize = 30)
fig.tight_layout()

fig.savefig(res_dir + "percent_celltype.png", bbox_inches = "tight")

# In[]
# 2. Plot the stack bar plot of cell type composition
cmap = cm.get_cmap("tab10")
plt.rcParams["font.size"] = 15
fig = plt.figure(figsize = (10,7))
ax = fig.subplots(nrows = 1, ncols = 1)

cum_ctrl_12wk = 100
cum_ctrl_18wk = 100
cum_ko_12wk = 100
cum_ko_18wk = 100
legend_ct = []

# the unique cell types
unique_ct = ["Luminal", "Basal", "Immune", "Luminal (Spink1+)", "Mesenchymal", "Club epithelia", "Other"]
for idx, ct in enumerate(unique_ct[::-1]):
    # percentage
    pct_ctrl_12wk = ct_composition.loc[ct, "CTRL (12wk)"]
    pct_ctrl_18wk = ct_composition.loc[ct, "CTRL (18wk)"]
    pct_ko_12wk = ct_composition.loc[ct, "CXCR7 KO (12wk)"]
    pct_ko_18wk = ct_composition.loc[ct, "CXCR7 KO (18wk)"]

    data = pd.DataFrame(columns = ["condition", "percentage"])
    data["Condition"] = ["CTRL (12wk)", "CXCR7 KO (12wk)", "CTRL (18wk)", "CXCR7 KO (18wk)"]
    data["Percentage"] = [cum_ctrl_12wk, cum_ko_12wk, cum_ctrl_18wk, cum_ko_18wk]
    # plot
    sns.barplot(data, x = "Condition", y = "Percentage", ax = ax, color = cmap(idx))

    cum_ctrl_12wk = cum_ctrl_12wk - pct_ctrl_12wk
    cum_ko_12wk = cum_ko_12wk - pct_ko_12wk
    cum_ctrl_18wk = cum_ctrl_18wk - pct_ctrl_18wk
    cum_ko_18wk = cum_ko_18wk - pct_ko_18wk
    legend_ct.append(mpatches.Patch(color = cmap(idx), label = ct))

leg = ax.legend(loc='upper left', prop={'size': 15}, frameon = False, handles = legend_ct, bbox_to_anchor=(1.04, 1))
fig.savefig(res_dir + "barplot_celltype.png", bbox_inches = "tight")


# In[]
# --------------------------------------------------------------------------------
#
# Calculate the percentage of immune subtypes, use consistent markers as scRNA-seq
#
# --------------------------------------------------------------------------------
markers_lymphoid = {}
markers_lymphoid["CD4_T"] = ["Cd4"]
markers_lymphoid["CD8_T"] = ["Cd8a"]
markers_lymphoid["cytotoxicity T"] = ["Gzmb", "Prf1"]
markers_lymphoid["exhausted T"] = ["Tigit", "Havcr2"]
markers_lymphoid["NK"] = ["Klrd1"]
markers_lymphoid["B"] = ["Ms4a1"]

markers_macrophage = {}
markers_macrophage["macrophage_m1"] = ["Cd80", "Cd86"]
markers_macrophage["macrophage_m2"] = ["Cd163", "Arg1", "Mrc1"]

markers_immune = {"Lymphoid": markers_lymphoid, "Macrophage": markers_macrophage}

# Consider each spot is a bulk of various cell types, the subtype percentage should be proportion to 
# proportion of immune cell * the expression abundance of the markers
prop_lymphoid_pp12 = deconv_score_pp12["Lymphoid"].values.squeeze()
prop_lymphoid_ppc12 = deconv_score_ppc12["Lymphoid"].values.squeeze()
prop_lymphoid_pp18 = deconv_score_pp18["Lymphoid"].values.squeeze()
prop_lymphoid_ppc18 = deconv_score_ppc18["Lymphoid"].values.squeeze()

prop_macrophage_pp12 = deconv_score_pp12["Macrophage"].values.squeeze()
prop_macrophage_ppc12 = deconv_score_ppc12["Macrophage"].values.squeeze()
prop_macrophage_pp18 = deconv_score_pp18["Macrophage"].values.squeeze()
prop_macrophage_ppc18 = deconv_score_ppc18["Macrophage"].values.squeeze()

prop_subtype = {"Lymphoid": {"pp12": prop_lymphoid_pp12, "ppc12": prop_lymphoid_ppc12, "pp18": prop_lymphoid_pp18, "ppc18": prop_lymphoid_ppc18},
                "Macrophage": {"pp12": prop_macrophage_pp12, "ppc12": prop_macrophage_ppc12, "pp18": prop_macrophage_pp18, "ppc18": prop_macrophage_ppc18}}

# calculate the average log-normalized marker expression
# log-normalize the data
import anndata
adata_merge = anndata.concat([x for x in adata_visium.values()], join = "inner", axis = 0, label = "sample", keys = [x for x in adata_visium.keys()])
sc.pp.normalize_total(adata_merge, 1e4)
sc.pp.log1p(adata_merge)

for subtype, markers in markers_immune.items():
    for ct in markers.keys():
        expr_ct_pp12 = adata_merge[adata_merge.obs["sample"] == "pp12", markers[ct]].X.toarray().mean(axis = 1)
        expr_ct_ppc12 = adata_merge[adata_merge.obs["sample"] == "ppc12", markers[ct]].X.toarray().mean(axis = 1).squeeze()
        expr_ct_pp18 = adata_merge[adata_merge.obs["sample"] == "pp18", markers[ct]].X.toarray().mean(axis = 1).squeeze()
        expr_ct_ppc18 = adata_merge[adata_merge.obs["sample"] == "ppc18", markers[ct]].X.toarray().mean(axis = 1).squeeze()

        prop_ct_pp12 = prop_subtype[subtype]["pp12"] * expr_ct_pp12
        prop_ct_ppc12 = prop_subtype[subtype]["ppc12"] * expr_ct_ppc12
        prop_ct_pp18 = prop_subtype[subtype]["pp18"] * expr_ct_pp18
        prop_ct_ppc18 = prop_subtype[subtype]["ppc18"] * expr_ct_ppc18

        # save the subtype proportion in adata_merge
        adata_merge.obs[ct] = np.concatenate([prop_ct_pp12, prop_ct_ppc12, prop_ct_pp18, prop_ct_ppc18], axis = 0)


# In[]
# Plot the immune subtype distribution in space
if not os.path.exists(res_dir + "immune_subtype/"):
    os.makedirs(res_dir + "immune_subtype/")

adata_pp12 = adata_merge[adata_merge.obs["sample"] == "pp12"]
adata_ppc12 = adata_merge[adata_merge.obs["sample"] == "ppc12"]
adata_pp18 = adata_merge[adata_merge.obs["sample"] == "pp18"]
adata_ppc18 = adata_merge[adata_merge.obs["sample"] == "ppc18"]
# Can barely see any lymphoid subtypes across genotypes
fig = utils.plot_embeds_continuous(embed = adata_pp12.obsm["spatial"], annos = adata_pp12.obs[[x for x in markers_lymphoid.keys()]], 
                                   figsize = (12, 7), s = 30, alpha = 1, markerscale = 10, ncols = 2, axis_label = "Spatial",
                                   vmax = np.max(adata_merge.obs[[x for x in markers_lymphoid.keys()]].values, axis = 0))
fig.suptitle("PP12")
fig.savefig(res_dir + "immune_subtype/" + "lymphoid_pp12.png", bbox_inches = "tight", dpi = 150)
fig = utils.plot_embeds_continuous(embed = adata_ppc12.obsm["spatial"], annos = adata_ppc12.obs[[x for x in markers_lymphoid.keys()]], 
                                   figsize = (12, 7), s = 30, alpha = 1, markerscale = 10, ncols = 2, axis_label = "Spatial",
                                   vmax = np.max(adata_merge.obs[[x for x in markers_lymphoid.keys()]].values, axis = 0))
fig.suptitle("PPC12")
fig.savefig(res_dir + "immune_subtype/" + "lymphoid_ppc12.png", bbox_inches = "tight", dpi = 150)
fig = utils.plot_embeds_continuous(embed = adata_pp18.obsm["spatial"], annos = adata_pp18.obs[[x for x in markers_lymphoid.keys()]],
                                   figsize = (12, 7), s = 30, alpha = 1, markerscale = 10, ncols = 2, axis_label = "Spatial",
                                   vmax = np.max(adata_merge.obs[[x for x in markers_lymphoid.keys()]].values, axis = 0))
fig.suptitle("PP18")
fig.savefig(res_dir + "immune_subtype/" + "lymphoid_pp18.png", bbox_inches = "tight", dpi = 150)
fig = utils.plot_embeds_continuous(embed = adata_ppc18.obsm["spatial"], annos = adata_ppc18.obs[[x for x in markers_lymphoid.keys()]],
                                   figsize = (12, 7), s = 30, alpha = 1, markerscale = 10, ncols = 2, axis_label = "Spatial",
                                   vmax = np.max(adata_merge.obs[[x for x in markers_lymphoid.keys()]].values, axis = 0))
fig.suptitle("PPC18")
fig.savefig(res_dir + "immune_subtype/" + "lymphoid_ppc18.png", bbox_inches = "tight", dpi = 150)


# Macrophage is more obvious
fig = utils.plot_embeds_continuous(embed = adata_pp12.obsm["spatial"], annos = adata_pp12.obs[[x for x in markers_macrophage.keys()]],
                                   figsize = (12, 7), s = 30, alpha = 1, markerscale = 10, ncols = 2, axis_label = "Spatial",
                                   vmax = np.max(adata_merge.obs[[x for x in markers_macrophage.keys()]].values, axis = 0))
fig.suptitle("PP12")
fig.savefig(res_dir + "immune_subtype/" + "macrophage_pp12.png", bbox_inches = "tight", dpi = 150)
fig = utils.plot_embeds_continuous(embed = adata_ppc12.obsm["spatial"], annos = adata_ppc12.obs[[x for x in markers_macrophage.keys()]],
                                   figsize = (12, 7), s = 30, alpha = 1, markerscale = 10, ncols = 2, axis_label = "Spatial",
                                   vmax = np.max(adata_merge.obs[[x for x in markers_macrophage.keys()]].values, axis = 0))
fig.suptitle("PPC12")
fig.savefig(res_dir + "immune_subtype/" + "macrophage_ppc12.png", bbox_inches = "tight", dpi = 150)
fig = utils.plot_embeds_continuous(embed = adata_pp18.obsm["spatial"], annos = adata_pp18.obs[[x for x in markers_macrophage.keys()]],
                                   figsize = (12, 7), s = 30, alpha = 1, markerscale = 10, ncols = 2, axis_label = "Spatial",
                                   vmax = np.max(adata_merge.obs[[x for x in markers_macrophage.keys()]].values, axis = 0))
fig.suptitle("PP18")
fig.savefig(res_dir + "immune_subtype/" + "macrophage_pp18.png", bbox_inches = "tight", dpi = 150)
fig = utils.plot_embeds_continuous(embed = adata_ppc18.obsm["spatial"], annos = adata_ppc18.obs[[x for x in markers_macrophage.keys()]],
                                   figsize = (12, 7), s = 30, alpha = 1, markerscale = 10, ncols = 2, axis_label = "Spatial",
                                   vmax = np.max(adata_merge.obs[[x for x in markers_macrophage.keys()]].values, axis = 0))
fig.suptitle("PPC18")
fig.savefig(res_dir + "immune_subtype/" + "macrophage_ppc18.png", bbox_inches = "tight", dpi = 150)

# In[]
# Calculate the abundance
props = []
for subtype, markers in markers_immune.items():
    for ct in markers.keys():
        prop = pd.DataFrame(columns = ["Abundance", "GT", "CT"])
        # sum then normalize by total number of spots
        prop_ct_pp12 = adata_merge[adata_merge.obs["sample"] == "pp12"].obs[ct].mean()
        prop_ct_ppc12 = adata_merge[adata_merge.obs["sample"] == "ppc12"].obs[ct].mean()
        prop_ct_pp18 = adata_merge[adata_merge.obs["sample"] == "pp18"].obs[ct].mean()
        prop_ct_ppc18 = adata_merge[adata_merge.obs["sample"] == "ppc18"].obs[ct].mean()
        prop["Abundance"] = np.array([prop_ct_pp12, prop_ct_ppc12, prop_ct_pp18, prop_ct_ppc18])
        prop["GT"] = np.array(["PP 12", "PPC 12", "PP 18", "PPC 18"])
        prop["CT"] = ct
        props.append(prop)
props = pd.concat(props, axis = 0, ignore_index = True)
props.to_csv(res_dir + "immune_subtype/" + "abundance.csv")

# NOTE: the abundance value is not the percentage, cannot use pie chart, use barplot for cross modality comparison instead
sns.set_theme(style = "darkgrid", font_scale = 1.5)
fig = plt.figure(figsize = (13, 14))
axs = fig.subplots(nrows = 2, ncols = 1)
props_lymphoid_12 = props[(props["CT"].isin([x for x in markers_lymphoid.keys()])) & (props["GT"].isin(["PP 12", "PPC 12"]))]
sns.barplot(data = props_lymphoid_12, x = "CT", hue = "GT", y = "Abundance", ax = axs[0])

props_lymphoid_18 = props[(props["CT"].isin([x for x in markers_lymphoid.keys()])) & (props["GT"].isin(["PP 18", "PPC 18"]))]
sns.barplot(data = props_lymphoid_18, x = "CT", hue = "GT", y = "Abundance", ax = axs[1])

handles, previous_labels = axs[0].get_legend_handles_labels()
leg = axs[0].legend(handles = handles, labels = ["PP", "PPC"], loc='upper left', prop={'size': 25}, frameon = False, bbox_to_anchor=(1.04, 1), fontsize = 20)
handles, previous_labels = axs[1].get_legend_handles_labels()
leg = axs[1].legend(handles = handles, labels = ["PP", "PPC"], loc='upper left', prop={'size': 25}, frameon = False, bbox_to_anchor=(1.04, 1), fontsize = 20)

axs[0].set_xlabel(None)
axs[1].set_xlabel(None)
axs[0].set_ylim([0, 1.1e-4])
axs[1].set_ylim([0, 1.1e-4])
axs[0].set_xticklabels(axs[0].get_xticklabels(), rotation=45, ha='right')
axs[1].set_xticklabels(axs[1].get_xticklabels(), rotation=45, ha='right')
axs[0].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
axs[1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
axs[0].set_title("12 wk", fontsize = 30)
axs[1].set_title("18 wk", fontsize = 30)
fig.tight_layout()
fig.savefig(res_dir + "immune_subtype/" + "lymphoid_subtypes.png", bbox_inches = "tight")

fig = plt.figure(figsize = (10, 14))
axs = fig.subplots(nrows = 2, ncols = 1)
props_macrophage_12 = props[(props["CT"].isin([x for x in markers_macrophage.keys()])) & (props["GT"].isin(["PP 12", "PPC 12"]))]
sns.barplot(data = props_macrophage_12, x = "CT", hue = "GT", y = "Abundance", ax = axs[0])

props_macrophage_18 = props[(props["CT"].isin([x for x in markers_macrophage.keys()])) & (props["GT"].isin(["PP 18", "PPC 18"]))]
sns.barplot(data = props_macrophage_18, x = "CT", hue = "GT", y = "Abundance", ax = axs[1])

handles, previous_labels = axs[0].get_legend_handles_labels()
leg = axs[0].legend(handles = handles, labels = ["PP", "PPC"], loc='upper left', prop={'size': 25}, frameon = False, bbox_to_anchor=(1.04, 1), fontsize = 20)
handles, previous_labels = axs[1].get_legend_handles_labels()
leg = axs[1].legend(handles = handles, labels = ["PP", "PPC"], loc='upper left', prop={'size': 25}, frameon = False, bbox_to_anchor=(1.04, 1), fontsize = 20)

axs[0].set_xlabel(None)
axs[1].set_xlabel(None)
# axs[0].set_ylim([0, 1.1e-4])
# axs[1].set_ylim([0, 1.1e-4])
# axs[0].set_xticklabels(axs[0].get_xticklabels(), rotation=45, ha='right')
# axs[1].set_xticklabels(axs[1].get_xticklabels(), rotation=45, ha='right')
axs[0].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
axs[1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
axs[0].set_title("12 wk", fontsize = 30)
axs[1].set_title("18 wk", fontsize = 30)
fig.tight_layout()
fig.savefig(res_dir + "immune_subtype/" + "macrophage_subtypes.png", bbox_inches = "tight")




# In[]
# # select the maximum scoring label of card
# adata_visium["pp12"].obs["label_card"] = deconv_score_pp12.columns.values[np.argmax(deconv_score_pp12.values, axis = 1)]
# adata_visium["pp12"].obs["label_card"] = pd.Categorical(adata_visium["pp12"].obs["label_card"], categories=celltypes_categories)

# adata_visium["pp18"].obs["label_card"] = deconv_score_pp18.columns.values[np.argmax(deconv_score_pp18.values, axis = 1)]
# adata_visium["pp18"].obs["label_card"] = pd.Categorical(adata_visium["pp18"].obs["label_card"], categories=celltypes_categories)

# adata_visium["ppc12"].obs["label_card"] = deconv_score_ppc12.columns.values[np.argmax(deconv_score_ppc12.values, axis = 1)]
# adata_visium["ppc12"].obs["label_card"] = pd.Categorical(adata_visium["ppc12"].obs["label_card"], categories=celltypes_categories)

# adata_visium["ppc18"].obs["label_card"] = deconv_score_ppc18.columns.values[np.argmax(deconv_score_ppc18.values, axis = 1)]
# adata_visium["ppc18"].obs["label_card"] = pd.Categorical(adata_visium["ppc18"].obs["label_card"], categories=celltypes_categories)

# # In[]
# # Load data
# sample = "ppc18"
# adata_vis = adata_visium[sample]
# result_dir = "./results_card/results_" + sample + "/"
# if not os.path.exists(result_dir):
#     os.makedirs(result_dir)

# fig = plt.figure(figsize = (12,7))
# ax = fig.add_subplot()
# sc.pl.spatial(adata_vis, img_key = "hires", color = "label_card", ax = ax)
# fig.savefig(result_dir + f"label_{sample}.pdf", bbox_inches = "tight")
# fig.savefig(result_dir + f"label_{sample}.png", bbox_inches = "tight", dpi = 300)

# fig = plt.figure(figsize = (12,7))
# ax = fig.add_subplot()
# sc.pl.spatial(adata_vis, img_key = "hires", color = "Spink1+ (Luminal)", ax = ax)
# fig.savefig(result_dir + f"label_{sample}_spink1.pdf", bbox_inches = "tight")
# fig.savefig(result_dir + f"label_{sample}_spink1.png", bbox_inches = "tight")

# fig = plt.figure(figsize = (12,7))
# ax = fig.add_subplot()
# sc.pl.spatial(adata_vis, img_key = "hires", color = "Luminal", ax = ax)
# fig.savefig(result_dir + f"label_{sample}_Luminal.pdf", bbox_inches = "tight")
# fig.savefig(result_dir + f"label_{sample}_Luminal.png", bbox_inches = "tight")

# fig = plt.figure(figsize = (12,7))
# ax = fig.add_subplot()
# sc.pl.spatial(adata_vis, img_key = "hires", color = "Basal", ax = ax)
# fig.savefig(result_dir + f"label_{sample}_Basal.pdf", bbox_inches = "tight")
# fig.savefig(result_dir + f"label_{sample}_Basal.png", bbox_inches = "tight")

# fig = plt.figure(figsize = (12,7))
# ax = fig.add_subplot()
# sc.pl.spatial(adata_vis, img_key = "hires", color = "Club epithelia (Luminal)", ax = ax)
# fig.savefig(result_dir + f"label_{sample}_Club_epithelia_(Luminal).pdf", bbox_inches = "tight")
# fig.savefig(result_dir + f"label_{sample}_Club_epithelia_(Luminal).png", bbox_inches = "tight")

# fig = plt.figure(figsize = (12,7))
# ax = fig.add_subplot()
# sc.pl.spatial(adata_vis, img_key = "hires", color = "Endothelial", ax = ax)
# fig.savefig(result_dir + f"label_{sample}_Endothelial.pdf", bbox_inches = "tight")
# fig.savefig(result_dir + f"label_{sample}_Endothelial.png", bbox_inches = "tight")

# fig = plt.figure(figsize = (12,7))
# ax = fig.add_subplot()
# sc.pl.spatial(adata_vis, img_key = "hires", color = "Hillock epithelia (Basal)", ax = ax)
# fig.savefig(result_dir + f"label_{sample}_Hillock_epithelia_(Basal).pdf", bbox_inches = "tight")
# fig.savefig(result_dir + f"label_{sample}_Hillock_epithelia_(Basal).png", bbox_inches = "tight")

# fig = plt.figure(figsize = (12,7))
# ax = fig.add_subplot()
# sc.pl.spatial(adata_vis, img_key = "hires", color = "Lymphoid (Total immune)", ax = ax)
# fig.savefig(result_dir + f"label_{sample}_Lymphoid_(Total immune).pdf", bbox_inches = "tight")
# fig.savefig(result_dir + f"label_{sample}_Lymphoid_(Total immune).png", bbox_inches = "tight")

# fig = plt.figure(figsize = (12,7))
# ax = fig.add_subplot()
# sc.pl.spatial(adata_vis, img_key = "hires", color = "Macrophage (Myeloid, Total immune)", ax = ax)
# fig.savefig(result_dir + f"label_{sample}_Macrophage_(Myeloid, Total immune).pdf", bbox_inches = "tight")
# fig.savefig(result_dir + f"label_{sample}_Macrophage_(Myeloid, Total immune).png", bbox_inches = "tight")

# fig = plt.figure(figsize = (12,7))
# ax = fig.add_subplot()
# sc.pl.spatial(adata_vis, img_key = "hires", color = "Mesenchymal", ax = ax)
# fig.savefig(result_dir + f"label_{sample}_Mesenchymal.pdf", bbox_inches = "tight")
# fig.savefig(result_dir + f"label_{sample}_Mesenchymal.png", bbox_inches = "tight")

# fig = plt.figure(figsize = (12,7))
# ax = fig.add_subplot()
# sc.pl.spatial(adata_vis, img_key = "hires", color = "Monocytes (Myeloid, Total immune)", ax = ax)
# fig.savefig(result_dir + f"label_{sample}_Monocytes (Myeloid, Total immune).pdf", bbox_inches = "tight")
# fig.savefig(result_dir + f"label_{sample}_Monocytes (Myeloid, Total immune).png", bbox_inches = "tight")

# fig = plt.figure(figsize = (12,7))
# ax = fig.add_subplot()
# sc.pl.spatial(adata_vis, img_key = "hires", color = "SV", ax = ax)
# fig.savefig(result_dir + f"label_{sample}_SV.pdf", bbox_inches = "tight")
# fig.savefig(result_dir + f"label_{sample}_SV.png", bbox_inches = "tight")

# %%
