# In[]
from scDisInFact import scdisinfact, create_scdisinfact_dataset, utils
import scanpy as sc
from anndata import AnnData
import torch
import numpy as np
import os
from umap import UMAP
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams["font.size"] = 15
import warnings
warnings.filterwarnings("ignore")
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import ttest_1samp

result_dir = "results_scdisinfact/"
use_luminal = False
# In[]
adata = sc.read_h5ad("Cell_Ranger_output/adata_seurat.h5ad")
# sc.pp.filter_cells(adata, min_genes = 200)
sc.pp.filter_genes(adata, min_cells = 10)
# save raw data
counts = adata.X.toarray()
meta_cells = adata.obs

# select highly-variable genes, or there will be too many genes
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes = 5000)

# if only use hvg
hvg = adata.var["highly_variable"].values.squeeze()
adata = adata[:,hvg]
counts = counts[:, hvg]
# if only use Luminal cells
if use_luminal:
    counts = counts[meta_cells["annot"] == "Luminal",:]
    adata = adata[meta_cells["annot"] == "Luminal",:]
    meta_cells = meta_cells[meta_cells["annot"] == "Luminal"]

meta_cells["genotype"] = meta_cells["genotype"].astype("str")
meta_cells["annot"] = meta_cells["annot"].astype("str")
meta_cells["age"] = meta_cells["age"].astype("str")
data_dict = create_scdisinfact_dataset(counts, meta_cells, condition_key = ["genotype"], batch_key = "age")

# In[]
# default setting of hyper-parameters
reg_mmd_comm = 1e-4
reg_mmd_diff = 1e-4
reg_kl_comm = 1e-5
reg_kl_diff = 1
reg_class = 1
reg_gl = 1

Ks = [8, 2]

batch_size = 64
nepochs = 100
interval = 10
lr = 5e-4
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

lambs = [reg_mmd_comm, reg_mmd_diff, reg_kl_comm, reg_kl_diff, reg_class, reg_gl]
model = scdisinfact(data_dict = data_dict, Ks = Ks, batch_size = batch_size, interval = interval, lr = lr, 
                    reg_mmd_comm = reg_mmd_comm, reg_mmd_diff = reg_mmd_diff, reg_gl = reg_gl, reg_class = reg_class, 
                    reg_kl_comm = reg_kl_comm, reg_kl_diff = reg_kl_diff, seed = 0, device = device)
model.train()
losses = model.train_model(nepochs = nepochs, recon_loss = "NB")
_ = model.eval()

if use_luminal:
    # torch.save(model.state_dict(), result_dir + f"model_luminal.pth")
    model.load_state_dict(torch.load(result_dir + f"model_luminal.pth", map_location = device))
else:    
    # torch.save(model.state_dict(), result_dir + f"model_full.pth")
    model.load_state_dict(torch.load(result_dir + f"model_full.pth", map_location = device))


# In[]
z_cs = []
z_ds = []
zs = []

with torch.no_grad():
    for dataset in data_dict["datasets"]:
        # pass through the encoders
        dict_inf = model.inference(counts = dataset.counts_norm.to(model.device), batch_ids = dataset.batch_id[:,None].to(model.device), print_stat = False)
        # pass through the decoder
        dict_gen = model.generative(z_c = dict_inf["mu_c"], z_d = dict_inf["mu_d"], batch_ids = dataset.batch_id[:,None].to(model.device))
        z_c = dict_inf["mu_c"]
        z_d = dict_inf["mu_d"]
        z = torch.cat([z_c] + z_d, dim = 1)
        mu = dict_gen["mu"]    
        z_ds.append([x.cpu().detach().numpy() for x in z_d])
        z_cs.append(z_c.cpu().detach().numpy())
        zs.append(np.concatenate([z_cs[-1]] + z_ds[-1], axis = 1))

# UMAP
umap_op = UMAP(min_dist = 0.1, random_state = 0)
z_cs_umap = umap_op.fit_transform(np.concatenate(z_cs, axis = 0))
z_ds_umap = []
# z_ds_umap.append(umap_op.fit_transform(np.concatenate([z_d[0] for z_d in z_ds], axis = 0)))
# z_ds_umap.append(umap_op.fit_transform(np.concatenate([z_d[1] for z_d in z_ds], axis = 0)))
z_ds_umap.append(np.concatenate([z_d[0] for z_d in z_ds], axis = 0))

if use_luminal:
    comment = f'figures_luminal/'
else:
    comment = f'figures_full/'
if not os.path.exists(result_dir + comment):
    os.makedirs(result_dir + comment)


utils.plot_latent(zs = z_cs_umap, annos = None, batches = np.concatenate([x["sample"].to_numpy().squeeze() for x in data_dict["meta_cells"]]), \
    mode = "batches", axis_label = "UMAP", figsize = (12,7), save = (result_dir + comment+"common_dims_batches.png".format()) if result_dir else None, markerscale = 9, s = 1, alpha = 0.5)
utils.plot_latent(zs = z_cs_umap, annos = None, batches = np.concatenate([x["annot"].to_numpy().squeeze() for x in data_dict["meta_cells"]]), \
    mode = "batches", axis_label = "UMAP", figsize = (12,7), save = (result_dir + comment+"common_dims_annot.png".format()) if result_dir else None, markerscale = 9, s = 1, alpha = 0.5)
utils.plot_latent(zs = z_ds_umap[0], annos = np.concatenate([x["genotype"].to_numpy().squeeze() for x in data_dict["meta_cells"]]), batches = np.concatenate([x["sample"].to_numpy().squeeze() for x in data_dict["meta_cells"]]), \
    mode = "annos", axis_label = "UMAP", figsize = (10,7), save = (result_dir + comment+"diff_dims1_cond1.png".format()) if result_dir else None, markerscale = 9, s = 1, alpha = 0.5)

# In[]
gene_scores = model.extract_gene_scores()[0]
keygenes = adata.var.index.values[np.argsort(gene_scores)[::-1]]
keygenes = pd.DataFrame(index = keygenes, data = np.sort(gene_scores)[::-1,None], columns = ["score"])
if use_luminal:
    keygenes.to_csv(result_dir + "keygenes_luminal.csv")
else:
    keygenes.to_csv(result_dir + "keygenes_full.csv")
keygenes = keygenes[keygenes["score"] > 0.8]
print(keygenes)

# In[]
# make prediction of Luminal cells under new condition
counts_ctr_12wk = counts[(meta_cells["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)") & (meta_cells["age"] == "12wk") & (meta_cells["annot"] == "Luminal"),:]
counts_ctr_18wk = counts[(meta_cells["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)") & (meta_cells["age"] == "18wk") & (meta_cells["annot"] == "Luminal"),:]
counts_ko_12wk = counts[(meta_cells["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") & (meta_cells["age"] == "12wk") & (meta_cells["annot"] == "Luminal"),:]
counts_ko_18wk = counts[(meta_cells["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") & (meta_cells["age"] == "18wk") & (meta_cells["annot"] == "Luminal"),:]

meta_ctr_12wk = meta_cells.loc[(meta_cells["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)") & (meta_cells["age"] == "12wk") & (meta_cells["annot"] == "Luminal"),:]
meta_ctr_18wk = meta_cells.loc[(meta_cells["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)") & (meta_cells["age"] == "18wk") & (meta_cells["annot"] == "Luminal"),:]
meta_ko_12wk = meta_cells.loc[(meta_cells["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") & (meta_cells["age"] == "12wk") & (meta_cells["annot"] == "Luminal"),:]
meta_ko_18wk = meta_cells.loc[(meta_cells["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") & (meta_cells["age"] == "18wk") & (meta_cells["annot"] == "Luminal"),:]


counts_pred_ko_12wk = model.predict_counts(input_counts = counts_ctr_12wk, meta_cells = meta_ctr_12wk, condition_keys = ["genotype"], 
                                              batch_key = "age", predict_conds = ["PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"], predict_batch = "12wk")
counts_pred_ctl_12wk = model.predict_counts(input_counts = counts_ctr_12wk, meta_cells = meta_ctr_12wk, condition_keys = ["genotype"], 
                                               batch_key = "age", predict_conds = None, predict_batch = None)

counts_pred_ko_18wk = model.predict_counts(input_counts = counts_ctr_18wk, meta_cells = meta_ctr_18wk, condition_keys = ["genotype"], 
                                              batch_key = "age", predict_conds = ["PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"], predict_batch = "18wk")
counts_pred_ctl_18wk = model.predict_counts(input_counts = counts_ctr_18wk, meta_cells = meta_ctr_18wk, condition_keys = ["genotype"], 
                                               batch_key = "age", predict_conds = None, predict_batch = None)

# normalize and log-1p transform
counts_pred_ko_12wk = counts_pred_ko_12wk/(np.sum(counts_pred_ko_12wk, axis = 1, keepdims = True) + 1e-6) * 100
counts_pred_ctl_12wk = counts_pred_ctl_12wk/(np.sum(counts_pred_ctl_12wk, axis = 1, keepdims = True) + 1e-6) * 100
counts_pred_ko_18wk = counts_pred_ko_18wk/(np.sum(counts_pred_ko_18wk, axis = 1, keepdims = True) + 1e-6) * 100
counts_pred_ctl_18wk = counts_pred_ctl_18wk/(np.sum(counts_pred_ctl_18wk, axis = 1, keepdims = True) + 1e-6) * 100
counts_pred_ko_12wk = np.log1p(counts_pred_ko_12wk)
counts_pred_ctl_12wk = np.log1p(counts_pred_ctl_12wk)
counts_pred_ko_18wk = np.log1p(counts_pred_ko_18wk)
counts_pred_ctl_18wk = np.log1p(counts_pred_ctl_18wk)

counts_pred_diff_12wk = counts_pred_ko_12wk - counts_pred_ctl_12wk
counts_pred_diff_18wk = counts_pred_ko_18wk - counts_pred_ctl_18wk
counts_pred_diff_12wk = pd.DataFrame(data = counts_pred_diff_12wk, index = meta_ctr_12wk.index.values, columns = adata.var.index.values)
counts_pred_diff_18wk = pd.DataFrame(data = counts_pred_diff_18wk, index = meta_ctr_18wk.index.values, columns = adata.var.index.values)


# In[]
# load interest genes 
interest_genes = pd.read_csv("spatial_data/interest_genes.csv", index_col = 0)
interest_genes2 = pd.read_excel("markers_info/Gene list for cell count in genotype.xlsx")
gene_names = np.array([x for x in interest_genes["Genes"]])
gene_names2 = np.array([x for x in interest_genes2["Genes"]])
gene_names = np.union1d(gene_names, gene_names2)

for idx, counts_diff in enumerate([counts_pred_diff_12wk, counts_pred_diff_18wk]):
    gene_names_overlap = np.intersect1d(counts_diff.columns.values, gene_names)
    mean_expr_diff = np.mean(counts_diff.loc[:, gene_names_overlap].values, axis = 0)
    gene_names_overlap = gene_names_overlap[np.argsort(mean_expr_diff)[::-1]]

    expr_diffs = pd.DataFrame(columns = ["Expr", "Gene"])
    for gene in gene_names_overlap:
        expr_diff = pd.DataFrame(columns = ["Expr", "Gene"])
        expr_diff["Expr"] = counts_diff.loc[:, gene].values
        expr_diff["Gene"] = gene
        expr_diffs = pd.concat([expr_diffs, expr_diff], axis = 0)

    fig = plt.figure(figsize = (35, 5))
    ax = fig.add_subplot()
    sns.boxplot(data = expr_diffs, x = "Gene", y = "Expr", ax = ax)
    
    mean_expr_diff = np.mean(counts_diff.values, axis = 0)
    scores = pd.DataFrame(columns = ["gene", "scores"])
    scores["gene"] = counts_diff.columns.values[np.argsort(mean_expr_diff)[::-1]]
    scores["scores"] = np.sort(mean_expr_diff)[::-1]
    
    if idx == 0:
        ax.set_title("12wk")
        if use_luminal:
            fig.savefig(result_dir + "expr_change_genotype_12wk_luminal.png", bbox_inches = "tight")
            scores.to_csv(result_dir + "expr_change_genotype_12wk_luminal.csv")
        else:    
            fig.savefig(result_dir + "expr_change_genotype_12wk_full.png", bbox_inches = "tight")
            scores.to_csv(result_dir + "expr_change_genotype_12wk_full.csv")
    else:
        ax.set_title("18wk")
        if use_luminal:
            fig.savefig(result_dir + "expr_change_genotype_18wk_luminal.png", bbox_inches = "tight")
            scores.to_csv(result_dir + "expr_change_genotype_18wk_luminal.csv")
        else:
            fig.savefig(result_dir + "expr_change_genotype_18wk_full.png", bbox_inches = "tight")
            scores.to_csv(result_dir + "expr_change_genotype_18wk_full.csv")

    
# In[]
result_dir = "results_scdisinfact/include_castrated/"
use_luminal = True

adata = sc.read_h5ad("Cell_Ranger_output/adata_seurat.h5ad")
adata_castrated = sc.read_h5ad("Cell_Ranger_output/adata_castrated_seurat.h5ad")
adata_castrated.obs = adata_castrated.obs.rename(columns = {"annot_transfer": "annot"})

# currently consider only the condition and genotype effect, assuming they are independent
adata_merge = AnnData.concatenate(adata[adata.obs["age"] == "18wk",:], adata_castrated, join = "inner")
adata_merge.obs = adata_merge.obs[["sample", "age", "condition", "genotype", "annot"]]
adata_merge.obs["age"] = "12wk after"

# sc.pp.filter_cells(adata, min_genes = 200)
sc.pp.filter_genes(adata_merge, min_cells = 10)
# save raw data
counts = adata_merge.X.toarray()
meta_cells = adata_merge.obs

# select highly-variable genes, or there will be too many genes
sc.pp.normalize_per_cell(adata_merge)
sc.pp.log1p(adata_merge)
sc.pp.highly_variable_genes(adata_merge, n_top_genes = 5000)

# if only use hvg
hvg = adata_merge.var["highly_variable"].values.squeeze()
adata_merge = adata_merge[:,hvg]
counts = counts[:, hvg]
# if only use Luminal cells
if use_luminal:
    counts = counts[meta_cells["annot"] == "Luminal",:]
    adata_merge = adata_merge[meta_cells["annot"] == "Luminal",:]
    meta_cells = meta_cells[meta_cells["annot"] == "Luminal"]

meta_cells["genotype"] = meta_cells["genotype"].astype("str")
meta_cells["annot"] = meta_cells["annot"].astype("str")
meta_cells["condition"] = meta_cells["condition"].astype("str")

# data_dict = create_scdisinfact_dataset(counts, meta_cells, condition_key = ["genotype", "condition"], batch_key = "age")
# comparing only the condition effect of genotype, treat condition as batch
data_dict = create_scdisinfact_dataset(counts, meta_cells, condition_key = ["genotype"], batch_key = "condition")

# In[]
# default setting of hyper-parameters
reg_mmd_comm = 1e-4
reg_mmd_diff = 1e-4
reg_kl_comm = 1e-5
reg_kl_diff = 1e-2
reg_class = 1
reg_gl = 1

Ks = [8, 2]

batch_size = 64
nepochs = 100
interval = 10
lr = 5e-4
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

lambs = [reg_mmd_comm, reg_mmd_diff, reg_kl_comm, reg_kl_diff, reg_class, reg_gl]
model = scdisinfact(data_dict = data_dict, Ks = Ks, batch_size = batch_size, interval = interval, lr = lr, 
                    reg_mmd_comm = reg_mmd_comm, reg_mmd_diff = reg_mmd_diff, reg_gl = reg_gl, reg_class = reg_class, 
                    reg_kl_comm = reg_kl_comm, reg_kl_diff = reg_kl_diff, seed = 0, device = device)
model.train()
losses = model.train_model(nepochs = nepochs, recon_loss = "NB")
_ = model.eval()

if use_luminal:
    torch.save(model.state_dict(), result_dir + f"model_luminal.pth")
    model.load_state_dict(torch.load(result_dir + f"model_luminal.pth", map_location = device))
else:    
    torch.save(model.state_dict(), result_dir + f"model_full.pth")
    model.load_state_dict(torch.load(result_dir + f"model_full.pth", map_location = device))

# In[]
z_cs = []
z_ds = []
zs = []

with torch.no_grad():
    for dataset in data_dict["datasets"]:
        # pass through the encoders
        dict_inf = model.inference(counts = dataset.counts_norm.to(model.device), batch_ids = dataset.batch_id[:,None].to(model.device), print_stat = False)
        # pass through the decoder
        dict_gen = model.generative(z_c = dict_inf["mu_c"], z_d = dict_inf["mu_d"], batch_ids = dataset.batch_id[:,None].to(model.device))
        z_c = dict_inf["mu_c"]
        z_d = dict_inf["mu_d"]
        z = torch.cat([z_c] + z_d, dim = 1)
        mu = dict_gen["mu"]    
        z_ds.append([x.cpu().detach().numpy() for x in z_d])
        z_cs.append(z_c.cpu().detach().numpy())
        zs.append(np.concatenate([z_cs[-1]] + z_ds[-1], axis = 1))

# UMAP
umap_op = UMAP(min_dist = 0.1, random_state = 0)
z_cs_umap = umap_op.fit_transform(np.concatenate(z_cs, axis = 0))
z_ds_umap = []
# z_ds_umap.append(umap_op.fit_transform(np.concatenate([z_d[0] for z_d in z_ds], axis = 0)))
# z_ds_umap.append(umap_op.fit_transform(np.concatenate([z_d[1] for z_d in z_ds], axis = 0)))
z_ds_umap.append(np.concatenate([z_d[0] for z_d in z_ds], axis = 0))
# z_ds_umap.append(np.concatenate([z_d[1] for z_d in z_ds], axis = 0))

if use_luminal:
    comment = f'figures_luminal/'
else:
    comment = f'figures_full/'
if not os.path.exists(result_dir + comment):
    os.makedirs(result_dir + comment)


utils.plot_latent(zs = z_cs_umap, annos = None, batches = np.concatenate([x["sample"].to_numpy().squeeze() for x in data_dict["meta_cells"]]), \
    mode = "batches", axis_label = "UMAP", figsize = (12,7), save = (result_dir + comment+"common_dims_batches.png".format()) if result_dir else None, markerscale = 9, s = 1, alpha = 0.5)
utils.plot_latent(zs = z_cs_umap, annos = None, batches = np.concatenate([x["annot"].to_numpy().squeeze() for x in data_dict["meta_cells"]]), \
    mode = "batches", axis_label = "UMAP", figsize = (12,7), save = (result_dir + comment+"common_dims_annot.png".format()) if result_dir else None, markerscale = 9, s = 1, alpha = 0.5)
utils.plot_latent(zs = z_ds_umap[0], annos = np.concatenate([x["genotype"].to_numpy().squeeze() for x in data_dict["meta_cells"]]), batches = np.concatenate([x["sample"].to_numpy().squeeze() for x in data_dict["meta_cells"]]), \
    mode = "annos", axis_label = "UMAP", figsize = (10,7), save = (result_dir + comment+"diff_dims1_cond1.png".format()) if result_dir else None, markerscale = 9, s = 1, alpha = 0.5)
# utils.plot_latent(zs = z_ds_umap[1], annos = np.concatenate([x["condition"].to_numpy().squeeze() for x in data_dict["meta_cells"]]), batches = np.concatenate([x["sample"].to_numpy().squeeze() for x in data_dict["meta_cells"]]), \
#     mode = "annos", axis_label = "UMAP", figsize = (10,7), save = (result_dir + comment+"diff_dims2_cond2.png".format()) if result_dir else None, markerscale = 9, s = 1, alpha = 0.5)

# In[]
# make prediction of Luminal cells under new condition
counts_intact_ctr = counts[(meta_cells["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)") & (meta_cells["condition"] == "intact") & (meta_cells["annot"] == "Luminal"),:]
counts_intact_ko = counts[(meta_cells["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") & (meta_cells["condition"] == "intact") & (meta_cells["annot"] == "Luminal"),:]
counts_cas_ctr = counts[(meta_cells["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)") & (meta_cells["condition"] == "castrated") & (meta_cells["annot"] == "Luminal"),:]
counts_cas_ko = counts[(meta_cells["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") & (meta_cells["condition"] == "castrated") & (meta_cells["annot"] == "Luminal"),:]

meta_intact_ctr = meta_cells.loc[(meta_cells["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)") & (meta_cells["condition"] == "intact") & (meta_cells["annot"] == "Luminal"),:]
meta_intact_ko = meta_cells.loc[(meta_cells["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") & (meta_cells["condition"] == "intact") & (meta_cells["annot"] == "Luminal"),:]
meta_cas_ctr = meta_cells.loc[(meta_cells["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)") & (meta_cells["condition"] == "castrated") & (meta_cells["annot"] == "Luminal"),:]
meta_cas_ko = meta_cells.loc[(meta_cells["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") & (meta_cells["condition"] == "castrated") & (meta_cells["annot"] == "Luminal"),:]

counts_pred_intact_ctr = model.predict_counts(input_counts = counts_intact_ctr, meta_cells = meta_intact_ctr, condition_keys = ["genotype"], 
                                              batch_key = "condition", predict_conds = None, predict_batch = None)
counts_pred_intact_ko = model.predict_counts(input_counts = counts_intact_ctr, meta_cells = meta_intact_ctr, condition_keys = ["genotype"], 
                                              batch_key = "condition", predict_conds = ["PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"], predict_batch = "intact")

counts_pred_cas_ctr = model.predict_counts(input_counts = counts_cas_ctr, meta_cells = meta_cas_ctr, condition_keys = ["genotype"], 
                                              batch_key = "condition", predict_conds = None, predict_batch = None)
counts_pred_cas_ko = model.predict_counts(input_counts = counts_cas_ctr, meta_cells = meta_cas_ctr, condition_keys = ["genotype"], 
                                              batch_key = "condition", predict_conds = ["PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)"], predict_batch = "castrated")


# normalize and log-1p transform
counts_pred_intact_ctr = counts_pred_intact_ctr/(np.sum(counts_pred_intact_ctr, axis = 1, keepdims = True) + 1e-6) * 100
counts_pred_intact_ko = counts_pred_intact_ko/(np.sum(counts_pred_intact_ko, axis = 1, keepdims = True) + 1e-6) * 100
counts_pred_cas_ctr = counts_pred_cas_ctr/(np.sum(counts_pred_cas_ctr, axis = 1, keepdims = True) + 1e-6) * 100
counts_pred_cas_ko = counts_pred_cas_ko/(np.sum(counts_pred_cas_ko, axis = 1, keepdims = True) + 1e-6) * 100
counts_pred_intact_ctr = np.log1p(counts_pred_intact_ctr)
counts_pred_intact_ko = np.log1p(counts_pred_intact_ko)
counts_pred_cas_ctr = np.log1p(counts_pred_cas_ctr)
counts_pred_cas_ko = np.log1p(counts_pred_cas_ko)

counts_pred_diff_intact = counts_pred_intact_ko - counts_pred_intact_ctr
counts_pred_diff_cas = counts_pred_cas_ko - counts_pred_cas_ctr
counts_pred_diff_intact = pd.DataFrame(data = counts_pred_diff_intact, index = meta_intact_ctr.index.values, columns = adata_merge.var.index.values)
counts_pred_diff_cas = pd.DataFrame(data = counts_pred_diff_cas, index = meta_cas_ctr.index.values, columns = adata_merge.var.index.values)

# In[]
ttest_result_intact = ttest_1samp(counts_pred_diff_intact.values.squeeze(), popmean = 0, axis = 0, alternative = "two-sided")
_, pval_fdr_intact = fdrcorrection(ttest_result_intact.pvalue, alpha = 0.05)
ttest_result_cas = ttest_1samp(counts_pred_diff_cas.values.squeeze(), popmean = 0, axis = 0, alternative = "two-sided")
_, pval_fdr_cas = fdrcorrection(ttest_result_cas.pvalue, alpha = 0.05)

pvals = pd.DataFrame(data = 0, index = counts_pred_diff_intact.columns.values, columns = ["intact", "castrated", "intact (fdr)", "castrated (fdr)", "intact (mean)", "castrated (mean)"])
pvals["intact"] = ttest_result_intact.pvalue
pvals["intact (fdr)"] = pval_fdr_intact
pvals["castrated"] = ttest_result_cas.pvalue
pvals["castrated (fdr)"] = pval_fdr_cas
pvals["intact (mean)"] = np.mean(counts_pred_diff_intact.values, axis = 0)
pvals["castrated (mean)"] = np.mean(counts_pred_diff_cas.values, axis = 0)
pvals["change (mean)"] = np.mean(counts_pred_diff_cas.values, axis = 0) - np.mean(counts_pred_diff_intact.values, axis = 0)

pvals = pvals.iloc[np.argsort(np.abs(pvals["change (mean)"].values.squeeze()))[::-1],:]

pvals.to_csv(result_dir + "change_pval.csv")


# In[]
# Check the expression change of CRPC related genes, castration-resistant prostate cancer (CRPC)




# %%
