# ------------------------------------------------------------ # 
# 
# Running RNA velocity on Luly's project data
#
# ------------------------------------------------------------ # 

# In[]
import numpy as np 
import scanpy as sc 
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd

# In[]
# Read the data
data_dir = "../dataset/data_scrnaseq/FOXA1/"

M1385_adata = scv.read(data_dir + "M1385.loom", cache=True)
M1386_adata = scv.read(data_dir + "M1386.loom", cache=True)
M1387_adata = scv.read(data_dir + "M1387.loom", cache=True)
M1415_adata = scv.read(data_dir + "M1415.loom", cache=True)
M1418_adata = scv.read(data_dir + "M1418.loom", cache=True)
M1475_adata = scv.read(data_dir + "M1475.loom", cache=True)

# In[]
M1385_adata.var_names_make_unique()
M1386_adata.var_names_make_unique()
M1387_adata.var_names_make_unique()
M1415_adata.var_names_make_unique()
M1418_adata.var_names_make_unique()
M1475_adata.var_names_make_unique()

# included cells 
meta_cells = pd.read_csv(data_dir + "included_cell_barcode.csv", index_col = 0)
# rename the index of meta cells
new_idx = []
for idx in meta_cells.index.values:
    if idx.split("-")[1] == "1_1":
        new_idx.append("M1385:" + idx.split("-")[0] + "x")
    elif idx.split("-")[1] == "1_2":
        new_idx.append("M1386:" + idx.split("-")[0] + "x")
    elif idx.split("-")[1] == "1_3":
        new_idx.append("M1387:" + idx.split("-")[0] + "x")
    elif idx.split("-")[1] == "1_4":
        new_idx.append("M1415:" + idx.split("-")[0] + "x")
    elif idx.split("-")[1] == "1_5":
        new_idx.append("M1418:" + idx.split("-")[0] + "x")
    elif idx.split("-")[1] == "1_6":
        new_idx.append("M1475:" + idx.split("-")[0] + "x")
new_idx = np.array(new_idx)

meta_cells.index = new_idx
# Merge adata
adata_merge = sc.concat([M1385_adata, M1386_adata, M1387_adata, M1415_adata, M1418_adata, M1475_adata], join = "inner")
# select according to the meta data
adata_merge = adata_merge[meta_cells.index.values,:]


# %%
