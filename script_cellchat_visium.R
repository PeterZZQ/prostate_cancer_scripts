rm(list=ls())
gc()

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)

setwd("/localscratch/ziqi/Prostate_Cancer/")

visium_1161_ppc <- Load10X_Spatial(data.dir = "spatial_data/Visium/1161_PPC", filename = "filtered_feature_bc_matrix.h5")
SpatialDimPlot(visium_1161_ppc, label = F)
data.input <- GetAssayData(visium_1161_ppc, slot = "data", assay = "Spatial") # normalized data matrix
# have to annotate cell types using Seurat, etc
# meta = data.frame(labels = Idents(visium.brain), row.names = names(Idents(visium.brain))) # manually create a dataframe consisting of the cell labels
