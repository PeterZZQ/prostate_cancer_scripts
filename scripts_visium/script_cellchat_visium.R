rm(list=ls())
gc()

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)

setwd("/localscratch/ziqi/Prostate_Cancer/")


# 12WK double knock-out
visium_pp12 <- Load10X_Spatial(data.dir = "dataset/data_visium/spatial_data/Visium/225_PP12", filename = "filtered_feature_bc_matrix.h5")
visium_pp18 <- Load10X_Spatial(data.dir = "dataset/data_visium/spatial_data/Visium/1687_PP18", filename = "filtered_feature_bc_matrix.h5")
visium_ppc12 <- Load10X_Spatial(data.dir = "dataset/data_visium/spatial_data/Visium/1161_PPC", filename = "filtered_feature_bc_matrix.h5")
visium_ppc18 <- Load10X_Spatial(data.dir = "dataset/data_visium/spatial_data/Visium/1660_PPC18", filename = "filtered_feature_bc_matrix.h5")

# card deconvolution result
deconv_score_pp12 <- read.table("results_visium/results_card/deconv_score_pp12.csv", sep = "\t")
deconv_score_pp18 <- read.table("results_visium/results_card/deconv_score_pp18.csv", sep = "\t")
deconv_score_ppc12 <- read.table("results_visium/results_card/deconv_score_ppc12.csv", sep = "\t")
deconv_score_ppc18 <- read.table("results_visium/results_card/deconv_score_ppc18.csv", sep = "\t")

max.label.pp12 <- names(deconv_score_pp12)[max.col(deconv_score_pp12)]
max.label.pp18 <- names(deconv_score_pp18)[max.col(deconv_score_pp18)]
max.label.ppc12 <- names(deconv_score_ppc12)[max.col(deconv_score_ppc12)]
max.label.ppc18 <- names(deconv_score_ppc18)[max.col(deconv_score_ppc18)]

visium_pp12@meta.data["predict.label"] <- max.label.pp12
visium_pp18@meta.data["predict.label"] <- max.label.pp18
visium_ppc12@meta.data["predict.label"] <- max.label.ppc12
visium_ppc18@meta.data["predict.label"] <- max.label.ppc18

visium_pp12 <- SCTransform(visium_pp12, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
visium_pp18 <- SCTransform(visium_pp18, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
visium_ppc12 <- SCTransform(visium_ppc12, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
visium_ppc18 <- SCTransform(visium_ppc18, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)

# Prepare input data for CelChat analysis
data.pp12 <- GetAssayData(visium_pp12, slot = "data", assay = "SCT") # normalized data matrix
data.pp18 <- GetAssayData(visium_pp18, slot = "data", assay = "SCT") # normalized data matrix
data.ppc12 <- GetAssayData(visium_ppc12, slot = "data", assay = "SCT") # normalized data matrix
data.ppc18 <- GetAssayData(visium_ppc18, slot = "data", assay = "SCT") # normalized data matrix

meta.pp12 <- data.frame(labels = visium_pp12@meta.data$predict.label, row.names = rownames(visium_pp12@meta.data)) # manually create a dataframe consisting of the cell labels
meta.pp18 <- data.frame(labels = visium_pp18@meta.data$predict.label, row.names = rownames(visium_pp18@meta.data)) # manually create a dataframe consisting of the cell labels
meta.ppc12 <- data.frame(labels = visium_ppc12@meta.data$predict.label, row.names = rownames(visium_ppc12@meta.data)) # manually create a dataframe consisting of the cell labels
meta.ppc18 <- data.frame(labels = visium_ppc18@meta.data$predict.label, row.names = rownames(visium_ppc18@meta.data)) # manually create a dataframe consisting of the cell labels

meta.pp12["samples"] <- "PP12"
meta.pp18["samples"] <- "PP18"
meta.ppc12["samples"] <- "PPC12"
meta.ppc18["samples"] <- "PPC18"

# ------------------------------------------------------------------------------
#
# Run CellChat
#
# ------------------------------------------------------------------------------
visium_list <- list("225_PP12" = visium_pp12, "1687_PP18" = visium_pp18, "1161_PPC" = visium_ppc12, "1660_PPC18" = visium_ppc18)
data_list <- list("225_PP12" = data.pp12, "1687_PP18" = data.pp18, "1161_PPC" = data.ppc12, "1660_PPC18" = data.ppc18)
meta_list <- list("225_PP12" = meta.pp12, "1687_PP18" = meta.pp18, "1161_PPC" = meta.ppc12, "1660_PPC18" = meta.ppc18)

# load database
CellChatDB <- CellChatDB.mouse
# interestDB <- read.table("dataset/data_visium/spatial_data/ligand_receptor_reformat.csv", sep = ",", row.names = 1)
# CellChatDB <- subsetDB(CellChatDB, search = interestDB$CellChat.interaction.name, key = "CellChat interaction name")
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

for(sample in c("225_PP12", "1687_PP18", "1161_PPC", "1660_PPC18")){
  visium_sample <- visium_list[[sample]]
  # load spatial imaging information
  # Spatial locations of spots from full (NOT high/low) resolution images are required
  spatial.locs.sample <- GetTissueCoordinates(visium_sample, scale = NULL, cols = c("imagerow", "imagecol")) 
  # Scale factors and spot diameters of the full resolution images 
  # USER can also extract scale factors from a Seurat object, but the `spot` value here is different from the one in Seurat. Thus, USER still needs to get the `spot` value from the json file. 
  scale.factors.sample <- jsonlite::fromJSON(txt = file.path(paste0("/localscratch/ziqi/Prostate_Cancer/dataset/data_visium/spatial_data/Visium/", sample, "/spatial"), 'scalefactors_json.json'))
  spot.size = 65 # the theoretical spot size (um) in 10X Visium
  conversion.factor.sample <- spot.size/scale.factors.sample$spot_diameter_fullres
  spatial.factors.sample = data.frame(ratio = conversion.factor.sample, tol = spot.size/2)
  d.spatial <- computeCellDistance(coordinates = spatial.locs.sample, ratio = spatial.factors.sample$ratio, tol = spatial.factors.sample$tol)
  print(min(d.spatial[d.spatial!=0]))

  # create cellchat object
  cellchat.sample <- createCellChat(object = data_list[[sample]], meta = meta_list[[sample]], group.by = "labels", datatype = "spatial", coordinates = spatial.locs.sample, spatial.factors = spatial.factors.sample)
  
  # failed to discover signaling pathway
  cellchat.sample@DB <- CellChatDB.use
  cellchat.sample <- subsetData(cellchat.sample)
  future::plan("multisession", workers = 4) 
  cellchat.sample <- identifyOverExpressedGenes(cellchat.sample)
  cellchat.sample <- identifyOverExpressedInteractions(cellchat.sample)
  cellchat.sample <- computeCommunProb(cellchat.sample, type = "truncatedMean", trim = 0.1, distance.use = TRUE, interaction.range = 250, scale.distance = 0.01, contact.dependent = TRUE, contact.range = 100)
  
  # filter communication
  cellchat.sample <- filterCommunication(cellchat.sample, min.cells = 10)
  cellchat.sample <- computeCommunProbPathway(cellchat.sample)
  cellchat.sample <- aggregateNet(cellchat.sample)
  cci.sample <- subsetCommunication(cellchat.sample)
  cciP.sample <- subsetCommunication(cellchat.sample, slot.name = "netP")
  saveRDS(cellchat.sample, file = paste0("results_visium/results_cellchat/cellchat_", sample, ".rds"))
  write.table(cci.sample, paste0("results_visium/results_cellchat/cci_", sample, ".csv"), sep = "\t", quote = F)
  write.table(cciP.sample, paste0("results_visium/results_cellchat/cciP_", sample, ".csv"), sep = "\t", quote = F)
  
}



