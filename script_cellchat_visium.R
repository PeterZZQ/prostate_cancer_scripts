rm(list=ls())
gc()

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)

setwd("/localscratch/ziqi/Prostate_Cancer/")

# 12WK double knock-out
visium_pp12 <- Load10X_Spatial(data.dir = "spatial_data/Visium/225_PP12", filename = "filtered_feature_bc_matrix.h5")
visium_pp18 <- Load10X_Spatial(data.dir = "spatial_data/Visium/1687_PP18", filename = "filtered_feature_bc_matrix.h5")
visium_ppc12 <- Load10X_Spatial(data.dir = "spatial_data/Visium/1161_PPC", filename = "filtered_feature_bc_matrix.h5")
visium_ppc18 <- Load10X_Spatial(data.dir = "spatial_data/Visium/1660_PPC18", filename = "filtered_feature_bc_matrix.h5")

meta.pp12 <- read.table("spatial_data/Visium/transfer_labels_225_pp12.csv")
meta.pp18 <- read.table("spatial_data/Visium/transfer_labels_1687_pp18.csv")
meta.ppc12 <- read.table("spatial_data/Visium/transfer_labels_1161_ppc.csv")
meta.ppc18 <- read.table("spatial_data/Visium/transfer_labels_1660_ppc18.csv")

visium_pp12@meta.data["predict.label"] <- meta.pp12["predict.label"]
visium_pp18@meta.data["predict.label"] <- meta.pp18["predict.label"]
visium_ppc12@meta.data["predict.label"] <- meta.ppc12["predict.label"]
visium_ppc18@meta.data["predict.label"] <- meta.ppc18["predict.label"]

visium_pp12 <- SCTransform(visium_pp12, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
visium_pp18 <- SCTransform(visium_pp18, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
visium_ppc12 <- SCTransform(visium_ppc12, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
visium_ppc18 <- SCTransform(visium_ppc18, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

# Prepare input data for CelChat analysis
data.pp12 = GetAssayData(visium_pp12, slot = "data", assay = "SCT") # normalized data matrix
data.pp18 = GetAssayData(visium_pp18, slot = "data", assay = "SCT") # normalized data matrix
data.ppc12 = GetAssayData(visium_ppc12, slot = "data", assay = "SCT") # normalized data matrix
data.ppc18 = GetAssayData(visium_ppc18, slot = "data", assay = "SCT") # normalized data matrix

meta.pp12 = data.frame(labels = visium_pp12@meta.data$predict.label, row.names = rownames(visium_pp12@meta.data)) # manually create a dataframe consisting of the cell labels
meta.pp18 = data.frame(labels = visium_pp18@meta.data$predict.label, row.names = rownames(visium_pp18@meta.data)) # manually create a dataframe consisting of the cell labels
meta.ppc12 = data.frame(labels = visium_ppc12@meta.data$predict.label, row.names = rownames(visium_ppc12@meta.data)) # manually create a dataframe consisting of the cell labels
meta.ppc18 = data.frame(labels = visium_ppc18@meta.data$predict.label, row.names = rownames(visium_ppc18@meta.data)) # manually create a dataframe consisting of the cell labels

# ------------------------------------------------------------------------------
#
# Run CellChat
#
# ------------------------------------------------------------------------------
# load spatial imaging information
# Spatial locations of spots from full (NOT high/low) resolution images are required
spatial.locs.pp12 = GetTissueCoordinates(visium_pp12, scale = NULL, cols = c("imagerow", "imagecol")) 
# Scale factors and spot diameters of the full resolution images 
# USER can also extract scale factors from a Seurat object, but the `spot` value here is different from the one in Seurat. Thus, USER still needs to get the `spot` value from the json file. 
scale.factors.pp12 = jsonlite::fromJSON(txt = file.path("/localscratch/ziqi/Prostate_Cancer/spatial_data/Visium/225_PP12/spatial/", 'scalefactors_json.json'))
scale.factors.pp12 = list(spot.diameter = 65, spot = scale.factors.pp12$spot_diameter_fullres, # these two information are required
                          fiducial = scale.factors.pp12$fiducial_diameter_fullres, hires = scale.factors.pp12$tissue_hires_scalef, lowres = scale.factors.pp12$tissue_lowres_scalef # these three information are not required
)
# create seurat object
cellchat.pp12 <- createCellChat(object = data.pp12, meta = meta.pp12, group.by = "labels", datatype = "spatial", coordinates = spatial.locs.pp12, scale.factors = scale.factors.pp12)
# load database
CellChatDB <- CellChatDB.mouse
interestDB <- read.table("spatial_data/ligand_receptor_reformat.csv", sep = ",", row.names = 1)
CellChatDB <- subsetDB(CellChatDB, search = interestDB$CellChat.interaction.name, key = "CellChat interaction name")
# failed to discover signaling pathway
# cellchat.pp12@DB <- CellChatDB
cellchat.pp12@DB <- CellChatDB.mouse
cellchat.pp12 <- subsetData(cellchat.pp12)
cellchat.pp12 <- identifyOverExpressedGenes(cellchat.pp12)
cellchat.pp12 <- identifyOverExpressedInteractions(cellchat.pp12)
cellchat.pp12 <- computeCommunProb(cellchat.pp12, type = "truncatedMean", trim = 0.1, distance.use = TRUE, interaction.length = 200, scale.distance = 0.01)
# filter communication
cellchat.pp12 <- filterCommunication(cellchat.pp12, min.cells = 10)
cellchat.pp12 <- computeCommunProbPathway(cellchat.pp12)
cellchat.pp12 <- aggregateNet(cellchat.pp12)
cci.pp12 <- subsetCommunication(cellchat.pp12)
cciP.pp12 <- subsetCommunication(cellchat.pp12, slot.name = "netP")
saveRDS(cellchat.pp12, file = "results_cellchat_visium/cellchat_pp12.rds")
write.table(cci.pp12, "results_cellchat_visium/cci_pp12.csv", sep = "\t", quote = F)
write.table(cciP.pp12, "results_cellchat_visium/cciP_pp12.csv", sep = "\t", quote = F)

netVisual_bubble(cellchat.pp12, sources.use = c( "Luminal 1", "Mesenchymal"), targets.use = c( "Luminal 1", "Mesenchymal"), signaling = c("CCL","CXCL"), remove.isolate = FALSE)


spatial.locs.pp18 = GetTissueCoordinates(visium_pp18, scale = NULL, cols = c("imagerow", "imagecol")) 
# Scale factors and spot diameters of the full resolution images 
# USER can also extract scale factors from a Seurat object, but the `spot` value here is different from the one in Seurat. Thus, USER still needs to get the `spot` value from the json file. 
scale.factors.pp18 = jsonlite::fromJSON(txt = file.path("/localscratch/ziqi/Prostate_Cancer/spatial_data/Visium/1687_PP18/spatial/", 'scalefactors_json.json'))
scale.factors.pp18 = list(spot.diameter = 65, spot = scale.factors.pp18$spot_diameter_fullres, # these two information are required
                          fiducial = scale.factors.pp18$fiducial_diameter_fullres, hires = scale.factors.pp18$tissue_hires_scalef, lowres = scale.factors.pp18$tissue_lowres_scalef # these three information are not required
)
# create seurat object
cellchat.pp18 <- createCellChat(object = data.pp18, meta = meta.pp18, group.by = "labels", datatype = "spatial", coordinates = spatial.locs.pp18, scale.factors = scale.factors.pp18)
# load database
CellChatDB <- CellChatDB.mouse
cellchat.pp18@DB <- CellChatDB
cellchat.pp18 <- subsetData(cellchat.pp18)
cellchat.pp18 <- identifyOverExpressedGenes(cellchat.pp18)
cellchat.pp18 <- identifyOverExpressedInteractions(cellchat.pp18)
cellchat.pp18 <- computeCommunProb(cellchat.pp18, type = "truncatedMean", trim = 0.1, distance.use = TRUE, interaction.length = 200, scale.distance = 0.01)
# filter communication
cellchat.pp18 <- filterCommunication(cellchat.pp18, min.cells = 10)
cellchat.pp18 <- computeCommunProbPathway(cellchat.pp18)
cellchat.pp18 <- aggregateNet(cellchat.pp18)
cci.pp18 <- subsetCommunication(cellchat.pp18)
cciP.pp18 <- subsetCommunication(cellchat.pp18, slot.name = "netP")
saveRDS(cellchat.pp18, file = "results_cellchat_visium/cellchat_pp18.rds")
write.table(cci.pp18, "results_cellchat_visium/cci_pp18.csv", sep = "\t", quote = F)
write.table(cciP.pp18, "results_cellchat_visium/cciP_pp18.csv", sep = "\t", quote = F)

# netVisual_bubble(cellchat.pp18, sources.use = c( "Luminal 1", "Mesenchymal", "Club epithelia (Luminal)", "SV, Spink1+"), targets.use = c( "Luminal 1", "Mesenchymal", "Club epithelia (Luminal)", "SV, Spink1+"), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
netVisual_chord_gene(cellchat.pp18, sources.use = "Luminal 1", targets.use = c( "Luminal 1", "Mesenchymal"), lab.cex = 0.5, legend.pos.y = 30)
netVisual_chord_gene(cellchat.pp18, sources.use = "Mesenchymal", targets.use = c( "Luminal 1", "Mesenchymal"), lab.cex = 0.5, legend.pos.y = 30)
netVisual_chord_gene(cellchat.pp18, sources.use = "SV, Spink1+", targets.use = c( "Luminal 1", "Mesenchymal", "SV, Spink1+"), lab.cex = 0.5, legend.pos.y = 30)


spatial.locs.ppc12 = GetTissueCoordinates(visium_ppc12, scale = NULL, cols = c("imagerow", "imagecol")) 
# Scale factors and spot diameters of the full resolution images 
# USER can also extract scale factors from a Seurat object, but the `spot` value here is different from the one in Seurat. Thus, USER still needs to get the `spot` value from the json file. 
scale.factors.ppc12 = jsonlite::fromJSON(txt = file.path("/localscratch/ziqi/Prostate_Cancer/spatial_data/Visium/1161_PPC/spatial/", 'scalefactors_json.json'))
scale.factors.ppc12 = list(spot.diameter = 65, spot = scale.factors.ppc12$spot_diameter_fullres, # these two information are required
                          fiducial = scale.factors.ppc12$fiducial_diameter_fullres, hires = scale.factors.ppc12$tissue_hires_scalef, lowres = scale.factors.ppc12$tissue_lowres_scalef # these three information are not required
)
# create seurat object
cellchat.ppc12 <- createCellChat(object = data.ppc12, meta = meta.ppc12, group.by = "labels", datatype = "spatial", coordinates = spatial.locs.ppc12, scale.factors = scale.factors.ppc12)
# load database
CellChatDB <- CellChatDB.mouse
cellchat.ppc12@DB <- CellChatDB
cellchat.ppc12 <- subsetData(cellchat.ppc12)
cellchat.ppc12 <- identifyOverExpressedGenes(cellchat.ppc12)
cellchat.ppc12 <- identifyOverExpressedInteractions(cellchat.ppc12)
cellchat.ppc12 <- computeCommunProb(cellchat.ppc12, type = "truncatedMean", trim = 0.1, distance.use = TRUE, interaction.length = 200, scale.distance = 0.01)
# filter communication
cellchat.ppc12 <- filterCommunication(cellchat.ppc12, min.cells = 10)
cellchat.ppc12 <- computeCommunProbPathway(cellchat.ppc12)
cellchat.ppc12 <- aggregateNet(cellchat.ppc12)
cci.ppc12 <- subsetCommunication(cellchat.ppc12)
cciP.ppc12 <- subsetCommunication(cellchat.ppc12, slot.name = "netP")
saveRDS(cellchat.ppc12, file = "results_cellchat_visium/cellchat_ppc12.rds")
write.table(cci.ppc12, "results_cellchat_visium/cci_ppc12.csv", sep = "\t", quote = F)
write.table(cciP.ppc12, "results_cellchat_visium/cciP_ppc12.csv", sep = "\t", quote = F)

spatial.locs.ppc18 = GetTissueCoordinates(visium_ppc18, scale = NULL, cols = c("imagerow", "imagecol")) 
# Scale factors and spot diameters of the full resolution images 
# USER can also extract scale factors from a Seurat object, but the `spot` value here is different from the one in Seurat. Thus, USER still needs to get the `spot` value from the json file. 
scale.factors.ppc18 = jsonlite::fromJSON(txt = file.path("/localscratch/ziqi/Prostate_Cancer/spatial_data/Visium/1660_PPC18/spatial/", 'scalefactors_json.json'))
scale.factors.ppc18 = list(spot.diameter = 65, spot = scale.factors.ppc18$spot_diameter_fullres, # these two information are required
                           fiducial = scale.factors.ppc18$fiducial_diameter_fullres, hires = scale.factors.ppc18$tissue_hires_scalef, lowres = scale.factors.ppc18$tissue_lowres_scalef # these three information are not required
)
# create seurat object
cellchat.ppc18 <- createCellChat(object = data.ppc18, meta = meta.ppc18, group.by = "labels", datatype = "spatial", coordinates = spatial.locs.ppc18, scale.factors = scale.factors.ppc18)
# load database
CellChatDB <- CellChatDB.mouse
cellchat.ppc18@DB <- CellChatDB
cellchat.ppc18 <- subsetData(cellchat.ppc18)
cellchat.ppc18 <- identifyOverExpressedGenes(cellchat.ppc18)
cellchat.ppc18 <- identifyOverExpressedInteractions(cellchat.ppc18)
cellchat.ppc18 <- computeCommunProb(cellchat.ppc18, type = "truncatedMean", trim = 0.1, distance.use = TRUE, interaction.length = 200, scale.distance = 0.01)
# filter communication
cellchat.ppc18 <- filterCommunication(cellchat.ppc18, min.cells = 10)
cellchat.ppc18 <- computeCommunProbPathway(cellchat.ppc18)
cellchat.ppc18 <- aggregateNet(cellchat.ppc18)
cci.ppc18 <- subsetCommunication(cellchat.ppc18)
cciP.ppc18 <- subsetCommunication(cellchat.ppc18, slot.name = "netP")
saveRDS(cellchat.ppc18, file = "results_cellchat_visium/cellchat_ppc18.rds")
write.table(cci.ppc18, "results_cellchat_visium/cci_ppc18.csv", sep = "\t", quote = F)
write.table(cciP.ppc18, "results_cellchat_visium/cciP_ppc18.csv", sep = "\t", quote = F)





