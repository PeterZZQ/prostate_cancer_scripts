rm(list=ls())
gc()

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

prostate.combined <- LoadH5Seurat("Cell_Ranger_output/seurat_integrate.h5Seurat")
# annotate cell types
prostate.combined <- RenameIdents(prostate.combined, `0` = "Luminal 1", `1` = "Luminal 1", `2` = "Hillock epithelia (Basal)",
                                  `3` = "Club epithelia (Luminal)", `4` = "Macrophage (Myeloid, Total immune)", `5` = "Monocytes (Myeloid, Total immune)", 
                                  `6` = "Mesenchymal", `7` = "Spink1+ (Luminal)", `8` = "Basal", `9` = "Lymphoid (Total immune)",
                                  `10` = "SV, Spink1+", `11` = "Endothelial", `12` = "Luminal 2", `13` = "Luminal 3")

# remove Luminal 2 & 3
prostate.combined <- prostate.combined[,(prostate.combined@active.ident != "Luminal 2")&(prostate.combined@active.ident != "Luminal 3")]

prostate.M1417_PP12 <- prostate.combined[, prostate.combined@meta.data["sample"] == "M1417-PP12"]
prostate.M1416_PP18 <- prostate.combined[, prostate.combined@meta.data["sample"] == "M1416-PP18"]
prostate.M1436_PPC12 <- prostate.combined[, prostate.combined@meta.data["sample"] == "M1436-PPC12"]
prostate.M1437_PPC18 <- prostate.combined[, prostate.combined@meta.data["sample"] == "M1437-PPC18"]

# log-normalized count
counts.norm <- prostate.M1417_PP12@assays$RNA@data
# meta-information
meta.cells <- prostate.M1417_PP12@meta.data
meta.cells["annotations"] <- prostate.M1417_PP12@active.ident
# create cellchat
cellchat.M1417_PP12 <- createCellChat(object = counts.norm, meta = meta.cells, group.by = "annotations")

counts.norm <- prostate.M1416_PP18@assays$RNA@data
meta.cells <- prostate.M1416_PP18@meta.data
meta.cells["annotations"] <- prostate.M1416_PP18@active.ident
cellchat.M1416_PP18 <- createCellChat(object = counts.norm, meta = meta.cells, group.by = "annotations")

counts.norm <- prostate.M1436_PPC12@assays$RNA@data
meta.cells <- prostate.M1436_PPC12@meta.data
meta.cells["annotations"] <- prostate.M1436_PPC12@active.ident
cellchat.M1436_PPC12 <- createCellChat(object = counts.norm, meta = meta.cells, group.by = "annotations")

counts.norm <- prostate.M1437_PPC18@assays$RNA@data
meta.cells <- prostate.M1437_PPC18@meta.data
meta.cells["annotations"] <- prostate.M1437_PPC18@active.ident
cellchat.M1437_PPC18 <- createCellChat(object = counts.norm, meta = meta.cells, group.by = "annotations")


# ellChatDB is a manually curated database of literature-supported ligand-receptor interactions in both human and mouse. 
# CellChatDB in mouse contains 2,021 validated molecular interactions, including 60% of secrete autocrine/paracrine signaling interactions, 
# 21% of extracellular matrix (ECM)-receptor interactions and 19% of cell-cell contact interactions.
CellChatDB <- CellChatDB.mouse
# showDatabaseCategory(CellChatDB)
# optional: select sub-set of CellChatDB
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat.M1417_PP12@DB <- CellChatDB
# subset the expression data of signaling genes for saving computation cost
cellchat.M1417_PP12 <- subsetData(cellchat.M1417_PP12) # This step is necessary even if using the whole database
# construct CCI using over-expression of ligand & target
cellchat.M1417_PP12 <- identifyOverExpressedGenes(cellchat.M1417_PP12)
cellchat.M1417_PP12 <- identifyOverExpressedInteractions(cellchat.M1417_PP12)
# Compute the communication probability and infer cellular communication network
# cellchat.M1417_PP12 <- computeCommunProb(cellchat.M1417_PP12, type =  "truncatedMean", trim = 0.1, population.size = TRUE)
cellchat.M1417_PP12 <- computeCommunProb(cellchat.M1417_PP12)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.M1417_PP12 <- filterCommunication(cellchat.M1417_PP12, min.cells = 10)
# compute signal pathway
cellchat.M1417_PP12 <- computeCommunProbPathway(cellchat.M1417_PP12)
# visualization
cellchat.M1417_PP12 <- aggregateNet(cellchat.M1417_PP12)
# save cellchat object
saveRDS(cellchat.M1417_PP12, file = "results_cellchat/cellchat_M1417_PP12.rds")
# Extract the inferred cellular communication network as a data frame
cci.M1417_PP12 <- subsetCommunication(cellchat.M1417_PP12)
# Extract the inferred interacting pathway
cciP.M1417_PP12 <- subsetCommunication(cellchat.M1417_PP12, slot.name = "netP")
write.table(cci.M1417_PP12, "results_cellchat/cci_M1417_PP12.csv", sep = ",", quote = F)
write.table(cciP.M1417_PP12, "results_cellchat/cciP_M1417_PP12.csv", sep = ",", quote = F)

# Running cellchat for sample M1416_PP18
cellchat.M1416_PP18@DB <- CellChatDB
cellchat.M1416_PP18 <- subsetData(cellchat.M1416_PP18) 
cellchat.M1416_PP18 <- identifyOverExpressedGenes(cellchat.M1416_PP18)
cellchat.M1416_PP18 <- identifyOverExpressedInteractions(cellchat.M1416_PP18)
cellchat.M1416_PP18 <- computeCommunProb(cellchat.M1416_PP18)
cellchat.M1416_PP18 <- filterCommunication(cellchat.M1416_PP18, min.cells = 10)
cellchat.M1416_PP18 <- computeCommunProbPathway(cellchat.M1416_PP18)
cellchat.M1416_PP18 <- aggregateNet(cellchat.M1416_PP18)
cci.M1416_PP18 <- subsetCommunication(cellchat.M1416_PP18)
cciP.M1416_PP18 <- subsetCommunication(cellchat.M1416_PP18, slot.name = "netP")
saveRDS(cellchat.M1416_PP18, file = "results_cellchat/cellchat_M1416_PP18.rds")
write.table(cci.M1416_PP18, "results_cellchat/cci_M1416_PP18.csv", sep = ",", quote = F)
write.table(cciP.M1416_PP18, "results_cellchat/cciP_M1416_PP18.csv", sep = ",", quote = F)

# Running cellchat for sample M1416_PP18
cellchat.M1436_PPC12@DB <- CellChatDB
cellchat.M1436_PPC12 <- subsetData(cellchat.M1436_PPC12) 
cellchat.M1436_PPC12 <- identifyOverExpressedGenes(cellchat.M1436_PPC12)
cellchat.M1436_PPC12 <- identifyOverExpressedInteractions(cellchat.M1436_PPC12)
cellchat.M1436_PPC12 <- computeCommunProb(cellchat.M1436_PPC12)
cellchat.M1436_PPC12 <- filterCommunication(cellchat.M1436_PPC12, min.cells = 10)
cellchat.M1436_PPC12 <- computeCommunProbPathway(cellchat.M1436_PPC12)
cellchat.M1436_PPC12 <- aggregateNet(cellchat.M1436_PPC12)
cci.M1436_PPC12 <- subsetCommunication(cellchat.M1436_PPC12)
cciP.M1436_PPC12 <- subsetCommunication(cellchat.M1436_PPC12, slot.name = "netP")
saveRDS(cellchat.M1436_PPC12, file = "results_cellchat/cellchat_M1436_PPC12.rds")
write.table(cci.M1436_PPC12, "results_cellchat/cci_M1436_PPC12.csv", sep = ",", quote = F)
write.table(cciP.M1436_PPC12, "results_cellchat/cciP_M1436_PPC12.csv", sep = ",", quote = F)


# Running cellchat for sample M1416_PP18
cellchat.M1437_PPC18@DB <- CellChatDB
cellchat.M1437_PPC18 <- subsetData(cellchat.M1437_PPC18) 
cellchat.M1437_PPC18 <- identifyOverExpressedGenes(cellchat.M1437_PPC18)
cellchat.M1437_PPC18 <- identifyOverExpressedInteractions(cellchat.M1437_PPC18)
cellchat.M1437_PPC18 <- computeCommunProb(cellchat.M1437_PPC18)
cellchat.M1437_PPC18 <- filterCommunication(cellchat.M1437_PPC18, min.cells = 10)
cellchat.M1437_PPC18 <- computeCommunProbPathway(cellchat.M1437_PPC18)
cellchat.M1437_PPC18 <- aggregateNet(cellchat.M1437_PPC18)
cci.M1437_PPC18 <- subsetCommunication(cellchat.M1437_PPC18)
cciP.M1437_PPC18 <- subsetCommunication(cellchat.M1437_PPC18, slot.name = "netP")
saveRDS(cellchat.M1437_PPC18, file = "results_cellchat/cellchat_M1437_PPC18.rds")
write.table(cci.M1437_PPC18, "results_cellchat/cci_M1437_PPC18.csv", sep = ",", quote = F)
write.table(cciP.M1437_PPC18, "results_cellchat/cciP_M1437_PPC18.csv", sep = ",", quote = F)


# Compare across conditions
cellchat.M1417_PP12 <- readRDS("results_cellchat/cellchat_M1417_PP12.rds")
cellchat.M1436_PPC12 <- readRDS("results_cellchat/cellchat_M1436_PPC12.rds")
cellchat.M1417_PP12 <- updateCellChat(cellchat.M1417_PP12)
cellchat.M1436_PPC12 <- updateCellChat(cellchat.M1436_PPC12)
object.list <- list(M1417_PP12 = cellchat.M1417_PP12, M1436_PPC12 = cellchat.M1436_PPC12)
cellchat.12wk <- mergeCellChat(object.list, add.names = names(object.list))

cellchat.M1416_PP18 <- readRDS("results_cellchat/cellchat_M1416_PP18.rds")
cellchat.M1437_PPC18 <- readRDS("results_cellchat/cellchat_M1437_PPC18.rds")
cellchat.M1416_PP18 <- updateCellChat(cellchat.M1416_PP18)
cellchat.M1437_PPC18 <- updateCellChat(cellchat.M1437_PPC18)
# object.list <- list(M1416_PP18 = cellchat.M1416_PP18, M1437_PPC18 = cellchat.M1437_PPC18)
object.list <- list(M1416_PP18 = cellchat.M1416_PP18, M1437_PPC18 = cellchat.M1437_PPC18)
cellchat.18wk <- mergeCellChat(object.list, add.names = names(object.list))


# Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat.12wk, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat.12wk, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

# Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat.18wk, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat.18wk, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


# Change of interaction across conditions for different pairs of cell-types
gg1 <- netVisual_heatmap(cellchat.12wk)
gg2 <- netVisual_heatmap(cellchat.12wk, measure = "weight")
gg1 + gg2

gg1 <- netVisual_heatmap(cellchat.18wk)
gg2 <- netVisual_heatmap(cellchat.18wk, measure = "weight")
gg1 + gg2

gg1 <- netVisual_heatmap(cellchat.18wk, slot.name = "netP")
gg2 <- netVisual_heatmap(cellchat.18wk, slot.name = "netP", measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat.18wk, weight.scale = T)
netVisual_diffInteraction(cellchat.18wk, weight.scale = T, measure = "weight")

