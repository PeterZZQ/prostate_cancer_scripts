rm(list=ls())
gc()

library(CellChat)
library(anndata)
library(patchwork)
options(stringsAsFactors = FALSE)

setwd("/localscratch/ziqi/Prostate_Cancer/")
# read in the annotated intact scRNA-seq dataset
adata.intact <- read_h5ad("dataset/data_scrnaseq/seurat_integration/adata_intact_seurat.h5ad")
# read in the raw count matrix
counts.intact <- Matrix::readMM("dataset/data_scrnaseq/data/qc_data/X_intact.mtx")
counts.intact <- as(t(as.matrix(counts.intact)), "dgCMatrix")
colnames(counts.intact) <- rownames(adata.intact$obs)
rownames(counts.intact) <- rownames(adata.intact$var)
meta.intact <- adata.intact$obs
# normalize the raw count data
library.size <- Matrix::colSums(counts.intact)
counts.intact.norm <- as(log1p(Matrix::t(Matrix::t(counts.intact)/library.size) * 10000), "dgCMatrix")
# separate the count matrix into different samples
counts.M1417_PP12 <- counts.intact.norm[, as.vector(meta.intact["sample"] == "M1417-PP12")]
counts.M1436_PPC12 <- counts.intact.norm[, as.vector(meta.intact["sample"] == "M1436-PPC12")]
counts.M1416_PP18 <- counts.intact.norm[, as.vector(meta.intact["sample"] == "M1416-PP18")]
counts.M1437_PPC18 <- counts.intact.norm[, as.vector(meta.intact["sample"] == "M1437-PPC18")]
# separate the meta data into different samples
meta.M1417_PP12 <- meta.intact[meta.intact["sample"] == "M1417-PP12",]
meta.M1436_PPC12 <- meta.intact[meta.intact["sample"] == "M1436-PPC12",]
meta.M1416_PP18 <- meta.intact[meta.intact["sample"] == "M1416-PP18",]
meta.M1437_PPC18 <- meta.intact[meta.intact["sample"] == "M1437-PPC18",]
names(meta.M1417_PP12)[1] <- "samples"
names(meta.M1436_PPC12)[1] <- "samples"
names(meta.M1416_PP18)[1] <- "samples"
names(meta.M1437_PPC18)[1] <- "samples"

# create cellchat
cellchat.M1417_PP12 <- createCellChat(object = counts.M1417_PP12, meta = meta.M1417_PP12, group.by = "annot")
cellchat.M1416_PP18 <- createCellChat(object = counts.M1416_PP18, meta = meta.M1416_PP18, group.by = "annot")
cellchat.M1436_PPC12 <- createCellChat(object = counts.M1436_PPC12, meta = meta.M1436_PPC12, group.by = "annot")
cellchat.M1437_PPC18 <- createCellChat(object = counts.M1437_PPC18, meta = meta.M1437_PPC18, group.by = "annot")

cellchat.list <- list("PP12" = cellchat.M1417_PP12, "PP18" = cellchat.M1416_PP18, "PPC12" = cellchat.M1436_PPC12, "PPC18" = cellchat.M1437_PPC18) 

# CellChatDB is a manually curated database of literature-supported ligand-receptor interactions in both human and mouse. 
# CellChatDB in mouse contains 2,021 validated molecular interactions, including 60% of secrete autocrine/paracrine signaling interactions, 
# 21% of extracellular matrix (ECM)-receptor interactions and 19% of cell-cell contact interactions.
CellChatDB <- CellChatDB.mouse
# showDatabaseCategory(CellChatDB)
# optional: select sub-set of CellChatDB
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling


for(sample in names(cellchat.list)){
  cellchat <- cellchat.list[[sample]]
  cellchat@DB <- CellChatDB.use
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat)
  # do parallel
  future::plan("multisession", workers = 4)
  # construct CCI using over-expression of ligand & target
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # Compute the communication probability and infer cellular communication network
  # cellchat <- computeCommunProb(cellchat, type =  "truncatedMean", trim = 0.1, population.size = TRUE)
  cellchat <- computeCommunProb(cellchat)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  # compute signal pathway
  cellchat <- computeCommunProbPathway(cellchat)
  # visualization
  cellchat <- aggregateNet(cellchat)
  # save cellchat object
  saveRDS(cellchat, file = paste0("results_scrnaseq/results_cellchat/cellchat_", sample, ".rds"))
  # Extract the inferred cellular communication network as a data frame
  cci <- subsetCommunication(cellchat)
  # Extract the inferred interacting pathway
  cciP <- subsetCommunication(cellchat, slot.name = "netP")
  write.table(cci, paste0("results_scrnaseq/results_cellchat/cci_", sample, ".csv"), sep = ",", quote = F)
  write.table(cciP, paste0("results_scrnaseq/results_cellchat/cciP_", sample, ".csv"), sep = ",", quote = F)
  
}


# # Compare across conditions
# cellchat.M1417_PP12 <- readRDS("results_cellchat/cellchat_M1417_PP12.rds")
# cellchat.M1436_PPC12 <- readRDS("results_cellchat/cellchat_M1436_PPC12.rds")
# cellchat.M1417_PP12 <- updateCellChat(cellchat.M1417_PP12)
# cellchat.M1436_PPC12 <- updateCellChat(cellchat.M1436_PPC12)
# object.list <- list(M1417_PP12 = cellchat.M1417_PP12, M1436_PPC12 = cellchat.M1436_PPC12)
# cellchat.12wk <- mergeCellChat(object.list, add.names = names(object.list))
# 
# cellchat.M1416_PP18 <- readRDS("results_cellchat/cellchat_M1416_PP18.rds")
# cellchat.M1437_PPC18 <- readRDS("results_cellchat/cellchat_M1437_PPC18.rds")
# cellchat.M1416_PP18 <- updateCellChat(cellchat.M1416_PP18)
# cellchat.M1437_PPC18 <- updateCellChat(cellchat.M1437_PPC18)
# # object.list <- list(M1416_PP18 = cellchat.M1416_PP18, M1437_PPC18 = cellchat.M1437_PPC18)
# object.list <- list(M1416_PP18 = cellchat.M1416_PP18, M1437_PPC18 = cellchat.M1437_PPC18)
# cellchat.18wk <- mergeCellChat(object.list, add.names = names(object.list))
# 
# 
# # Compare the total number of interactions and interaction strength
# gg1 <- compareInteractions(cellchat.12wk, show.legend = F, group = c(1,2))
# gg2 <- compareInteractions(cellchat.12wk, show.legend = F, group = c(1,2), measure = "weight")
# gg1 + gg2
# 
# # Compare the total number of interactions and interaction strength
# gg1 <- compareInteractions(cellchat.18wk, show.legend = F, group = c(1,2))
# gg2 <- compareInteractions(cellchat.18wk, show.legend = F, group = c(1,2), measure = "weight")
# gg1 + gg2
# 
# 
# # Change of interaction across conditions for different pairs of cell-types
# gg1 <- netVisual_heatmap(cellchat.12wk)
# gg2 <- netVisual_heatmap(cellchat.12wk, measure = "weight")
# gg1 + gg2
# 
# gg1 <- netVisual_heatmap(cellchat.18wk)
# gg2 <- netVisual_heatmap(cellchat.18wk, measure = "weight")
# gg1 + gg2
# 
# gg1 <- netVisual_heatmap(cellchat.18wk, slot.name = "netP")
# gg2 <- netVisual_heatmap(cellchat.18wk, slot.name = "netP", measure = "weight")
# gg1 + gg2
# 
# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_diffInteraction(cellchat.18wk, weight.scale = T)
# netVisual_diffInteraction(cellchat.18wk, weight.scale = T, measure = "weight")
# 
