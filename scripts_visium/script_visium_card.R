rm(list=ls())
gc()
setwd("/localscratch/ziqi/Prostate_Cancer/")

# ----------------Pre-process the input data for CARD---------------------------
library(Seurat)
library(SeuratDisk)
# Spatial (VISIUM)
visium_pp12 <- Load10X_Spatial(data.dir = "dataset/data_visium/spatial_data/Visium/225_PP12", filename = "filtered_feature_bc_matrix.h5")
visium_pp18 <- Load10X_Spatial(data.dir = "dataset/data_visium/spatial_data/Visium/1687_PP18", filename = "filtered_feature_bc_matrix.h5")
visium_ppc12 <- Load10X_Spatial(data.dir = "dataset/data_visium/spatial_data/Visium/1161_PPC", filename = "filtered_feature_bc_matrix.h5")
visium_ppc18 <- Load10X_Spatial(data.dir = "dataset/data_visium/spatial_data/Visium/1660_PPC18", filename = "filtered_feature_bc_matrix.h5")
visium_list <- list("pp12" = visium_pp12, "pp18" = visium_pp18, "ppc12" = visium_ppc12, "ppc18" = visium_ppc18)

# reference scRNA-seq
# Load the merged anndata object (obtained from preprocess.py)
adata_intact <- read_h5ad("dataset/data_scrnaseq/seurat_integration/adata_intact_seurat.h5ad")
# annotation
meta_ref <- adata_intact$obs
# use raw count
counts_ref <- Matrix::readMM("dataset/data_scrnaseq/data/qc_data/X_intact.mtx")
counts_ref <- as(t(as.matrix(counts_ref)), "dgCMatrix")
colnames(counts_ref) <- rownames(adata_intact$obs)
rownames(counts_ref) <- rownames(adata_intact$var)

# count spatial
counts_pp12 <- Assays(visium_pp12, slot = "Spatial")@layers$counts
rownames(counts_pp12) <- rownames(Assays(visium_pp12, slot = "Spatial"))
colnames(counts_pp12) <- colnames(Assays(visium_pp12, slot = "Spatial"))
counts_pp18 <- Assays(visium_pp18, slot = "Spatial")@layers$counts
rownames(counts_pp18) <- rownames(Assays(visium_pp18, slot = "Spatial"))
colnames(counts_pp18) <- colnames(Assays(visium_pp18, slot = "Spatial"))
counts_ppc12 <- Assays(visium_ppc12, slot = "Spatial")@layers$counts
rownames(counts_ppc12) <- rownames(Assays(visium_ppc12, slot = "Spatial"))
colnames(counts_ppc12) <- colnames(Assays(visium_ppc12, slot = "Spatial"))
counts_ppc18 <- Assays(visium_ppc18, slot = "Spatial")@layers$counts
rownames(counts_ppc18) <- rownames(Assays(visium_ppc18, slot = "Spatial"))
colnames(counts_ppc18) <- colnames(Assays(visium_ppc18, slot = "Spatial"))

# row image coordinate (in pixel), not transformed (GetTissueCoordinates(...))
coord_pp12 <- visium_pp12@images$slice1@coordinates[c("imagerow", "imagecol")]
coord_pp18 <- visium_pp18@images$slice1@coordinates[c("imagerow", "imagecol")]
coord_ppc12 <- visium_ppc12@images$slice1@coordinates[c("imagerow", "imagecol")]
coord_ppc18 <- visium_ppc18@images$slice1@coordinates[c("imagerow", "imagecol")]

# --------------------------------RUN CARD--------------------------------------
library(CARD)

colnames(coord_pp12) <- c("x", "y")
colnames(coord_pp18) <- c("x", "y")
colnames(coord_ppc12) <- c("x", "y")
colnames(coord_ppc18) <- c("x", "y")

counts_visium <- list("pp12" = counts_pp12, "pp18" = counts_pp18, "ppc12" = counts_ppc12, "ppc18" = counts_ppc18)
coords_visium <- list("pp12" = coord_pp12, "pp18" = coord_pp18, "ppc12" = coord_ppc12, "ppc18" = coord_ppc18)
deconv_scores <- list()
for(sample in names(counts_visium)){
  # CARD object
  CARD_sample <- createCARDObject(
    sc_count = counts_ref,
    sc_meta = meta_ref,
    spatial_count = counts_visium[[sample]],
    spatial_location = coords_visium[[sample]],
    ct.varname = "annot",
    ct.select = unique(meta_ref$predict_label),
    sample.varname = "sample",
    minCountGene = 0,
    minCountSpot = 5) 
  
  # Deconvolution using CARD
  CARD_sample <- CARD_deconvolution(CARD_object = CARD_sample)
  scores.deconv <- CARD_sample@Proportion_CARD
  # load scores of the other methods

  # set the colors. Here, I just use the colors in the manuscript, if the color is not provided, the function will use default color in the package. 
  # colors <- c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
  #     "#FDC086","#FFFF99","#386CB0")
  # # colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
  # #            "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
  # #            "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
  # # TOO MANY CELLS
  # p1 <- CARD.visualize.pie(
  #   proportion = scores.deconv,
  #   spatial_location = CARD_sample@spatial_location,
  #   colors = colors,
  #     radius = 0.52) ### You can choose radius = NULL or your own radius number
  # print(p1)
  deconv_scores[[sample]] <- CARD_sample@Proportion_CARD
  write.table(deconv_scores[[sample]], file = paste0("results_visium/results_card/deconv_score_", sample, ".csv"), sep = "\t", quote = FALSE)
}







