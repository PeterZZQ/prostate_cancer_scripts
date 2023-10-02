rm(list=ls())
gc()

# ----------------Pre-process the input data for CARD---------------------------
# library(Seurat)
# library(SeuratDisk)
# setwd("/localscratch/ziqi/Prostate_Cancer/")
# # Spatial (VISIUM)
# visium_pp12 <- Load10X_Spatial(data.dir = "spatial_data/Visium/225_PP12", filename = "filtered_feature_bc_matrix.h5")
# visium_pp18 <- Load10X_Spatial(data.dir = "spatial_data/Visium/1687_PP18", filename = "filtered_feature_bc_matrix.h5")
# visium_ppc12 <- Load10X_Spatial(data.dir = "spatial_data/Visium/1161_PPC", filename = "filtered_feature_bc_matrix.h5")
# visium_ppc18 <- Load10X_Spatial(data.dir = "spatial_data/Visium/1660_PPC18", filename = "filtered_feature_bc_matrix.h5")
# visium_list <- list("pp12" = visium_pp12, "pp18" = visium_pp18, "ppc12" = visium_ppc12, "ppc18" = visium_ppc18)
# 
# # reference scRNA-seq
# reference.obj <- LoadH5Seurat("Cell_Ranger_output/seurat_integrate.h5Seurat")
# reference.obj <- RenameIdents(reference.obj, `0` = "Luminal", `1` = "Luminal", `2` = "Hillock epithelia (Basal)",
#                                   `3` = "Club epithelia (Luminal)", `4` = "Macrophage (Myeloid, Total immune)", `5` = "Monocytes (Myeloid, Total immune)", 
#                                   `6` = "Mesenchymal", `7` = "Spink1+ (Luminal)", `8` = "Basal", `9` = "Lymphoid (Total immune)",
#                                   `10` = "SV", `11` = "Endothelial", `12` = "Luminal 2", `13` = "Luminal 3")
# reference.obj <- reference.obj[,(reference.obj@active.ident!= "Luminal 3")&(reference.obj@active.ident!= "Luminal 2")]
# 
# # count spatial
# counts_pp12 <- Assays(visium_pp12, slot = "Spatial")@counts
# counts_pp18 <- Assays(visium_pp18, slot = "Spatial")@counts
# counts_ppc12 <- Assays(visium_ppc12, slot = "Spatial")@counts
# counts_ppc18 <- Assays(visium_ppc18, slot = "Spatial")@counts
# # row image coordinate (in pixel), not transformed (GetTissueCoordinates(...))
# coord_pp12 <- visium_pp12@images$slice1@coordinates[c("imagerow", "imagecol")]
# coord_pp18 <- visium_pp18@images$slice1@coordinates[c("imagerow", "imagecol")]
# coord_ppc12 <- visium_ppc12@images$slice1@coordinates[c("imagerow", "imagecol")]
# coord_ppc18 <- visium_ppc18@images$slice1@coordinates[c("imagerow", "imagecol")]
# 
# # count & meta reference
# counts_ref <- Assays(reference.obj, slot = "RNA")@counts
# meta_ref <- reference.obj@meta.data
# meta_ref["predict_label"] <- reference.obj@active.ident
# 
# # save the processed input of CARD
# input_dir <- "results_card/inputdata/"
# dir.create(input_dir, showWarnings = F)
# save(counts_pp12, file = paste0(input_dir, "counts_pp12.RData"))
# save(counts_pp18, file = paste0(input_dir, "counts_pp18.RData"))
# save(counts_ppc12, file = paste0(input_dir, "counts_ppc12.RData"))
# save(counts_ppc18, file = paste0(input_dir, "counts_ppc18.RData"))
# 
# save(coord_pp12, file = paste0(input_dir, "coord_pp12.RData"))
# save(coord_pp18, file = paste0(input_dir, "coord_pp18.RData"))
# save(coord_ppc12, file = paste0(input_dir, "coord_ppc12.RData"))
# save(coord_ppc18, file = paste0(input_dir, "coord_ppc18.RData"))
# 
# save(counts_ref, file = paste0(input_dir, "counts_ref.RData"))
# save(meta_ref, file = paste0(input_dir, "meta_ref.RData"))

# --------------------------------RUN CARD--------------------------------------
library(CARD)
setwd("/localscratch/ziqi/Prostate_Cancer/")
#### load the example spatial transcriptomics count data
input_dir <- "results_card/inputdata/"
load(paste0(input_dir, "counts_pp12.RData"))
load(paste0(input_dir, "counts_pp18.RData"))
load(paste0(input_dir, "counts_ppc12.RData"))
load(paste0(input_dir, "counts_ppc18.RData"))
load(paste0(input_dir, "coord_pp12.RData"))
load(paste0(input_dir, "coord_pp18.RData"))
load(paste0(input_dir, "coord_ppc12.RData"))
load(paste0(input_dir, "coord_ppc18.RData"))

colnames(coord_pp12) <- c("x", "y")
colnames(coord_pp18) <- c("x", "y")
colnames(coord_ppc12) <- c("x", "y")
colnames(coord_ppc18) <- c("x", "y")

load(paste0(input_dir, "counts_ref.RData"))
load(paste0(input_dir, "meta_ref.RData"))

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
    ct.varname = "predict_label",
    ct.select = unique(meta_ref$predict_label),
    sample.varname = "sample",
    minCountGene = 0,
    minCountSpot = 5) 
  
  # Deconvolution using CARD
  CARD_sample <- CARD_deconvolution(CARD_object = CARD_sample)
  deconv_scores[[sample]] <- CARD_sample@Proportion_CARD
  write.table(deconv_scores[[sample]], file = paste0("results_card/deconv_score_", sample, ".csv"), sep = "\t", quote = FALSE)
}







