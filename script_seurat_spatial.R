rm(list=ls())
gc()

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(anndata)


setwd("/localscratch/ziqi/Prostate_Cancer/")


# ------------------------------------------------------------------------------
# 
# 1. cluster cells by gene expression data (https://satijalab.org/seurat/articles/spatial_vignette.html)
#
# ------------------------------------------------------------------------------
visium_pp12 <- Load10X_Spatial(data.dir = "spatial_data/Visium/225_PP12", filename = "filtered_feature_bc_matrix.h5")
visium_pp18 <- Load10X_Spatial(data.dir = "spatial_data/Visium/1687_PP18", filename = "filtered_feature_bc_matrix.h5")
visium_ppc12 <- Load10X_Spatial(data.dir = "spatial_data/Visium/1161_PPC", filename = "filtered_feature_bc_matrix.h5")
visium_ppc18 <- Load10X_Spatial(data.dir = "spatial_data/Visium/1660_PPC18", filename = "filtered_feature_bc_matrix.h5")

visium_list <- list("pp12" = visium_pp12, "pp18" = visium_pp18, "ppc12" = visium_ppc12, "ppc18" = visium_ppc18)

marker.known <- read.csv("Dataset/Cell annotation gene markers.csv")
marker.luminal <- marker.known[marker.known["List"] == "Luminal markers", "Name"]
marker.basal <- marker.known[marker.known["List"] == "Basal markers", "Name"]
marker.mesenchymal <- marker.known[marker.known["List"] == "mesenchymal", "Name"]
marker.sv <- marker.known[marker.known["List"] == "SV", "Name"]
marker.total_immune <- marker.known[marker.known["List"] == "total immune cells", "Name"]
marker.myeloid <- marker.known[marker.known["List"] == "Myeloid cells", "Name"]
marker.lymphoid <- marker.known[marker.known["List"] == "Lymphoid cells", "Name"]
marker.hillock_epithelia <- marker.known[marker.known["List"] == "Hillock epithelia", "Name"]
marker.endothelial <- marker.known[marker.known["List"] == "endothelial", "Name"]
marker.club_epithelia <- marker.known[marker.known["List"] == "club epithelia", "Name"]
marker.spink1 <- marker.known[marker.known["List"] == "Spink1+", "Name"]

for(key in names(visium_list)){
  visium_data <- visium_list[[key]]
  # normalize the input data (sctransform)
  visium_data <- SCTransform(visium_data, assay = "Spatial", verbose = FALSE)
  # cluster by gene expression data
  visium_data <- RunPCA(visium_data, assay = "SCT", verbose = FALSE)
  visium_data <- FindNeighbors(visium_data, reduction = "pca", dims = 1:30)
  visium_data <- FindClusters(visium_data, verbose = FALSE, resolution = 0.8)
  visium_data <- RunUMAP(visium_data, reduction = "pca", dims = 1:30)
  
  result_markers_dir = paste0("results_seurat_visium/markers_", key, "/")
  dir.create(result_markers_dir, showWarnings = FALSE)
  p1 <- DimPlot(visium_data, reduction = "umap", label = TRUE)
  ggsave(filename = paste0(result_markers_dir, "clusters.pdf"), device = "pdf", plot = p1, width = 12, height = 10)
  p2 <- SpatialDimPlot(visium_data, label = TRUE, label.size = 3)
  ggsave(filename = paste0(result_markers_dir, "clusters_spatial.pdf"), device = "pdf", plot = p2, width = 12, height = 10)
  
  p1 <- FeaturePlot(visium_data, features = marker.luminal[1:3], min.cutoff = "q9")
  ggsave(filename = paste0(result_markers_dir, "marker_luminal.pdf"), device = "pdf", plot = p1, width = 12, height = 10)
  ggsave(filename = paste0(result_markers_dir, "marker_luminal.png"), device = "png", plot = p1, width = 12, height = 10, dpi = 100)
  p1 <- FeaturePlot(visium_data, features = marker.basal, min.cutoff = "q9")
  ggsave(filename = paste0(result_markers_dir, "marker_basal.pdf"), device = "pdf", plot = p1, width = 12, height = 10)
  ggsave(filename = paste0(result_markers_dir, "marker_basal.png"), device = "png", plot = p1, width = 12, height = 10, dpi = 100)
  FeaturePlot(visium_data, features = marker.mesenchymal, min.cutoff = "q9")
  ggsave(filename = paste0(result_markers_dir, "marker_mesenchymal.pdf"), device = "pdf", plot = p1, width = 12, height = 8)
  ggsave(filename = paste0(result_markers_dir, "marker_mesenchymal.png"), device = "png", plot = p1, width = 12, height = 8, dpi = 100)
  p1 <- FeaturePlot(visium_data, features = marker.sv, min.cutoff = "q9")
  ggsave(filename = paste0(result_markers_dir, "marker_sv.pdf"), device = "pdf", plot = p1, width = 12, height = 10)
  ggsave(filename = paste0(result_markers_dir, "marker_sv.png"), device = "png", plot = p1, width = 12, height = 10, dpi = 100)
  p1 <- FeaturePlot(visium_data, features = marker.total_immune, min.cutoff = "q9")
  ggsave(filename = paste0(result_markers_dir, "marker_total_immune.pdf"), device = "pdf", plot = p1, width = 6, height = 5)
  ggsave(filename = paste0(result_markers_dir, "marker_total_immune.png"), device = "png", plot = p1, width = 6, height = 5, dpi = 100)
  p1 <- FeaturePlot(visium_data, features = marker.myeloid, min.cutoff = "q9")
  ggsave(filename = paste0(result_markers_dir, "marker_myeloid.pdf"), device = "pdf", plot = p1, width = 12, height = 10)
  ggsave(filename = paste0(result_markers_dir, "marker_myeloid.png"), device = "png", plot = p1, width = 12, height = 10, dpi = 100)
  p1 <- FeaturePlot(visium_data, features = marker.lymphoid, min.cutoff = "q9")
  ggsave(filename = paste0(result_markers_dir, "marker_lymphoid.pdf"), device = "pdf", plot = p1, width = 12, height = 10)
  ggsave(filename = paste0(result_markers_dir, "marker_lymphoid.png"), device = "png", plot = p1, width = 12, height = 10, dpi = 100)
  p1 <- FeaturePlot(visium_data, features = marker.hillock_epithelia, min.cutoff = "q9")
  ggsave(filename = paste0(result_markers_dir, "marker_epithelia.pdf"), device = "pdf", plot = p1, width = 12, height = 10)
  ggsave(filename = paste0(result_markers_dir, "marker_epithelia.png"), device = "png", plot = p1, width = 12, height = 10, dpi = 100)
  p1 <- FeaturePlot(visium_data, features = marker.endothelial, min.cutoff = "q9")
  ggsave(filename = paste0(result_markers_dir, "marker_endothelial.pdf"), device = "pdf", plot = p1, width = 12, height = 5)
  ggsave(filename = paste0(result_markers_dir, "marker_endothelial.png"), device = "png", plot = p1, width = 12, height = 5, dpi = 100)
  p1 <- FeaturePlot(visium_data, features = marker.club_epithelia, min.cutoff = "q9")
  ggsave(filename = paste0(result_markers_dir, "marker_club_epithelia.pdf"), device = "pdf", plot = p1, width = 12, height = 5)
  ggsave(filename = paste0(result_markers_dir, "marker_club_epithelia.png"), device = "png", plot = p1, width = 12, height = 5, dpi = 100)
  p1 <- FeaturePlot(visium_data, features = marker.spink1, min.cutoff = "q9")
  ggsave(filename = paste0(result_markers_dir, "marker_spink1.pdf"), device = "pdf", plot = p1, width = 6, height = 5)
  ggsave(filename = paste0(result_markers_dir, "marker_spink1.png"), device = "png", plot = p1, width = 6, height = 5, dpi = 100)
  
}




# ------------------------------------------------------------------------------
# 
# 2. Transfer label from scRNA-seq data (https://satijalab.org/seurat/articles/spatial_vignette.html)
#
# ------------------------------------------------------------------------------  
visium_pp12 <- Load10X_Spatial(data.dir = "spatial_data/Visium/225_PP12", filename = "filtered_feature_bc_matrix.h5")
visium_pp18 <- Load10X_Spatial(data.dir = "spatial_data/Visium/1687_PP18", filename = "filtered_feature_bc_matrix.h5")
visium_ppc12 <- Load10X_Spatial(data.dir = "spatial_data/Visium/1161_PPC", filename = "filtered_feature_bc_matrix.h5")
visium_ppc18 <- Load10X_Spatial(data.dir = "spatial_data/Visium/1660_PPC18", filename = "filtered_feature_bc_matrix.h5")

# reference.adata <- read_h5ad("Cell_Ranger_output/adata_merge.h5ad")
# reference.adata.pp12 <- reference.adata[(reference.adata$obs["age"] == "12wk") & (reference.adata$obs["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)"),]
# reference.pp12 <- CreateSeuratObject(counts = t(reference.adata.pp12$X), meta.data = reference.adata.pp12$obs)
# assign annotation
reference.obj <- LoadH5Seurat("Cell_Ranger_output/seurat_integrate.h5Seurat")
reference.obj <- RenameIdents(reference.obj, `0` = "Luminal 1", `1` = "Luminal 1", `2` = "Hillock epithelia (Basal)",
                              `3` = "Club epithelia (Luminal)", `4` = "Macrophage (Myeloid, Total immune)", `5` = "Monocytes (Myeloid, Total immune)", 
                              `6` = "Mesenchymal", `7` = "Spink1+ (Luminal)", `8` = "Basal", `9` = "Lymphoid (Total immune)",
                              `10` = "SV, Spink1+", `11` = "Endothelial", `12` = "Luminal 2", `13` = "Luminal 3")

# 1. BATCH PP12
reference.pp12 <- reference.obj[,(reference.obj@meta.data["age"] == "12wk") & (reference.obj@meta.data["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")]
counts <- reference.pp12@assays$RNA@counts
meta.data <- reference.pp12@meta.data
meta.data$labels <- reference.pp12@active.ident
reference.pp12 <- CreateSeuratObject(counts = counts, meta.data = meta.data)
# use sctransform for reference
reference.pp12 <- SCTransform(reference.pp12, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
# use sctransform for spatial
visium_pp12 <- SCTransform(visium_pp12, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
anchors <- FindTransferAnchors(reference = reference.pp12, query = visium_pp12, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = reference.pp12$labels, prediction.assay = TRUE,
                                  weight.reduction = visium_pp12[["pca"]], dims = 1:30)
visium_pp12[["predictions"]] <- predictions.assay
DefaultAssay(visium_pp12) <- "predictions"  
SpatialFeaturePlot(visium_pp12, features = c("Luminal 1", "Hillock epithelia (Basal)", "Club epithelia (Luminal)", "Macrophage (Myeloid, Total immune)",
                                             "Monocytes (Myeloid, Total immune)", "Mesenchymal", "Spink1+ (Luminal)", "Basal", "Lymphoid (Total immune)", 
                                             "SV, Spink1+", "Endothelial", "Luminal 2", "Luminal 3"), pt.size.factor = 1.6, ncol = 4, crop = TRUE)
transfer.labels <- rownames(predictions.assay@data)[max.col(t(predictions.assay@data), ties.method = "first")]
visium_pp12@meta.data["predict.label"] <- transfer.labels

# 2. BATCH PP18
reference.pp18 <- reference.obj[,(reference.obj@meta.data["age"] == "18wk") & (reference.obj@meta.data["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-)")]
counts <- reference.pp18@assays$RNA@counts
meta.data <- reference.pp18@meta.data
meta.data$labels <- reference.pp18@active.ident
reference.pp18 <- CreateSeuratObject(counts = counts, meta.data = meta.data)
# use sctransform for reference
reference.pp18 <- SCTransform(reference.pp18, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
# use sctransform for spatial
visium_pp18 <- SCTransform(visium_pp18, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
anchors <- FindTransferAnchors(reference = reference.pp18, query = visium_pp18, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = reference.pp18$labels, prediction.assay = TRUE,
                                  weight.reduction = visium_pp18[["pca"]], dims = 1:30)
visium_pp18[["predictions"]] <- predictions.assay
DefaultAssay(visium_pp18) <- "predictions"  
SpatialFeaturePlot(visium_pp18, features = c("Luminal 1", "Hillock epithelia (Basal)", "Club epithelia (Luminal)", "Macrophage (Myeloid, Total immune)",
                                             "Monocytes (Myeloid, Total immune)", "Mesenchymal", "Spink1+ (Luminal)", "Basal", "Lymphoid (Total immune)", 
                                             "SV, Spink1+", "Endothelial", "Luminal 2", "Luminal 3"), pt.size.factor = 1.6, ncol = 4, crop = TRUE)
transfer.labels <- rownames(predictions.assay@data)[max.col(t(predictions.assay@data), ties.method = "first")]
visium_pp18@meta.data["predict.label"] <- transfer.labels


# 3. BATCH PPC12
reference.ppc12 <- reference.obj[,(reference.obj@meta.data["age"] == "18wk") & (reference.obj@meta.data["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)")]
counts <- reference.ppc12@assays$RNA@counts
meta.data <- reference.ppc12@meta.data
meta.data$labels <- reference.ppc12@active.ident
reference.ppc12 <- CreateSeuratObject(counts = counts, meta.data = meta.data)
# use sctransform for reference
reference.ppc12 <- SCTransform(reference.ppc12, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
# use sctransform for spatial
visium_ppc12 <- SCTransform(visium_ppc12, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
anchors <- FindTransferAnchors(reference = reference.ppc12, query = visium_ppc12, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = reference.ppc12$labels, prediction.assay = TRUE,
                                  weight.reduction = visium_ppc12[["pca"]], dims = 1:30)
visium_ppc12[["predictions"]] <- predictions.assay
DefaultAssay(visium_ppc12) <- "predictions"  
SpatialFeaturePlot(visium_ppc12, features = c("Luminal 1", "Hillock epithelia (Basal)", "Club epithelia (Luminal)", "Macrophage (Myeloid, Total immune)",
                                              "Monocytes (Myeloid, Total immune)", "Mesenchymal", "Spink1+ (Luminal)", "Basal", "Lymphoid (Total immune)", 
                                              "SV, Spink1+", "Endothelial", "Luminal 2", "Luminal 3"), pt.size.factor = 1.6, ncol = 4, crop = TRUE)
transfer.labels <- rownames(predictions.assay@data)[max.col(t(predictions.assay@data), ties.method = "first")]
visium_ppc12@meta.data["predict.label"] <- transfer.labels


# 4. BATCH PPC18
reference.ppc18 <- reference.obj[,(reference.obj@meta.data["age"] == "18wk") & (reference.obj@meta.data["genotype"] == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)")]
counts <- reference.ppc18@assays$RNA@counts
meta.data <- reference.ppc18@meta.data
meta.data$labels <- reference.ppc18@active.ident
reference.ppc18 <- CreateSeuratObject(counts = counts, meta.data = meta.data)
# use sctransform for reference
reference.ppc18 <- SCTransform(reference.ppc18, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
# use sctransform for spatial
visium_ppc18 <- SCTransform(visium_ppc18, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
anchors <- FindTransferAnchors(reference = reference.ppc18, query = visium_ppc18, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = reference.ppc18$labels, prediction.assay = TRUE,
                                  weight.reduction = visium_ppc18[["pca"]], dims = 1:30)
visium_ppc18[["predictions"]] <- predictions.assay
DefaultAssay(visium_ppc18) <- "predictions"  
SpatialFeaturePlot(visium_ppc18, features = c("Luminal 1", "Hillock epithelia (Basal)", "Club epithelia (Luminal)", "Macrophage (Myeloid, Total immune)",
                                             "Monocytes (Myeloid, Total immune)", "Mesenchymal", "Spink1+ (Luminal)", "Basal", "Lymphoid (Total immune)", 
                                             "SV, Spink1+", "Endothelial", "Luminal 2", "Luminal 3"), pt.size.factor = 1.6, ncol = 4, crop = TRUE)
transfer.labels <- rownames(predictions.assay@data)[max.col(t(predictions.assay@data), ties.method = "first")]
visium_ppc18@meta.data["predict.label"] <- transfer.labels


# Save the transferred labels
write.table(file = "spatial_data/Visium/transfer_labels_225_pp12.csv", x = visium_pp12@meta.data)
write.table(file = "spatial_data/Visium/transfer_labels_1687_pp18.csv", x = visium_pp18@meta.data)
write.table(file = "spatial_data/Visium/transfer_labels_1161_ppc.csv", x = visium_ppc12@meta.data)
write.table(file = "spatial_data/Visium/transfer_labels_1660_ppc18.csv", x = visium_ppc18@meta.data)

