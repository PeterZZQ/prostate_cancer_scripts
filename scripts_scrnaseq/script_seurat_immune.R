rm(list=ls())
gc()

setwd("/localscratch/ziqi/Prostate_Cancer")
library("Seurat")
library("anndata")
library("SeuratDisk")
library("patchwork")

raw_dir <- "dataset/data_scrnaseq/data/qc_data/"
seurat_dir <- "dataset/data_scrnaseq/seurat_integration/"

# 1. Loading only data of immune cell ----------------------------------------------------
# extract annotation from anndata, generated from `script_plots_markers.py`
adata <- read_h5ad(paste0(seurat_dir, "adata_intact_seurat.h5ad"))
X_intact <- Matrix::readMM(paste0(raw_dir, "X_intact.mtx"))
X_intact <- as(t(as.matrix(X_intact)), "dgCMatrix")
meta_cells_intact <- adata$obs
meta_genes_intact <- read.csv(paste0(raw_dir, "gene_intact.csv"), sep = ",", row.names = 1)
colnames(X_intact) <- rownames(meta_cells_intact)
rownames(X_intact) <- rownames(meta_genes_intact)
# Create the Seurat object
prostate.intact <- CreateSeuratObject(counts = X_intact, meta.data = meta_cells_intact, assay = "RNA")
prostate.intact[["RNA"]] <- AddMetaData(prostate.intact[["RNA"]], meta_genes_intact)
Idents(prostate.intact) <- "annot"
# subset only immune cell from seurat object
prostate.immune <- subset(prostate.intact, idents = c("Lymphoid", "Macrophage", "Monocytes"))
# split the dataset into layers corresponding to different samples
prostate.immune[["RNA"]] <- split(prostate.immune[["RNA"]], f = prostate.immune$sample)

# 2. Seurat integration ----------------------------------------------------
# normalize and identify variable features for each dataset independently
prostate.immune <- NormalizeData(prostate.immune)
prostate.immune <- FindVariableFeatures(prostate.immune, selection.method = "vst", nfeatures = 2000)
prostate.immune <- ScaleData(prostate.immune)
# calculate the PCA of the data
prostate.immune <- RunPCA(prostate.immune)
# run Seurat integration, original CCA method
prostate.immune <- IntegrateLayers(object = prostate.immune, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
# Obtain cluster result
prostate.immune <- FindNeighbors(prostate.immune, reduction = "integrated.cca", dims = 1:30)
prostate.immune <- FindClusters(prostate.immune, resolution = 0.3, cluster.name = "cca_clusters")
# visualize the clusters
prostate.immune <- RunUMAP(prostate.immune, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
p1 <- DimPlot(prostate.immune, reduction = "umap.cca", group.by = c("sample", "cca_clusters"), combine = FALSE, label.size = 2)
p1

# 2. Save files ------------------------------------------------------
seurat_dir <- "dataset/data_scrnaseq/seurat_integration/immune_subset/"
SaveH5Seurat(prostate.immune, filename = paste0(seurat_dir, "prostate_immune.h5Seurat"), overwrite = TRUE)
meta.seurat <- prostate.immune@meta.data
pca_cca.seurat <- prostate.immune@reductions$integrated.cca@cell.embeddings
umap_cca.seurat <- prostate.immune@reductions$umap.cca@cell.embeddings
write.table(meta.seurat, file = paste0(seurat_dir, "meta_immune.csv"), sep = "\t", quote = FALSE)
write.table(pca_cca.seurat, file = paste0(seurat_dir, "cca_pca_immune.csv"), sep = "\t", quote = FALSE)
write.table(umap_cca.seurat, file = paste0(seurat_dir, "cca_umap_immune.csv"), sep = "\t", quote = FALSE)




# # 2. After annotated cells, provide high-resolution annotation on immune cells -----------------------------------------
# # Select the immune cells and re-run the integration methods
# immune_idx <- (prostate.combined@meta.data$seurat_clusters == 4)|(prostate.combined@meta.data$seurat_clusters == 5)|(prostate.combined@meta.data$seurat_clusters == 9)
# immune_barcodes <- rownames(prostate.combined@meta.data[immune_idx,])
# immune.h5ad <- prostate.h5ad[immune_barcodes,]
# 
# immune.obj <- CreateSeuratObject(counts = t(immune.h5ad$X), meta.data = immune.h5ad$obs)
# 
# # split the dataset into a list of two seurat objects (stim and CTRL)
# immune.obj.list <- SplitObject(immune.obj, split.by = "sample")
# 
# # normalize and identify variable features for each dataset independently
# immune.obj.list <- lapply(X = immune.obj.list, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
# 
# # select features that are repeatedly variable across datasets for integration
# features <- SelectIntegrationFeatures(object.list = immune.obj.list)
# # all immune cells
# immune.anchors <- FindIntegrationAnchors(object.list = immune.obj.list, anchor.features = features)
# immune.combined <- IntegrateData(anchorset = immune.anchors)
# DefaultAssay(immune.combined) <- "integrated"
# 
# immune.combined <- ScaleData(immune.combined, verbose = FALSE)
# immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30, n.neighbors = 50L, min.dist = 0.3)
# immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30, k.param = 50)
# immune.combined <- FindClusters(immune.combined, resolution = 0.5)
# 
# SaveH5Seurat(immune.combined, filename = "Cell_Ranger_output/seurat_integrate_immune.h5Seurat", overwrite = TRUE)
# 
# # --------------------------------------------------------
# # load data
# immune.combined <- LoadH5Seurat("Cell_Ranger_output/seurat_integrate_immune.h5Seurat")
# 
# # Visualization
# p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "sample")
# p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
# p1 + p2
# 
# DefaultAssay(immune.combined) <- "RNA"
# FeaturePlot(immune.combined, features = total.immune, min.cutoff = "q9")
# 
# # Myeloid
# FeaturePlot(immune.combined, features = total.myeloid, min.cutoff = "q9")
# FeaturePlot(immune.combined, features = monocytes, min.cutoff = "q9")
# FeaturePlot(immune.combined, features = total.macrophages, min.cutoff = "q9")
# FeaturePlot(immune.combined, features = m2.macrophages, min.cutoff = "q9")
# FeaturePlot(immune.combined, features = m1.macrophages, min.cutoff = "q9")
# FeaturePlot(immune.combined, features = neutrophils, min.cutoff = "q9")
# FeaturePlot(immune.combined, features = cdcs, min.cutoff = "q9")
# 
# # Lymphoid
# FeaturePlot(immune.combined, features = total.cd4.T, min.cutoff = "q9")
# FeaturePlot(immune.combined, features = total.cd8.T, min.cutoff = "q9")
# FeaturePlot(immune.combined, features = tregs, min.cutoff = "q9")
# FeaturePlot(immune.combined, features = pd1.cd4.T, min.cutoff = "q9")
# FeaturePlot(immune.combined, features = pd1.cd8.T, min.cutoff = "q9")
# FeaturePlot(immune.combined, features = exhausted.cd4.T, min.cutoff = "q9")
# FeaturePlot(immune.combined, features = exhausted.cd8.T, min.cutoff = "q9")
# FeaturePlot(immune.combined, features = cytotoxic.cd4.T, min.cutoff = "q9")
# FeaturePlot(immune.combined, features = cytotoxic.cd8.T, min.cutoff = "q9")
# FeaturePlot(immune.combined, features = nk, min.cutoff = "q9")
# FeaturePlot(immune.combined, features = b, min.cutoff = "q9")
# 
# 
# plots <- VlnPlot(immune.combined, features = cdcs, group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
# wrap_plots(plots = plots, ncol = 1)
# 
# plots <- VlnPlot(immune.combined, features = m1.macrophages, group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
# wrap_plots(plots = plots, ncol = 1)
# 
# plots <- VlnPlot(immune.combined, features = m2.macrophages, group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
# wrap_plots(plots = plots, ncol = 1)
# 
# # differential expression analysis
# immune.combined <- RenameIdents(immune.combined, `0` = "Monocytes", `1` = "Macrophage", `2` = "Monocytes", `3` = "Lymphoid", `4` = "Macrophage",
#                                 `5` = "Macrophage", `6` = "Monocytes", `7` = "Dendritic 1", `8` = "Macrophage", `9` = "Dendritic 2", `10` = "Dendritic 3")
# DimPlot(immune.combined, label = TRUE)
# 
# # dendritic.1.markers <- FindMarkers(immune.combined, ident.1 = "Dendritic 1", ident.2 = "Macrophage")
# # dendritic.2.markers <- FindMarkers(immune.combined, ident.1 = "Dendritic 2", ident.2 = "Macrophage")
# # dendritic.3.markers <- FindMarkers(immune.combined, ident.1 = "Dendritic 3", ident.2 = "Macrophage")
# dendritic.1.markers <- FindMarkers(immune.combined, ident.1 = "Dendritic 1", ident.2 = NULL, features = immune.combined@assays$integrated@var.features)
# dendritic.2.markers <- FindMarkers(immune.combined, ident.1 = "Dendritic 2", ident.2 = NULL, features = immune.combined@assays$integrated@var.features)
# dendritic.3.markers <- FindMarkers(immune.combined, ident.1 = "Dendritic 3", ident.2 = NULL, features = immune.combined@assays$integrated@var.features)
# macrophage.markers <- FindMarkers(immune.combined, ident.1 = "Macrophage", ident.2 = NULL, features = immune.combined@assays$integrated@var.features)
# dendritic.1.markers <- dendritic.1.markers[dendritic.1.markers$p_val_adj < 0.05,]
# dendritic.2.markers <- dendritic.2.markers[dendritic.2.markers$p_val_adj < 0.05,]
# dendritic.3.markers <- dendritic.3.markers[dendritic.3.markers$p_val_adj < 0.05,]
# macrophage.markers <- macrophage.markers[macrophage.markers$p_val_adj < 0.05,]
# dendritic.1.markers.enriched <- dendritic.1.markers[dendritic.1.markers$avg_log2FC >= 0,]
# dendritic.1.markers.depleted <- dendritic.1.markers[dendritic.1.markers$avg_log2FC < 0,]
# dendritic.2.markers.enriched <- dendritic.2.markers[dendritic.2.markers$avg_log2FC >= 0,]
# dendritic.2.markers.depleted <- dendritic.2.markers[dendritic.2.markers$avg_log2FC < 0,]
# dendritic.3.markers.enriched <- dendritic.3.markers[dendritic.3.markers$avg_log2FC >= 0,]
# dendritic.3.markers.depleted <- dendritic.3.markers[dendritic.3.markers$avg_log2FC < 0,]
# macrophage.markers.enriched <- macrophage.markers[macrophage.markers$avg_log2FC >= 0,]
# macrophage.markers.depleted <- macrophage.markers[macrophage.markers$avg_log2FC < 0,]
# 
# dendritic.1.markers.enriched <- dendritic.1.markers.enriched[order(x = dendritic.1.markers.enriched$p_val_adj, y = -dendritic.1.markers.enriched$avg_log2FC),]
# dendritic.1.markers.depleted <- dendritic.1.markers.depleted[order(x = dendritic.1.markers.depleted$p_val_adj, y = dendritic.1.markers.depleted$avg_log2FC),]
# dendritic.2.markers.enriched <- dendritic.2.markers.enriched[order(x = dendritic.2.markers.enriched$p_val_adj, y = -dendritic.2.markers.enriched$avg_log2FC),]
# dendritic.2.markers.depleted <- dendritic.2.markers.depleted[order(x = dendritic.2.markers.depleted$p_val_adj, y = dendritic.2.markers.depleted$avg_log2FC),]
# dendritic.3.markers.enriched <- dendritic.3.markers.enriched[order(x = dendritic.3.markers.enriched$p_val_adj, y = -dendritic.3.markers.enriched$avg_log2FC),]
# dendritic.3.markers.depleted <- dendritic.3.markers.depleted[order(x = dendritic.3.markers.depleted$p_val_adj, y = dendritic.3.markers.depleted$avg_log2FC),]
# macrophage.markers.enriched <- macrophage.markers.enriched[order(x = macrophage.markers.enriched$p_val_adj, y = -macrophage.markers.enriched$avg_log2FC),]
# macrophage.markers.depleted <- macrophage.markers.depleted[order(x = macrophage.markers.depleted$p_val_adj, y = macrophage.markers.depleted$avg_log2FC),]
# 
# write.table(dendritic.1.markers.enriched, file = "immune_integration/DE_dendritic_enriched_(cluster7).csv", sep = "\t", quote = FALSE)
# write.table(dendritic.1.markers.depleted, file = "immune_integration/DE_dendritic_depleted_(cluster7).csv", sep = "\t", quote = FALSE)
# write.table(dendritic.2.markers.enriched, file = "immune_integration/DE_dendritic_enriched_(cluster9).csv", sep = "\t", quote = FALSE)
# write.table(dendritic.2.markers.depleted, file = "immune_integration/DE_dendritic_depleted_(cluster9).csv", sep = "\t", quote = FALSE)
# write.table(dendritic.3.markers.enriched, file = "immune_integration/DE_dendritic_enriched_(cluster10).csv", sep = "\t", quote = FALSE)
# write.table(dendritic.3.markers.depleted, file = "immune_integration/DE_dendritic_depleted_(cluster10).csv", sep = "\t", quote = FALSE)
# write.table(macrophage.markers.enriched, file = "immune_integration/DE_macrophage_enriched.csv", sep = "\t", quote = FALSE)
# write.table(macrophage.markers.depleted, file = "immune_integration/DE_macrophage_depleted.csv", sep = "\t", quote = FALSE)
# 
# dendritic.1.2.markers <- FindMarkers(immune.combined, ident.1 = "Dendritic 1", ident.2 = "Dendritic 2")

