rm(list=ls())
gc()

setwd("/localscratch/ziqi/Prostate_Cancer")
library("Seurat")
library("anndata")
library("SeuratDisk")
library("patchwork")

# 1. Run Seurat integration ----------------------------------------------------
prostate.h5ad <- read_h5ad("Cell_Ranger_output/adata_merge.h5ad")
prostate.obj <- CreateSeuratObject(counts = t(prostate.h5ad$X), meta.data = prostate.h5ad$obs)

# split the dataset into a list of two seurat objects (stim and CTRL)
prostate.obj.list <- SplitObject(prostate.obj, split.by = "sample")

# normalize and identify variable features for each dataset independently
prostate.obj.list <- lapply(X = prostate.obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = prostate.obj.list)

# run seurat integration
prostate.anchors <- FindIntegrationAnchors(object.list = prostate.obj.list, anchor.features = features)
prostate.combined <- IntegrateData(anchorset = prostate.anchors)
DefaultAssay(prostate.combined) <- "integrated"
prostate.combined <- ScaleData(prostate.combined, verbose = FALSE)
prostate.combined <- RunPCA(prostate.combined, npcs = 30, verbose = FALSE)
prostate.combined <- RunUMAP(prostate.combined, reduction = "pca", dims = 1:30, n.neighbors = 50L, min.dist = 0.3)
prostate.combined <- FindNeighbors(prostate.combined, reduction = "pca", dims = 1:30, k.param = 50)
prostate.combined <- FindClusters(prostate.combined, resolution = 0.5)

# save files
SaveH5Seurat(prostate.combined, filename = "Cell_Ranger_output/seurat_integrate.h5Seurat", overwrite = TRUE)


# 2. Load Seurat object and visualize ------------------------------------------
# Load Seurat object
prostate.combined <- LoadH5Seurat("Cell_Ranger_output/seurat_integrate.h5Seurat")
# x_pca <- prostate.combined@reductions$pca@cell.embeddings
# write.table(x_pca, "Cell_Ranger_output/seurat_pca.txt", sep = "\t", quote = F)
# write.table(prostate.combined@meta.data, "Cell_Ranger_output/seurat_meta.txt", sep = "\t", quote = F)

# Visualization
p1 <- DimPlot(prostate.combined, reduction = "umap", group.by = "sample")
p2 <- DimPlot(prostate.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
ggsave(filename = "integration/seurat_cluster.pdf", device = "pdf", plot = p1 + p2, width = 12, height = 5)
ggsave(filename = "integration/seurat_cluster.png", device = "png", plot = p1 + p2, width = 12, height = 5, dpi = 150)

prostate.m1416_pp18 <- prostate.combined[,prostate.combined$sample == "M1416-PP18"]
prostate.m1417_pp12 <- prostate.combined[,prostate.combined$sample == "M1417-PP12"]
prostate.m1436_ppc12 <- prostate.combined[,prostate.combined$sample == "M1436-PPC12"]
prostate.m1437_ppc18 <- prostate.combined[,prostate.combined$sample == "M1437-PPC18"]

p1 <- DimPlot(prostate.m1416_pp18, reduction = "umap", group.by = "seurat_clusters")
p2 <- DimPlot(prostate.m1417_pp12, reduction = "umap", group.by = "seurat_clusters")
p3 <- DimPlot(prostate.m1436_ppc12, reduction = "umap", group.by = "seurat_clusters")
p4 <- DimPlot(prostate.m1437_ppc18, reduction = "umap", group.by = "seurat_clusters")
p1 + p2 + p3 + p4
ggsave(filename = "integration/seurat_cluster_separate.pdf", device = "pdf", plot = p1 + p2 + p3 + p4, width = 12, height = 10)
ggsave(filename = "integration/seurat_cluster_separate.png", device = "png", plot = p1 + p2 + p3 + p4, width = 12, height = 10, dpi = 150)


# 3. Plot marker genes ---------------------------------------------------------

# markers, columns: List, Name, ID, FeatureType
marker.known <- read.csv("Dataset/Cell annotation gene markers.csv")
DefaultAssay(prostate.combined) <- "RNA"
marker.luminal <- marker.known[marker.known["List"] == "Luminal markers", "Name"]
marker.basal <- marker.known[marker.known["List"] == "Basal markers", "Name"]
marker.ne <- marker.known[marker.known["List"] == "NE markers", "Name"]
marker.early_squamous <- marker.known[marker.known["List"] == "early squamous diff", "Name"]
marker.late_squamous <- marker.known[marker.known["List"] == "late squamous diff", "Name"]
marker.fibroblasts <- marker.known[marker.known["List"] == "fibroblasts", "Name"]
marker.mesenchymal <- marker.known[marker.known["List"] == "mesenchymal", "Name"]
marker.sv <- marker.known[marker.known["List"] == "SV", "Name"]
marker.total_immune <- marker.known[marker.known["List"] == "total immune cells", "Name"]
marker.myeloid <- marker.known[marker.known["List"] == "Myeloid cells", "Name"]
marker.lymphoid <- marker.known[marker.known["List"] == "Lymphoid cells", "Name"]
marker.l1 <- marker.known[marker.known["List"] == "L1", "Name"]
marker.l3 <- marker.known[marker.known["List"] == "L3", "Name"]
marker.l2 <- marker.known[marker.known["List"] == "L2 stem cell like", "Name"]
marker.myofibroblasts <- marker.known[marker.known["List"] == "myofibroblasts and smooth muscle populations", "Name"]
marker.hillock_epithelia <- marker.known[marker.known["List"] == "Hillock epithelia", "Name"]
marker.endothelial <- marker.known[marker.known["List"] == "endothelial", "Name"]
marker.club_epithelia <- marker.known[marker.known["List"] == "club epithelia", "Name"]
marker.spink1 <- marker.known[marker.known["List"] == "Spink1+", "Name"]


p1 <- FeaturePlot(prostate.combined, features = marker.luminal, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_luminal.pdf", device = "pdf", plot = p1, width = 12, height = 10)
ggsave(filename = "integration/markers/marker_luminal.png", device = "png", plot = p1, width = 12, height = 10, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.basal, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_basal.pdf", device = "pdf", plot = p1, width = 12, height = 10)
ggsave(filename = "integration/markers/marker_basal.png", device = "png", plot = p1, width = 12, height = 10, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.ne, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_ne.pdf", device = "pdf", plot = p1, width = 12, height = 8)
ggsave(filename = "integration/markers/marker_ne.png", device = "png", plot = p1, width = 12, height = 8, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.early_squamous, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_early_squamous.pdf", device = "pdf", plot = p1, width = 12, height = 5)
ggsave(filename = "integration/markers/marker_early_squamous.png", device = "png", plot = p1, width = 12, height = 5, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.late_squamous, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_late_squamous.pdf", device = "pdf", plot = p1, width = 12, height = 5)
ggsave(filename = "integration/markers/marker_late_squamous.png", device = "png", plot = p1, width = 12, height = 5, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.fibroblasts, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_fibroblasts.pdf", device = "pdf", plot = p1, width = 12, height = 10)
ggsave(filename = "integration/markers/marker_fibroblasts.png", device = "png", plot = p1, width = 12, height = 10, dpi = 100)
FeaturePlot(prostate.combined, features = marker.mesenchymal, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_mesenchymal.pdf", device = "pdf", plot = p1, width = 12, height = 8)
ggsave(filename = "integration/markers/marker_mesenchymal.png", device = "png", plot = p1, width = 12, height = 8, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.sv, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_sv.pdf", device = "pdf", plot = p1, width = 12, height = 10)
ggsave(filename = "integration/markers/marker_sv.png", device = "png", plot = p1, width = 12, height = 10, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.total_immune, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_total_immune.pdf", device = "pdf", plot = p1, width = 6, height = 5)
ggsave(filename = "integration/markers/marker_total_immune.png", device = "png", plot = p1, width = 6, height = 5, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.myeloid, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_myeloid.pdf", device = "pdf", plot = p1, width = 12, height = 10)
ggsave(filename = "integration/markers/marker_myeloid.png", device = "png", plot = p1, width = 12, height = 10, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.lymphoid, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_lymphoid.pdf", device = "pdf", plot = p1, width = 12, height = 10)
ggsave(filename = "integration/markers/marker_lymphoid.png", device = "png", plot = p1, width = 12, height = 10, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.l1, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_l1.pdf", device = "pdf", plot = p1, width = 12, height = 15)
ggsave(filename = "integration/markers/marker_l1.png", device = "png", plot = p1, width = 12, height = 15, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.l3, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_l3.pdf", device = "pdf", plot = p1, width = 12, height = 10)
ggsave(filename = "integration/markers/marker_l3.png", device = "png", plot = p1, width = 12, height = 10, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.l2, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_l2.pdf", device = "pdf", plot = p1, width = 12, height = 15)
ggsave(filename = "integration/markers/marker_l2.png", device = "png", plot = p1, width = 12, height = 15, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.myofibroblasts, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_myofibroblasts.pdf", device = "pdf", plot = p1, width = 12, height = 10)
ggsave(filename = "integration/markers/marker_myofibroblasts.png", device = "png", plot = p1, width = 12, height = 10, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.hillock_epithelia, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_epithelia.pdf", device = "pdf", plot = p1, width = 12, height = 10)
ggsave(filename = "integration/markers/marker_epithelia.png", device = "png", plot = p1, width = 12, height = 10, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.endothelial, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_endothelial.pdf", device = "pdf", plot = p1, width = 12, height = 5)
ggsave(filename = "integration/markers/marker_endothelial.png", device = "png", plot = p1, width = 12, height = 5, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.club_epithelia, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_club_epithelia.pdf", device = "pdf", plot = p1, width = 12, height = 5)
ggsave(filename = "integration/markers/marker_club_epithelia.png", device = "png", plot = p1, width = 12, height = 5, dpi = 100)
p1 <- FeaturePlot(prostate.combined, features = marker.spink1, min.cutoff = "q9")
ggsave(filename = "integration/markers/marker_spink1.pdf", device = "pdf", plot = p1, width = 6, height = 5)
ggsave(filename = "integration/markers/marker_spink1.png", device = "png", plot = p1, width = 6, height = 5, dpi = 100)


# Plot interest genes
library(readxl)
interest.genes <- read_excel("Dataset/Interest genes.xlsx")
immune.markers <- read_excel("Dataset/immune cell markers.xlsx")
# No Tp53, should be Trp53?
interest.genes <- c("Trp53", "Pten", "Ackr3", "Ar", "Krt8", "Krt5", "Trp63", "Cd4", "Cd8a", "Cd68", "Foxp3", "Ctla4", 
                    "Itgam", "Itgax", "Il2", "Il4", "Il6", "Cxcl2", "Cxcl12", "Spink1", "Tmprss4", "Syp", "Chga", "Sox2", "Mycn", 
                    "Dcn", "Pate4", "Nkx3-1", "Pbsn", "Lgr5", "Fn1", "Ccl21a", "Ccl21b", "Ccl21c", "Fgf10", "Cxcr5")
FeaturePlot(prostate.combined, features = interest.genes[1:12], min.cutoff = "q9")
FeaturePlot(prostate.combined, features = interest.genes[12:24], min.cutoff = "q9")
FeaturePlot(prostate.combined, features = interest.genes[24:36], min.cutoff = "q9")


# Plot immune cell marker genes 
total.immune <- c("Ptprc")
total.myeloid <- c("Ptprc", "Itgam")
monocytes <- c("Ptprc", "Itgam", "Cd14", "S100a8", "S100a9")
total.macrophages <- c("Ptprc", "Itgam", "Adgre1")
m2.macrophages <- c("Ptprc", "Itgam", "Adgre1", "Mrc1", "Arg1")
m1.macrophages <- c("Ptprc", "Itgam", "Adgre1", "Cd68", "Nos2")
neutrophils <- c("Ptprc", "Itgam", "Adgre1", "Ly6g")
cdcs <- c("Ptprc", "Itgam", "Itgax", "H2-Eb1", "H2-Aa", "H2-Ab1")
total.cd4.T <- c("Ptprc", "Cd3e", "Cd4")
total.cd8.T <- c("Ptprc", "Cd3e", "Cd8a")
tregs <- c("Ptprc", "Cd3e", "Cd4", "Foxp3")
pd1.cd4.T <- c("Ptprc", "Cd3e",  "Cd4", "Pdcd1")
pd1.cd8.T <- c("Ptprc", "Cd3e", "Cd8a", "Pdcd1")
exhausted.cd4.T <- c("Ptprc", "Cd3e", "Cd4", "Havcr2")
exhausted.cd8.T <- c("Ptprc", "Cd3e", "Cd8a", "Havcr2")
cytotoxic.cd4.T <- c("Ptprc", "Cd3e", "Cd4", "Gzmb")
cytotoxic.cd8.T <- c("Ptprc", "Cd3e", "Cd8a", "Gzmb")
nk <- c("Ptprc", "Klrd1")
b <- c("Ptprc", "Ms4a1")

FeaturePlot(prostate.combined, features = total.immune, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = total.myeloid, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = monocytes, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = total.macrophages, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = m2.macrophages, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = m1.macrophages, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = neutrophils, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = cdcs, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = total.cd4.T, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = total.cd8.T, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = tregs, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = pd1.cd4.T, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = pd1.cd8.T, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = exhausted.cd4.T, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = exhausted.cd8.T, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = cytotoxic.cd4.T, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = cytotoxic.cd8.T, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = nk, min.cutoff = "q9")
FeaturePlot(prostate.combined, features = b, min.cutoff = "q9")



# 4. differential expression analysis across cell types ------------------------
# differential expression analysis between clusters 0, 1, and 7
# cluster0.markers <- FindMarkers(prostate.combined, ident.1 = 0, ident.2 = NULL, features = prostate.combined@assays$integrated@var.features)
# cluster1.markers <- FindMarkers(prostate.combined, ident.1 = 1, ident.2 = NULL, features = prostate.combined@assays$integrated@var.features)
# cluster7.markers <- FindMarkers(prostate.combined, ident.1 = 7, ident.2 = NULL, features = prostate.combined@assays$integrated@var.features)
# cluster0.1.markers <- FindMarkers(prostate.combined, ident.1 = 0, ident.2 = 1, features = prostate.combined@assays$integrated@var.features)
# cluster0.7.markers <- FindMarkers(prostate.combined, ident.1 = 0, ident.2 = 7, features = prostate.combined@assays$integrated@var.features)
# 
# cluster0.markers <- cluster0.markers[cluster0.markers$p_val_adj < 0.05,]
# cluster1.markers <- cluster1.markers[cluster1.markers$p_val_adj < 0.05,]
# cluster7.markers <- cluster7.markers[cluster7.markers$p_val_adj < 0.05,]
# cluster0.1.markers <- cluster0.1.markers[cluster0.1.markers$p_val_adj < 0.05,]
# cluster0.7.markers <- cluster0.7.markers[cluster0.7.markers$p_val_adj < 0.05,]
# 
# cluster0.markers.enriched <- cluster0.markers[cluster0.markers$avg_log2FC >= 0,]
# cluster0.markers.depleted <- cluster0.markers[cluster0.markers$avg_log2FC < 0,]
# cluster1.markers.enriched <- cluster1.markers[cluster1.markers$avg_log2FC >= 0,]
# cluster1.markers.depleted <- cluster1.markers[cluster1.markers$avg_log2FC < 0,]
# cluster7.markers.enriched <- cluster7.markers[cluster7.markers$avg_log2FC >= 0,]
# cluster7.markers.depleted <- cluster7.markers[cluster7.markers$avg_log2FC < 0,]
# cluster0.1.markers.enriched <- cluster0.1.markers[cluster0.1.markers$avg_log2FC >= 0,]
# cluster0.1.markers.depleted <- cluster0.7.markers[cluster0.1.markers$avg_log2FC < 0,]
# cluster0.7.markers.enriched <- cluster0.7.markers[cluster0.7.markers$avg_log2FC >= 0,]
# cluster0.7.markers.depleted <- cluster0.7.markers[cluster0.7.markers$avg_log2FC < 0,]
# 
# cluster0.markers.enriched <- cluster0.markers.enriched[order(x = cluster0.markers.enriched$p_val_adj, y = -cluster0.markers.enriched$avg_log2FC),]
# cluster0.markers.depleted <- cluster0.markers.depleted[order(x = cluster0.markers.depleted$p_val_adj, y = cluster0.markers.depleted$avg_log2FC),]
# cluster1.markers.enriched <- cluster1.markers.enriched[order(x = cluster1.markers.enriched$p_val_adj, y = -cluster1.markers.enriched$avg_log2FC),]
# cluster1.markers.depleted <- cluster1.markers.depleted[order(x = cluster1.markers.depleted$p_val_adj, y = cluster1.markers.depleted$avg_log2FC),]
# cluster7.markers.enriched <- cluster7.markers.enriched[order(x = cluster7.markers.enriched$p_val_adj, y = -cluster7.markers.enriched$avg_log2FC),]
# cluster7.markers.depleted <- cluster7.markers.depleted[order(x = cluster7.markers.depleted$p_val_adj, y = cluster7.markers.depleted$avg_log2FC),]
# cluster0.1.markers.enriched <- cluster0.1.markers.enriched[order(x = cluster0.1.markers.enriched$p_val_adj, y = -cluster0.1.markers.enriched$avg_log2FC),]
# cluster0.1.markers.depleted <- cluster0.1.markers.depleted[order(x = cluster0.1.markers.depleted$p_val_adj, y = cluster0.1.markers.depleted$avg_log2FC),]
# cluster0.7.markers.enriched <- cluster0.7.markers.enriched[order(x = cluster0.7.markers.enriched$p_val_adj, y = -cluster0.7.markers.enriched$avg_log2FC),]
# cluster0.7.markers.depleted <- cluster0.7.markers.depleted[order(x = cluster0.7.markers.depleted$p_val_adj, y = cluster0.7.markers.depleted$avg_log2FC),]
# 
# write.table(cluster0.markers.enriched, file = "integration/DE_check_0&1&7/DE_cluster0_enriched.csv", sep = ",", quote = FALSE)
# write.table(cluster0.markers.depleted, file = "integration/DE_check_0&1&7/DE_cluster0_depleted.csv", sep = ",", quote = FALSE)
# write.table(cluster1.markers.enriched, file = "integration/DE_check_0&1&7/DE_cluster1_enriched.csv", sep = ",", quote = FALSE)
# write.table(cluster1.markers.depleted, file = "integration/DE_check_0&1&7/DE_cluster1_depleted.csv", sep = ",", quote = FALSE)
# write.table(cluster7.markers.enriched, file = "integration/DE_check_0&1&7/DE_cluster7_enriched.csv", sep = ",", quote = FALSE)
# write.table(cluster7.markers.depleted, file = "integration/DE_check_0&1&7/DE_cluster7_depleted.csv", sep = ",", quote = FALSE)
# write.table(cluster0.1.markers.enriched, file = "integration/DE_check_0&1&7/DE_cluster0&1_enriched.csv", sep = ",", quote = FALSE)
# write.table(cluster0.1.markers.depleted, file = "integration/DE_check_0&1&7/DE_cluster0&1_depleted.csv", sep = ",", quote = FALSE)
# write.table(cluster0.7.markers.enriched, file = "integration/DE_check_0&1&7/DE_cluster0&7_enriched.csv", sep = ",", quote = FALSE)
# write.table(cluster0.7.markers.depleted, file = "integration/DE_check_0&1&7/DE_cluster0&7_depleted.csv", sep = ",", quote = FALSE)

# differential expression analysis of all clusters
# for(cluster.id in seq(0,13)){
#   cluster.marker <- FindMarkers(prostate.combined, ident.1 = cluster.id, ident.2 = NULL, features = prostate.combined@assays$integrated@var.features)
#   cluster.marker <- cluster.marker[cluster.marker$p_val_adj < 0.05,]
#   cluster.marker.enriched <- cluster.marker[cluster.marker$avg_log2FC >= 0,]
#   cluster.marker.depleted <- cluster.marker[cluster.marker$avg_log2FC < 0,]
#   cluster.marker.enriched <- cluster.marker.enriched[order(x = cluster.marker.enriched$p_val_adj, y = -cluster.marker.enriched$avg_log2FC),]
#   cluster.marker.depleted <- cluster.marker.depleted[order(x = cluster.marker.depleted$p_val_adj, y = cluster.marker.depleted$avg_log2FC),]
#   write.table(cluster.marker.enriched, file = paste0("integration/DE/DEGeneList/DE_cluster", cluster.id, "_enriched.csv"), sep = ",", quote = FALSE)
#   write.table(cluster.marker.depleted, file = paste0("integration/DE/DEGeneList/DE_cluster", cluster.id, "_depleted.csv"), sep = ",", quote = FALSE)
# }

# annotate cell types
prostate.combined <- RenameIdents(prostate.combined, `0` = "Luminal 1", `1` = "Luminal 1", `2` = "Hillock epithelia (Basal)",
                                `3` = "Club epithelia (Luminal)", `4` = "Macrophage (Myeloid, Total immune)", `5` = "Monocytes (Myeloid, Total immune)", 
                                `6` = "Mesenchymal", `7` = "Spink1+ (Luminal)", `8` = "Basal", `9` = "Lymphoid (Total immune)",
                                `10` = "SV, Spink1+", `11` = "Endothelial", `12` = "Luminal 2", `13` = "Luminal 3")
DimPlot(prostate.combined, label = TRUE)

print("calculate de gene list...")
for(cluster.id in c("Luminal 1", "Hillock epithelia (Basal)", "Club epithelia (Luminal)", "Macrophage (Myeloid, Total immune)", "Monocytes (Myeloid, Total immune)", 
                 "Mesenchymal", "Spink1+ (Luminal)", "Basal", "Lymphoid (Total immune)", "SV, Spink1+", "Endothelial", "Luminal 2", "Luminal 3")){
  print(cluster.id)
  # de genes selected from top-2000 highly-variable genes
  cluster.marker <- FindMarkers(prostate.combined, ident.1 = cluster.id, ident.2 = NULL, features = prostate.combined@assays$integrated@var.features)
  cluster.marker <- cluster.marker[cluster.marker$p_val_adj < 0.05,]
  cluster.marker.enriched <- cluster.marker[cluster.marker$avg_log2FC >= 0,]
  cluster.marker.depleted <- cluster.marker[cluster.marker$avg_log2FC < 0,]
  cluster.marker.enriched <- cluster.marker.enriched[order(x = cluster.marker.enriched$p_val_adj, y = -cluster.marker.enriched$avg_log2FC),]
  cluster.marker.depleted <- cluster.marker.depleted[order(x = cluster.marker.depleted$p_val_adj, y = cluster.marker.depleted$avg_log2FC),]
  write.table(cluster.marker.enriched, file = paste0("integration/DE/DEGeneList/DE_", cluster.id, "_enriched.csv"), sep = ",", quote = FALSE)
  write.table(cluster.marker.depleted, file = paste0("integration/DE/DEGeneList/DE_", cluster.id, "_depleted.csv"), sep = ",", quote = FALSE)
}

# filter background genes
min.cells <- 200
num.cells <- rowSums(prostate.combined@assays$RNA@counts > 0)
genes.filtered <- names(num.cells[which(num.cells >= min.cells)])
# ranked gene list: log fold change
print("calculated ranked gene list...")
for(cluster.id in c("Luminal 1", "Hillock epithelia (Basal)", "Club epithelia (Luminal)", "Macrophage (Myeloid, Total immune)", "Monocytes (Myeloid, Total immune)", 
                    "Mesenchymal", "Spink1+ (Luminal)", "Basal", "Lymphoid (Total immune)", "SV, Spink1+", "Endothelial", "Luminal 2", "Luminal 3")){
  print(cluster.id)
  # selected from filtered background gene
  expr.fold.change <- FoldChange(prostate.combined, ident.1 = cluster.id, ident.2 = NULL, features = genes.filtered)
  expr.fold.change["ENSEMBL"] <- prostate.h5ad$var[rownames(expr.fold.change), "gene_ids"]
  expr.fold.change <- expr.fold.change[order(x = -expr.fold.change$avg_log2FC),]
  write.table(expr.fold.change, file = paste0("integration/DE/RankedGeneList/rgl_", cluster.id, ".csv"), sep = ",", quote = FALSE)
}



# 5. Cluster-specific differential-expression genes across genotypes------------
# for (cluster.id in seq(0,13)){
#   print(cluster.id)
#   prostate.cluster <- prostate.combined[,prostate.combined$seurat_clusters == cluster.id]
#   prostate.12wk <- prostate.cluster[,prostate.cluster$age == "12wk"]
#   prostate.18wk <- prostate.cluster[,prostate.cluster$age == "18wk"]
#   print(dim(prostate.12wk)[2])
#   print(dim(prostate.18wk)[2])
#   if((sum(prostate.12wk$genotype == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") > 10) & (sum(prostate.12wk$genotype == "PbCre(+/-),Pten(-/-),P53(-/-)") > 10)){
#     de.12wk <- FindMarkers(prostate.12wk, ident.1 = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ident.2 = "PbCre(+/-),Pten(-/-),P53(-/-)", group.by = "genotype", 
#                            features = prostate.combined@assays$integrated@var.features)
#     de.12wk <- de.12wk[de.12wk$p_val_adj < 0.05,]
#     de.12wk.enriched <- de.12wk[de.12wk$avg_log2FC >= 0,]
#     de.12wk.depleted <- de.12wk[de.12wk$avg_log2FC < 0,]
#     de.12wk.enriched <- de.12wk.enriched[order(x = de.12wk.enriched$p_val_adj, y = -de.12wk.enriched$avg_log2FC),]
#     de.12wk.depleted <- de.12wk.depleted[order(x = de.12wk.depleted$p_val_adj, y = de.12wk.depleted$avg_log2FC),]
#     write.table(de.12wk.enriched, file = paste0("integration/DE_genotype/DE_cluster", cluster.id, "_12wk_enriched.csv"), sep = ",", quote = FALSE)
#     write.table(de.12wk.depleted, file = paste0("integration/DE_genotype/DE_cluster", cluster.id, "_12wk_depleted.csv"), sep = ",", quote = FALSE)
#   }
#   if((sum(prostate.18wk$genotype == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") > 10) & (sum(prostate.18wk$genotype == "PbCre(+/-),Pten(-/-),P53(-/-)") > 10)){
#     de.18wk <- FindMarkers(prostate.18wk, ident.1 = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ident.2 = "PbCre(+/-),Pten(-/-),P53(-/-)", group.by = "genotype", 
#                            features = prostate.combined@assays$integrated@var.features)
#     de.18wk <- de.18wk[de.18wk$p_val_adj < 0.05,]
#     de.18wk.enriched <- de.18wk[de.18wk$avg_log2FC >= 0,]
#     de.18wk.depleted <- de.18wk[de.18wk$avg_log2FC < 0,]
#     de.18wk.enriched <- de.18wk.enriched[order(x = de.18wk.enriched$p_val_adj, y = -de.18wk.enriched$avg_log2FC),]
#     de.18wk.depleted <- de.18wk.depleted[order(x = de.18wk.depleted$p_val_adj, y = de.18wk.depleted$avg_log2FC),]
#     write.table(de.18wk.enriched, file = paste0("integration/DE_genotype/DE_cluster", cluster.id, "_18wk_enriched.csv"), sep = ",", quote = FALSE)
#     write.table(de.18wk.depleted, file = paste0("integration/DE_genotype/DE_cluster", cluster.id, "_18wk_depleted.csv"), sep = ",", quote = FALSE)
#   }
# }


# filter background gene
min.cells <- 200
num.cells <- rowSums(prostate.combined@assays$RNA@counts > 0)
genes.filtered <- names(num.cells[which(num.cells >= min.cells)])
for(cluster.id in c("Luminal 1", "Hillock epithelia (Basal)", "Club epithelia (Luminal)", "Macrophage (Myeloid, Total immune)", "Monocytes (Myeloid, Total immune)", 
                    "Mesenchymal", "Spink1+ (Luminal)", "Basal", "Lymphoid (Total immune)", "SV, Spink1+", "Endothelial", "Luminal 2", "Luminal 3")){

  print(cluster.id)
  prostate.cluster <- prostate.combined[,prostate.combined@active.ident == cluster.id]
  prostate.12wk <- prostate.cluster[,prostate.cluster$age == "12wk"]
  prostate.18wk <- prostate.cluster[,prostate.cluster$age == "18wk"]
  print(dim(prostate.12wk)[2])
  print(dim(prostate.18wk)[2])
  
  if((sum(prostate.12wk$genotype == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") > 10) & (sum(prostate.12wk$genotype == "PbCre(+/-),Pten(-/-),P53(-/-)") > 10)){
    de.12wk <- FindMarkers(prostate.12wk, ident.1 = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ident.2 = "PbCre(+/-),Pten(-/-),P53(-/-)", group.by = "genotype", 
                           features = prostate.combined@assays$integrated@var.features)
    de.12wk <- de.12wk[de.12wk$p_val_adj < 0.05,]
    de.12wk.enriched <- de.12wk[de.12wk$avg_log2FC >= 0,]
    de.12wk.depleted <- de.12wk[de.12wk$avg_log2FC < 0,]
    de.12wk.enriched <- de.12wk.enriched[order(x = de.12wk.enriched$p_val_adj, y = -de.12wk.enriched$avg_log2FC),]
    de.12wk.depleted <- de.12wk.depleted[order(x = de.12wk.depleted$p_val_adj, y = de.12wk.depleted$avg_log2FC),]
    write.table(de.12wk.enriched, file = paste0("integration/DE_genotype/DEGeneList/DE_", cluster.id, "_12wk_enriched.csv"), sep = ",", quote = FALSE)
    write.table(de.12wk.depleted, file = paste0("integration/DE_genotype/DEGeneList/DE_", cluster.id, "_12wk_depleted.csv"), sep = ",", quote = FALSE)
    
    expr.fold.change <- FoldChange(prostate.12wk, ident.1 = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ident.2 = "PbCre(+/-),Pten(-/-),P53(-/-)", group.by = "genotype", features = genes.filtered)
    expr.fold.change["ENSEMBL"] <- prostate.h5ad$var[rownames(expr.fold.change), "gene_ids"]
    expr.fold.change <- expr.fold.change[order(x = -expr.fold.change$avg_log2FC),]
    write.table(expr.fold.change, file = paste0("integration/DE_genotype/RankedGeneList/rgl_", cluster.id, "_12wk.csv"), sep = ",", quote = FALSE)
  }
  
  if((sum(prostate.18wk$genotype == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") > 10) & (sum(prostate.18wk$genotype == "PbCre(+/-),Pten(-/-),P53(-/-)") > 10)){
    de.18wk <- FindMarkers(prostate.18wk, ident.1 = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ident.2 = "PbCre(+/-),Pten(-/-),P53(-/-)", group.by = "genotype", 
                           features = prostate.combined@assays$integrated@var.features)
    de.18wk <- de.18wk[de.18wk$p_val_adj < 0.05,]
    de.18wk.enriched <- de.18wk[de.18wk$avg_log2FC >= 0,]
    de.18wk.depleted <- de.18wk[de.18wk$avg_log2FC < 0,]
    de.18wk.enriched <- de.18wk.enriched[order(x = de.18wk.enriched$p_val_adj, y = -de.18wk.enriched$avg_log2FC),]
    de.18wk.depleted <- de.18wk.depleted[order(x = de.18wk.depleted$p_val_adj, y = de.18wk.depleted$avg_log2FC),]
    write.table(de.18wk.enriched, file = paste0("integration/DE_genotype/DEGeneList/DE_", cluster.id, "_18wk_enriched.csv"), sep = ",", quote = FALSE)
    write.table(de.18wk.depleted, file = paste0("integration/DE_genotype/DEGeneList/DE_", cluster.id, "_18wk_depleted.csv"), sep = ",", quote = FALSE)
    
    expr.fold.change <- FoldChange(prostate.18wk, ident.1 = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ident.2 = "PbCre(+/-),Pten(-/-),P53(-/-)", group.by = "genotype", features = genes.filtered)
    expr.fold.change["ENSEMBL"] <- prostate.h5ad$var[rownames(expr.fold.change), "gene_ids"]
    expr.fold.change <- expr.fold.change[order(x = -expr.fold.change$avg_log2FC),]
    write.table(expr.fold.change, file = paste0("integration/DE_genotype/RankedGeneList/rgl_", cluster.id, "_18wk.csv"), sep = ",", quote = FALSE)
    
  }
}








# ------------------------------------------------------
# Select the immune cells and re-run the integration methods
immune_idx <- (prostate.combined@meta.data$seurat_clusters == 4)|(prostate.combined@meta.data$seurat_clusters == 5)|(prostate.combined@meta.data$seurat_clusters == 9)
immune_barcodes <- rownames(prostate.combined@meta.data[immune_idx,])
immune.h5ad <- prostate.h5ad[immune_barcodes,]

immune.obj <- CreateSeuratObject(counts = t(immune.h5ad$X), meta.data = immune.h5ad$obs)

# split the dataset into a list of two seurat objects (stim and CTRL)
immune.obj.list <- SplitObject(immune.obj, split.by = "sample")

# normalize and identify variable features for each dataset independently
immune.obj.list <- lapply(X = immune.obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = immune.obj.list)
# all immune cells
immune.anchors <- FindIntegrationAnchors(object.list = immune.obj.list, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"

immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30, n.neighbors = 50L, min.dist = 0.3)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30, k.param = 50)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

SaveH5Seurat(immune.combined, filename = "Cell_Ranger_output/seurat_integrate_immune.h5Seurat", overwrite = TRUE)

# --------------------------------------------------------
# load data
immune.combined <- LoadH5Seurat("Cell_Ranger_output/seurat_integrate_immune.h5Seurat")

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "sample")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DefaultAssay(immune.combined) <- "RNA"
FeaturePlot(immune.combined, features = total.immune, min.cutoff = "q9")

# Myeloid
FeaturePlot(immune.combined, features = total.myeloid, min.cutoff = "q9")
FeaturePlot(immune.combined, features = monocytes, min.cutoff = "q9")
FeaturePlot(immune.combined, features = total.macrophages, min.cutoff = "q9")
FeaturePlot(immune.combined, features = m2.macrophages, min.cutoff = "q9")
FeaturePlot(immune.combined, features = m1.macrophages, min.cutoff = "q9")
FeaturePlot(immune.combined, features = neutrophils, min.cutoff = "q9")
FeaturePlot(immune.combined, features = cdcs, min.cutoff = "q9")

# Lymphoid
FeaturePlot(immune.combined, features = total.cd4.T, min.cutoff = "q9")
FeaturePlot(immune.combined, features = total.cd8.T, min.cutoff = "q9")
FeaturePlot(immune.combined, features = tregs, min.cutoff = "q9")
FeaturePlot(immune.combined, features = pd1.cd4.T, min.cutoff = "q9")
FeaturePlot(immune.combined, features = pd1.cd8.T, min.cutoff = "q9")
FeaturePlot(immune.combined, features = exhausted.cd4.T, min.cutoff = "q9")
FeaturePlot(immune.combined, features = exhausted.cd8.T, min.cutoff = "q9")
FeaturePlot(immune.combined, features = cytotoxic.cd4.T, min.cutoff = "q9")
FeaturePlot(immune.combined, features = cytotoxic.cd8.T, min.cutoff = "q9")
FeaturePlot(immune.combined, features = nk, min.cutoff = "q9")
FeaturePlot(immune.combined, features = b, min.cutoff = "q9")


plots <- VlnPlot(immune.combined, features = cdcs, group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

plots <- VlnPlot(immune.combined, features = m1.macrophages, group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

plots <- VlnPlot(immune.combined, features = m2.macrophages, group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

# differential expression analysis
immune.combined <- RenameIdents(immune.combined, `0` = "Monocytes", `1` = "Macrophage", `2` = "Monocytes", `3` = "Lymphoid", `4` = "Macrophage",
                                `5` = "Macrophage", `6` = "Monocytes", `7` = "Dendritic 1", `8` = "Macrophage", `9` = "Dendritic 2", `10` = "Dendritic 3")
DimPlot(immune.combined, label = TRUE)

# dendritic.1.markers <- FindMarkers(immune.combined, ident.1 = "Dendritic 1", ident.2 = "Macrophage")
# dendritic.2.markers <- FindMarkers(immune.combined, ident.1 = "Dendritic 2", ident.2 = "Macrophage")
# dendritic.3.markers <- FindMarkers(immune.combined, ident.1 = "Dendritic 3", ident.2 = "Macrophage")
dendritic.1.markers <- FindMarkers(immune.combined, ident.1 = "Dendritic 1", ident.2 = NULL, features = immune.combined@assays$integrated@var.features)
dendritic.2.markers <- FindMarkers(immune.combined, ident.1 = "Dendritic 2", ident.2 = NULL, features = immune.combined@assays$integrated@var.features)
dendritic.3.markers <- FindMarkers(immune.combined, ident.1 = "Dendritic 3", ident.2 = NULL, features = immune.combined@assays$integrated@var.features)
macrophage.markers <- FindMarkers(immune.combined, ident.1 = "Macrophage", ident.2 = NULL, features = immune.combined@assays$integrated@var.features)
dendritic.1.markers <- dendritic.1.markers[dendritic.1.markers$p_val_adj < 0.05,]
dendritic.2.markers <- dendritic.2.markers[dendritic.2.markers$p_val_adj < 0.05,]
dendritic.3.markers <- dendritic.3.markers[dendritic.3.markers$p_val_adj < 0.05,]
macrophage.markers <- macrophage.markers[macrophage.markers$p_val_adj < 0.05,]
dendritic.1.markers.enriched <- dendritic.1.markers[dendritic.1.markers$avg_log2FC >= 0,]
dendritic.1.markers.depleted <- dendritic.1.markers[dendritic.1.markers$avg_log2FC < 0,]
dendritic.2.markers.enriched <- dendritic.2.markers[dendritic.2.markers$avg_log2FC >= 0,]
dendritic.2.markers.depleted <- dendritic.2.markers[dendritic.2.markers$avg_log2FC < 0,]
dendritic.3.markers.enriched <- dendritic.3.markers[dendritic.3.markers$avg_log2FC >= 0,]
dendritic.3.markers.depleted <- dendritic.3.markers[dendritic.3.markers$avg_log2FC < 0,]
macrophage.markers.enriched <- macrophage.markers[macrophage.markers$avg_log2FC >= 0,]
macrophage.markers.depleted <- macrophage.markers[macrophage.markers$avg_log2FC < 0,]

dendritic.1.markers.enriched <- dendritic.1.markers.enriched[order(x = dendritic.1.markers.enriched$p_val_adj, y = -dendritic.1.markers.enriched$avg_log2FC),]
dendritic.1.markers.depleted <- dendritic.1.markers.depleted[order(x = dendritic.1.markers.depleted$p_val_adj, y = dendritic.1.markers.depleted$avg_log2FC),]
dendritic.2.markers.enriched <- dendritic.2.markers.enriched[order(x = dendritic.2.markers.enriched$p_val_adj, y = -dendritic.2.markers.enriched$avg_log2FC),]
dendritic.2.markers.depleted <- dendritic.2.markers.depleted[order(x = dendritic.2.markers.depleted$p_val_adj, y = dendritic.2.markers.depleted$avg_log2FC),]
dendritic.3.markers.enriched <- dendritic.3.markers.enriched[order(x = dendritic.3.markers.enriched$p_val_adj, y = -dendritic.3.markers.enriched$avg_log2FC),]
dendritic.3.markers.depleted <- dendritic.3.markers.depleted[order(x = dendritic.3.markers.depleted$p_val_adj, y = dendritic.3.markers.depleted$avg_log2FC),]
macrophage.markers.enriched <- macrophage.markers.enriched[order(x = macrophage.markers.enriched$p_val_adj, y = -macrophage.markers.enriched$avg_log2FC),]
macrophage.markers.depleted <- macrophage.markers.depleted[order(x = macrophage.markers.depleted$p_val_adj, y = macrophage.markers.depleted$avg_log2FC),]

write.table(dendritic.1.markers.enriched, file = "immune_integration/DE_dendritic_enriched_(cluster7).csv", sep = "\t", quote = FALSE)
write.table(dendritic.1.markers.depleted, file = "immune_integration/DE_dendritic_depleted_(cluster7).csv", sep = "\t", quote = FALSE)
write.table(dendritic.2.markers.enriched, file = "immune_integration/DE_dendritic_enriched_(cluster9).csv", sep = "\t", quote = FALSE)
write.table(dendritic.2.markers.depleted, file = "immune_integration/DE_dendritic_depleted_(cluster9).csv", sep = "\t", quote = FALSE)
write.table(dendritic.3.markers.enriched, file = "immune_integration/DE_dendritic_enriched_(cluster10).csv", sep = "\t", quote = FALSE)
write.table(dendritic.3.markers.depleted, file = "immune_integration/DE_dendritic_depleted_(cluster10).csv", sep = "\t", quote = FALSE)
write.table(macrophage.markers.enriched, file = "immune_integration/DE_macrophage_enriched.csv", sep = "\t", quote = FALSE)
write.table(macrophage.markers.depleted, file = "immune_integration/DE_macrophage_depleted.csv", sep = "\t", quote = FALSE)

dendritic.1.2.markers <- FindMarkers(immune.combined, ident.1 = "Dendritic 1", ident.2 = "Dendritic 2")




