rm(list=ls())
gc()

setwd("/localscratch/ziqi/Prostate_Cancer")
library("Seurat")
library("anndata")
library("SeuratDisk")
library("patchwork")

raw_dir <- "dataset/data_scrnaseq/data/qc_data/"

# 1. Run Seurat integration ----------------------------------------------------
# Load the merged anndata object (obtained from preprocess.py)
X_intact <- Matrix::readMM(paste0(raw_dir, "X_intact.mtx"))
X_intact <- as(t(as.matrix(X_intact)), "dgCMatrix")
meta_cells_intact <- read.csv(paste0(raw_dir, "meta_intact.csv"), sep = ",", row.names = 1)
meta_genes_intact <- read.csv(paste0(raw_dir, "gene_intact.csv"), sep = ",", row.names = 1)
colnames(X_intact) <- rownames(meta_cells_intact)
rownames(X_intact) <- rownames(meta_genes_intact)
# Create the Seurat object
prostate.intact <- CreateSeuratObject(counts = X_intact, meta.data = meta_cells_intact, assay = "RNA")
prostate.intact[["RNA"]] <- AddMetaData(prostate.intact[["RNA"]], meta_genes_intact)

# split the dataset into layers corresponding to different samples
prostate.intact[["RNA"]] <- split(prostate.intact[["RNA"]], f = prostate.intact$sample)

# normalize and identify variable features for each dataset independently
prostate.intact <- NormalizeData(prostate.intact)
prostate.intact <- FindVariableFeatures(prostate.intact, selection.method = "vst", nfeatures = 2000)
prostate.intact <- ScaleData(prostate.intact)

# calculate the PCA of the data
prostate.intact <- RunPCA(prostate.intact)

# run Seurat integration, original CCA method
prostate.intact <- IntegrateLayers(object = prostate.intact, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)

# Obtain cluster result
prostate.intact <- FindNeighbors(prostate.intact, reduction = "integrated.cca", dims = 1:30)
prostate.intact <- FindClusters(prostate.intact, resolution = 0.3, cluster.name = "cca_clusters")

# visualize the clusters
prostate.intact <- RunUMAP(prostate.intact, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
p1 <- DimPlot(prostate.intact, reduction = "umap.cca", group.by = c("sample", "cca_clusters"), combine = FALSE, label.size = 2)

# save files
seurat_dir <- "dataset/data_scrnaseq/"
SaveH5Seurat(prostate.intact, filename = paste0(seurat_dir, "prostate_intact.h5Seurat"), overwrite = TRUE)
meta.seurat <- prostate.intact@meta.data
pca_cca.seurat <- prostate.intact@reductions$integrated.cca@cell.embeddings
umap_cca.seurat <- prostate.intact@reductions$umap.cca@cell.embeddings
write.table(meta.seurat, file = paste0(seurat_dir, "meta_intact.csv"), sep = "\t", quote = FALSE)
write.table(pca_cca.seurat, file = paste0(seurat_dir, "cca_pca_intact.csv"), sep = "\t", quote = FALSE)
write.table(umap_cca.seurat, file = paste0(seurat_dir, "cca_umap_intact.csv"), sep = "\t", quote = FALSE)


# # 2. Transfer labels to the castrated dataset ----------------------------------
# # Load the query dataset: castrated dataset
# X_castrated <- Matrix::readMM(paste0(raw_dir, "X_castrated.mtx"))
# X_castrated <- as(t(as.matrix(X_castrated)), "dgCMatrix")
# meta_cells_castrated <- read.csv(paste0(raw_dir, "meta_castrated.csv"), sep = ",", row.names = 1)
# meta_genes_castrated <- read.csv(paste0(raw_dir, "gene_castrated.csv"), sep = ",", row.names = 1)
# colnames(X_castrated) <- rownames(meta_cells_castrated)
# rownames(X_castrated) <- rownames(meta_genes_castrated)
# # Create the Seurat object, two samples
# prostate.castrated <- CreateSeuratObject(counts = X_castrated, meta.data = meta_cells_castrated, assay = "RNA")
# prostate.castrated.pp <- prostate.castrated[,prostate.castrated$sample == "M1476-PP"]
# prostate.castrated.ppc <- prostate.castrated[,prostate.castrated$sample == "M1477-PPC"]

# # Load the reference dataset: intact dataset at 18wk
# prostate.intact.pp18 <- prostate.intact[,prostate.intact$sample == "M1416-PP18"]
# prostate.intact.ppc18 <- prostate.intact[,prostate.intact$sample == "M1437-PPC18"]

# # run label transfer: pp
# prostate.castrated.pp <- NormalizeData(prostate.castrated.pp)
# prostate.anchors <- FindTransferAnchors(reference = prostate.intact.pp18, query = prostate.castrated.pp, dims = 1:30, reference.reduction = "pca")
# predictions <- TransferData(anchorset = prostate.anchors, refdata = prostate.intact.pp18$seurat_clusters, dims = 1:30)
# prostate.castrated.pp <- AddMetaData(prostate.castrated.pp, metadata = predictions)

# # run label transfer: ppc
# prostate.castrated.ppc <- NormalizeData(prostate.castrated.ppc)
# prostate.anchors <- FindTransferAnchors(reference = prostate.intact.ppc18, query = prostate.castrated.ppc, dims = 1:30, reference.reduction = "pca")
# predictions <- TransferData(anchorset = prostate.anchors, refdata = prostate.intact.ppc18$seurat_clusters, dims = 1:30)
# prostate.castrated.ppc <- AddMetaData(prostate.castrated.ppc, metadata = predictions)

# meta_cells_castrated[colnames(prostate.castrated.pp), "seurat_cluster"] <- prostate.castrated.pp$predicted.id
# meta_cells_castrated[colnames(prostate.castrated.ppc), "seurat_cluster"] <- prostate.castrated.ppc$predicted.id
# write.table(meta_cells_castrated, file = paste0(seurat_dir, "meta_castrated_transfer.csv"), sep = "\t", quote = FALSE)



# 5. differential expression analysis across cell types ------------------------
# differential expression analysis of all clusters
# for(cluster.id in seq(0,13)){
#   cluster.marker <- FindMarkers(prostate.combined, ident.1 = cluster.id, ident.2 = NULL, features = prostate.combined@assays$integrated@var.features)
#   cluster.marker <- cluster.marker[cluster.marker$p_val_adj < 0.05,]
#   cluster.marker.enriched <- cluster.marker[cluster.marker$avg_log2FC >= 0,]
#   cluster.marker.depleted <- cluster.marker[cluster.marker$avg_log2FC < 0,]
#   cluster.marker.enriched <- cluster.marker.enriched[order(x = cluster.marker.enriched$p_val_adj, y = -cluster.marker.enriched$avg_log2FC),]
#   cluster.marker.depleted <- cluster.marker.depleted[order(x = cluster.marker.depleted$p_val_adj, y = cluster.marker.depleted$avg_log2FC),]
#   write.table(cluster.marker.enriched, file = paste0("results_seurat_scrnaseq/DE/DEGeneList/DE_cluster", cluster.id, "_enriched.csv"), sep = ",", quote = FALSE)
#   write.table(cluster.marker.depleted, file = paste0("results_seurat_scrnaseq/DE/DEGeneList/DE_cluster", cluster.id, "_depleted.csv"), sep = ",", quote = FALSE)
# }

prostate.combined <- LoadH5Seurat("Cell_Ranger_output/seurat_integrate.h5Seurat")
prostate.combined <- SetIdent(prostate.combined, value = "annot")
DimPlot(prostate.combined, label = TRUE)

print("calculate de gene list...")
for(cluster.id in c("Luminal", "Basal", "Club epithelia", "Macrophage", "Monocytes", 
                 "Mesenchymal", "Luminal (Spink1+)", "Basal", "Lymphoid", "SV", "Endothelial")){
  print(cluster.id)
  # de genes selected from top-2000 highly-variable genes
  cluster.marker <- FindMarkers(prostate.combined, ident.1 = cluster.id, ident.2 = NULL, features = prostate.combined@assays$integrated@var.features)
  cluster.marker <- cluster.marker[cluster.marker$p_val_adj < 0.05,]
  cluster.marker.enriched <- cluster.marker[cluster.marker$avg_log2FC >= 0,]
  cluster.marker.depleted <- cluster.marker[cluster.marker$avg_log2FC < 0,]
  cluster.marker.enriched <- cluster.marker.enriched[order(x = cluster.marker.enriched$p_val_adj, y = -cluster.marker.enriched$avg_log2FC),]
  cluster.marker.depleted <- cluster.marker.depleted[order(x = cluster.marker.depleted$p_val_adj, y = cluster.marker.depleted$avg_log2FC),]
  write.table(cluster.marker.enriched, file = paste0("results_seurat_scrnaseq/DE/DEGeneList/DE_", cluster.id, "_enriched.csv"), sep = ",", quote = FALSE)
  write.table(cluster.marker.depleted, file = paste0("results_seurat_scrnaseq/DE/DEGeneList/DE_", cluster.id, "_depleted.csv"), sep = ",", quote = FALSE)
}

# filter background genes
min.cells <- 200
num.cells <- rowSums(prostate.combined@assays$RNA@counts > 0)
genes.filtered <- names(num.cells[which(num.cells >= min.cells)])
# ranked gene list: log fold change
print("calculated ranked gene list...")
for(cluster.id in c("Luminal", "Basal", "Club epithelia", "Macrophage", "Monocytes", 
                    "Mesenchymal", "Luminal (Spink1+)", "Lymphoid", "SV", "Endothelial")){
  print(cluster.id)
  # selected from filtered background gene
  expr.fold.change <- FoldChange(prostate.combined, ident.1 = cluster.id, ident.2 = NULL, features = genes.filtered)
  expr.fold.change["ENSEMBL"] <- prostate.h5ad$var[rownames(expr.fold.change), "gene_ids"]
  expr.fold.change <- expr.fold.change[order(x = -expr.fold.change$avg_log2FC),]
  write.table(expr.fold.change, file = paste0("results_seurat_scrnaseq/DE/RankedGeneList/rgl_", cluster.id, ".csv"), sep = ",", quote = FALSE)
}



# 6. differential-expression genes across genotypes-----------------------------
# intact
prostate.intact <- LoadH5Seurat("Cell_Ranger_output/seurat_integrate.h5Seurat")
prostate.intact <- SetIdent(prostate.intact, value = "annot")
DimPlot(prostate.intact, label = TRUE)
# castrated
X_castrated <- Matrix::readMM("Cell_Ranger_output/X_castrated.mtx")
meta_castrated <- read.csv("Cell_Ranger_output/meta_castrated.csv", row.names = 1, sep = "\t")
gene_castrated <- read.csv("Cell_Ranger_output/gene_castrated.csv", row.names = 1, sep = ",")
rownames(X_castrated) <- rownames(meta_castrated)
colnames(X_castrated) <- rownames(gene_castrated)
prostate.castrated <- CreateSeuratObject(counts = t(X_castrated), meta.data = meta_castrated)
prostate.castrated <- SetIdent(prostate.castrated, value = "annot.transfer")
# find variable features
prostate.castrated <- NormalizeData(prostate.castrated)
prostate.castrated <- FindVariableFeatures(prostate.castrated, selection.method = "vst", nfeatures = 2000)

# INTACT, filter background gene
min.cells <- 200
num.cells <- rowSums(prostate.intact@assays$RNA@counts > 0)
genes.filtered <- names(num.cells[which(num.cells >= min.cells)])
for(cluster.id in c("Luminal", "Basal", "Club epithelia", "Macrophage", "Monocytes", 
                    "Mesenchymal", "Luminal (Spink1+)", "Lymphoid", "SV", "Endothelial")){
  
  print(cluster.id)
  prostate.cluster <- prostate.intact[,prostate.intact@active.ident == cluster.id]
  prostate.12wk <- prostate.cluster[,prostate.cluster$age == "12wk"]
  prostate.18wk <- prostate.cluster[,prostate.cluster$age == "18wk"]
  print(dim(prostate.12wk)[2])
  print(dim(prostate.18wk)[2])
  
  if((sum(prostate.12wk$genotype == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") > 10) & (sum(prostate.12wk$genotype == "PbCre(+/-),Pten(-/-),P53(-/-)") > 10)){
    de.12wk <- FindMarkers(prostate.12wk, ident.1 = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ident.2 = "PbCre(+/-),Pten(-/-),P53(-/-)", group.by = "genotype", 
                           features = prostate.intact@assays$integrated@var.features)
    de.12wk <- de.12wk[de.12wk$p_val_adj < 0.05,]
    de.12wk.enriched <- de.12wk[de.12wk$avg_log2FC >= 0,]
    de.12wk.depleted <- de.12wk[de.12wk$avg_log2FC < 0,]
    de.12wk.enriched <- de.12wk.enriched[order(x = de.12wk.enriched$p_val_adj, y = -de.12wk.enriched$avg_log2FC),]
    de.12wk.depleted <- de.12wk.depleted[order(x = de.12wk.depleted$p_val_adj, y = de.12wk.depleted$avg_log2FC),]
    write.table(de.12wk.enriched, file = paste0("results_seurat_scrnaseq/DE_genotype/DEGeneList/DE_", cluster.id, "_12wk_enriched.csv"), sep = ",", quote = FALSE)
    write.table(de.12wk.depleted, file = paste0("results_seurat_scrnaseq/DE_genotype/DEGeneList/DE_", cluster.id, "_12wk_depleted.csv"), sep = ",", quote = FALSE)
    
    expr.fold.change <- FoldChange(prostate.12wk, ident.1 = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ident.2 = "PbCre(+/-),Pten(-/-),P53(-/-)", group.by = "genotype", features = genes.filtered)
    expr.fold.change["ENSEMBL"] <- prostate.h5ad$var[rownames(expr.fold.change), "gene_ids"]
    expr.fold.change <- expr.fold.change[order(x = -expr.fold.change$avg_log2FC),]
    write.table(expr.fold.change, file = paste0("results_seurat_scrnaseq/DE_genotype/RankedGeneList/rgl_", cluster.id, "_12wk.csv"), sep = ",", quote = FALSE)
  }
  
  if((sum(prostate.18wk$genotype == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") > 10) & (sum(prostate.18wk$genotype == "PbCre(+/-),Pten(-/-),P53(-/-)") > 10)){
    de.18wk <- FindMarkers(prostate.18wk, ident.1 = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ident.2 = "PbCre(+/-),Pten(-/-),P53(-/-)", group.by = "genotype", 
                           features = prostate.intact@assays$integrated@var.features)
    de.18wk <- de.18wk[de.18wk$p_val_adj < 0.05,]
    de.18wk.enriched <- de.18wk[de.18wk$avg_log2FC >= 0,]
    de.18wk.depleted <- de.18wk[de.18wk$avg_log2FC < 0,]
    de.18wk.enriched <- de.18wk.enriched[order(x = de.18wk.enriched$p_val_adj, y = -de.18wk.enriched$avg_log2FC),]
    de.18wk.depleted <- de.18wk.depleted[order(x = de.18wk.depleted$p_val_adj, y = de.18wk.depleted$avg_log2FC),]
    write.table(de.18wk.enriched, file = paste0("results_seurat_scrnaseq/DE_genotype/DEGeneList/DE_", cluster.id, "_18wk_enriched.csv"), sep = ",", quote = FALSE)
    write.table(de.18wk.depleted, file = paste0("results_seurat_scrnaseq/DE_genotype/DEGeneList/DE_", cluster.id, "_18wk_depleted.csv"), sep = ",", quote = FALSE)
    
    expr.fold.change <- FoldChange(prostate.18wk, ident.1 = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ident.2 = "PbCre(+/-),Pten(-/-),P53(-/-)", group.by = "genotype", features = genes.filtered)
    expr.fold.change["ENSEMBL"] <- prostate.h5ad$var[rownames(expr.fold.change), "gene_ids"]
    expr.fold.change <- expr.fold.change[order(x = -expr.fold.change$avg_log2FC),]
    write.table(expr.fold.change, file = paste0("results_seurat_scrnaseq/DE_genotype/RankedGeneList/rgl_", cluster.id, "_18wk.csv"), sep = ",", quote = FALSE)
    
  }
}


# CASTRATED, filter background gene
min.cells <- 200
num.cells <- rowSums(prostate.castrated@assays$RNA@counts > 0)
genes.filtered <- names(num.cells[which(num.cells >= min.cells)])
for(cluster.id in c("Luminal", "Basal", "Club epithelia", "Macrophage", "Monocytes", 
                    "Mesenchymal", "Luminal (Spink1+)", "Lymphoid", "SV", "Endothelial")){
  
  print(cluster.id)
  prostate.cluster <- prostate.castrated[,prostate.castrated@active.ident == cluster.id]
  print(dim(prostate.cluster)[1])

  if((sum(prostate.cluster$genotype == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") > 10) & (sum(prostate.cluster$genotype == "PbCre(+/-),Pten(-/-),P53(-/-)") > 10)){
    # ident.2: reference, the control
    de <- FindMarkers(prostate.cluster, ident.1 = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ident.2 = "PbCre(+/-),Pten(-/-),P53(-/-)", group.by = "genotype", 
                           features = prostate.castrated@assays$RNA@var.features)
    de <- de[de$p_val_adj < 0.05,]
    de.enriched <- de[de$avg_log2FC >= 0,]
    de.depleted <- de[de$avg_log2FC < 0,]
    de.enriched <- de.enriched[order(x = de.enriched$p_val_adj, y = -de.enriched$avg_log2FC),]
    de.depleted <- de.depleted[order(x = de.depleted$p_val_adj, y = de.depleted$avg_log2FC),]
    write.table(de.enriched, file = paste0("results_seurat_scrnaseq/DE_genotype/castrated/DEGeneList/DE_", cluster.id, "_12wk_enriched.csv"), sep = ",", quote = FALSE)
    write.table(de.depleted, file = paste0("results_seurat_scrnaseq/DE_genotype/castrated/DEGeneList/DE_", cluster.id, "_12wk_depleted.csv"), sep = ",", quote = FALSE)
    
    expr.fold.change <- FoldChange(prostate.cluster, ident.1 = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ident.2 = "PbCre(+/-),Pten(-/-),P53(-/-)", group.by = "genotype", features = genes.filtered)
    expr.fold.change["ENSEMBL"] <- prostate.h5ad$var[rownames(expr.fold.change), "gene_ids"]
    expr.fold.change <- expr.fold.change[order(x = -expr.fold.change$avg_log2FC),]
    write.table(expr.fold.change, file = paste0("results_seurat_scrnaseq/DE_genotype/castrated/RankedGeneList/rgl_", cluster.id, "_12wk.csv"), sep = ",", quote = FALSE)
  }
}


# 7. transfer labels from the reference scRNA-seq dataset ----------------------
prostate.reference.h5ad <- read_h5ad("reference_data/reference_raw.h5ad")
# select batch T00_intact_1
prostate.reference.h5ad <- prostate.reference.h5ad[prostate.reference.h5ad$obs$batchID == "T00_intact_1",]
prostate.reference.obj <- CreateSeuratObject(counts = t(prostate.reference.h5ad$X), meta.data = prostate.reference.h5ad$obs)
prostate.reference.obj <- NormalizeData(prostate.reference.obj)
prostate.reference.obj <- FindVariableFeatures(prostate.reference.obj, selection.method = "vst", nfeatures = 2000)
# IntegrateData
DefaultAssay(prostate.reference.obj) <- "RNA"
# Run the standard workflow for visualization and clustering
prostate.reference.obj <- ScaleData(prostate.reference.obj, verbose = FALSE)
prostate.reference.obj <- RunPCA(prostate.reference.obj, npcs = 30, verbose = FALSE)
prostate.reference.obj <- RunUMAP(prostate.reference.obj, reduction = "pca", dims = 1:30, verbose = FALSE)
# p1 <- DimPlot(prostate.reference.obj, reduction = "umap", group.by = "IntType")
# p2 <- DimPlot(prostate.reference.obj, reduction = "umap", group.by = "batchID")
# p1 + p2

prostate.h5ad <- read_h5ad("Cell_Ranger_output/adata_merge.h5ad")
prostate.obj <- CreateSeuratObject(counts = t(prostate.h5ad$X), meta.data = prostate.h5ad$obs)
# split the dataset into a list of seurat objects by batches
prostate.obj.list <- SplitObject(prostate.obj, split.by = "sample")

# normalize and identify variable features for each dataset independently
prostate.obj.list <- lapply(X = prostate.obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x.anchors <- FindTransferAnchors(reference = prostate.reference.obj, query = x, dims = 1:30, reference.reduction = "pca")
  predictions <- TransferData(anchorset = x.anchors, refdata = prostate.reference.obj$IntType, dims = 1:30)
  x <- AddMetaData(x, metadata = predictions)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = prostate.obj.list)
prostate.anchors <- FindIntegrationAnchors(object.list = prostate.obj.list, anchor.features = features)
prostate.combined <- IntegrateData(anchorset = prostate.anchors)
DefaultAssay(prostate.combined) <- "integrated"

prostate.combined <- ScaleData(prostate.combined, verbose = FALSE)
prostate.combined <- RunPCA(prostate.combined, npcs = 30, verbose = FALSE)
prostate.combined <- RunUMAP(prostate.combined, reduction = "pca", dims = 1:30, n.neighbors = 50L, min.dist = 0.3)

prostate.combined2 <- LoadH5Seurat("Cell_Ranger_output/seurat_integrate.h5Seurat")
meta.data <- prostate.combined2@meta.data
meta.data$transfer.label <- prostate.combined@meta.data$predicted.id 
prostate.h5ad$obs$transfer_label <- meta.data[rownames(prostate.h5ad$obs), "transfer.label"] 
prostate.h5ad$obs$seurat_clusters <- meta.data[rownames(prostate.h5ad$obs), "seurat_clusters"]
prostate.h5ad$obs$annot <- meta.data[rownames(prostate.h5ad$obs), "annot"]
prostate.h5ad$obs$highres_annot <- meta.data[rownames(prostate.h5ad$obs), "highres.annot"]
prostate.combined2@meta.data$transfer.label <- meta.data$transfer.label
write_h5ad(prostate.h5ad, filename = "Cell_Ranger_output/adata_seurat.h5ad")
SaveH5Seurat(prostate.combined2, filename = "Cell_Ranger_output/seurat_integrate.h5Seurat", overwrite = TRUE)
write.table(meta.data, file = "Cell_Ranger_output/seurat_meta.txt", sep = "\t", quote = FALSE)

# Visualization
p1 <- DimPlot(prostate.combined, reduction = "umap", group.by = "sample")
p2 <- DimPlot(prostate.combined, reduction = "umap", group.by = "predicted.id")
p1 + p2
ggsave(filename = "integration/transfer_label.pdf", device = "pdf", plot = p1 + p2, width = 12, height = 5)
ggsave(filename = "integration/transfer_label.png", device = "png", plot = p1 + p2, width = 12, height = 5, dpi = 150)


prostate.m1416_pp18 <- prostate.combined[,prostate.combined$sample == "M1416-PP18"]
prostate.m1417_pp12 <- prostate.combined[,prostate.combined$sample == "M1417-PP12"]
prostate.m1436_ppc12 <- prostate.combined[,prostate.combined$sample == "M1436-PPC12"]
prostate.m1437_ppc18 <- prostate.combined[,prostate.combined$sample == "M1437-PPC18"]

p1 <- DimPlot(prostate.m1416_pp18, reduction = "umap", group.by = "predicted.id")
p2 <- DimPlot(prostate.m1417_pp12, reduction = "umap", group.by = "predicted.id")
p3 <- DimPlot(prostate.m1436_ppc12, reduction = "umap", group.by = "predicted.id")
p4 <- DimPlot(prostate.m1437_ppc18, reduction = "umap", group.by = "predicted.id")
p1 + p2 + p3 + p4
ggsave(filename = "integration/transfer_label_separate.pdf", device = "pdf", plot = p1 + p2 + p3 + p4, width = 12, height = 8)
ggsave(filename = "integration/transfer_label_separate.png", device = "png", plot = p1 + p2 + p3 + p4, width = 12, height = 8, dpi = 150)





# # ------------------------------------------------------
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




