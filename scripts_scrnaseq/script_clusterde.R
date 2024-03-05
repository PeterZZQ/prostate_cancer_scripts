rm(list=ls())
gc()

setwd("/localscratch/ziqi/Prostate_Cancer")
library("Seurat")
library("anndata")
library("SeuratDisk")
library("patchwork")
library("ClusterDE")

prostate.combined <- LoadH5Seurat("Cell_Ranger_output/seurat_integrate.h5Seurat")
# filtered out not-expressing genes
counts <- GetAssayData(prostate.combined, slot="counts", assay="RNA")  
genes.percent.expressed <- rowMeans(counts>0 )*100  
genes.filter <- names(genes.percent.expressed[genes.percent.expressed>1])  #select genes expressed in at least 1% of cells
prostate.combined <- prostate.combined[genes.filter,]
use.hvg <- FALSE
# DE gene analysis across genotypes
# for(cluster.id in c("Luminal", "Basal", "Club epithelia", "Macrophage", "Monocytes", 
#                     "Mesenchymal", "Luminal (Spink1+)", "Lymphoid", "SV", "Endothelial")){
for(cluster.id in c("Luminal", "Basal")){
  print(cluster.id)
  prostate.cluster <- prostate.combined[,prostate.combined$annot == cluster.id]
  for(age in c("12wk", "18wk")){
    prostate.cluster.age <- prostate.cluster[,prostate.cluster$age == age]
    print(dim(prostate.cluster.age)[2])
    
    # check if there are enough cells
    if((sum(prostate.cluster.age$genotype == "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)") > 10) & (sum(prostate.cluster.age$genotype == "PbCre(+/-),Pten(-/-),P53(-/-)") > 10)){
      # extract the count matrix (UMI)
      counts.cluster.age <- GetAssay(prostate.cluster.age, assay = "RNA")@counts
      # find DE in the original real dataset
      if(use.hvg){
        original_markers <- FindMarkers(prostate.cluster.age, ident.1 = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ident.2 = "PbCre(+/-),Pten(-/-),P53(-/-)", group.by = "genotype", 
                                        features = prostate.combined@assays$integrated@var.features)
      }else{
        original_markers <- FindMarkers(prostate.cluster.age, ident.1 = "PbCre(+/-),Pten(-/-),P53(-/-),CXCR7(-/-)", ident.2 = "PbCre(+/-),Pten(-/-),P53(-/-)", group.by = "genotype")
      }
      # generate the synthetic null dataset
      synthetic_null <- constructNull(mat = counts.cluster.age, nCores = 8, fastVersion = FALSE, corrCut = 0.2)
      saveRDS(synthetic_null, file = paste0("results_seurat_scrnaseq/DE/clusterDE/synthetic_null_", cluster.id, "_", age, ".RDS"))
      synthetic_null <- readRDS(paste0("results_seurat_scrnaseq/DE/clusterDE/synthetic_null_", cluster.id, "_", age, ".RDS"))
      # synthetic_null <- constructNull(counts.cluster.age, nCores = 8, fastVersion = TRUE)
      # following the same clustering pipeline as real dataset
      synthetic_null_seurat <- CreateSeuratObject(counts = synthetic_null)
      synthetic_null_seurat <- NormalizeData(synthetic_null_seurat)
      synthetic_null_seurat <- FindVariableFeatures(synthetic_null_seurat, selection.method = "vst", nfeatures = 2000)
      synthetic_null_seurat <- ScaleData(object = synthetic_null_seurat)
      #> Centering and scaling data matrix
      synthetic_null_seurat <- RunPCA(object = synthetic_null_seurat, npcs = 30, verbose = FALSE)
      synthetic_null_seurat <- RunUMAP(synthetic_null_seurat, reduction = "pca", dims = 1:30, n.neighbors = 50L, min.dist = 0.3)
      synthetic_null_seurat <- FindNeighbors(synthetic_null_seurat, reduction = "pca", dims = 1:30, k.param = 50)
      # Obtain cluster result, make sure only two clusters
      synthetic_null_seurat <- FindClusters(synthetic_null_seurat, resolution = 0.3)
      # find DE in the synthetic null dataset
      if(use.hvg){
        null_markers <- FindMarkers(synthetic_null_seurat, ident.1 = 0, ident.2 = 1, min.pct = 0, logfc.threshold = 0, features = prostate.combined@assays$integrated@var.features)
      }else{
        null_markers <- FindMarkers(synthetic_null_seurat, ident.1 = 0, ident.2 = 1, min.pct = 0, logfc.threshold = 0)
      }
      
      # adjusted p-values
      original_pval <- original_markers$p_val
      names(original_pval) <- rownames(original_markers)
      null_pval <- null_markers$p_val
      names(null_pval) <- rownames(null_markers)
      res <- ClusterDE::callDE(original_pval, null_pval, nlogTrans = TRUE)
      head(res$summaryTable)
      if(use.hvg){
        write.table(res$summaryTable, file = paste0("results_seurat_scrnaseq/DE/clusterDE/DE_", cluster.id, "_", age, ".csv"), sep = ",", quote = FALSE)
      }else{
        write.table(res$summaryTable, file = paste0("results_seurat_scrnaseq/DE/clusterDE/DE_", cluster.id, "_", age, "_full.csv"), sep = ",", quote = FALSE)
      }
      
    }
  }
}


