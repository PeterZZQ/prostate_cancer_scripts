rm(list = ls())
gc()

library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db)
library(Seurat)
library(SeuratDisk)
library(enrichplot)

setwd("/localscratch/ziqi/Prostate_Cancer")

clusters <- c("Luminal 1", "Hillock epithelia (Basal)", "Club epithelia (Luminal)", "Macrophage (Myeloid, Total immune)", "Monocytes (Myeloid, Total immune)", 
              "Mesenchymal", "Spink1+ (Luminal)", "Basal", "Lymphoid (Total immune)", "SV, Spink1+", "Endothelial", "Luminal 2", "Luminal 3")

# First create the gene universe. 
prostate.combined <- LoadH5Seurat("Cell_Ranger_output/seurat_integrate.h5Seurat")
# set background gene set
min.cells <- 200
num.cells <- rowSums(prostate.combined@assays$RNA@counts > 0)
background_gene <- names(num.cells[which(num.cells >= min.cells)])

for(cluster.id in clusters){
  print(cluster.id)
  print("GO ORA analysis...")
  # 1. GO ORA analysis on enriched de genes
  if(file.exists(paste0("integration/DE_genotype/DEGeneList/DE_", cluster.id, "_18wk_enriched.csv"))){
    # read de gene list
    de_gene_enriched <- rownames(read.csv(paste0("integration/DE_genotype/DEGeneList/DE_", cluster.id, "_18wk_enriched.csv"), sep = ","))
    # GO over-representation analysis
    ego_enriched <- enrichGO(gene = de_gene_enriched, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", universe = background_gene, 
                             ont = "ALL", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.2, readable = TRUE)
    # save results
    if(!is.null(ego_enriched)){
      write.table(ego_enriched@result, file = paste0("integration/DE_genotype/enrichment_analysis/GO_ORA_", cluster.id, "_18wk_enriched.csv"), sep = ",", quote = FALSE)
    }
  }
  
  # 2. GO ORA analysis on depleted de genes
  if(file.exists(paste0("integration/DE_genotype/DEGeneList/DE_", cluster.id, "_18wk_depleted.csv"))){
    # read de gene list
    de_gene_depleted <- rownames(read.csv(paste0("integration/DE_genotype/DEGeneList/DE_", cluster.id, "_18wk_depleted.csv"), sep = ","))
    # GO over-representation analysis
    ego_depleted <- enrichGO(gene = de_gene_depleted, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", universe = background_gene, 
                             ont = "ALL", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.2, readable = TRUE)
    # save results
    if(!is.null(ego_depleted)){
      write.table(ego_depleted@result, file = paste0("integration/DE_genotype/enrichment_analysis/GO_ORA_", cluster.id, "_18wk_depleted.csv"), sep = ",", quote = FALSE)  
    }
  }
  
  # 3. GO GSEA analysis on ranked gene
  print("GO GSEA analysis...")
  if(file.exists(paste0("integration/DE_genotype/RankedGeneList/rgl_", cluster.id, "_18wk.csv"))){
    # read ranked gene list  
    ranked_gene <- read.csv(paste0("integration/DE_genotype/RankedGeneList/rgl_", cluster.id, "_18wk.csv"), sep = ",")
    ranked_gene <- setNames(ranked_gene$avg_log2FC, rownames(ranked_gene))
    # GO GSEA, ranked by log2FC (not p-value): https://www.biostars.org/p/375584/
    ggo <- gseGO(geneList = ranked_gene, ont = "ALL", OrgDb = org.Mm.eg.db, keyType = "SYMBOL", pvalueCutoff = 0.05, eps = 0)
    # save reults
    if(!is.null(ggo)){
      write.table(ggo@result, file = paste0("integration/DE_genotype/enrichment_analysis/GO_GSEA_", cluster.id, "_18wk.csv"), sep = ",", quote = FALSE)
      geneSetID <- ggo@result[1:5,"ID"]
      # plot gsea enrichment score
      p1 <- gseaplot2(ggo, geneSetID = geneSetID)
      pdf(paste0("integration/DE_genotype/enrichment_analysis/GO_GSEA_", cluster.id, "_18wk.pdf"))
      print(p1)
      dev.off()
    }
  }
}

# 4. HM GSEA analysis on ranked gene
print("Hull mark gene GSEA analysis...")
# retrive hull mark gene set
hull_mark <- msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select(gs_name, entrez_gene)
for(cluster.id in clusters){
  # read ranked gene list
  if(file.exists(paste0("integration/DE_genotype/RankedGeneList/rgl_", cluster.id, "_18wk.csv"))){
    print(cluster.id)
    ranked_gene <- read.csv(paste0("integration/DE_genotype/RankedGeneList/rgl_", cluster.id, "_18wk.csv"), sep = ",")
    rownames(ranked_gene) <- ranked_gene$ENSEMBL
    ranked_gene.df <- bitr(ranked_gene$ENSEMBL, fromType = "ENSEMBL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
    ranked_gene <- setNames(ranked_gene[ranked_gene.df$ENSEMBL,"avg_log2FC"], ranked_gene.df$ENTREZID)
    # GSEA analysis
    ghullmark <- GSEA(geneList = ranked_gene, TERM2GENE = hull_mark)
    if(!is.null(ghullmark)){
      write.table(ghullmark@result, file = paste0("integration/DE_genotype/enrichment_analysis/HM_GSEA_", cluster.id, "_18wk.csv"), sep = ",", quote = FALSE) 
      geneSetID <- ghullmark@result[1:5,"ID"]
      # plot gsea enrichment score
      p1 <- gseaplot2(ghullmark, geneSetID = geneSetID)
      pdf(paste0("integration/DE_genotype/enrichment_analysis/HM_GSEA_", cluster.id, "_18wk.pdf"))
      print(p1)
      dev.off()
      
    }
  }
}