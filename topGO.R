rm(list = ls())
gc()

# https://jdblischak.github.io/nw/analysis/mouse/go.html
library(topGO)
library(org.Mm.eg.db)

setwd("/localscratch/ziqi/Prostate_Cancer")

# First create the gene universe. 
# This is all the genes tested for differential expression assigned a 1 for differentially expressed and 0 if not.
prostate.combined <- LoadH5Seurat("Cell_Ranger_output/seurat_integrate.h5Seurat")
background_gene <- prostate.combined@assays$integrated@var.features
# immune.combined <- LoadH5Seurat("Cell_Ranger_output/seurat_integrate_immune.h5Seurat")
# background_gene <- immune.combined@assays$integrated@var.features

de_gene <- rownames(read.csv("integration/DE/DE_cluster13_depleted.csv", sep = ","))[1:150]
gene_universe <- as.integer(background_gene%in%de_gene)
gene_universe <- factor(gene_universe)
names(gene_universe) <- background_gene

# Create the topGO data object. 
# Only consider “Biological Process” categories(which is mostly the case) and use the Mouse Ensembl database for annotation.
go_data <- new("topGOdata",
               ontology = "BP",
               allGenes = gene_universe,
               nodeSize = 5,
               annotationFun = annFUN.org,
               mapping = "org.Mm.eg.db",
               ID = "symbol")

# performing enrichment test
# Use the classic algorithm and score the tests with Fisher’s exact test.

# https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf see all algorithm
# default is weight01
go_test <- runTest(go_data, algorithm = "weight01", statistic = "fisher")

# analysis of the result
# Keep the results with a Fisher’s exact test p-value < 0.05.
go_table <- GenTable(go_data, weightFisher = go_test,
                     orderBy = "weightFisher", ranksOf = "weightFisher",
                     topNodes = sum(score(go_test) < .05),
                     numChar=1000)
head(go_table)

write.table(go_table, file = paste0("integration/DE/GO_cluster13_depleted.csv"), sep = ",", quote = FALSE)
