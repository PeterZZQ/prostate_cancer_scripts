rm(list=ls())
gc()

library("Seurat")
setwd("/localscratch/ziqi/Prostate_Cancer")
data_dir <- "spatial_data/Visium/"

SpatialExperiment::read10xVisium()