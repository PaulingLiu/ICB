library(tidyverse)
library(Seurat)

matr <- load10x("/raid1/pauling/projects/01_data/01_lung_immune_therapy/43.PXX.ZXZ.post/filtered_feature_bc_matrix.h5", pr = "P40.tr.1")

var.genes <- readr::read_rds("/raid1/pauling/projects/04_lung_immune_therapy/02_clustering/version1/12.features.rds.gz")

meta <- tibble(cellid = colnames(matr)) %>%
  dplyr::mutate(name = cellid) %>%
  tidyr::separate(name, c("patient","treatment","num","barcodes"), sep = "\\.") %>%
  dplyr::mutate(barcode = ifelse(is.na(barcodes), num, barcodes)) %>%
  dplyr::mutate(num = ifelse(is.na(barcodes), NA, num)) %>%
  dplyr::select(-barcodes)

meta <- as.data.frame(meta)
rownames(meta) <- meta$cellid
 
sce.pre <- CreateSeuratObject(counts = matr, meta.data = meta)
sce.pre <- NormalizeData(sce.pre, normalization.method = "LogNormalize", scale.factor = 10000)
sce.pre <- FindVariableFeatures(sce.pre, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sce.pre)
sce.pre <- ScaleData(sce.pre, features = all.genes, do.center = F, do.scale = F)
sce.pre <- RunPCA(sce.pre, features = var.genes)
sce.pre <- RunUMAP(sce.pre, dims = 1:10)
sce.pre <- FindNeighbors(sce.pre, dims = 1:10)
sce.pre <- FindClusters(sce.pre, resolution = 2)

DimPlot(sce.pre, label = T)
markers <- FindMarkers(sce.pre, ident.1 = 14, only.pos = T)
FeaturePlot(sce.filt, features = "LYZ")

sce.filt <- sce.pre[,!(Idents(sce.pre) %in% c(27,33,1,4,20,19,29,18,26,28,32))]
DimPlot(sce.filt, label = T) + geom_vline(xintercept = 0.1)

sce.filt <- sce.filt[,!(Idents(sce.filt) == 11 & sce.filt@reductions$umap@cell.embeddings[,1] > 0.1)]

FeaturePlot(sce.filt, features = "CLDN5")
p43.sce %>% readr::write_rds("projects/01_data/01_lung_immune_therapy/43.PXX.ZXZ.post/P43.sce.rds.gz", compress = "gz")

DimPlot(sce.filt)

p43.sce <- readr::read_rds("projects/01_data/01_lung_immune_therapy/43.PXX.ZXZ.post/P43.sce.rds.gz")
DimPlot(p43.sce)
FeaturePlot(p43.sce, features = "MMP1")

cellnames <- colnames(p43.sce)
cellnames.new <- stringr::str_replace(cellnames, "P40","P41")
p43.sce <- RenameCells(p43.sce, new.names = cellnames.new, old.names = cellnames)

p43.sce@meta.data$cellid <- rownames(p43.sce@meta.data)
p43.sce@meta.data$patient <- "P41"
