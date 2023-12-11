
#--- library

library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(ggsci)
library(data.table)
library(reticulate)
use_python("/home/pauling/anaconda3/bin/python")

#--- load data

sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/01.all.sce.rds.gz")
mda <- readr::read_tsv("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/06.tumor.classification/04.NM.DL.tsv")

#DimPlot(sce, label = T)
#FeaturePlot(sce, features = "CD4")
new.cluster.ids.sort <- c("Cancer cell","Myeloid","CD8 T","CD4 T","Prolif. T","CD4 T","NK","MAST","CD4 T","CD8 T","Myeloid","Myeloid","Myeloid","B cell","CAF","Cancer cell","Plasma B","Cancer cell","Plasma B","Plasma B","Plasma B", "Myeloid","pDC","Endo.")
names(new.cluster.ids) <- levels(sce)
all.sce.pro <- RenameIdents(sce, new.cluster.ids)
DimPlot(all.sce.pro, label = T, pt.size = 0.001) + scale_color_d3(palette = "category20c") #+ NoLegend() + NoAxes()

#--- T cells

DimPlot(sce, label = T) + NoLegend() + NoAxes()
VlnPlot(sce, features = "CXCL13")

#--- CD4 T cells

cd4.sce <- sce[,Idents(sce) %in% c(3,10,11)]
cd4.sce <- NormalizeData(cd4.sce)
cd4.sce <- ScaleData(cd4.sce, features = rownames(cd4.sce), do.center = F, do.scale = F)
cd4.sce <- FindVariableFeatures(cd4.sce, selection.method = "vst", nfeatures = 2000)
cd4.sce <- RunPCA(cd4.sce, features = VariableFeatures(cd4.sce))
cd4.sce <- bbknn.batch(cd4.sce, r = 1, n.pc = 15)

DimPlot(sce, label = T) + geom_abline(slope = -1.5, intercept = 8.5)
FeaturePlot(cd4.sce, features = "IFNG")
treg.mda <- tibble(
  cellid = colnames(cd4.sce[,Idents(cd4.sce) %in% c(0,8)]),
  cluster = "30"
)

new.idents <- tibble(
  cellid = colnames(sce),
  cluster = as.character(Idents(sce)),
  UMAP1 = unlist(sce@reductions$umap@cell.embeddings[,1]),
  UMAP2 = unlist(sce@reductions$umap@cell.embeddings[,2])
) %>%
  dplyr::mutate(cluster = ifelse(cellid %in% treg.mda$cellid, "30", cluster)) %>%
  dplyr::mutate(cluster = ifelse(cluster == 17 & (UMAP1*-1.5+8.5 > UMAP2), "31", cluster))

sce.pro <- sce
Idents(sce.pro) <- new.idents$cluster

new.idents %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/02.add.Tcell.sub.annotation.rds.gz", compress = "gz")
new.idents <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/02.add.Tcell.sub.annotation.rds.gz")
all.sce.pro <- sce
Idents(all.sce.pro) <- new.idents$cluster
DimPlot(all.sce.pro, label = T)

FeaturePlot(all.sce.pro, features = "LILRA4")

new.cluster.ids <- levels(all.sce.pro)
new.cluster.ids[new.cluster.ids %in% c(0,8,5)] <- "CD8 T"
new.cluster.ids[new.cluster.ids %in% c(7)] <- "NK"
new.cluster.ids[new.cluster.ids %in% c(3,10,30,11,13)] <- "CD4 T"
new.cluster.ids[new.cluster.ids %in% c(17,31)] <- "Prolif. T"
new.cluster.ids[new.cluster.ids %in% c(2,6,9,21,24)] <- "Cancer cell"
new.cluster.ids[new.cluster.ids %in% c(1,12,16)] <- "Myeloid"
new.cluster.ids[new.cluster.ids %in% c(4)] <- "B cell"
new.cluster.ids[new.cluster.ids %in% c(15,19,20,23)] <- "Plasma B"
new.cluster.ids[new.cluster.ids %in% c(22)] <- "MAST"
new.cluster.ids[new.cluster.ids %in% c(14)] <- "CAF"
new.cluster.ids[new.cluster.ids %in% c(18)] <- "Endo."

names(new.cluster.ids) <- levels(all.sce.pro)
all.sce.pro <- RenameIdents(all.sce.pro, new.cluster.ids)
levels(all.sce.pro) <- unique(new.cluster.ids.sort)
DimPlot(all.sce.pro, label = T)

#--- macrophages

DimPlot(sce, label = T)
mye.sce <- sce[,Idents(sce) %in% c(1,12,16)]
DimPlot(mye.sce, label = T, pt.size = 0.001) + scale_color_d3(palette = "category20") #+ NoLegend() + NoAxes()
FeaturePlot(mye.sce, features = "CLEC10A")

mac.sce <- mye.sce
mac.sce <- NormalizeData(mac.sce)
mac.sce <- ScaleData(mac.sce, features = rownames(mac.sce), do.center = F, do.scale = F)
mac.sce <- FindVariableFeatures(mac.sce, selection.method = "vst", nfeatures = 2000)
mac.sce <- RunPCA(mac.sce, features = VariableFeatures(mac.sce))
mac.sce <- bbknn.batch(mac.sce, r = 2.5, n.pc = 30)

DimPlot(mac.sce, label = T)
FeaturePlot(mac.sce, features = "FCN1")

new.cluster.ids <- levels(mac.sce)
new.cluster.ids[new.cluster.ids %in% c(1,5,18)] <- "FCN1"
new.cluster.ids[new.cluster.ids %in% c(0,11,19,6,4,25,13,16)] <- "SPP1"
new.cluster.ids[new.cluster.ids %in% c(2,10,14,9,21)] <- "FOLR2"
new.cluster.ids[new.cluster.ids %in% c(3,20,24,12,17)] <- "C1QA"
new.cluster.ids[new.cluster.ids %in% c(7)] <- "DC2.CD1A"
new.cluster.ids[new.cluster.ids %in% c(15)] <- "DC2.CD1C"
new.cluster.ids[new.cluster.ids %in% c(22)] <- "DC1.LAMP3.DC"
new.cluster.ids[new.cluster.ids %in% c(23)] <- "pDC"
new.cluster.ids[new.cluster.ids %in% c(8)] <- "MARCO"

mac.pro <- mac.sce
names(new.cluster.ids) <- levels(mac.sce)
mac.pro <- RenameIdents(mac.sce, new.cluster.ids)
DimPlot(mac.pro, label = T)

#DC FCN1 SPP1 FOLR2 C1QA MARCO

mac.pro %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/03.myeloid.rds.gz", compress = "gz")
mac.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/03.myeloid.rds.gz")

DimPlot(mac.sce, label = T)
FeaturePlot(mac.sce, features = "CD1A") + geom_abline(slope = -1, intercept = -8.2)

mac.idents <- tibble(
  UMAP1 = unlist(mac.sce@reductions$umap@cell.embeddings[,1]),
  UMAP2 = unlist(mac.sce@reductions$umap@cell.embeddings[,2]),
  label = as.character(Idents(mac.sce))
) %>%
  dplyr::mutate(label = ifelse(label == "DC1.LAMP3.DC" & (-1*UMAP1-8.2 >= UMAP2), "DC1", label)) %>%
  dplyr::mutate(label = ifelse(label == "DC1.LAMP3.DC" & (-1*UMAP1-8.2 < UMAP2), "LAMP3.DC", label))

Idents(mac.sce) <- mac.idents$label
DimPlot(mac.sce, label = T)

#--- endothelials

DimPlot(sce, label = T)
endo.sce <- sce[,Idents(sce) %in% c(18)]

endo.sce <- NormalizeData(endo.sce)
endo.sce <- ScaleData(endo.sce, features = rownames(endo.sce), do.center = F, do.scale = F)
endo.sce <- FindVariableFeatures(endo.sce, selection.method = "vst", nfeatures = 2000)
endo.sce <- RunPCA(endo.sce, features = VariableFeatures(endo.sce))
endo.sce <- bbknn.batch(endo.sce, r = 1, n.pc = 30)

DimPlot(endo.sce, label = T)

markers <- FindMarkers(endo.sce, ident.1 = 4, only.pos = T)
FeaturePlot(endo.sce, features = "INSR")
genes <- c("ACKR1","CCL14","ADIRF","MT-CO1","MT-CO3","CCL21","PDPN","PROX1","IGFBP5","PTGDS","FLT1","INSR","KDR","PDXN","PGF","APLN","STMN1","MKI67","TYMS")

new.cluster.ids <- levels(endo.sce)
new.cluster.ids[new.cluster.ids %in% c(0,4)] <- "ACKR1"
new.cluster.ids[new.cluster.ids %in% c(5)] <- "lymp"
new.cluster.ids[new.cluster.ids %in% c(1)] <- "INSR"
new.cluster.ids[new.cluster.ids %in% c(3,2,6,7)] <- "other"


names(new.cluster.ids) <- levels(endo.sce)
endo.pro <- RenameIdents(endo.sce, new.cluster.ids)
DimPlot(endo.pro, label = T)

endo.pro %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/04.endothelial.rds.gz", compress = "gz")
endo.pro <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/04.endothelial.rds.gz")

levels(endo.pro)
new.cluster.ids <- c("Endo.other","Endo.lymp","Endo.INST","Endo.ACKR1")
names(new.cluster.ids) <- levels(endo.pro)
endo.pro <- RenameIdents(endo.pro, new.cluster.ids)

#--- fibroblast

DimPlot(fb.sce, label = T)
fb.sce <- sce[,Idents(sce) %in% c(14)]

fb.sce <- NormalizeData(fb.sce)
fb.sce <- ScaleData(fb.sce, features = rownames(fb.sce), do.center = F, do.scale = F)
fb.sce <- FindVariableFeatures(fb.sce, selection.method = "vst", nfeatures = 2000)
fb.sce <- RunPCA(fb.sce, features = VariableFeatures(fb.sce))
fb.sce <- bbknn.batch(fb.sce, r = 1, n.pc = 30)

DimPlot(fb.sce, label = T)
FeaturePlot(fb.sce, features = "COL18A1")
markers <- FindAllMarkers(fb.sce, only.pos = T, test.use = "roc")

fb.sce %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/05.fibroblast.rds.gz", compress = "gz")
fb.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/05.fibroblast.rds.gz")
