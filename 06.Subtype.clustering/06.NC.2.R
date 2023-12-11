
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

sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/10.external.data/04.NC2/01.all.sce.rds.gz")
DimPlot(sce, label = T)
FeaturePlot(sce, features = c("MKI67"))
mda <- readr::read_tsv("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/06.tumor.classification/04.NM.DL.tsv")

#--- T cells

DimPlot(sce, label = T) + NoLegend() + NoAxes()
FeaturePlot(sce, features = "MMP2")

#--- CD4 T cells

cd4.sce <- sce[,Idents(sce) %in% c(13,31)]
cd4.sce <- NormalizeData(cd4.sce)
cd4.sce <- ScaleData(cd4.sce, features = rownames(cd4.sce), do.center = F, do.scale = F)
cd4.sce <- FindVariableFeatures(cd4.sce, selection.method = "vst", nfeatures = 2000)
cd4.sce <- RunPCA(cd4.sce, features = VariableFeatures(cd4.sce))
cd4.sce <- cd4.sce[,cd4.sce$patient %in% names(table(cd4.sce$patient)[table(cd4.sce$patient) > 4])]
cd4.sce <- bbknn.batch(cd4.sce, r = 2, n.pc = 30)

DimPlot(cd4.sce, label = T) #+ geom_abline(slope = -1.5, intercept = 8.5)
FeaturePlot(cd4.sce, features = "CD8A")

new.cluster.ids <- levels(cd4.sce)
new.cluster.ids[new.cluster.ids %in% c(0)] <- "CD8.CXCL13"
new.cluster.ids[new.cluster.ids %in% c(6,13,12,7,10)] <- "CD8.other"
new.cluster.ids[new.cluster.ids %in% c(16)] <- "CD4.CXCL13"
new.cluster.ids[new.cluster.ids %in% c(8,4,9,14,15,11)] <- "CD4.other"
new.cluster.ids[new.cluster.ids %in% c(1)] <- "CD4.Treg"
new.cluster.ids[new.cluster.ids %in% c(5)] <- "CD8.Prolif"
new.cluster.ids[new.cluster.ids %in% c(2)] <- "CD4.Prolif"
new.cluster.ids[new.cluster.ids %in% c(3)] <- "NK"


names(new.cluster.ids) <- levels(cd4.sce)
tcell.sce <- RenameIdents(cd4.sce, new.cluster.ids)
DimPlot(tcell.sce)

tcell.sce %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/04.NC2/02.tcell.sce.rds.gz", compress = "gz")
tcell.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/04.NC2/02.tcell.sce.rds.gz")

levels(tcell.sce)
new.cluster.ids <- c("CD8-c01_CXCL13-","CD8-c03_Prolif","CD4-c03_Prolif","CD4-c01_CXCL13","CD8-c02_CXCL13+","CD4-c02_Treg","NK","CD4-c04_other")
names(new.cluster.ids) <- levels(tcell.sce)
tcell.sce <- RenameIdents(tcell.sce, new.cluster.ids)
DimPlot(tcell.sce)

#--- macrophages

DimPlot(sce, label = T)
new.idents <- tibble(
  idents = as.character(Idents(sce)),
  UMAP1 = unlist(sce@reductions$umap@cell.embeddings[,1]),
  UMAP2 = unlist(sce@reductions$umap@cell.embeddings[,2])
) %>%
  dplyr::mutate(idents = ifelse(idents == 25 & (UMAP1*1.3-7.5 > UMAP2), 38, idents))

Idents(sce) <- new.idents$idents
FeaturePlot(sce, features = "LYZ")
mye.sce <- sce[,Idents(sce) %in% c(5,1,11,29,38)]
DimPlot(mye.sce, label = T, pt.size = 0.001) + scale_color_d3(palette = "category20") #+ NoLegend() + NoAxes()
FeaturePlot(mye.sce, features = "MS4A1") + geom_abline(slope = 1.3, intercept = -7.5)


mac.sce <- mye.sce
mac.sce <- NormalizeData(mac.sce)
mac.sce <- ScaleData(mac.sce, features = rownames(mac.sce), do.center = F, do.scale = F)
mac.sce <- FindVariableFeatures(mac.sce, selection.method = "vst", nfeatures = 2000)
mac.sce <- RunPCA(mac.sce, features = VariableFeatures(mac.sce))
mac.sce <- mac.sce[,mac.sce$patient %in% names(table(mac.sce$patient)[table(mac.sce$patient) > 4])]
mac.sce <- bbknn.batch(mac.sce, r = 1.5, n.pc = 30)

DimPlot(mac.sce, label = T)
markers <- FindMarkers(mac.sce, ident.1 = 8, only.pos = T)
FeaturePlot(mac.sce, features = c("SPP1"))
mac.sce <- mac.sce[,Idents(mac.sce) != 11]

new.cluster.ids <- levels(mac.sce)
new.cluster.ids[new.cluster.ids %in% c(1)] <- "FCN1"
new.cluster.ids[new.cluster.ids %in% c(2,0,12,5,3,9)] <- "SPP1"
new.cluster.ids[new.cluster.ids %in% c(7)] <- "FOLR2"
new.cluster.ids[new.cluster.ids %in% c(10)] <- "C1QA"
new.cluster.ids[new.cluster.ids %in% c(8,6)] <- "DC"
new.cluster.ids[new.cluster.ids %in% c(4)] <- "MARCO"

names(new.cluster.ids) <- levels(mac.sce)
mac.pro <- RenameIdents(mac.sce, new.cluster.ids)
DimPlot(mac.pro, label = T)




#DC FCN1 SPP1 FOLR2 C1QA MARCO

mac.pro %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/04.NC2/03.myeloid.rds.gz", compress = "gz")
mac.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/04.NC2/03.myeloid.rds.gz")
DimPlot(mac.sce)
mac.sce <- mac.sce[,Idents(mac.sce) != "DC"]

#--- DC analysis

dc.sce <- mac.pro[,Idents(mac.pro) == "DC"]
dc.sce <- FindVariableFeatures(dc.sce, selection.method = "vst", nfeatures = 2000)
dc.sce <- RunPCA(dc.sce, features = VariableFeatures(dc.sce))
dc.sce <- dc.sce[,dc.sce$patient %in% names(table(dc.sce$patient)[table(dc.sce$patient) > 4])]
dc.sce <- bbknn.batch(dc.sce, r = 1.5, n.pc = 30)

DimPlot(dc.sce, label = T)
FeaturePlot(dc.sce, features = "CD1A")

dc.sce <- dc.sce[,Idents(dc.sce) %in% c(1,5,9,11,6,10)]

new.cluster.ids <- levels(dc.sce)
new.cluster.ids[new.cluster.ids %in% c(6)] <- "LAMP3"
new.cluster.ids[new.cluster.ids %in% c(10)] <- "pDC"
new.cluster.ids[new.cluster.ids %in% c(11)] <- "DC1"
new.cluster.ids[new.cluster.ids %in% c(5)] <- "DC2.CD1A"
new.cluster.ids[new.cluster.ids %in% c(1,9)] <- "DC2.CD1C"

names(new.cluster.ids) <- levels(dc.sce)
dc.pro <- RenameIdents(dc.sce, new.cluster.ids)
DimPlot(dc.pro)

dc.pro <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/04.NC2/03.dc.rds.gz")
levels(dc.pro)
new.cluster.ids <- c("DC2.CD1C","LAMP3.DC","pDC","DC1","DC2.CD1A")
names(new.cluster.ids) <- levels(dc.pro)
dc.pro <- RenameIdents(dc.pro, new.cluster.ids)

#--- endothelials

DimPlot(sce, label = T)
FeaturePlot(sce, features = "CLDN5")
endo.sce <- sce[,Idents(sce) %in% c(27)]

endo.sce <- NormalizeData(endo.sce)
endo.sce <- ScaleData(endo.sce, features = rownames(endo.sce), do.center = F, do.scale = F)
endo.sce <- FindVariableFeatures(endo.sce, selection.method = "vst", nfeatures = 2000)
endo.sce <- RunPCA(endo.sce, features = VariableFeatures(endo.sce))
endo.sce <- endo.sce[,endo.sce$patient %in% names(table(endo.sce$patient)[table(endo.sce$patient) > 4])]
endo.sce <- bbknn.batch(endo.sce, r = 1, n.pc = 30)

DimPlot(endo.sce, label = T)
FeaturePlot(endo.sce, features = "INSR")

markers <- FindMarkers(endo.sce, ident.1 = 4, only.pos = T)
FeaturePlot(endo.sce, features = "INSR")
genes <- c("ACKR1","CCL14","ADIRF","MT-CO1","MT-CO3","CCL21","PDPN","PROX1","IGFBP5","PTGDS","FLT1","INSR","KDR","PDXN","PGF","APLN","STMN1","MKI67","TYMS")

new.cluster.ids <- levels(endo.sce)
new.cluster.ids[new.cluster.ids %in% c(0)] <- "ACKR1"
#new.cluster.ids[new.cluster.ids %in% c(3)] <- "lymp"
new.cluster.ids[new.cluster.ids %in% c(1,3,2)] <- "INSR"
new.cluster.ids[new.cluster.ids %in% c(4,5,6,7)] <- "other"


names(new.cluster.ids) <- levels(endo.sce)
endo.pro <- RenameIdents(endo.sce, new.cluster.ids)
DimPlot(endo.pro, label = T)

endo.pro %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/04.NC2/04.endothelial.rds.gz", compress = "gz")
endo.pro <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/04.NC2/04.endothelial.rds.gz")
levels(endo.pro)
new.cluster.ids <- c("Endo.INSR","Endo.other","Endo.ACKR1")
names(new.cluster.ids) <- levels(endo.pro)
endo.pro <- RenameIdents(endo.pro, new.cluster.ids)

#--- fibroblast

DimPlot(fb.sce, label = T)
fb.sce <- sce[,Idents(sce) %in% c(33,24,12,34)]

fb.sce <- NormalizeData(fb.sce)
fb.sce <- ScaleData(fb.sce, features = rownames(fb.sce), do.center = F, do.scale = F)
fb.sce <- FindVariableFeatures(fb.sce, selection.method = "vst", nfeatures = 2000)
fb.sce <- RunPCA(fb.sce, features = VariableFeatures(fb.sce))

fb.counts <- table(fb.sce$patient)
fb.sce <- fb.sce[,fb.sce$patient %in% names(fb.counts[fb.counts > 9])]
fb.sce <- bbknn.batch(fb.sce, r = 1, n.pc = 30)

DimPlot(fb.sce, label = T)
FeaturePlot(fb.sce, features = "PDGFRA")
markers <- FindAllMarkers(fb.sce, only.pos = T, test.use = "roc")

fb.sce %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/04.NC2/05.fibroblast.rds.gz", compress = "gz")

