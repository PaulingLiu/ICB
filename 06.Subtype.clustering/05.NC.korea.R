
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

sce <- readr::read_rds("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/10.external.data/03.NC/01.all.sce.rds.gz")
DimPlot(sce, label = T)
FeaturePlot(sce, features = c("COL1A1"))
mda <- readr::read_tsv("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/06.tumor.classification/04.NM.DL.tsv")

DimPlot(sce, label = T)
FeaturePlot(sce, features = "CD4")
new.cluster.ids <- c("Cancer cell","Myeloid","CD8 T","CD4 T","Prolif. T","CD4 T","NK","MAST","CD4 T","CD8 T","Myeloid","Myeloid","Myeloid","B cell","CAF","Cancer cell","Plasma B","Cancer cell","Plasma B","Plasma B","Plasma B", "Myeloid","pDC","Endo.")
names(new.cluster.ids) <- levels(sce)
all.sce.pro <- RenameIdents(sce, new.cluster.ids)
DimPlot(sce, label = T, pt.size = 0.001) + scale_color_d3(palette = "category20") #+ NoLegend() + NoAxes()

#--- T cells

DimPlot(sce, label = T) + NoLegend() + NoAxes()
VlnPlot(sce, features = "CXCL13")

#--- CD4 T cells

cd4.sce <- sce[,Idents(sce) %in% c(21,6,4,19,10,3,5,8)]
cd4.sce <- NormalizeData(cd4.sce)
cd4.sce <- ScaleData(cd4.sce, features = rownames(cd4.sce), do.center = F, do.scale = F)
cd4.sce <- FindVariableFeatures(cd4.sce, selection.method = "vst", nfeatures = 2000)
cd4.sce <- RunPCA(cd4.sce, features = VariableFeatures(cd4.sce))
cd4.sce <- bbknn.batch(cd4.sce, r = 3, n.pc = 30)

DimPlot(cd4.sce, label = T) #+ geom_abline(slope = -1.5, intercept = 8.5)
FeaturePlot(cd4.sce, features = "FOXP3")

new.cluster.ids <- levels(cd4.sce)
new.cluster.ids[new.cluster.ids %in% c(0,24)] <- "CD8.CXCL13"
new.cluster.ids[new.cluster.ids %in% c(14,11,23,3,4,8,5)] <- "CD8.other"
new.cluster.ids[new.cluster.ids %in% c(15)] <- "CD4.CXCL13"
new.cluster.ids[new.cluster.ids %in% c(6,9,1,20,22,12,19,13,2,18,21,10)] <- "CD4.other"
new.cluster.ids[new.cluster.ids %in% c(7,16,17)] <- "CD4.Treg"

mac.pro <- mac.sce
names(new.cluster.ids) <- levels(cd4.sce)
tcell.sce <- RenameIdents(cd4.sce, new.cluster.ids)
DimPlot(tcell.sce)
FeaturePlot(tcell.sce, features = "FOXP3")

tcell.sce %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/03.NC/02.tcell.sce.rds.gz", compress = "gz")
tcell.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/03.NC/02.tcell.sce.rds.gz")
levels(tcell.sce)
new.cluster.ids <- c("CD8-c01_CXCL13-","CD4-c04_other","CD8-c02_CXCL13+","CD4-c02_Treg","CD4-c01_CXCL13")
names(new.cluster.ids) <- levels(tcell.sce)
tcell.sce <- RenameIdents(tcell.sce, new.cluster.ids)
DimPlot(tcell.sce)


#--- macrophages

DimPlot(sce, label = T)
mye.sce <- sce[,Idents(sce) %in% c(24,16,20,13,9,27,18,23)]
DimPlot(mye.sce, label = T, pt.size = 0.001) + scale_color_d3(palette = "category20") #+ NoLegend() + NoAxes()
FeaturePlot(mye.sce, features = "CLEC10A")

mac.sce <- mye.sce
mac.sce <- NormalizeData(mac.sce)
mac.sce <- ScaleData(mac.sce, features = rownames(mac.sce), do.center = F, do.scale = F)
mac.sce <- FindVariableFeatures(mac.sce, selection.method = "vst", nfeatures = 2000)
mac.sce <- RunPCA(mac.sce, features = VariableFeatures(mac.sce))
mac.sce <- bbknn.batch(mac.sce, r = 2, n.pc = 30)

DimPlot(mac.sce, label = T)
FeaturePlot(mac.sce, features = "APOE")

new.cluster.ids <- levels(mac.sce)
new.cluster.ids[new.cluster.ids %in% c(4,10,16)] <- "FCN1"
new.cluster.ids[new.cluster.ids %in% c(8,11,9,5)] <- "SPP1"
new.cluster.ids[new.cluster.ids %in% c(1,6,7)] <- "FOLR2"
new.cluster.ids[new.cluster.ids %in% c(12,13,14)] <- "C1QA"
new.cluster.ids[new.cluster.ids %in% c(2,15,17,18)] <- "DC"
new.cluster.ids[new.cluster.ids %in% c(0,3)] <- "MARCO"

names(new.cluster.ids) <- levels(mac.sce)
mac.pro <- RenameIdents(mac.sce, new.cluster.ids)
DimPlot(mac.pro, label = T)




#DC FCN1 SPP1 FOLR2 C1QA MARCO

mac.pro %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/03.NC/03.myeloid.rds.gz", compress = "gz")
mac.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/03.NC/03.myeloid.rds.gz")
mac.sce <- mac.sce[,Idents(mac.sce) != "DC"]
DimPlot(mac.sce)

#--- DC analysis

dc.sce <- mac.pro[,Idents(mac.pro) == "DC"]
dc.sce <- FindVariableFeatures(dc.sce, selection.method = "vst", nfeatures = 2000)
dc.sce <- RunPCA(dc.sce, features = VariableFeatures(dc.sce))
dc.sce <- bbknn.batch(dc.sce, r = 1, n.pc = 30)

DimPlot(dc.sce, label = T)
FeaturePlot(dc.sce, features = "LILRA4")

new.cluster.ids <- levels(dc.sce)
new.cluster.ids[new.cluster.ids %in% c(3)] <- "LAMP3"
new.cluster.ids[new.cluster.ids %in% c(4)] <- "pDC"
new.cluster.ids[new.cluster.ids %in% c(7)] <- "DC1"
new.cluster.ids[new.cluster.ids %in% c(6)] <- "DC2.CD1A"
new.cluster.ids[new.cluster.ids %in% c(2,5,0,1,8)] <- "DC2.CD1C"

names(new.cluster.ids) <- levels(dc.sce)
dc.pro <- RenameIdents(dc.sce, new.cluster.ids)
DimPlot(dc.pro)

dc.pro %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/03.NC/03.dc.rds.gz", compress = "gz")
dc.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/03.NC/03.dc.rds.gz")


#--- endothelials

DimPlot(sce, label = T)
FeaturePlot(sce, features = "CLDN5")
endo.sce <- sce[,Idents(sce) %in% c(22)]

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
new.cluster.ids[new.cluster.ids %in% c(0,5)] <- "ACKR1"
new.cluster.ids[new.cluster.ids %in% c(3)] <- "lymp"
new.cluster.ids[new.cluster.ids %in% c(1)] <- "INSR"
new.cluster.ids[new.cluster.ids %in% c(2,4,6)] <- "other"


names(new.cluster.ids) <- levels(endo.sce)
endo.pro <- RenameIdents(endo.sce, new.cluster.ids)
DimPlot(endo.pro, label = T)

endo.pro %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/03.NC/04.endothelial.rds.gz", compress = "gz")
endo.pro <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/03.NC/04.endothelial.rds.gz")
DimPlot(endo.pro)
levels(endo.pro)
new.cluster.ids <- c("Endo.other","Endo.ACKR1","Endo.lymp","Endo.INSR")
names(new.cluster.ids) <- levels(endo.pro)
endo.pro <- RenameIdents(endo.pro, new.cluster.ids)
DimPlot(endo.pro)

#--- fibroblast

DimPlot(fb.sce, label = T)
fb.sce <- sce[,Idents(sce) %in% c(15)]

fb.sce <- NormalizeData(fb.sce)
fb.sce <- ScaleData(fb.sce, features = rownames(fb.sce), do.center = F, do.scale = F)
fb.sce <- FindVariableFeatures(fb.sce, selection.method = "vst", nfeatures = 2000)
fb.sce <- RunPCA(fb.sce, features = VariableFeatures(fb.sce))
fb.sce <- bbknn.batch(fb.sce, r = 1, n.pc = 30)

DimPlot(fb.sce, label = T)
FeaturePlot(fb.sce, features = "PDGFRA")
markers <- FindAllMarkers(fb.sce, only.pos = T, test.use = "roc")

fb.sce %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/03.NC/05.fibroblast.rds.gz", compress = "gz")

