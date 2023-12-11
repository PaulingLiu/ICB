
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

sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/01.sce.rds.gz")
mda <- readr::read_tsv("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/06.tumor.classification/03.CR.DL.tsv")

new.cluster.ids <- c("Cancer cell","Myeloid","CD8 T","CD4 T","Prolif. T","CD4 T","NK","MAST","CD4 T","CD8 T","Myeloid","Myeloid","Myeloid","B cell","CAF","Cancer cell","Plasma B","Cancer cell","Plasma B","Plasma B","Plasma B", "Myeloid","pDC","Endo.")
names(new.cluster.ids) <- levels(sce)
all.sce.pro <- RenameIdents(sce, new.cluster.ids)
DimPlot(all.sce.pro, label = T, pt.size = 0.001) + scale_color_d3(palette = "category20c") #+ NoLegend() + NoAxes()
DimPlot(all.sce.pro, label = T, pt.size = 0.001) + scale_color_manual(values = comb.d3[c(1,4,2,3,5:20)]) #+ NoLegend() + NoAxes()


#--- T cells

DimPlot(sce, label = T) + scale_color_d3(palette = "category20c") + NoLegend() + NoAxes()
VlnPlot(sce, features = "CXCL13")

#--- CD4 T cells

cd4.sce <- sce[,Idents(sce) %in% c(7,0,14)]
cd4.sce <- NormalizeData(cd4.sce)
cd4.sce <- ScaleData(cd4.sce, features = rownames(cd4.sce), do.center = F, do.scale = F)
cd4.sce <- FindVariableFeatures(cd4.sce, selection.method = "vst", nfeatures = 2000)
cd4.sce <- RunPCA(cd4.sce, features = VariableFeatures(cd4.sce))
cd4.sce <- bbknn.batch(cd4.sce, r = 1, n.pc = 15)

DimPlot(cd4.sce, label = T)
treg.mda <- tibble(
  cellid = colnames(cd4.sce[,Idents(cd4.sce) == 0]),
  cluster = "24"
)

new.idents <- tibble(
  cellid = colnames(sce),
  cluster = as.character(Idents(sce))
) %>%
  dplyr::mutate(cluster = ifelse(cellid %in% treg.mda$cellid, "24", cluster))

sce.pro <- sce
Idents(sce.pro) <- new.idents$cluster

new.idents <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/02.add.Tcell.sub.annotation.rds.gz")

#--- macrophages

mye.sce <- sce[,Idents(sce) %in% c(13,3,5,12,10,23)]
DimPlot(mye.sce, label = T, pt.size = 0.001) + scale_color_d3(palette = "category20") #+ NoLegend() + NoAxes()
FeaturePlot(mye.sce, features = "C1QA")

mac.sce <- mye.sce
mac.sce <- NormalizeData(mac.sce)
mac.sce <- ScaleData(mac.sce, features = rownames(mac.sce), do.center = F, do.scale = F)
mac.sce <- FindVariableFeatures(mac.sce, selection.method = "vst", nfeatures = 2000)
mac.sce <- RunPCA(mac.sce, features = VariableFeatures(mac.sce))
mac.sce <- bbknn.batch(mac.sce, r = 2, n.pc = 30)

DimPlot(mac.pro, label = T)
FeaturePlot(mac.sce, features = "CXCL9")
markers.7 <- FindMarkers(mac.sce[,Idents(mac.sce) %in% c(10,11,12)], ident.1 = 12, only.pos = T, test.use = "roc")

#DC FCN1 SPP1 FOLR2 C1QA MARCO

mac.sce %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/03.myeloid.rds.gz", compress = "gz")
mac.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/03.myeloid.rds.gz")
mac.pro <- mac.sce[,!(Idents(mac.sce) %in% c(18,3,19))]


levels(mac.pro)
new.cluster.ids <- as.character(levels(mac.pro))
new.cluster.ids[new.cluster.ids %in% c(1,6,4,20,5,14,9)] <- "SPP1+"
new.cluster.ids[new.cluster.ids %in% c(12,10,11)] <- "FOLR2+"
new.cluster.ids[new.cluster.ids %in% c(0)] <- "MARCO+"
new.cluster.ids[new.cluster.ids %in% c(17,7,13)] <- "FCN1+"
new.cluster.ids[new.cluster.ids %in% c(2,8,15)] <- "C1QA+"
new.cluster.ids[new.cluster.ids %in% c(16)] <- "RPL+"

names(new.cluster.ids) <- levels(mac.pro)
mac.pro <- RenameIdents(mac.pro, new.cluster.ids)

DimPlot(mac.pro, label = T)
FeaturePlot(mac.pro, features = "SPP1")

#--- DC clustering

dc.sce <- mac.sce[,Idents(mac.sce) %in% c(18,3,19)]
FeaturePlot(dc.sce, features = "LAMP3")
DimPlot(dc.sce, label = T) + 
  geom_abline(slope = -2, intercept = -12) +
  geom_hline(yintercept = 1.25)

dc.idents <- tibble(
  UMAP1 = unlist(dc.sce@reductions$umap@cell.embeddings[,1]),
  UMAP2 = unlist(dc.sce@reductions$umap@cell.embeddings[,2]),
  idents = as.character(Idents(dc.sce))
) %>%
  dplyr::mutate(idents = ifelse((-2*UMAP1-12 > UMAP2) & idents == 18, "39", idents)) %>%
  dplyr::mutate(idents = ifelse(UMAP2 > 1.25 & idents == 3, "30", idents))

Idents(dc.sce) <- dc.idents$idents
levels(dc.sce)
new.cluster.ids <- c("LAMP3+ DC","CLEC10A+ DC2","CD1A+ DC2","pDC","CLEC9A+ DC1")
names(new.cluster.ids) <- levels(dc.sce)
dc.sce <- RenameIdents(dc.sce, new.cluster.ids)
DimPlot(dc.sce)

#--- endothelials

DimPlot(sce, label = T)
endo.sce <- sce[,Idents(sce) %in% c(17)]

endo.sce <- NormalizeData(endo.sce)
endo.sce <- ScaleData(endo.sce, features = rownames(endo.sce), do.center = F, do.scale = F)
endo.sce <- FindVariableFeatures(endo.sce, selection.method = "vst", nfeatures = 2000)
endo.sce <- RunPCA(endo.sce, features = VariableFeatures(endo.sce))
endo.sce <- bbknn.batch(endo.sce, r = 1, n.pc = 30)

DimPlot(endo.sce, label = T)

markers <- FindMarkers(endo.sce, ident.1 = 4, only.pos = T)
FeaturePlot(endo.sce, features = "APLN")
genes <- c("ACKR1","CCL14","ADIRF","MT-CO1","MT-CO3","CCL21","PDPN","PROX1","IGFBP5","PTGDS","FLT1","INSR","KDR","PDXN","PGF","APLN","STMN1","MKI67","TYMS")

endo.sce %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/04.endothelial.rds.gz", compress = "gz")
endo.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/04.endothelial.rds.gz")

FeaturePlot(endo.sce, features = "INSR")

levels(endo.sce)
new.cluster.ids <- c("Endo_c04-other","Endo_c04-other","Endo_c04-other","Endo_c01-lymphatic","Endo_c03-INSR","Endo_c02-ACKR1","Endo_c04-other")
names(new.cluster.ids) <- levels(endo.sce)
endo.sce <- RenameIdents(endo.sce, new.cluster.ids)
DimPlot(endo.sce)

#--- fibroblast

DimPlot(fb.sce, label = T)
fb.sce <- sce[,Idents(sce) %in% c(11)]

fb.sce <- NormalizeData(fb.sce)
fb.sce <- ScaleData(fb.sce, features = rownames(fb.sce), do.center = F, do.scale = F)
fb.sce <- FindVariableFeatures(fb.sce, selection.method = "vst", nfeatures = 2000)
fb.sce <- RunPCA(fb.sce, features = VariableFeatures(fb.sce))
fb.sce <- bbknn.batch(fb.sce, r = 1, n.pc = 30)

DimPlot(fb.sce, label = T)
FeaturePlot(fb.sce, features = "PDGFRA")
markers <- FindAllMarkers(fb.sce, only.pos = T, test.use = "roc")

fb.sce %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/05.fibroblast.rds.gz", compress = "gz")
fb.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/05.fibroblast.rds.gz")

DimPlot(fb.sce)
FeaturePlot(fb.sce, features = "MMP1")

