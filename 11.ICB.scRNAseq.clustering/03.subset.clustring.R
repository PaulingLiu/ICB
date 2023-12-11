
#-- library

library(Seurat)
library(tidyverse)
library(ggplot2)
library(reticulate)
library(ggsci)
library(metacell)
library(ks)
library(RColorBrewer)

#-- load data

all.sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/15.ICB.sc/00.all.icb.pro.rds.gz")

DimPlot(all.sce, label = T)
FeaturePlot(all.sce, feature = "CXCL13")

#-- response information

rda <- tibble(
  cellid = all.sce$cellid,
  patient = all.sce$patient,
  num = all.sce$num
) %>%
  dplyr::mutate(sample = paste0(patient, ".", num))

rda.pro <- rda %>%
  dplyr::distinct(patient, num, sample) %>%
  as.data.frame() %>%
  dplyr::mutate(timepoint = ifelse(num == 0, "pre", "post")) %>%
  dplyr::mutate(response = ifelse(num == 2 | patient %in% c("P36","P37","P38","P39","P40","P41"), "NR", "R")) %>%
  dplyr::mutate(group = paste0(timepoint,".",response))

rda <- rda %>% dplyr::inner_join(rda.pro, by = "sample")

all.sce$response = rda$response
all.sce$timepoint = rda$timepoint
all.sce$group = rda$group

#--- T cell clustering

tcell.sce <- all.sce[,Idents(all.sce) %in% c("CD4 T","CD8 T","NK","Prolif. T")]
tcell.sce <- NormalizeData(tcell.sce)
tcell.sce <- ScaleData(tcell.sce, features = rownames(tcell.sce), do.center = F, do.scale = F)
tcell.sce <- FindVariableFeatures(tcell.sce, selection.method = "vst", nfeatures = 2000)
tcell.sce <- RunPCA(tcell.sce, features = VariableFeatures(tcell.sce))
tcell.sce <- tcell.sce[,tcell.sce$patient %in% names(table(tcell.sce$patient)[table(tcell.sce$patient) > 4])]
tcell.sce <- bbknn.batch(tcell.sce, r = 2, n.pc = 30)

DimPlot(tcell.sce, label = T) + geom_abline(slope = -2.3, intercept = -6)
FeaturePlot(tcell.sce, features = "CD8A")

new.cluster.ids <- levels(tcell.sce)
new.cluster.ids[new.cluster.ids %in% c(0,6,18)] <- "CD8.CXCL13"
new.cluster.ids[new.cluster.ids %in% c(23,3,2,5,17,9,15)] <- "CD8.other"
new.cluster.ids[new.cluster.ids %in% c(16,20)] <- "CD4.CXCL13"
new.cluster.ids[new.cluster.ids %in% c(4,7,22,10,24,14,19,11)] <- "CD4.other"
new.cluster.ids[new.cluster.ids %in% c(1,8)] <- "CD4.Treg"
new.cluster.ids[new.cluster.ids %in% c(21)] <- "Prolif"
new.cluster.ids[new.cluster.ids %in% c(13,12)] <- "NK"
tcell.sce$seurat_clusters <- Idents(tcell.sce)

names(new.cluster.ids) <- levels(tcell.sce)
tcell.sce <- RenameIdents(tcell.sce, new.cluster.ids)
DimPlot(tcell.sce, label = T) 

clusters <- as.character(Idents(tcell.sce))
clusters[clusters == "Prolif" & tcell.sce@reductions$umap@cell.embeddings[,1]*(-2.3)-6 < tcell.sce@reductions$umap@cell.embeddings[,2]] <- "CD8.Prolif"
clusters[clusters == "Prolif" & tcell.sce@reductions$umap@cell.embeddings[,1]*(-2.3)-6 > tcell.sce@reductions$umap@cell.embeddings[,2]] <- "CD4.Prolif"

Idents(tcell.sce) <- clusters
levels(tcell.sce)
new.cluster.ids <- c("NK","CD4-c02_Treg","CD4-c03_Prolif","CD4-c01_CXCL13",
                     "CD8-c01_CXCL13-","CD4-c04_other","CD8-c03_Prolif","CD8-c02_CXCL13+")
names(new.cluster.ids) <- levels(tcell.sce)
tcell.sce <- RenameIdents(tcell.sce, new.cluster.ids)

DimPlot(tcell.sce, label = T)

tcell.sce %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/15.ICB.sc/03.tcell.rds.gz", compress = "gz")
tcell.sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/15.ICB.sc/03.tcell.rds.gz")

table(tcell.sce$patient,tcell.sce$num)
tcell.sce$response = "R.pre"
tcell.sce$response[tcell.sce$patient %in% c("P36","P37","P38","P39","P40","P41") & tcell.sce$num == 1] = "NR.post"
tcell.sce$response[tcell.sce$response == "R.pre" & tcell.sce$num != 0] = "R.post"
tcell.sce$response[tcell.sce$num == 2] = "NR.post"

cd8.sce <- tcell.sce[,stringr::str_detect(Idents(tcell.sce),"CD8")]
cd4.sce <- tcell.sce[,stringr::str_detect(Idents(tcell.sce),"CD4")]
Idents(cd8.sce) <- cd8.sce$seurat_clusters

FeaturePlot(cd4.sce[,cd4.sce$patient == "P39" & cd4.sce$num == 1], features = "CXCL13")
DimPlot(cd8.sce)
FeaturePlot(cd8.sce[,cd8.sce$patient == "P33" & cd8.sce$num == 1], features = "CXCL13")

tibble(
  cluster = as.character(Idents(cd8.sce)),
  patient = cd8.sce$patient,
  num = cd8.sce$num,
  response = cd8.sce$response
) %>%
  dplyr::group_by(patient, num, response) %>%
  dplyr::count(cluster) %>%
  dplyr::mutate(freq = n/sum(n)) %>%
  dplyr::filter(cluster %in%  c("CD8-c02_CXCL13+")) %>%
  dplyr::mutate(f = sum(freq)) %>%
  ggplot(aes(response, f)) +
  geom_boxplot(outlier.size = -1) +
  geom_jitter(width = 0.1)

tibble(
  expr = unlist(cd8.sce@assays$RNA@scale.data["CXCL13",]),
  patient = cd8.sce$patient,
  num = cd8.sce$num,
  response = cd8.sce$response
) %>%
  dplyr::mutate(cluster = ifelse(expr > 2, "pos", "neg")) %>%
  dplyr::group_by(patient, num, response) %>%
  dplyr::count(cluster) %>%
  dplyr::mutate(freq = n/sum(n)) %>%
  dplyr::filter(cluster == "pos") %>%
  ggplot(aes(response, freq)) +
  geom_boxplot(outlier.size = -1) +
  geom_jitter(width = 0.1)

tibble(
  cluster = as.character(Idents(cd4.sce)),
  patient = cd4.sce$patient,
  num = cd4.sce$num,
  response = cd4.sce$response
) %>%
  dplyr::group_by(patient, num, response) %>%
  dplyr::count(cluster) %>%
  dplyr::mutate(freq = n/sum(n)) %>%
  dplyr::filter(cluster == "CD4-c01_CXCL13") %>%
  ggplot(aes(response, freq)) +
  geom_boxplot(outlier.size = -1) +
  geom_jitter(width = 0.1)

#--- macrophages

mye.sce <- all.sce[,Idents(all.sce) %in% c("Myeloid")]
DimPlot(mye.sce, label = T, pt.size = 0.001) + scale_color_d3(palette = "category20") #+ NoLegend() + NoAxes()


mac.sce <- mye.sce
mac.sce <- NormalizeData(mac.sce)
mac.sce <- ScaleData(mac.sce, features = rownames(mac.sce), do.center = F, do.scale = F)
mac.sce <- FindVariableFeatures(mac.pro, selection.method = "vst", nfeatures = 2000)
mac.sce <- RunPCA(mac.sce, features = VariableFeatures(mac.sce))
mac.sce <- bbknn.batch(mac.sce, r = 1.5, n.pc = 30)

DimPlot(mac.sce[,!(Idents(mac.sce) %in% c(7,10,13,14))], label = T)
FeaturePlot(mac.sce, features = c("APOE","FCN1","CLEC9A","CD1A","CD1C","LAMP3","CLEC10A"))
FeaturePlot(mac.sce, features = "C1QC")
FeaturePlot(mac.sce[,!(Idents(mac.sce) %in% c(7,10,13,14))], features = c("APOE","FCN1","SPP1","FOLR2","C1QA","MARCO"), ncol = 3)
markers.7 <- FindMarkers(mac.sce[,!(Idents(mac.sce) %in% c(7,10,13,14))], ident.1 = 3, only.pos = T, test.use = "roc")

#DC FCN1 SPP1 FOLR2 C1QA MARCO

mac.sce %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/15.ICB.sc/05.myeloid.rds.gz", compress = "gz")
mac.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/03.myeloid.rds.gz")
mac.pro <- mac.sce[,!(Idents(mac.sce) %in% c(14,13,7,10,6))]

mac.pro <- FindVariableFeatures(mac.pro, selection.method = "vst", nfeatures = 1000)
mac.pro <- RunPCA(mac.pro, features = VariableFeatures(mac.pro))
mac.pro <- bbknn.batch(mac.pro, r = 2, n.pc = 30)

DimPlot(mac.pro, label = T)
FeaturePlot(mac.pro, features = c("SPP1","FOLR2","MARCO","MT-CO1","FCN1"), ncol = 3)
markers.7 <- FindMarkers(mac.pro, ident.1 = 8, only.pos = T, test.use = "roc")

VlnPlot(mac.pro, features = "SPP1")

#mac.pro <- mac.pro[,Idents(mac.pro) != "9"]
levels(mac.pro)
new.cluster.ids <- as.character(levels(mac.pro))
new.cluster.ids[new.cluster.ids %in% c(13,15,6,12,14)] <- "SPP1+"
new.cluster.ids[new.cluster.ids %in% c(0,5)] <- "FOLR2+"
new.cluster.ids[new.cluster.ids %in% c(4)] <- "MARCO+"
new.cluster.ids[new.cluster.ids %in% c(1,8,16)] <- "FCN1+"
new.cluster.ids[new.cluster.ids %in% c(17,7,2,10,3)] <- "C1QA+"
new.cluster.ids[new.cluster.ids %in% c(11,9)] <- "MT-CO1"

names(new.cluster.ids) <- levels(mac.pro)
mac.pro.new <- RenameIdents(mac.pro, new.cluster.ids)

DimPlot(mac.pro.new, label = T)
FeaturePlot(mac.pro.new, features = "SPP1")

mac.pro.new %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/15.ICB.sc/05.macrophage.rds.gz", compress = "gz")


#--- DC clustering

dc.sce <- mac.sce[,Idents(mac.sce) %in% c(7,10,13,14,6)]
dc.pro <- all.sce[,Idents(all.sce) == "pDC" | all.sce$cellid %in% dc.sce$cellid]
dc.pro <- FindVariableFeatures(dc.pro, selection.method = "vst", nfeatures = 1000)
dc.pro <- RunPCA(dc.pro, features = VariableFeatures(dc.pro))
dc.pro <- bbknn.batch(dc.pro, r = 2, n.pc = 30)

DimPlot(dc.pro, label = T)

FeaturePlot(dc.pro, features = "LAMP3")

dc.pro <- dc.pro[,Idents(dc.pro) %in% c(0,11,5,17,13,4,12,15,6,16,8)]

levels(dc.pro)
dc.pro$idents <- Idents(dc.pro)
new.cluster.ids <- c("CD1A+ DC2","LAMP3+ DC","CLEC10A+ DC2","CD1A+ DC2","LAMP3+ DC","CLEC10A+ DC2","CLEC10A+ DC2","pDC","CLEC9A+ DC1","pDC","CLEC10A+ DC2")
names(new.cluster.ids) <- levels(dc.pro)
dc.pro <- RenameIdents(dc.pro, new.cluster.ids)
DimPlot(dc.pro, group.by = "idents", label = T)

table(dc.pro$idents, paste0(dc.pro$patient,"." ,dc.pro$num))
tibble(idents = Idents(dc.pro), patient = dc.pro$patient, group = dc.pro$group) %>%
  dplyr::count(patient, group, idents) %>%
  dplyr::group_by(patient, group) %>%
  dplyr::mutate(f = n/sum(n)) %>%
  dplyr::filter(idents == "LAMP3+ DC") %>%
  ggplot(aes(group, f)) +
  geom_boxplot(outlier.size = -1) +
  geom_jitter(width = 0.1)

dc.pro %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/15.ICB.sc/05.dc.rds.gz", compress = "gz")

#--- endothelials

DimPlot(sce, label = T)
endo.sce <- all.sce[,Idents(all.sce) == "Endo."]

endo.sce <- NormalizeData(endo.sce)
endo.sce <- ScaleData(endo.sce, features = rownames(endo.sce), do.center = F, do.scale = F)
endo.sce <- FindVariableFeatures(endo.sce, selection.method = "vst", nfeatures = 2000)
endo.sce <- RunPCA(endo.sce, features = VariableFeatures(endo.sce))
endo.sce <- endo.sce[,endo.sce$patient != "P41"]
endo.sce <- bbknn.batch(endo.sce, r = 1, n.pc = 30)

DimPlot(endo.sce, label = T)

markers <- FindMarkers(endo.sce, ident.1 = 4, only.pos = T)
FeaturePlot(endo.sce, features = "ACKR1")
genes <- c("ACKR1","CCL14","ADIRF","MT-CO1","MT-CO3","CCL21","PDPN","PROX1","IGFBP5","PTGDS","FLT1","INSR","KDR","PDXN","PGF","APLN","STMN1","MKI67","TYMS")

FeaturePlot(endo.sce, features = "INSR")

levels(endo.sce)
new.cluster.ids <- c("Endo_c01-lymphatic","Endo_c04-other","Endo_c04-other","Endo_c02-ACKR1","Endo_c03-INSR","Endo_c04-other","Endo_c04-other","Endo_c04-other","Endo_c02-ACKR1","Endo_c04-other","Endo_c04-other")
names(new.cluster.ids) <- levels(endo.sce)
endo.sce <- RenameIdents(endo.sce, new.cluster.ids)
DimPlot(endo.sce)

endo.sce %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/15.ICB.sc/06.endo.rds.gz", compress = "gz")
endo.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/15.ICB.sc/06.endo.rds.gz")

DimPlot(endo.sce, label = T)

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

