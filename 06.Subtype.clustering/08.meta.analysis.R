
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(ggsci)
library(data.table)
library(reticulate)
use_python("/home/pauling/anaconda3/bin/python")

#--- load data 1

all.sce <- readr::read_rds("projects/04_lung_immune_therapy/02_clustering/revision/00.re.clustering/01.all.sce.rds.gz")
new.cluster.ids <- c("Cancer cell","Myeloid","Plasma B","NK", "B cell","SMC","CD4 T","Prolif. T", "Cancer cell", "MAST", "CD4 T", "CAF", "CD8 T", "CD8 T", "CD4 T", "Cancer cell", "Cancer cell", "Endo.", "pDC", "Cancer cell")
names(new.cluster.ids) <- levels(all.sce)
all.sce.pro <- RenameIdents(all.sce, new.cluster.ids)

all.sce <- all.sce.pro[,!(all.sce.pro$patient %in% c("P1","P10","P19","P13","P29","P30","P33","P35","P36","P37","P38"))]
endo.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/01.clustering/05.endo.clustering.TN.sce.rds.gz")
mac.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/01.clustering/06.mac.tn.rds.gz")
dc.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/01.clustering/01.dc.sce.rds.gz")
dc.sce <- dc.sce[,!(dc.sce$patient %in% c("P1","P10","P19","P13","P29","P30","P33","P35","P36","P37","P38"))]
cd8.sce <- readr::read_rds("projects/04_lung_immune_therapy/02_clustering/revision/00.re.clustering/06.CD8.rds.gz")
cd8.sce  <- cd8.sce[,!(cd8.sce$patient %in% c("P1","P10","P19","P13","P29","P30","P33","P35","P36","P37","P38"))]

levels(endo.sce)
endo.names <- levels(endo.sce)
new.cluster.ids <- c(endo.names[1],"Endo_c04-other",endo.names[3],"Endo_c04-other","Endo_c03-INSR","Endo_c04-other","Endo_c04-other")
names(new.cluster.ids) <- levels(endo.sce)
endo.sce <- RenameIdents(endo.sce, new.cluster.ids)
DimPlot(endo.sce)

#--- load data 2 CR Lambrechts

sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/01.sce.rds.gz")
mda <- readr::read_tsv("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/06.tumor.classification/03.CR.DL.tsv")

new.cluster.ids <- c("Cancer cell","Myeloid","CD8 T","CD4 T","Prolif. T","CD4 T","NK","MAST","CD4 T","CD8 T","Myeloid","Myeloid","Myeloid","B cell","CAF","Cancer cell","Plasma B","Cancer cell","Plasma B","Plasma B","Plasma B", "Myeloid","pDC","Endo.")
names(new.cluster.ids) <- levels(sce)
all.sce.pro <- RenameIdents(sce, new.cluster.ids)
new.idents <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/02.add.Tcell.sub.annotation.rds.gz")
Idents(sce) <- new.idents$cluster
DimPlot(tcell.sce, label = T)
FeaturePlot(sce, features = "CD3D")

tcell.sce <- sce[,Idents(sce) %in% c(15,2,14,4,7,0,24)]
tcell.idents <- tibble(
  UMAP1 = unlist(tcell.sce@reductions$umap@cell.embeddings[,1]),
  UMAP2 = unlist(tcell.sce@reductions$umap@cell.embeddings[,2]),
  idents = as.character(Idents(tcell.sce))
) %>%
  dplyr::mutate(idents = ifelse(UMAP2>-5 & idents == 15, "65", idents))

Idents(tcell.sce) <- tcell.idents$idents
levels(tcell.sce)
new.cluster.ids <- c("CD8-c01_CXCL13-","CD4-c04_other","CD4-c03_Prolif","CD4-c04_other","CD4-c01_CXCL13","CD8-c02_CXCL13+","CD8-c03_Prolif","CD4-c02_Treg")
names(new.cluster.ids) <- levels(tcell.sce)
tcell.sce <- RenameIdents(tcell.sce, new.cluster.ids)
DimPlot(tcell.sce, label = T)
tcell.sce %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/02.tcell.sce.rds.gz", compress = "gz")

#--- load data 3 NM Lambrechts

sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/01.all.sce.rds.gz")
FeaturePlot(sce, features = c("CXCL13"))

DimPlot(all.sce.pro, label = T)
tcell.sce <- all.sce.pro[,Idents(all.sce.pro) %in% c(10,11,3,0,13,8,5,17,31,30)]
new.cluster.ids <- levels(tcell.sce)
new.cluster.ids[new.cluster.ids %in% c(30)] <- "CD4-c02_Treg"
new.cluster.ids[new.cluster.ids %in% c(3,10,11)] <- "CD4-c04_other"
new.cluster.ids[new.cluster.ids %in% c(0)] <- "CD8-c01_CXCL13-"
new.cluster.ids[new.cluster.ids %in% c(5,8)] <- "CD8-c02_CXCL13+"
new.cluster.ids[new.cluster.ids %in% c(13)] <- "CD4-c01_CXCL13"
new.cluster.ids[new.cluster.ids %in% c(17)] <- "CD8-c03_Prolif"
new.cluster.ids[new.cluster.ids %in% c(31)] <- "CD4-c03_Prolif"

names(new.cluster.ids) <- levels(tcell.sce)
tcell.sce <- RenameIdents(tcell.sce, new.cluster.ids)
DimPlot(tcell.sce, label = T)

tcell.sce %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/02.Tcell.sub.rds.gz", compress = "gz")

#--- functions to get cell clusters

get.major <- function(.x){
  tibble(
    patient = .x$patient,
    cellid = colnames(.x),
    ct1 = as.character(Idents(.x))
  )
}
get.sub <- function(.meta, .x){
  tmp.meta <- get.major(.x) %>% 
    dplyr::select(-patient) %>%
    dplyr::rename(ct3 = ct1)
  
  .meta %>%
    dplyr::left_join(tmp.meta, by = "cellid") %>%
    dplyr::mutate(ct2 = ifelse(!is.na(ct3), ct3, NA)) %>%
    dplyr::select(-ct3)
}
get.other <- function(.meta, .x){
  tmp.meta <- get.major(.x) %>% 
    dplyr::select(-patient) %>%
    dplyr::rename(ct3 = ct1)
  
  .meta %>%
    dplyr::left_join(tmp.meta, by = "cellid") %>%
    dplyr::mutate(ct2 = ifelse(!is.na(ct3), ct3, ct2)) %>%
    dplyr::select(-ct3)
}

#--- data 1

d1.meta <- get.major(all.sce)
d1.meta <- get.sub(d1.meta, endo.sce)
d1.meta <- get.other(d1.meta, mac.sce)
d1.meta <- get.other(d1.meta, dc.sce)
d1.meta <- get.other(d1.meta, cd8.filt)
d1.meta <- get.other(d1.meta, cd4.filt)

d1.meta <- d1.meta %>% dplyr::mutate(ct2 = ifelse(is.na(ct2), "ND", ct2))
d1.meta %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/01.our.data.TN.rds.gz", compress = "gz")
d1.meta <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/01.our.data.TN.rds.gz")

#--- data 2

d2.meta <- get.major(all.sce.pro)
d2.meta <- get.sub(d2.meta, tcell.sce)
d2.meta <- get.other(d2.meta, mac.pro)
d2.meta <- get.other(d2.meta, dc.sce)
d2.meta <- get.other(d2.meta, endo.sce)

d2.meta <- d2.meta %>% dplyr::mutate(ct2 = ifelse(is.na(ct2), "ND", ct2))
d2.meta %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/02.CR.DL.rds.gz", compress = "gz")


#--- data 3

d3.meta <- get.major(all.sce.pro)
d3.meta <- get.sub(d3.meta, tcell.sce)
d3.meta <- get.other(d3.meta, mac.sce)
d3.meta <- get.other(d3.meta, endo.pro)

d3.meta <- d3.meta %>% dplyr::mutate(ct2 = ifelse(is.na(ct2), "ND", ct2))
d3.meta %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/05.NM.rds.gz", compress = "gz")


#--- data 4

d3.meta <- get.major(all.sce.pro)
d3.meta <- get.sub(d3.meta, tcell.sce)
d3.meta <- get.other(d3.meta, mac.sce)
d3.meta <- get.other(d3.meta, dc.sce)
d3.meta <- get.other(d3.meta, endo.pro)

d3.meta <- d3.meta %>% dplyr::mutate(ct2 = ifelse(is.na(ct2), "ND", ct2))
d3.meta %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/03.NC.rds.gz", compress = "gz")

#--- data 5

d4.meta <- get.major(all.sce.pro)
d4.meta <- get.sub(d4.meta, tcell.sce)
d4.meta <- d4.meta %>%
  dplyr::mutate(ct1 = ifelse(ct2 == "NK" & !is.na(ct2), "NK", ct1))
d4.meta <- get.other(d4.meta, mac.sce)
d4.meta <- get.other(d4.meta, dc.pro)
d4.meta <- get.other(d4.meta, endo.pro)

d4.meta <- d4.meta %>% dplyr::mutate(ct2 = ifelse(is.na(ct2), "ND", ct2))
d4.meta %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/04.NC2.rds.gz")
d4.meta <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/04.NC2.rds.gz")
