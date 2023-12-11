#--- library

library(Seurat)
library(tidyverse)
library(reticulate)
library(ggsci)

#--- load data 2

new.cluster.ids.sort <- c("Cancer cell","Myeloid","CD8 T","CD4 T","Prolif. T","CD4 T","NK","MAST","CD4 T","CD8 T","Myeloid","Myeloid","Myeloid","B cell","CAF","SMC","Cancer cell","Plasma B","Cancer cell","Plasma B","Plasma B","Plasma B", "Myeloid","pDC","Endo.","Neutrophil")

d2.sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/01.all.sce.rds.gz")
new.idents <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/02.add.Tcell.sub.annotation.rds.gz")
Idents(d2.sce) <- new.idents$cluster
DimPlot(d2.sce, label = T)

#--- load macrophage

mac.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/03.myeloid.rds.gz")
fb.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/05.fibroblast.rds.gz")

new.cluster.ids <- levels(d2.sce)
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

names(new.cluster.ids) <- levels(d2.sce)
d2.sce <- RenameIdents(d2.sce, new.cluster.ids)

Idents(d2.sce) <- ifelse(d2.sce$CellID %in% mac.sce$CellID[Idents(mac.sce) == 'pDC'], 'pDC', as.character(Idents(d2.sce)))
Idents(d2.sce) <- ifelse(d2.sce$CellID %in% fb.sce$CellID[Idents(fb.sce) %in% c("SMC", "perycite")], 'SMC', as.character(Idents(d2.sce)))

levels(d2.sce) <- unique(new.cluster.ids.sort)

DimPlot(d2.sce, label = F, pt.size = 0.001) + scale_color_d3(palette = "category20c") + NoLegend() + NoAxes()
as.character(Idents(d2.sce)) %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/11.major.cell.prop/12.d2.idents.rds.gz", compress = "gz")
d2.idents <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/11.major.cell.prop/12.d2.idents.rds.gz")
Idents(d2.sce) <- d2.idents

#--- load data 3

d3.sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/01.sce.rds.gz")

new.cluster.ids <- c("Cancer cell","Myeloid","CD8 T","CD4 T","Prolif. T","CD4 T","NK","MAST","CD4 T","CD8 T","Myeloid","Myeloid","Myeloid","B cell","CAF","Cancer cell","Plasma B","Cancer cell","Plasma B","Plasma B","Plasma B", "Myeloid","pDC","Endo.")
names(new.cluster.ids) <- levels(d3.sce)
d3.sce <- RenameIdents(d3.sce, new.cluster.ids)

fb.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/05.fibroblast.rds.gz")
DimPlot(fb.sce, label = T)

Idents(d3.sce) <- ifelse(d3.sce$Cell %in% fb.sce$Cell[Idents(fb.sce) %in% c(2,4,7)], 'SMC', as.character(Idents(d3.sce)))
levels(d3.sce) <- unique(new.cluster.ids.sort)
DimPlot(d3.sce, label = F, pt.size = 0.001) + scale_color_d3(palette = "category20c") + NoLegend() + NoAxes()
as.character(Idents(d3.sce)) %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/11.major.cell.prop/13.d3.idents.rds.gz", compress = "gz")
d3.idents <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/11.major.cell.prop/13.d3.idents.rds.gz")
Idents(d3.sce) <- d3.idents

#--- load data 4

d4.sce <- readr::read_rds("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/10.external.data/03.NC/01.all.sce.rds.gz")
DimPlot(d4.sce, label = T)

d4.sce <- d4.sce[,Idents(d4.sce) != "25"]

new.ids <- c(rep("CD4 T",4),rep("CD8 T",4),"NK",rep("Cancer cell",4), rep("Myeloid",7), rep("B cell",2), "pDC","MAST","Plasma B","Endo.","CAF")
ids <- c(3,5,8,10,4,6,19,21,11,1,0,7,12,20,16,13,27,9,23,18,2,26,24,14,17,22,15)

ida <- tibble(new.ids = new.ids, ids = as.character(ids))
ida <- tibble(ids = levels(d4.sce)) %>%
  dplyr::left_join(ida)

new.cluster.ids <- ida$new.ids
names(new.cluster.ids) <- levels(d4.sce)
d4.sce <- RenameIdents(d4.sce, new.cluster.ids)
fb.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/03.NC/05.fibroblast.rds.gz")

Idents(d4.sce) <- ifelse(colnames(d4.sce) %in% colnames(fb.sce)[Idents(fb.sce) %in% c("SMC","perycite")], 'SMC', as.character(Idents(d4.sce)))
levels(d4.sce) <- unique(new.cluster.ids.sort)
DimPlot(d4.sce, label = F, pt.size = 0.001) + scale_color_manual(values = unique(unname(d3_20c))[-5]) + NoLegend() + NoAxes()
as.character(Idents(d4.sce)) %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/11.major.cell.prop/14.d4.idents.rds.gz", compress = "gz")

d4.idents <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/11.major.cell.prop/14.d4.idents.rds.gz")
Idents(d4.sce) <- d4.idents


#--- load data 5

d5.sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/10.external.data/04.NC2/01.all.sce.rds.gz")

d5.sce <- d5.sce[,Idents(d5.sce) != "34"]
d5.sce <- d5.sce[,Idents(d5.sce) != "32"]

new.ids <- c(rep("Cancer cell",19),rep("T cell",2), rep("Myeloid",4), "B cell","MAST","Plasma B","Endo.",rep("CAF",2),"SMC","Neutrophil")
ids <- c(3,10,2,7,8,0,4,17,21,22,15,16,20,28,6,18,9,23,26,13,31,11,29,1,5,25,30,14,27,12,33,24,19)

ida <- tibble(new.ids = new.ids, ids = as.character(ids))
ida <- tibble(ids = levels(d5.sce)) %>%
  dplyr::left_join(ida)

new.cluster.ids <- ida$new.ids
names(new.cluster.ids) <- levels(d5.sce)
d5.sce <- RenameIdents(d5.sce, new.cluster.ids)

new.idents <- tibble(
  idents = as.character(Idents(d5.sce)),
  UMAP1 = unlist(d5.sce@reductions$umap@cell.embeddings[,1]),
  UMAP2 = unlist(d5.sce@reductions$umap@cell.embeddings[,2])
) %>%
  dplyr::mutate(idents = ifelse(idents %in% c("Myeloid","B cell") & (UMAP1*1.3-7.5 > UMAP2), "Myeloid", idents))

Idents(d5.sce) <- new.idents$idents

tcell.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/04.NC2/02.tcell.sce.rds.gz")

Idents(d5.sce) <- ifelse(d5.sce$cellnames %in% tcell.sce$cellnames[stringr::str_detect(Idents(tcell.sce),"CD8")], 'CD8 T', as.character(Idents(d5.sce)))
Idents(d5.sce) <- ifelse(d5.sce$cellnames %in% tcell.sce$cellnames[stringr::str_detect(Idents(tcell.sce),"CD4")], 'CD4 T', as.character(Idents(d5.sce)))
Idents(d5.sce) <- ifelse(d5.sce$cellnames %in% tcell.sce$cellnames[stringr::str_detect(Idents(tcell.sce),"NK")], 'NK', as.character(Idents(d5.sce)))
d5.sce <- d5.sce[,Idents(d5.sce) != "T cell"]

levels(d5.sce) <- unique(new.cluster.ids.sort)

DimPlot(d5.sce, label = F, pt.size = 0.001) + scale_color_manual(values = unique(unname(d3_20c))[-c(5,6,12)]) + NoLegend() + NoAxes()
as.character(Idents(d5.sce)) %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/11.major.cell.prop/15.d5.idents.rds.gz", compress = "gz")

d5.idents <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/11.major.cell.prop/15.d5.idents.rds.gz")
Idents(d5.sce) <- d5.idents

