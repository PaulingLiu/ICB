
#--- major cell type plot (our data treatment-naive tumors)

#--- library

library(Seurat)
library(tidyverse)
library(reticulate)
library(ggsci)

#--- load data

all.sce <- readr::read_rds("projects/04_lung_immune_therapy/02_clustering/revision/00.re.clustering/01.all.sce.rds.gz")
filt.patients <- c(c("P1","P10","P19","P29","P30","P13","P33","P35"),c("P36","P37","P38"),"P21","P25") # ICB-treated or liver mets

all.sce.filt <- all.sce[,!(all.sce$patient %in% filt.patients)]
#sum(table(all.sce.filt$patient))

all.sce.filt$num <- NA
all.sce.filt$treatment <- NA

all.sce.filt$cellid <- stringr::str_replace_all(all.sce.filt$cellid, "tr.1", "ut")
colnames(all.sce.filt) <- all.sce.filt$cellid
aa = Seurat::RenameCells(all.sce.filt, all.sce.filt$cellid)

aa %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/discovery_cohort_1.rds.gz", compress = "gz")

#DimPlot(all.sce, label = T)
levels(all.sce)

new.cluster.ids <- c("Cancer cell","Myeloid","Plasma B","NK", "B cell","pericyte & SMC","CD4 T","Prolif. T", "Cancer cell", "MAST", "CD4 T", "fibroblast", "CD8 T", "CD8 T", "CD4 T", "Cancer cell", "Cancer cell", "Endo.", "pDC", "Cancer cell")
names(new.cluster.ids) <- levels(all.sce.filt)
all.sce.filt <- RenameIdents(all.sce.filt, new.cluster.ids)

#--- plot 

cols <- c("#C1CDCD","#EEA9B8","#C6B596", "#EED2EE", "#5497C0", "#00729A", "#62A49B", "#3CB371", "#9BCD9B", "#E66764","#7B68EE","#EEDFCC","#D31771")
cols <- c("#C1CDCD","#EEA9B8","#00729A", "#EEB4B4", "#5497C0", "#62A49B", "#62A49B", "#3CB371", "#9BCD9B", "#E66764","#62A49B","#EEDFCC","#D31771")

DimPlot(all.sce.filt, label = T, pt.size = 0.3) + scale_color_manual(values = cols) + NoLegend() + NoAxes()
DimPlot(all.sce.filt, label = F, pt.size = 0.1, group.by = "patient") + scale_color_manual(values = comb.d3[-c(1:3)])

all.sce.filt@reductions$umap@cell.embeddings %>%
  as.tibble() %>%
  dplyr::mutate(patient = all.sce.filt$patient) %>%
  dplyr::sample_n(nrow(.), replace = F) %>%
  ggplot(aes(UMAP_1,UMAP_2)) +
  geom_point(aes(color = patient), size = 0.1) +
  theme_classic() +
  scale_color_manual(values = comb.d3)


tibble(
  expr1 = all.sce.filt@assays$RNA@scale.data["LILRA4",],
  label = Idents(all.sce.filt) 
) %>%
  ggplot(aes(1, expr1)) +
  geom_jitter(aes(color = label), size = 0.2) +
  theme_classic() +
  scale_color_manual(values = cols)

#--- DC plot

dc.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/01.clustering/01.dc.sce.rds.gz")
filt.patients <- c(c("P1","P10","P19","P29","P30","P13","P33","P35"),c("P36","P37","P38"),"P21","P25") # ICB-treated or liver mets
dc.sce.filt <- dc.sce[,!(dc.sce$patient %in% filt.patients)]
cols <- c("#C1CDCD","#EEA9B8","#00729A", "#EEB4B4", "#5497C0", "#62A49B", "#62A49B", "#3CB371", "#9BCD9B", "#E66764","#62A49B","#EEDFCC","#D31771")

DimPlot(dc.sce.filt, label = T, pt.size = 1.5) + scale_color_manual(values = cols[c(2,3,8,5,10)]) + NoAxes() + NoLegend()

#--- endothelial cell plot

endo <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/01.clustering/05.endo.clustering.TN.sce.rds.gz")
filt.patients <- c(c("P1","P10","P19","P29","P30","P13","P33","P35"),c("P36","P37","P38"),"P21","P25") # ICB-treated or liver mets
endo.filt <- endo[,!(endo$patient %in% filt.patients)]
DimPlot(d1.endo, label = F, pt.size = 1) + scale_color_manual(values = cols[c(5,2,3,1,4,7,8)]) + NoAxes() + NoLegend()

#--- CD8 T cell plot

cd8 <- readr::read_rds("projects/04_lung_immune_therapy/02_clustering/revision/00.re.clustering/06.CD8.rds.gz")
filt.patients <- c(c("P1","P10","P19","P29","P30","P13","P33","P35"),c("P36","P37","P38"),"P21","P25") # ICB-treated or liver mets
cd8.filt <- cd8[,!(cd8$patient %in% filt.patients)]
levels(cd8.filt)

new.cluster.ids <- c("CD8_c01-CXCL13-","CD8_c01-CXCL13-","CD8_c01-CXCL13-","CD8_c02-CXCL13+","CD8_c01-CXCL13-","CD8_c01-CXCL13-","CD8_c03-Prolif")
names(new.cluster.ids) <- levels(cd8.filt)
cd8.filt <- RenameIdents(cd8.filt, new.cluster.ids)

DimPlot(cd8.filt, label = F, pt.size = 0.8) #+ scale_color_manual(values = cols[c(5,2,3,1,4,7,8)]) + NoAxes() + NoLegend()

#--- CD4 T cell plot

cd4 <- readr::read_rds("projects/04_lung_immune_therapy/02_clustering/revision/00.re.clustering/05.CD4.sct.new.rds.gz")
filt.patients <- c(c("P1","P10","P19","P29","P30","P13","P33","P35"),c("P36","P37","P38"),"P21","P25") # ICB-treated or liver mets
cd4.filt <- cd4[,!(cd4$patient %in% filt.patients)]
levels(cd4)

cd4.filt <- cd4.filt[,Idents(cd4.filt) != "XCL1"]
new.cluster.ids <- c("CD4-c02_Treg","CD4-c01_CXCL13","CD4-c04_other","CD4-c04_other","CD4-c04_other","CD4-c04_other","CD4-c04_other","CD4-c04_other","CD4-c03_Prolif")
names(new.cluster.ids) <- levels(cd4.filt)
cd4.filt <- RenameIdents(cd4.filt, new.cluster.ids)

DimPlot(cd4.filt, label = F, pt.size = 0.6, reduction = "umap") + scale_color_manual(values = cols[c(7,2,1,13)]) + NoAxes() + NoLegend()

#--- macrophage plot

mac.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/01.clustering/02.mac.rds.gz")

DimPlot(mac.sce.pro, label = F, pt.size = 1) + scale_color_manual(values = cols[c(4,5,7,6,11)]) + NoAxes() + NoLegend()
