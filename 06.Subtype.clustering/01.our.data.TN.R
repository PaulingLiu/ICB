
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(ggsci)
library(data.table)
library(reticulate)
use_python("/home/pauling/anaconda3/bin/python")

#--- load data

all.sce <- readr::read_rds("projects/04_lung_immune_therapy/02_clustering/revision/00.re.clustering/01.all.sce.rds.gz")
all.sce <- readr::read_rds("projects/04_lung_immune_therapy/02_clustering/revision/00.re.clustering/01.all.sce.rds.gz")
new.cluster.ids <- c("Cancer cell","Myeloid","Plasma B","NK", "B cell","SMC","CD4 T","Prolif. T", "Cancer cell", "MAST", "CD4 T", "CAF", "CD8 T", "CD8 T", "CD4 T", "Cancer cell", "Cancer cell", "Endo.", "pDC", "Cancer cell")
names(new.cluster.ids) <- levels(all.sce)
all.sce.pro <- RenameIdents(all.sce, new.cluster.ids)
DimPlot(all.sce.pro)

#--- clustering endothelial

endo.sce <- all.sce.pro[,Idents(all.sce.pro) == "Endo."]
endo.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/01.clustering/05.endo.sce.rds.gz")

endo.sce <- endo.sce[,!(endo.sce$patient %in% c("P1","P10","P19","P13","P29","P30","P33","P35","P36","P37","P38"))]
sample.count <- table(endo.sce$patient)
qualified.samples <- names(sample.count[sample.count > 14])

endo.sce <- endo.sce[,endo.sce$patient %in% qualified.samples]

endo.sce <- NormalizeData(endo.sce)
endo.sce <- ScaleData(endo.sce, features = rownames(endo.sce), do.center = F, do.scale = F)
endo.sce <- FindVariableFeatures(endo.sce, selection.method = "vst", nfeatures = 2000)
endo.sce <- RunPCA(endo.sce, features = VariableFeatures(endo.sce))
endo.sce <- RunUMAP(endo.sce, dims = 1:10)
endo.sce <- bbknn.batch(endo.sce, r = 1, n.pc = 30)

DimPlot(endo.sce, label = T)
FeaturePlot(endo.sce, features = "SELP")
markers <- FindAllMarkers(endo.sce[,], only.pos = T, test.use = "roc")

markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(endo.sce, features = top10$gene) + NoLegend()
endo.sce <- endo.sce[,Idents(endo.sce) != 9]

levels(endo.sce)
new.cluster.ids <- c("Endo_c02-ACKR1","Endo_c02-ACKR1","Endo_c07-MT-CO1","Endo_c01-lymphatic","Endo_c07-MT-CO1","Endo_c03-IGFBP5","Endo_c04-INSR","Endo_c05-PGF","Endo_c06-MKI67")
names(new.cluster.ids) <- levels(endo.sce)
endo.sce <- RenameIdents(endo.sce, new.cluster.ids)
DimPlot(endo.sce)

genes <- c("ACKR1","CCL14","ADIRF","MT-CO1","MT-CO3","CCL21","PDPN","PROX1","IGFBP5","PTGDS","FLT1","INSR","KDR","PDXN","PGF","APLN","STMN1","MKI67","TYMS")

endo.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/01.clustering/05.endo.clustering.TN.sce.rds.gz")

levels(endo.sce)
new.cluster.ids <- c("Endo_c02-ACKR1","Endo_c04-other","Endo_c01-lymphatic","Endo_c04-other","Endo_c03-INSR","Endo_c04-other","Endo_c04-other")
names(new.cluster.ids) <- levels(endo.sce)
d1.endo <- RenameIdents(endo.sce, new.cluster.ids)

#--- comparison

R.patients <- c("P4","P8","P24","P7","P11","P28")

table(paste0(endo.sce$patient, ".",endo.sce$num), Idents(endo.sce)) %>%
  as.data.frame.array() -> count.matr

count.matr <- count.matr/rowSums(count.matr)

count.matr %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::gather(key = "celltype", value = "prop", -sample) %>%
  tidyr::separate(sample, c("patient","num"), remove = F) %>%
  dplyr::left_join(omda[,c(1,7)], by = "patient") %>%
  dplyr::mutate(response = ifelse(patient %in% c("P1","P10","P19","P30","P13","P33","P35"), ifelse(num == 0, "R.real.pre","R.real.post"), response)) %>%
  dplyr::mutate(response = ifelse(patient %in% c("P36","P37","P38"), "NR.real", response)) %>%
  dplyr::mutate(response = ifelse(num %in% c(2), "NR.real", response)) %>%
  dplyr::left_join(ctda, by = "patient") %>%
  dplyr::mutate(response = ifelse(response %in% c("R","NR"), ifelse(patient %in% R.patients, "R", "NR"), response)) %>%
  #dplyr::filter(response %in% c("R.real.pre","R.real.post","NR.real")) %>%
  #dplyr::mutate(response = ifelse(response == "R.real.pre","R.real.post",response)) %>%
  dplyr::filter(response %in% c("R","NR")) %>%
  #dplyr::filter(celltype == "Prolif. T") %>%
  #dplyr::filter(type == "LUSC") %>%
  #ggplot(aes(factor(response, levels = c("R.real.pre","R.real.post","NR.real")), prop, color = response)) +
  #ggplot(aes(factor(response, levels = c("R","NR")), prop, color = response)) +
  ggplot(aes(celltype, prop)) +
  geom_boxplot(outlier.size = -1, aes(color = response)) +
  geom_point(aes(color = response), position = position_jitterdodge(jitter.width = 0.1)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  labs(
    x = "",
    y = "Frequency in stromal cells"
  )



#--- macrophage clustering

mac.sce <- mac.sce[,!(mac.sce$patient %in% c("P1","P10","P19","P13","P29","P30","P33","P35","P36","P37","P38"))]
sample.count <- table(mac.sce$patient)
qualified.samples <- names(sample.count[sample.count > 14])

#mac.sce <- mac.sce[,mac.sce$patient %in% qualified.samples]

mac.sce <- NormalizeData(mac.sce)
mac.sce <- ScaleData(mac.sce, features = rownames(mac.sce), do.center = F, do.scale = F)
mac.sce <- FindVariableFeatures(mac.sce, selection.method = "vst", nfeatures = 2000)
mac.sce <- RunPCA(mac.sce, features = VariableFeatures(mac.sce))
mac.sce <- RunUMAP(mac.sce, dims = 1:10)
mac.sce <- FindNeighbors(mac.sce, dims = 1:10)
mac.sce <- FindClusters(mac.sce, resolution = 0.5)

#--- BBKNN

mac.sce.filt <- bbknn.batch(mac.sce.filt, r = 2, n.pc = 30)
DimPlot(mac.sce.filt, label = T)
levels(mac.sce.filt)
FeaturePlot(mac.sce.filt, feature = "FCN1")
markers.9 <- FindMarkers(mac.sce.filt, ident.1 = 13)
#mac.sce.filt <- mac.sce.bbknn[,!(Idents(mac.sce.bbknn) %in% c(12,9,19,16,20))]
#mac.sce.filt <- mac.sce.filt[,!(Idents(mac.sce.filt) %in% c(13,11))]

new.cluster.ids <- c("SPP1+","SPP1+","SPP1+","C1QA+","SPP1+","FCN1+","FCN1+","FCN1+","FCN1+","SPP1+","MARCO+","FOLR2+","SPP1+","FCN1+","SPP1+","SPP1+")
names(new.cluster.ids) <- levels(mac.sce.filt)
mac.sce.pro <- RenameIdents(mac.sce.filt, new.cluster.ids)
DimPlot(mac.sce.pro)
mac.sce.pro %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/01.clustering/06.mac.tn.rds.gz", compress = "gz")
