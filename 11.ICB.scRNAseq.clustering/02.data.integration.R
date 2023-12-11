
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

all.sce <- readr::read_rds("projects/04_lung_immune_therapy/02_clustering/revision/00.re.clustering/01.all.sce.rds.gz")
icb.sce <- all.sce[,all.sce$patient %in% c("P1","P10","P13","P19","P30","P33","P35","P36","P37","P38")]
icb.sce %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/15.ICB.sc/01.icb.data1.rds.gz", compress = "gz")
icb.sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/15.ICB.sc/01.icb.data1.rds.gz")

#-- load new data

matr39 <- load10x("/raid1/pauling/projects/01_data/01_lung_immune_therapy/39.P39.post/filtered_feature_bc_matrix.h5", pr = "P39.tr.1")
matr40 <- load10x("/raid1/pauling/projects/01_data/01_lung_immune_therapy/42.Pxx.WFB.post/filtered_feature_bc_matrix.h5", pr = "P40.tr.1")

#-- Seurat object

var.genes <- readr::read_rds("/raid1/pauling/projects/04_lung_immune_therapy/02_clustering/version1/12.features.rds.gz")

all.matr <- list(matr39, matr40)

matr <- Reduce(cbind, all.matr)
meta <- tibble(cellid = colnames(matr)) %>%
  dplyr::mutate(name = cellid) %>%
  tidyr::separate(name, c("patient","treatment","num","barcodes"), sep = "\\.") %>%
  dplyr::mutate(barcode = ifelse(is.na(barcodes), num, barcodes)) %>%
  dplyr::mutate(num = ifelse(is.na(barcodes), 0, num)) %>%
  dplyr::select(-barcodes)

meta <- as.data.frame(meta)
rownames(meta) <- meta$cellid

sce <- CreateSeuratObject(counts = matr, meta.data = meta)
sce %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/15.ICB.sc/02.p39.40.rds.gz", compress = "gz")
sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/15.ICB.sc/02.p39.40.rds.gz")
 
#-- combine all data

all.sce <- merge(icb.sce, y = list(sce, p43.sce))

var.genes <- readr::read_rds("/raid1/pauling/projects/04_lung_immune_therapy/02_clustering/version1/12.features.rds.gz")

all.sce <- NormalizeData(all.sce, normalization.method = "LogNormalize", scale.factor = 10000)
all.sce <- ScaleData(all.sce, features = rownames(all.sce), do.center = F, do.scale = F)
all.sce <- RunPCA(all.sce, features = var.genes)
all.sce <- RunUMAP(all.sce, dims = 1:10)

pca <- all.sce@reductions$pca@cell.embeddings
anndata = import("anndata",convert=FALSE)
sc = import("scanpy",convert=FALSE)
np = import("numpy",convert=FALSE)
bbknn = import("bbknn", convert=FALSE)
adata = anndata$AnnData(X=pca, obs=all.sce$patient)
sc$tl$pca(adata)
adata$obsm$X_pca = pca
bbknn$bbknn(adata, batch_key=0)
sc$tl$umap(adata)
sc$tl$leiden(adata, resolution = 1)
umap = py_to_r(adata$obsm['X_umap'])

tibble::tibble(
  leiden.res.1.bbknn = py_to_r(np$asarray(adata$obs$leiden)), 
  UMAP1 = umap[,1], 
  UMAP2 = umap[,2]) %>%
  dplyr::mutate(Donor = all.sce$patient) -> bb.res

all.sce <- loadumap(all.sce, umap, bb.res)

all.sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/15.ICB.sc/00.all.icb.rds.gz")

#-- annotation

sce.pro <- all.sce[,!(all.sce$patient == "P13" & all.sce$num %in% c(0,1))]
DimPlot(sce.pro, label = T)
FeaturePlot(sce.pro, features = "MS4A1")

levels(sce.pro)

new.cluster.ids <- c("Cancer cell","Myeloid","Plasma B","NK",
                     "B cell","pericyte & SMC","CD4 T","Prolif. T",
                     "MAST", "CD4 T", "fibroblast","CD8 T",
                     "CD8 T", "CD4 T", "Endo.", "CD8 T", "Myeloid","pDC")

new.cluster.ids <- c("Cancer cell","Myeloid","Plasma B","NK", 
                     "B cell","pericyte & SMC","CD4 T","Prolif. T",
                     "Cancer cell", "MAST", "CD4 T", "fibroblast",
                     "CD8 T", "CD8 T", "CD4 T", "Cancer cell", "Cancer cell", "Endo.", "pDC", "Cancer cell")


names(new.cluster.ids) <- levels(sce.pro)
sce.pro <- RenameIdents(sce.pro, new.cluster.ids)
DimPlot(sce.pro)

sce.pro <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/15.ICB.sc/00.all.icb.pro.rds.gz")

#--- plot

cols <- c("#C1CDCD","#EEA9B8","#00729A", "#EEB4B4", "#5497C0", "#62A49B", "#62A49B", "#3CB371", "#9BCD9B", "#E66764","#62A49B","#EEDFCC","#D31771")

DimPlot(sce.pro, label = T, pt.size = 0.3) + 
  scale_color_manual(values = cols[c(1:8,7,9,10,12,13:20)]) + NoLegend() + NoAxes()

DimPlot(sce.pro, label = F, pt.size = 0.1, group.by = "patient") + scale_color_manual(values = comb.d3[-c(1:3)])

all.sce.filt@reductions$umap@cell.embeddings %>%
  as.tibble() %>%
  dplyr::mutate(patient = all.sce.filt$patient) %>%
  dplyr::sample_n(nrow(.), replace = F) %>%
  ggplot(aes(UMAP_1,UMAP_2)) +
  geom_point(aes(color = patient), size = 0.1) +
  theme_classic() +
  scale_color_manual(values = comb.d3)
