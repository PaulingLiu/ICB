
#--- library

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)

#--- load data

get.expr <- function(object, data){
  object <- ScaleData(object, do.center = T, do.scale = T, features = rownames(object), scale.max = 5)
  tmp.expr <- as.tibble(t(object@assays$RNA@scale.data))
  tmp.expr$cluster = as.character(Idents(object))
  tmp.expr %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise_if(is.double, ~ mean(.x, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cluster = paste(cluster, data, sep = "."))
}

re.idents <- function(.x){
  tmp.idents <- as.character(Idents(.x))
  tmp.idents[tmp.idents %in% c("CD4 T", "CD8 T")] = "T cell"
  tmp.idents[tmp.idents %in% c("pericyte & SMC")] = "SMC"
  Idents(.x) <- tmp.idents
  return(.x)
}

d1.sce <- re.idents(d1.sce)
d2.sce <- re.idents(d2.sce)
d3.sce <- re.idents(d3.sce)
d4.sce <- re.idents(d4.sce)
d5.sce <- re.idents(d5.sce)

genes <- c("CXCL13","PDCD1","BHLHE40","TOX","GNG4","IFNG","IL21","CD82","TPI1","HAVCR2",
           "TIGIT","ENTPD1","TNFRSF4","BATF","FOXP3","CCR8","IL2RA","LAYN",
           "ANXA1","LMNA","GZMA","GZMK","CD69","DUSP1",
           "TCF7","SELL","IL7R","CCR7")

expr.D1 <- get.expr(d1.sce, "D01")
expr.D2 <- get.expr(d2.sce, "D02")
expr.D3 <- get.expr(d3.sce, "D03")
expr.D4 <- get.expr(d4.sce, "D04")
expr.D5 <- get.expr(d5.sce, "D05")

expr.D1 %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/05.figure.data/02.marker.heatmap/01.major.celltype/01.D1.expr.mean.rds.gz", compress = "gz")
expr.D2 %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/05.figure.data/02.marker.heatmap/01.major.celltype/02.D2.expr.mean.rds.gz", compress = "gz")
expr.D3 %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/05.figure.data/02.marker.heatmap/01.major.celltype/03.D3.expr.mean.rds.gz", compress = "gz")
expr.D4 %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/05.figure.data/02.marker.heatmap/01.major.celltype/04.D4.expr.mean.rds.gz", compress = "gz")
expr.D5 %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/05.figure.data/02.marker.heatmap/01.major.celltype/05.D5.expr.mean.rds.gz", compress = "gz")

#
tibble(
  cellname = d5.sce$cellnames,
  cluster = as.character(Idents(d5.sce))
) %>%
  dplyr::group_by(cluster) %>%
  dplyr::sample_n(500, replace = T) %>%
  dplyr::distinct(cellname, cluster) %>%
  dplyr::ungroup() -> sample.data

tmp.sce <- d5.sce[,d5.sce$cellnames %in% sample.data$cellname]
markers <- FindAllMarkers(tmp.sce, only.pos = T, logfc.threshold = 0.1)
markers %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/05.figure.data/02.marker.heatmap/02.marker.list.rds.gz", compress = "gz")
markers %>%
  dplyr::group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  as.data.frame()

#--- heatmap

genes <- c("CD3D","CD3E","CD3G","CCL5",
           "STMN1","MKI67","CDK1","TUBB",
           "NKG7","GNLY","GZMA","PRF1",
           "MS4A1","BANK1","CD19","CD79A","JCHAIN","DERL3","MZB1",
           "LYZ","C1QA","CD14","CD74","APOE",
           "LILRA4","LILRB4","IRF7",
           "TPSB2","CLU",
           "VWF","PECAM1","CLDN5","FLT1",
           "KRT19","KRT17","FXYD3","PERP","EPCAM",
           "POSTN","PDGFRA","COL1A1","DCN",
           "ACTA2","RGS5","MYH11",
           "CXCL8","IL1R2","IL1RN","CSF3R")

expr.D1 <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/05.figure.data/02.marker.heatmap/01.major.celltype/01.D1.expr.mean.rds.gz")
expr.D2 <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/05.figure.data/02.marker.heatmap/01.major.celltype/02.D2.expr.mean.rds.gz")
expr.D3 <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/05.figure.data/02.marker.heatmap/01.major.celltype/03.D3.expr.mean.rds.gz")
expr.D4 <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/05.figure.data/02.marker.heatmap/01.major.celltype/04.D4.expr.mean.rds.gz")
expr.D5 <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/05.figure.data/02.marker.heatmap/01.major.celltype/05.D5.expr.mean.rds.gz")

process.expr <- function(.x){
  sub.tmp <- .x[,c("cluster",genes)]
  sub.tmp$cluster[stringr::str_detect(sub.tmp$cluster, "fibroblast")] <- "CAF.D01"
  sub.tmp %>%
    tidyr::gather(key = "gene", value = "expr", -cluster)
}

expr.d1 <- process.expr(expr.D1)
expr.d2 <- process.expr(expr.D2)
expr.d3 <- process.expr(expr.D3)
expr.d4 <- process.expr(expr.D4)
expr.d5 <- process.expr(expr.D5)

comb.expr <- rbind(expr.d1, expr.d2, expr.d3, expr.d4, expr.d5) %>%
  dplyr::arrange(cluster)

clusters <- as.character(unique(comb.expr$cluster))
cluster.order <- clusters[c(which(stringr::str_detect(clusters,"T cell")),
                            which(stringr::str_detect(clusters,"Prolif. T")),
                            which(stringr::str_detect(clusters,"NK")),
                            which(stringr::str_detect(clusters,"B cell")),
                            which(stringr::str_detect(clusters,"Plasma B")),
                            which(stringr::str_detect(clusters,"Myeloid")),
                            which(stringr::str_detect(clusters,"pDC")),
                            which(stringr::str_detect(clusters,"MAST")),
                            which(stringr::str_detect(clusters,"Endo")),
                            which(stringr::str_detect(clusters,"Cancer")),
                            which(stringr::str_detect(clusters,"CAF")),
                            which(stringr::str_detect(clusters,"SMC")),
                            which(stringr::str_detect(clusters,"Neutrophil")))]

comb.expr %>%
  dplyr::mutate(expr = ifelse(expr < -1, -1, expr)) %>%
  dplyr::mutate(expr = ifelse(expr > 1, 1, expr)) %>%
  ggplot(aes(factor(cluster, levels = cluster.order), factor(gene, levels = rev(genes)))) +
  geom_tile(aes(fill = expr)) +
  scale_fill_distiller(palette = "RdBu") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 10, colour = "black", angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(x = "", y = "")  

