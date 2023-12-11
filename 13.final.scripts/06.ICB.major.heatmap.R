
#--- library

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)

#--- load data

icb.sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/15.ICB.sc/00.all.icb.pro.rds.gz")

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

d1.sce <- re.idents(icb.sce)
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
           "ACTA2","RGS5","MYH11")

expr.D1 <- get.expr(d1.sce[genes,])

#--- heatmap 

process.expr <- function(.x){
  sub.tmp <- .x[,c("cluster",genes)]
  sub.tmp$cluster[stringr::str_detect(sub.tmp$cluster, "fibroblast")] <- "CAF.D01"
  sub.tmp %>%
    tidyr::gather(key = "gene", value = "expr", -cluster)
}
expr.D1 <- process.expr(expr.D1)

clusters <- as.character(unique(expr.D1$cluster))
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
                            which(stringr::str_detect(clusters,"SMC")))]

expr.D1 %>%
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
