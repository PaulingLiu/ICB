
#--- library

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)

#---

get.expr <- function(object, .genes, data){
  object <- ScaleData(object, do.center = T, do.scale = T, features = .genes, scale.max = 5)
  tmp.expr <- as.tibble(t(object@assays$RNA@scale.data))
  tmp.expr$cluster = as.character(Idents(object))
  tmp.expr %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise_if(is.double, ~ mean(.x, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cluster = paste(cluster, data, sep = "."))
}


#--- CD8 T cell subset

d1.cd8 <- cd8.filt
d1.cd4 <- cd4.filt
d2.tcell <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/02.Tcell.sub.rds.gz")
d3.tcell <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/02.tcell.sce.rds.gz")
d4.tcell <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/03.NC/02.tcell.sce.rds.gz")
d5.tcell <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/04.NC2/02.tcell.sce.rds.gz")

cd8.genes <- c("IL7R","GZMK","CXCR4","CD28","GPR183","GNLY",
               "CXCL13","PDCD1","HAVCR2","CTLA4","TIGIT","LAG3","TNFRSF9",
               "STMN1","MKI67","TUBB","CDK1","TYMS")

cd8.d1 <- get.expr(d1.cd8, cd8.genes, "D01")
cd8.d1$cluster <- c("CD8-c01_CXCL13-.D01","CD8-c02_CXCL13+.D01","CD8-c03_Prolif.D01")
cd8.d2 <- get.expr(d2.tcell[,stringr::str_detect(Idents(d2.tcell), "CD8")], cd8.genes, "D02")
cd8.d3 <- get.expr(d3.tcell[,stringr::str_detect(Idents(d3.tcell), "CD8")], cd8.genes, "D03")
cd8.d4 <- get.expr(d4.tcell[,stringr::str_detect(Idents(d4.tcell), "CD8")], cd8.genes, "D04")
cd8.d5 <- get.expr(d5.tcell[,stringr::str_detect(Idents(d5.tcell), "CD8")], cd8.genes, "D05")

comb.expr <- rbind(cd8.d1, cd8.d2, cd8.d3, cd8.d4, cd8.d5) %>%
  dplyr::arrange(cluster)

comb.expr %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/05.figure.data/02.marker.heatmap/02.sub.celltype/01.cd8.rds.gz", compress = "gz")

clusters <- as.character(unique(comb.expr$cluster))
comb.expr %>%
  tidyr::gather(key = "gene", value = "expr", -cluster) %>%
  dplyr::mutate(expr = ifelse(expr < -1, -1, expr)) %>%
  dplyr::mutate(expr = ifelse(expr > 1, 1, expr)) %>%
  ggplot(aes(factor(cluster, levels = clusters), factor(gene, levels = rev(cd8.genes)))) +
  geom_tile(aes(fill = expr)) +
  theme_classic() +
  theme(strip.text.x = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5, size = 10), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.grid.major = element_blank()) + 
  #facet_grid(. ~ group, scales = "free", space = "free") + 
  scale_fill_gradient2(low = "white", mid = "white", high = pal_material("teal")(10)[c(1,5,10)]) +
  labs(x = "", y = "")


#--- CD4 subset

cd4.genes <- c("IL7R","GZMK","CD69","ANXA1",
               "PDCD1","CXCL13","BHLHE40","TOX","IFNG","IL21",
               "TNFRSF4","FOXP3","IL2RA","LAYN","ENTPD1",
               "STMN1","MKI67","TUBB","CDK1","TYMS")

cd4.d1 <- get.expr(d1.cd4, cd4.genes, "D01")
cd4.d2 <- get.expr(d2.tcell[,stringr::str_detect(Idents(d2.tcell), "CD4")], cd4.genes, "D02")
cd4.d3 <- get.expr(d3.tcell[,stringr::str_detect(Idents(d3.tcell), "CD4")], cd4.genes, "D03")
cd4.d4 <- get.expr(d4.tcell[,stringr::str_detect(Idents(d4.tcell), "CD4")], cd4.genes, "D04")
cd4.d5 <- get.expr(d5.tcell[,stringr::str_detect(Idents(d5.tcell), "CD4")], cd4.genes, "D05")

comb.expr <- rbind(cd4.d1, cd4.d2, cd4.d3, cd4.d4, cd4.d5) %>%
  dplyr::arrange(cluster)

comb.expr %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/05.figure.data/02.marker.heatmap/02.sub.celltype/02.cd4.rds.gz", compress = "gz")

clusters <- as.character(unique(comb.expr$cluster))
cluster.order <- clusters[c(which(stringr::str_detect(clusters,"c04")),
                            which(stringr::str_detect(clusters,"c01")),
                            which(stringr::str_detect(clusters,"c02")),
                            which(stringr::str_detect(clusters,"c03")))]

comb.expr %>%
  tidyr::gather(key = "gene", value = "expr", -cluster) %>%
  dplyr::mutate(expr = ifelse(expr < -1.5, -1.5, expr)) %>%
  dplyr::mutate(expr = ifelse(expr > 1.5, 1.5, expr)) %>%
  ggplot(aes(factor(cluster, levels = cluster.order), factor(gene, levels = rev(cd4.genes)))) +
  geom_tile(aes(fill = expr)) +
  theme_classic() +
  theme(strip.text.x = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5, size = 10), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.grid.major = element_blank()) + 
  #facet_grid(. ~ group, scales = "free", space = "free") + 
  scale_fill_gradient2(low = "white", mid = "white", high = pal_material("teal")(10)[c(1,5,10)]) +
  labs(x = "", y = "")


#--- DC subset
mac.d2 <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/03.myeloid.rds.gz")



d1.dc <- readr::read_rds("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/01.clustering/01.dc.sce.rds.gz")
d2.dc <- mac.d2[,stringr::str_detect(Idents(mac.d2), "DC")]
d3.dc <- dc.sce
d4.dc <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/03.NC/03.dc.rds.gz")
d5.dc <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/10.external.data/04.NC2/03.dc.rds.gz")


markers <- FindAllMarkers(d1.dc, logfc.threshold = 0.1, only.pos = T)
markers %>%
  dplyr::group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  as.data.frame()

dc.genes <- c("CLEC9A","RAB7B","C1orf54","CPNE3",
               "CD1A","FCER1A","CD207",
               "CLEC10A","FCGR2B","ANXA1",
               "LAMP3","CCR7","CCL19","CCL22",
              "LILRA4","IRF7","ITM2C","IL3RA","NR3C1")

dc.d1 <- get.expr(d1.dc, dc.genes, "D01")
dc.d2 <- get.expr(d2.dc, dc.genes, "D02")
dc.d2$cluster <- c("LAMP3+ DC.D02","CD1A+ DC2.D02","CLEC10A+ DC2.D02","CLEC9A+ DC1.D02",'pDC.D02')
dc.d3 <- get.expr(d3.dc, dc.genes, "D03")
dc.d4 <- get.expr(d4.dc, dc.genes, "D04")
dc.d4$cluster <- c("CLEC9A+ DC1.D04","CD1A+ DC2.D04","CLEC10A+ DC2.D04","LAMP3+ DC.D04","pDC.D04")
dc.d5 <- get.expr(d5.dc, dc.genes, "D05")
dc.d5$cluster <- c("CLEC9A+ DC1.D05","CD1A+ DC2.D05","CLEC10A+ DC2.D05","LAMP3+ DC.D05","pDC.D05")

comb.expr <- rbind(dc.d1, dc.d2, dc.d3, dc.d4, dc.d5) %>%
  dplyr::arrange(cluster)

comb.expr %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/05.figure.data/02.marker.heatmap/02.sub.celltype/03.dc.rds.gz", compress = "gz")

clusters <- as.character(unique(comb.expr$cluster))
cluster.order <- clusters[c(which(stringr::str_detect(clusters,"DC1")),
                            which(stringr::str_detect(clusters,"CD1A")),
                            which(stringr::str_detect(clusters,"CLEC10A")),
                            which(stringr::str_detect(clusters,"LAMP3")),
                            which(stringr::str_detect(clusters,"pDC")))]

comb.expr %>%
  tidyr::gather(key = "gene", value = "expr", -cluster) %>%
  dplyr::mutate(expr = ifelse(expr < -1, -1, expr)) %>%
  dplyr::mutate(expr = ifelse(expr > 1, 1, expr)) %>%
  ggplot(aes(factor(cluster, levels = cluster.order), factor(gene, levels = rev(dc.genes)))) +
  geom_tile(aes(fill = expr)) +
  theme_classic() +
  theme(strip.text.x = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5, size = 10), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.grid.major = element_blank()) + 
  #facet_grid(. ~ group, scales = "free", space = "free") + 
  scale_fill_gradient2(low = "white", mid = "white", high = pal_material("teal")(10)[c(1,5,10)]) +
  labs(x = "", y = "")

#--- macrophage
d1.mac <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/01.clustering/06.mac.tn.rds.gz")
d2.mac <- mac.d2[,!(stringr::str_detect(Idents(mac.d2), "DC"))]
d3.mac <- mac.pro
d4.mac <- mac.sce
d5.mac <- mac.sce

markers4 <- FindAllMarkers(d4.mac, logfc.threshold = 0.1, only.pos = T)
markers2 %>%
  dplyr::group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  as.data.frame()

mac.genes <- c("FCN1","S100A8","VCAN","TIMP1","CD55","FPR1",
              "SPP1","FABP5","GAPDH","TPI1",
              "FOLR2","CCL18","F13A1","CCL13",
              "APOE","CXCL9","C1QA","C1QB","C1QC",
              "MARCO","MCEMP1","PPARG","FABP4","TREM1","LPL")


mac.d1 <- get.expr(d1.mac, mac.genes, "D01")
mac.d2 <- get.expr(d2.mac, mac.genes, "D02")
mac.d3 <- get.expr(d3.mac, mac.genes, "D03")
mac.d4 <- get.expr(d4.mac, mac.genes, "D04")
mac.d5 <- get.expr(d5.mac, mac.genes, "D05")


comb.expr <- rbind(mac.d1, mac.d2, mac.d3, mac.d4, mac.d5) %>%
  dplyr::mutate(cluster = stringr::str_remove(cluster, "\\+")) %>%
  dplyr::arrange(cluster) 

comb.expr %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/05.figure.data/02.marker.heatmap/02.sub.celltype/05.mac.rds.gz", compress = "gz")

clusters <- as.character(unique(comb.expr$cluster))
clusters <- clusters[-21]
cluster.order <- clusters[c(which(stringr::str_detect(clusters,"FCN1")),
                            which(stringr::str_detect(clusters,"SPP1")),
                            which(stringr::str_detect(clusters,"FOLR2")),
                            which(stringr::str_detect(clusters,"C1Q")),
                            which(stringr::str_detect(clusters,"MARCO")))]

comb.expr %>%
  tidyr::gather(key = "gene", value = "expr", -cluster) %>%
  dplyr::mutate(expr = ifelse(expr < -1, -1, expr)) %>%
  dplyr::mutate(expr = ifelse(expr > 1.2, 1.2, expr)) %>%
  ggplot(aes(factor(cluster, levels = cluster.order), factor(gene, levels = rev(mac.genes)))) +
  geom_tile(aes(fill = expr)) +
  theme_classic() +
  theme(strip.text.x = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5, size = 10), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.grid.major = element_blank()) + 
  #facet_grid(. ~ group, scales = "free", space = "free") + 
  scale_fill_gradient2(low = "white", mid = "white", high = pal_material("teal")(10)[c(1,5,10)]) +
  labs(x = "", y = "")

#--- endothelial

d1.endo <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/01.clustering/05.endo.clustering.TN.sce.rds.gz")
d2.endo <- endo.pro
d3.endo <- endo.sce
d4.endo <- endo.pro
d5.endo <- endo.pro

markers <- FindAllMarkers(d1.endo, logfc.threshold = 0.1, only.pos = T)
markers %>%
  dplyr::group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  as.data.frame()

endo.genes <- c("CCL21","PDPN","TFF3","FABP4",
               "ACKR1","CLU","IFITM1","IL33",
               "INSR","FLT1","IGFBP3","KDR")


endo.d1 <- get.expr(d1.endo, endo.genes, "D01")
endo.d1$cluster <- c("Endo.lymp.D01","Endo.ACKR1.D01","Endo.INSR.D01","Endo.other.D01")
endo.d2 <- get.expr(d2.endo, endo.genes, "D02")
endo.d2$cluster[2] <- "Endo.INSR.D02"
endo.d3 <- get.expr(d3.endo, endo.genes, "D03")
endo.d3$cluster <- c("Endo.lymp.D03","Endo.ACKR1.D03","Endo.INSR.D03","Endo.other.D03")
endo.d4 <- get.expr(d4.endo, endo.genes, "D04")
endo.d5 <- get.expr(d5.endo, endo.genes, "D05")


comb.expr <- rbind(endo.d1, endo.d2, endo.d3, endo.d4, endo.d5) %>%
  dplyr::mutate(cluster = stringr::str_remove(cluster, "\\+")) %>%
  dplyr::arrange(cluster) 

comb.expr %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/05.figure.data/02.marker.heatmap/02.sub.celltype/06.endo.rds.gz", compress = "gz")

clusters <- as.character(unique(comb.expr$cluster))
cluster.order <- clusters[c(which(stringr::str_detect(clusters,"lymp")),
                            which(stringr::str_detect(clusters,"ACKR1")),
                            which(stringr::str_detect(clusters,"INSR")),
                            which(stringr::str_detect(clusters,"other")))]

comb.expr %>%
  tidyr::gather(key = "gene", value = "expr", -cluster) %>%
  dplyr::mutate(expr = ifelse(expr < -1, -1, expr)) %>%
  dplyr::mutate(expr = ifelse(expr > 1.5, 1.5, expr)) %>%
  ggplot(aes(factor(cluster, levels = cluster.order), factor(gene, levels = rev(endo.genes)))) +
  geom_tile(aes(fill = expr)) +
  theme_classic() +
  theme(strip.text.x = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5, size = 10), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.grid.major = element_blank()) + 
  #facet_grid(. ~ group, scales = "free", space = "free") + 
  scale_fill_gradient2(low = "white", mid = "white", high = pal_material("teal")(10)[c(1,5,10)]) +
  labs(x = "", y = "")
