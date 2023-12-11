
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(ggsci)
library(data.table)
library(reticulate)

#--- load all data

all.sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/15.ICB.sc/00.all.icb.pro.rds.gz")
cell.cluster <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/14.fibroblast.metacell/08.metacell.cluster/05.icb.rds.gz")


tn.sce <- readr::read_rds("projects/04_lung_immune_therapy/02_clustering/revision/00.re.clustering/01.all.sce.rds.gz")
filt.patients <- c(c("P1","P10","P19","P29","P30","P13","P33","P35"),c("P36","P37","P38"),"P21","P25")
tn.filt <- tn.sce[,!(tn.sce$patient %in% filt.patients)]
new.cluster.ids <- c("Cancer cell","Myeloid","Plasma B","NK", "B cell","pericyte & SMC","CD4 T","Prolif. T", "Cancer cell", "MAST", "CD4 T", "fibroblast", "CD8 T", "CD8 T", "CD4 T", "Cancer cell", "Cancer cell", "Endo.", "pDC", "Cancer cell")
names(new.cluster.ids) <- levels(tn.filt)
tn.filt <- RenameIdents(tn.filt, new.cluster.ids)
cell.cluster.tn <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/14.fibroblast.metacell/08.metacell.cluster/01.our.data.metacell.cluster.rds.gz")

gene1 <- "F2R"
gene2 <- "WNT10A"

d1 <- tibble(
  sample = paste(all.sce$patient, all.sce$num, sep = "."),
  cellid = all.sce$cellid,
  cluster = as.character(Idents(all.sce)),
  MMP1 = all.sce@assays$RNA@scale.data[gene1,],
  CRABP1 = all.sce@assays$RNA@scale.data[gene2,]
)

fb.sub <- cell.cluster %>%
  dplyr::filter(cluster != "NA") %>%
  dplyr::select(cellid, cluster) %>%
  dplyr::rename(sub = cluster)


d2 <- tibble(
  sample = paste(tn.filt$patient, tn.filt$num, sep = "."),
  cellid = tn.filt$cellid,
  cluster = as.character(Idents(tn.filt)),
  MMP1 = tn.filt@assays$RNA@scale.data[gene1,],
  CRABP1 = tn.filt@assays$RNA@scale.data[gene2,]
)

fb.sub.tn <- cell.cluster.tn %>%
  dplyr::filter(cluster != "NA") %>%
  dplyr::select(cellid, cluster) %>%
  dplyr::rename(sub = cluster)


d2 <- d2 %>%
  dplyr::filter(cluster != "fibroblast" | (cluster == "fibroblast" & cellid %in% fb.sub.tn$cellid)) %>%
  dplyr::left_join(fb.sub.tn, by = "cellid") %>%
  dplyr::mutate(sub = ifelse(is.na(sub), cluster, sub))

d1 <- d1 %>%
  dplyr::filter(cluster != "fibroblast" | (cluster == "fibroblast" & cellid %in% fb.sub$cellid)) %>%
  dplyr::left_join(fb.sub, by = "cellid") %>%
  dplyr::mutate(sub = ifelse(is.na(sub), cluster, sub))

d1.sort <- rbind(d2) %>%
  #dplyr::mutate(MMP1 = MMP1 + rnorm(nrow(.),0,1)/1000) %>%
  dplyr::filter(sub != 'other') %>%
  dplyr::group_by(sub) %>%
  dplyr::mutate(mean.1 = mean(MMP1)) %>%
  dplyr::distinct(sub, mean.1) %>%
  dplyr::arrange(desc(mean.1))

rbind(d1, d2) %>%
  #dplyr::mutate(MMP1 = MMP1 + rnorm(nrow(.),0,1)/1000) %>%
  dplyr::filter(sub != 'other') %>%
  dplyr::group_by(sub) %>%
  dplyr::mutate(mean.1 = mean(MMP1)) %>%
  ggplot(aes(factor(sub, levels = d1.sort$sub[c(1:4,9,5:8,10:17)]), MMP1)) +
  geom_violin(aes(fill = mean.1), scale = "width", trim = T) +
  scale_fill_distiller(palette = "Oranges", direction = 1) +
  theme(
    #legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", angle = 45, hjust = 1)
  ) +
  labs(
    x = "",
    y = ""
  )

d2.sort <- rbind(d2) %>%
  #dplyr::mutate(MMP1 = MMP1 + rnorm(nrow(.),0,1)/1000) %>%
  dplyr::filter(sub != 'other') %>%
  dplyr::group_by(sub) %>%
  dplyr::mutate(mean.1 = mean(CRABP1)) %>%
  dplyr::distinct(sub, mean.1) %>%
  dplyr::arrange(desc(mean.1))

rbind(d2) %>%
  dplyr::group_by(sub) %>%
  dplyr::filter(sub != "other") %>%
  dplyr::mutate(mean.1 = mean(CRABP1)) %>%
  ggplot(aes(factor(sub, levels = d2.sort$sub[c(1:3,5,7,4,6,8:17)]), CRABP1)) +
  geom_violin(aes(fill = mean.1), scale = "width", trim = T) +
  scale_fill_distiller(palette = "Oranges", direction = 1) +
  theme(
    #legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", angle = 45, hjust = 1)
  ) +
  labs(
    x = "",
    y = ""
  )

comb.da <- rbind(d1,d2) %>%
  dplyr::filter(sub != 'other') %>%
  dplyr::group_by(sub) %>%
  dplyr::summarize(F2R = mean(MMP1), WNT10A = mean(CRABP1)) %>%
  dplyr::arrange(desc(F2R))

comb.da %>%
  ggplot(aes(factor(sub, levels = comb.da$sub), 1)) +
  geom_tile(aes(fill = F2R)) +
  scale_fill_distiller(palette = "RdBu") +
  theme(
    #legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", angle = 45, hjust = 1)
  ) +
  labs(
    x = "",
    y = ""
  )

comb.da$F2R <- (comb.da$F2R-mean(comb.da$F2R))/sd(comb.da$F2R)
comb.da$WNT10A <- (comb.da$WNT10A-mean(comb.da$WNT10A))/sd(comb.da$WNT10A)
