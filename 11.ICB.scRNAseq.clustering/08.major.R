
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(ggsci)
library(data.table)
library(reticulate)
use_python("/home/pauling/anaconda3/bin/python")

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

all.sce %>% readr::write_rds("/home/pauling/projects/99_other/00.icb.rds.gz", compress = "gz")

#-- load data

cd8.prop <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/01.cd8.sub.prop.rds.gz")
cd4.prop <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/02.cd4.sub.prop.new.rds.gz")
dc.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/15.ICB.sc/05.dc.rds.gz")
mac.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/15.ICB.sc/05.macrophage.rds.gz")
endo.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/15.ICB.sc/06.endo.rds.gz")

#-- freq calculation

get.major <- function(.x){
  tibble(
    cellid = colnames(.x),
    patient = .x$patient,
    num = .x$num,
    group = .x$group,
    ct1 = as.character(Idents(.x))
  )
}

get.sub <- function(.meta, .x){
  tmp.meta <- get.major(.x) %>% 
    dplyr::select(-patient, -num, -group) %>%
    dplyr::rename(ct3 = ct1) 
  
  .meta %>%
    dplyr::left_join(tmp.meta, by = "cellid") %>%
    dplyr::mutate(ct2 = ifelse(!is.na(ct3), ct3, ct2)) %>%
    dplyr::select(-ct3)
}

d1.meta <- get.major(all.sce)
d1.meta <- d1.meta %>% dplyr::mutate(ct1 = ifelse(ct1 %in% c("CD4 T", "CD8 T", "Prolif. T"), "T cell", ct1))
d1.meta <- d1.meta %>% dplyr::mutate(ct1 = ifelse(ct1 %in% c("pDC"), "Myeloid", ct1))
d1.meta <- d1.meta %>% dplyr::mutate(ct2 = ct1)
d1.meta <- d1.meta %>% dplyr::mutate(ct1 = ifelse(ct1 == "pericyte & SMC", "fibroblast", ct1))
d1.meta$ct2 <- NA

d1.meta <- get.sub(d1.meta, dc.sce)
d1.meta <- get.sub(d1.meta, mac.sce)
d1.meta <- get.sub(d1.meta, endo.sce)

d1.meta %>%
  dplyr::count(group, patient, ct1) %>%
  dplyr::group_by(group, patient) %>%
  dplyr::mutate(f = n/sum(n)) %>%
  ggplot(aes(factor(ct1, levels = c("T cell","NK","B cell","Plasma B","MAST","Myeloid","Cancer cell","Endo.","fibroblast")), f)) +
  geom_boxplot(aes(color = factor(group, levels = c("pre.R","post.R","post.NR"))), outlier.size = -1) +
  geom_jitter(aes(color =  factor(group, levels = c("pre.R","post.R","post.NR"))), position = position_dodge(width = 0.75)) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", angle = 45, hjust = 1)
  ) +
  scale_color_manual(values = c("#EEB4B4", "#CD2626", "#00868B")) +
  labs(
    x = "",
    y = ""
  )

t.test(a$f[a$group == "post.NR" & a$ct1 == "T cell"], a$f[a$group == "pre.R" & a$ct1 == "T cell"])

d1.meta %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/00.major.prop.rds.gz", compress = "gz")
