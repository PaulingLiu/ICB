
#--- library

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)
library(esc)
library(meta)

#--- load cell cluster data

d1.meta <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/01.our.data.TN.rds.gz")
d2.meta <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/02.CR.DL.rds.gz")
d3.meta <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/05.NM.rds.gz")
d4.meta <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/03.NC.rds.gz")
d5.meta <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/04.NC2.rds.gz")

#--- load response data

R.patients <- c("P4","P8","P24","P7","P11","P28","P17")

resp1 <- tibble(patient = unique(d1.meta$patient)) %>%
  dplyr::filter(!(patient %in% c("P21","P25","P29"))) %>%
  dplyr::mutate(response = ifelse(patient %in% R.patients, "R", "NR")) %>%
  dplyr::mutate(data = "d1")

resp2 <- readr::read_tsv("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/06.tumor.classification/03.CR.DL.tsv")
resp3 <- readr::read_tsv("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/06.tumor.classification/04.NM.DL.tsv")
resp4 <- readr::read_tsv("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/06.tumor.classification/05.NC.tsv")
resp5 <- readr::read_tsv("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/06.tumor.classification/06.NC2.tsv")

process.rda <- function(.x, .y){
  .x %>%
    dplyr::rename(response = group) %>%
    dplyr::select(patient, response) %>%
    dplyr::mutate(patient = as.character(patient)) %>%
    dplyr::mutate(data = .y)
}

resp2 <- process.rda(resp2, "d2")
resp3 <- process.rda(resp3, "d3")
resp4 <- process.rda(resp4, "d4")
resp5 <- process.rda(resp5, "d5")

d1.meta <- d1.meta %>% dplyr::mutate(patient = as.character(patient)) %>% dplyr::inner_join(resp1, by = "patient")
d2.meta <- d2.meta %>% dplyr::mutate(patient = as.character(patient)) %>% dplyr::inner_join(resp2, by = "patient")
d3.meta <- d3.meta %>% dplyr::mutate(patient = as.character(patient)) %>% dplyr::inner_join(resp3, by = "patient")
d4.meta <- d4.meta %>% dplyr::mutate(patient = as.character(patient)) %>% dplyr::inner_join(resp4, by = "patient")
d5.meta <- d5.meta %>% dplyr::mutate(patient = as.character(patient)) %>% dplyr::inner_join(resp5, by = "patient")
d5.meta <- d5.meta %>%
  dplyr::mutate(ct1 = ifelse(stringr::str_detect(ct2,"CD8"), "CD8 T", ct1)) %>%
  dplyr::mutate(ct1 = ifelse(stringr::str_detect(ct2,"CD4"), "CD4 T", ct1)) %>%
  dplyr::mutate(ct1 = ifelse(stringr::str_detect(ct2,"Prolif"), "Prolif. T", ct1))

all.meta <- rbind(d1.meta, d2.meta, d3.meta, d4.meta, d5.meta)
all.meta <- all.meta %>% dplyr::mutate(ct2 = ifelse(is.na(ct2), "ND", ct2))

#--- label process

table(all.meta$ct1, all.meta$data)

all.meta <- all.meta %>%
  dplyr::mutate(ct1 = ifelse(ct1 %in% c("CAF","Fibroblast","SMC"), "CAF", ct1)) %>%
  dplyr::mutate(ct1 = ifelse(ct1 %in% c("Endo."), "Endo", ct1)) %>%
  dplyr::mutate(ct1 = ifelse(ct1 %in% c("pDC"), "Myeloid", ct1)) %>%
  dplyr::mutate(ct = ct1) %>%
  dplyr::mutate(ct = ifelse(ct %in% c("CD4 T","CD8 T", "Prolif. T"), "T cell", ct))

all.meta <- all.meta %>%
  dplyr::mutate(ct2 = stringr::str_replace(ct2, "CD8_","CD8-")) %>%
  dplyr::mutate(ct2 = stringr::str_replace(ct2, "CD4_","CD4-")) %>%
  dplyr::mutate(ct2 = stringr::str_replace(ct2, "_CXCL","-CXCL")) %>%
  dplyr::mutate(ct2 = stringr::str_replace(ct2, "-Prolif","_Prolif"))

all.meta <- all.meta %>%
  dplyr::mutate(ct2 = ifelse(ct1 == "pDC", "pDC", ct2)) %>%
  dplyr::mutate(ct2 = stringr::str_replace(ct2, "\\+ ","\\.")) %>%
  dplyr::mutate(ct2 = stringr::str_replace(ct2, "\\+","")) %>%
  dplyr::mutate(ct2 = stringr::str_replace(ct2, "CD1A.DC2","DC2.CD1A")) %>%
  dplyr::mutate(ct2 = stringr::str_replace(ct2, "CLEC10A.DC2","DC2.CD1C")) %>%
  dplyr::mutate(ct2 = stringr::str_replace(ct2, "CLEC9A.DC1","DC1")) %>%
  dplyr::mutate(ct2 = ifelse(stringr::str_detect(ct2, "LAMP3"),"LAMP3.DC", ct2))

all.meta <- all.meta %>%
  dplyr::mutate(ct2 = stringr::str_replace(ct2, "\\.other","_c04-other")) %>%
  dplyr::mutate(ct2 = stringr::str_replace(ct2, "\\.lymp","_c01-lymphatic")) %>%
  dplyr::mutate(ct2 = stringr::str_replace(ct2, "\\.ACKR1","_c02-ACKR1")) %>%
  dplyr::mutate(ct2 = stringr::str_replace(ct2, "\\.INS[RT]","_c03-INSR"))

table(all.meta$ct2, all.meta$data)

#--- filter datasets that did not detect certain celltypes

filt.fun <- function(.x){
  .x %>%
    dplyr::filter(n > 9) %>%
    dplyr::count(celltype, data) %>%
    dplyr::filter(celltype %in% celltypes) %>%
    dplyr::filter(n > 2) %>%
    dplyr::mutate(tag = paste(data, celltype, sep = "."))
}

#--- effect size calculation --- stage 1 

all.meta <- all.meta %>%
  dplyr::mutate(patient = ifelse(data == "d2", paste0("CR.", patient), patient)) %>%
  dplyr::mutate(patient = ifelse(data == "d3", paste0("NM.", patient), patient)) %>%
  dplyr::mutate(data = ifelse(data == "d3", "d2", data))

celltypes = unique(all.meta$ct)
celltypes = celltypes[-10]
datasets = unique(all.meta$data)

s1.input <- all.meta %>%
  dplyr::count(data, patient, response, ct) %>%
  dplyr::group_by(data, patient) %>%
  dplyr::mutate(prop = n/sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(celltype = ct)

cal.es <- function(input){
  tmp.data <- list()
  for (j in 1:length(datasets)) {
    tmp.res <- list()
    for (i in 1:length(celltypes)) {
      sda <- input %>%
        dplyr::filter(data == datasets[j]) %>%
        dplyr::filter(celltype == celltypes[i]) %>%
        dplyr::filter(response %in% c("R","NR"))
      
      if(nrow(sda) > 0){
        a <- t.test(sda$prop[sda$response == 'R'], sda$prop[sda$response == "NR"])
        
        b <- esc_t(t = a$statistic,     
                   grp1n = length(sda$prop[sda$response == 'R']), 
                   grp2n = length(sda$prop[sda$response == 'NR']), 
                   es.type="d")
        
        tmp.res[[i]] <- tibble(effect = b$es, SE = b$se, celltype = celltypes[i], pvalue = a$p.value)
      }
      if(nrow(sda) == 0){
        tmp.res[[i]] <- tibble(effect = NA, SE = NA, celltype = NA, pvalue = NA)
      }
      
    }
    tmp.data[[j]] <- Reduce(rbind, tmp.res) %>% dplyr::mutate(data = datasets[j])
  }
  Reduce(rbind, tmp.data)
}

s1 <- cal.es(s1.input)

qualify.ct <- filt.fun(s1.input) %>% dplyr::pull(tag)
s1 <- s1 %>% 
  dplyr::filter(!is.na(effect)) %>%
  dplyr::mutate(tag = paste0(data,".",celltype)) %>%
  dplyr::filter(tag %in% qualify.ct)

s1.op <- list()

for (i in 1:length(celltypes)) {
  tmp.data <- s1 %>% dplyr::filter(celltype == celltypes[i])
  a <- metagen(effect, SE, studlab = data, sm = "SMD", data = tmp.data)
  s1.op[[i]] <- tibble(
    SMD = a$TE.random,
    low = a$lower.random,
    up = a$upper.random,
    p = a$pval.random,
    celltype = celltypes[i]
  )
}

s1.op <- Reduce(rbind, s1.op)

#--- effect size calculation --- stage 2

celltypes <- c("CD4 T", "CD8 T", "Prolif. T")

s2.input <- all.meta %>%
  dplyr::count(data, patient, response, ct1) %>%
  dplyr::group_by(data, patient) %>%
  dplyr::mutate(prop = n/sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(celltype = ct1) %>%
  dplyr::filter(celltype %in% celltypes)

s2 <- cal.es(s2.input)
qualify.ct <- filt.fun(s2.input) %>% dplyr::pull(tag)
s2 <- s2 %>% 
  dplyr::filter(!is.na(effect)) %>%
  dplyr::mutate(tag = paste0(data,".",celltype)) %>%
  dplyr::filter(tag %in% qualify.ct)

s2.op <- list()

for (i in 1:length(celltypes)) {
  tmp.data <- s2 %>% dplyr::filter(celltype == celltypes[i])
  a <- metagen(effect, SE, studlab = data, sm = "SMD", data = tmp.data)
  s2.op[[i]] <- tibble(
    SMD = a$TE.random,
    low = a$lower.random,
    up = a$upper.random,
    p = a$pval.random,
    celltype = celltypes[i]
  )
}

s2.op <- Reduce(rbind, s2.op)

#--- effect size calculation --- stage 3 CD4 subsets

celltypes <- unique(all.meta$ct2) %>% sort()
celltypes <- celltypes[2:5] #celltypes[4:7]

s3.input <- all.meta %>%
  dplyr::filter(ct2 %in% celltypes) %>%
  dplyr::filter(ct2 != "ND") %>%
  dplyr::count(data, patient, response, ct2) %>%
  dplyr::group_by(data, patient) %>%
  dplyr::mutate(prop = n/sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(celltype = ct2) %>%
  dplyr::filter(celltype %in% celltypes)

s3 <- cal.es(s3.input)
qualify.ct <- filt.fun(s3.input) %>% dplyr::pull(tag)
s3 <- s3 %>% 
  dplyr::filter(!is.na(effect)) %>%
  dplyr::mutate(tag = paste0(data,".",celltype)) %>%
  dplyr::filter(tag %in% qualify.ct)

s3.op <- list()

for (i in 1:length(celltypes)) {
  tmp.data <- s3 %>% dplyr::filter(celltype == celltypes[i])
  a <- metagen(effect, SE, studlab = data, sm = "SMD", data = tmp.data)
  s3.op[[i]] <- tibble(
    SMD = a$TE.random,
    low = a$lower.random,
    up = a$upper.random,
    p = a$pval.random,
    celltype = celltypes[i]
  )
}

s3.op <- Reduce(rbind, s3.op)

#--- effect size calculation --- stage 3 CD8 subsets

celltypes <- unique(all.meta$ct2) %>% sort()
celltypes <- celltypes[6:8] #celltypes[8:10]

s4.input <- all.meta %>%
  dplyr::filter(ct2 %in% celltypes) %>%
  dplyr::filter(ct2 != "ND") %>%
  dplyr::count(data, patient, response, ct2) %>%
  dplyr::group_by(data, patient) %>%
  dplyr::mutate(prop = n/sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(celltype = ct2) %>%
  dplyr::filter(celltype %in% celltypes)

s4 <- cal.es(s4.input)
qualify.ct <- filt.fun(s4.input) %>% dplyr::pull(tag)
s4 <- s4 %>% 
  dplyr::filter(!is.na(effect)) %>%
  dplyr::mutate(tag = paste0(data,".",celltype)) %>%
  dplyr::filter(tag %in% qualify.ct)

s4.op <- list()

for (i in 1:length(celltypes)) {
  tmp.data <- s4 %>% dplyr::filter(celltype == celltypes[i])
  a <- metagen(effect, SE, studlab = data, sm = "SMD", data = tmp.data)
  s4.op[[i]] <- tibble(
    SMD = a$TE.random,
    low = a$lower.random,
    up = a$upper.random,
    p = a$pval.random,
    celltype = celltypes[i]
  )
}

s4.op <- Reduce(rbind, s4.op)

#--- effect size calculation --- stage 3 macrophage subsets

celltypes <- all.meta %>%
  dplyr::filter(ct == "Myeloid") %>%
  dplyr::distinct(ct2) %>%
  dplyr::pull(ct2)

celltypes <- celltypes[-2]

s5.input <- all.meta %>%
  dplyr::filter(ct2 %in% celltypes) %>%
  dplyr::filter(ct2 != "ND") %>%
  dplyr::count(data, patient, response, ct2) %>%
  dplyr::group_by(data, patient) %>%
  dplyr::mutate(prop = n/sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(celltype = ct2) %>%
  dplyr::filter(celltype %in% celltypes)

s5 <- cal.es(s5.input)
qualify.ct <- filt.fun(s5.input) %>% dplyr::pull(tag)
s5 <- s5 %>% 
  dplyr::filter(!is.na(effect)) %>%
  dplyr::mutate(tag = paste0(data,".",celltype)) %>%
  dplyr::filter(tag %in% qualify.ct)

s5.op <- list()

for (i in 1:length(celltypes)) {
  tmp.data <- s5 %>% dplyr::filter(celltype == celltypes[i])
  a <- metagen(effect, SE, studlab = data, sm = "SMD", data = tmp.data)
  s5.op[[i]] <- tibble(
    SMD = a$TE.random,
    low = a$lower.random,
    up = a$upper.random,
    p = a$pval.random,
    celltype = celltypes[i]
  )
}

s5.op <- Reduce(rbind, s5.op)

#--- effect size calculation --- stage 3 endothelial subsets

all.meta <- all.meta %>% 
  dplyr::mutate(ct2 = ifelse(ct1 == "CAF", "CAF", ct2))

celltypes <- all.meta %>%
  dplyr::filter(ct1 %in% c("Endo","CAF")) %>%
  dplyr::distinct(ct2) %>%
  dplyr::pull(ct2)

celltypes <- celltypes[-c(2,6)]

s6.input <- all.meta %>%
  #dplyr::filter(ct2 %in% celltypes) %>%
  #dplyr::filter(ct2 != "ND") %>%
  dplyr::count(data, patient, response, ct2) %>%
  dplyr::group_by(data, patient) %>%
  dplyr::mutate(prop = n/sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(celltype = ct2) %>%
  dplyr::filter(celltype %in% celltypes) %>%
  dplyr::filter(data != "d5")

s6 <- cal.es(s6.input)
qualify.ct <- filt.fun(s6.input) %>% dplyr::pull(tag)
s6 <- s6 %>% 
  dplyr::filter(!is.na(effect)) %>%
  dplyr::mutate(tag = paste0(data,".",celltype)) %>%
  dplyr::filter(tag %in% qualify.ct)

s6.op <- list()

for (i in 1:length(celltypes)) {
  tmp.data <- s6 %>% dplyr::filter(celltype == celltypes[i])
  a <- metagen(effect, SE, studlab = data, sm = "SMD", data = tmp.data)
  s6.op[[i]] <- tibble(
    SMD = a$TE.random,
    low = a$lower.random,
    up = a$upper.random,
    p = a$pval.random,
    celltype = celltypes[i]
  )
}

s6.op <- Reduce(rbind, s6.op)


#--- plot
wda <- rbind(s1.op, s2.op, s3.op, s4.op, s5.op, s6.op, fibro.op)
wda %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/07.wda.rds.gz", compress = "gz")
wda <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/07.wda.rds.gz")
   
cts <- wda$celltype
rank.cts <- tibble(
  celltype = cts[c(1,12,11,18,17,19,10,13,16,14,15,4,3,7,9,6,20,21,26,28,29,22,25,24,23,27,2,5,33,31,34,32,8,35:40)],
  group = 1:39
)

wda <- wda %>%
  dplyr::inner_join(rank.cts, by = "celltype") %>%
  dplyr::arrange(group)

S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
IS_sqrt <- function(x){x^2*sign(x)}
S_sqrt_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt)

tcellsubsets <- wda[c(4:6,8:11),]$celltype
tcellsubgroups <- celltypes
majortypes <- celltypes

wda.sub <- wda %>%
  dplyr::filter(celltype %in% c(tcellsubsets))

wda.sub %>%
  dplyr::mutate(group = 1:nrow(.)) %>%
  dplyr::mutate(label = ifelse(p < 0.05 & SMD > 0.5, "sig.up", "not")) %>%
  dplyr::mutate(label = ifelse(p < 0.05 & SMD < -0.5, "sig.down", label)) %>%
  #dplyr::mutate(SMD = ifelse(SMD > 0, sqrt(SMD), -sqrt(-SMD))) %>%
  #dplyr::mutate(low = ifelse(low > 0, sqrt(low), -sqrt(-low))) %>%
  #dplyr::mutate(up = ifelse(up > 0, sqrt(up), -sqrt(-up))) %>%
  ggplot() +
  geom_hline(yintercept = 0, color = "black", lwd = 0.6, linetype = "dashed") +
  geom_hline(yintercept = 0.8, color = "grey", lwd = 0.6, linetype = "dashed") +
  geom_hline(yintercept = -0.8, color = "grey", lwd = 0.6, linetype = "dashed") +
  geom_segment(aes(x = group, xend = group, y = low, yend = up, color = label), lwd = 0.9) +
  geom_point(aes(group, SMD, color = label), size = 2.6) +
  theme_classic() +
  scale_x_continuous(
    breaks = c(1:length(wda.sub$celltype)),
    label = wda.sub$celltype
  ) +
  scale_y_continuous(trans="S_sqrt",breaks=c(-1,0,1,4,-4)) +
  #coord_flip() +
  labs(
    x = "",
    y = "Effect size (SMD)"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", angle = 90, hjust = 1)
  ) +
  scale_color_manual(values = c("grey70", "#00868B", "#CD2626"))# +
  ylim(-2.5,2)

#--- plot 2

fibro.s1 <- fibro.s1 %>%
  dplyr::mutate(data = ifelse(data == "d3","d4",data)) %>%
  dplyr::mutate(data = ifelse(data == "d4","d5",data))

ida <- rbind(s1,s2,s3,s4,s5,s6,fibro.s1)
ida %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/08.ida.rds.gz", compress = "gz")
ida <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/08.ida.rds.gz")

rank.cts <- tibble(
  celltype = cts[c(1,12,11,18,17,19,10,13,16,14,15,4,3,7,9,6,20,21,26,28,29,22,25,24,23,27,2,5,33,31,34,32,8,35:40)]
)

rank1 <- rank.cts %>%
  dplyr::filter(celltype %in% tcellsubsets)

rank2 <- rank.cts %>%
  dplyr::filter(celltype %in% majortypes)

rank3 <- rank.cts %>%
  dplyr::filter(celltype %in% tcellsubgroups)

rank4 <- rank.cts %>%
  dplyr::filter(!(celltype %in% c(tcellsubsets, tcellsubgroups, majortypes)))

rank.cts <- rbind(rank1, rank2, rank3, rank4) %>% dplyr::mutate(group = 1:39)

ida <- ida %>%
  dplyr::inner_join(rank.cts, by = "celltype") %>%
  dplyr::arrange(group)

ida %>%
  dplyr::mutate(effect = ifelse(effect > 3, 3, effect)) %>%
  dplyr::mutate(effect = ifelse(effect < -3, -3, effect)) %>%
  dplyr::mutate(name = ifelse(pvalue < 0.05,"*",NA)) %>%
  ggplot(aes(factor(celltype, levels = rank.cts$celltype), factor(data, levels = paste0("d",5:1)))) +
  geom_tile(lwd = 0.3,color = "white",aes(fill = ifelse(effect >= 0, sqrt(effect), -sqrt(-effect)))) +
  scale_fill_distiller(palette = "RdBu") +
  geom_text(aes(label = name)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", angle = 90, hjust = 1, vjust = 0.5)
  )
