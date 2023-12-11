
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(ggsci)
library(data.table)
library(reticulate)
library(esc)
library(meta)

#--- major cell prop

major.prop <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/00.major.prop.rds.gz")
major.prop <- d1.meta
major.prop %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/00.major.prop.new.rds.gz", compress = "gz")
mp <- major.prop %>%
  dplyr::count(group, patient, ct2) %>%
  dplyr::group_by(group, patient) %>%
  dplyr::mutate(f = n/sum(n)) %>%
  dplyr::ungroup()

tmp.pre <- list()
tmp.postR <- list()

cts <- unique(mp$ct2)
for (i in 1:length(cts)) {
  tmp.mp <- mp %>%
    dplyr::filter(ct2 == cts[i])
  
  a <- t.test(tmp.mp$f[tmp.mp$group == 'pre.R'], tmp.mp$f[tmp.mp$group == "post.NR"])
  
  b <- esc_t(t = a$statistic,     
             grp1n = length(tmp.mp$f[tmp.mp$group == 'pre.R']), 
             grp2n = length(tmp.mp$f[tmp.mp$group == 'post.NR']), 
             es.type="d")
  
  tmp.pre[[i]] <- tibble(effect = b$es, SE = b$se, celltype = cts[i], pvalue = a$p.value, group = "preR")
}

for (i in 1:length(cts)) {
  tmp.mp <- mp %>%
    dplyr::filter(ct2 == cts[i])
  
  a <- t.test(tmp.mp$f[tmp.mp$group == 'post.R'], tmp.mp$f[tmp.mp$group == "post.NR"])
  
  b <- esc_t(t = a$statistic,     
             grp1n = length(tmp.mp$f[tmp.mp$group == 'post.R']), 
             grp2n = length(tmp.mp$f[tmp.mp$group == 'post.NR']), 
             es.type="d")
  
  tmp.postR[[i]] <- tibble(effect = b$es, SE = b$se, celltype = cts[i], pvalue = a$p.value, group = "postR")
}

tmp.pre <- Reduce(rbind, tmp.pre)
tmp.postR <- Reduce(rbind, tmp.postR)

rbind(tmp.pre, tmp.postR) %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/06.prop.effect.size/01.major.effect.size.rds.gz", compress = "gz")

#--- cell subtype proportion

get.es <- function(.x){
  cts <- unique(.x$cluster)
  
  res.pre <- list()
  res.postR <- list()
  
  for (i in 1:length(cts)) {
    tmp.mp <- .x %>%
      dplyr::filter(cluster == cts[i])
    
    a <- t.test(tmp.mp$prop[tmp.mp$group == 'pre.R'], tmp.mp$prop[tmp.mp$group == "post.NR"])
    
    b <- esc_t(t = a$statistic,     
               grp1n = length(tmp.mp$prop[tmp.mp$group == 'pre.R']), 
               grp2n = length(tmp.mp$prop[tmp.mp$group == 'post.NR']), 
               es.type="d")
    
    res.pre[[i]] <- tibble(effect = b$es, SE = b$se, celltype = cts[i], pvalue = a$p.value, group = "preR")
  }
  
  for (i in 1:length(cts)) {
    tmp.mp <- .x %>%
      dplyr::filter(cluster == cts[i])
    
    a <- t.test(tmp.mp$prop[tmp.mp$group == 'post.R'], tmp.mp$prop[tmp.mp$group == "post.NR"])
    
    b <- esc_t(t = a$statistic,     
               grp1n = length(tmp.mp$prop[tmp.mp$group == 'post.R']), 
               grp2n = length(tmp.mp$prop[tmp.mp$group == 'post.NR']), 
               es.type="d")
    
    res.postR[[i]] <- tibble(effect = b$es, SE = b$se, celltype = cts[i], pvalue = a$p.value, group = "postR")
  }
  res.pre <- Reduce(rbind, res.pre)
  res.postR <- Reduce(rbind, res.postR)
  rbind(res.pre, res.postR)
}

get.es(rbind(dc.prop, mac.prop, fb.prop, endo.prop, ps.prop)) %>%
  readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/06.prop.effect.size/02.nonT.sub.effect.size.rds.gz", compress = "gz")

#--- CD8 and CD4 subsets

cd8.prop.pro <- cd8.prop %>%
  dplyr::mutate(group = response) %>%
  dplyr::mutate(group = ifelse(response %in% c("pre"),"pre.R", group)) %>%
  dplyr::mutate(group = ifelse(response %in% c("PR"),"post.R", group)) %>%
  dplyr::mutate(group = ifelse(response %in% c("PD"),"post.NR", group)) %>%
  dplyr::mutate(`CD8_c02-CXCL13-` = 1-freq) %>%
  dplyr::rename(`CD8_c01-CXCL13+` = freq) %>%
  dplyr::select(patient, `CD8_c01-CXCL13+`, `CD8_c02-CXCL13-`, group) %>%
  tidyr::gather(key = cluster, value = prop, -patient, -group)

cd4.prop.pro <- cd4.prop %>%
  dplyr::mutate(group = response) %>%
  dplyr::mutate(group = ifelse(response %in% c("pre"),"pre.R", group)) %>%
  dplyr::mutate(group = ifelse(response %in% c("PR"),"post.R", group)) %>%
  dplyr::mutate(group = ifelse(response %in% c("PD"),"post.NR", group)) %>%
  dplyr::rename(`CD4_c01-CXCL13+` = `CXCL13+`) %>%
  dplyr::rename(`CD4_c02-CXCL13-` = `CXCL13-`) %>%
  dplyr::rename(`CD4_c03-Treg` = Treg) %>%
  dplyr::select(-response) %>%
  tidyr::gather(key = cluster, value = prop, -patient, -group)

get.es(rbind(cd8.prop.pro, cd4.prop.pro)) %>%
  readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/06.prop.effect.size/03.T.sub.effect.size.rds.gz", compress = "gz")

#--- validation data

tcell <- get.es(rbind(cd8.prop.pro, cd4.prop.pro))
other <- get.es(rbind(dc.prop, mac.prop, fb.prop, endo.prop, ps.prop))
major.p <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/06.prop.effect.size/01.major.effect.size.rds.gz")


crn <- rbind(tcell, other, major.p) %>%
  dplyr::distinct(celltype) %>%
  dplyr::rename(cellname = celltype) %>%
  dplyr::mutate(celltype = cellname)

dm <- c(3.346,-3.0169,-1.028,1.065,0.75,-0.316605858,-0.875425821,-0.209816237,-0.365295434,0.181864318,0.203609721,-0.01922208,0.1081663,-0.512979,NA,0.2978347,NA,-0.3568433,-0.2800554,-0.123673,-1.263161,NA,NA,-0.30768663,-0.11837877,-0.13046263,-0.07786908,0.2442592,0.3611994,-0.84263,	-0.05626146,0.11027024,-0.7364435,-0.2491697,-0.2290711,-0.02787833,0.7291312,-0.6150711,NA)

crn$celltype <- c(
  crn$celltype[1:5],"DC2-CD1A","DC2-CD1C","DC1",crn$celltype[9:10],"Marco_c04-C1QA","Marco_c01-FCN1",
  "Marco_c03-FOLR2","Marco_c05-MARCO","MT-CO1","Marco_c02-SPP1","Fibro_c01-PI16","Fibro_c02-ADH1B",
  "Fibro_c03-ACTA2","Fibro_c04-COL18A1","Fibro_c05-MMP1","NA","other",crn$celltype[24:37],"CAF","pericyte & SMC"
)

crn <- crn %>% dplyr::mutate(vd.meta = as.double(dm))

rbind(tcell, other, major.p) %>%
  dplyr::rename(cellname = celltype) %>%
  dplyr::inner_join(crn, by = "cellname") %>%
  dplyr::filter(group == "postR") %>%
  dplyr::filter(!is.na(vd.meta)) %>%
  dplyr::filter(!stringr::str_detect(celltype, "CD8")) %>%
  dplyr::filter(!stringr::str_detect(celltype, "CD4")) %>%
  ggplot(aes(vd.meta, effect)) +
  geom_point(aes(color = celltype)) +
  theme_light() +
  xlim(-1.4,1.4) +
  ylim(-1.6,1.6) +
  theme_linedraw() +
  geom_text_repel(aes(label = celltype)) +
  scale_color_manual(values = comb.d3) +
  theme(
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )

#--- discovery cohort

wda <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/12.sub.cell.prop.meta/07.wda.rds.gz")
