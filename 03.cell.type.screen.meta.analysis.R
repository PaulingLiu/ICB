
#--- library

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)
library(esc)
library(meta)

#--- load data

d1 <- readr::read_rds("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/11.major.cell.prop/01.NSCLC.rds.gz")
d2 <- readr::read_tsv("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/11.major.cell.prop/02.CR.DL.tsv")
d3 <- readr::read_tsv("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/11.major.cell.prop/03.NM.DL.tsv")
d4 <- readr::read_tsv("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/11.major.cell.prop/04.NC.tsv")
d5 <- readr::read_tsv("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/11.major.cell.prop/05.NC2.tsv")

#--- effect size calculation

d2.pro <- rbind(d2, d3)
celltypes <- intersect(unique(d1$celltype), unique(d5$celltype))

cal.es <- function(input, data){
  tmp.res <- list()
  
  for (i in 1:length(celltypes)) {
    sda <- input %>%
      dplyr::filter(celltype == celltypes[i]) %>%
      dplyr::filter(response %in% c("R","NR"))
    
    a <- t.test(sda$prop[sda$response == 'R'], sda$prop[sda$response == "NR"])
    
    b <- esc_t(t = a$statistic,     
               grp1n = length(sda$prop[sda$response == 'R']), 
               grp2n = length(sda$prop[sda$response == 'NR']), 
               es.type="d")
    
    tmp.res[[i]] <- tibble(effect = b$es, SE = b$se, celltype = celltypes[i], data = data, pvalue = a$p.value)
  }
  
  Reduce(rbind, tmp.res)
}

res1 <- cal.es(d1, "d1") #4 PD 7 PR
res2 <- cal.es(d2.pro, "d2") #6 PD 4 PR
res3 <- cal.es(d4, "d3") #5 PD 6 PR
res4 <- cal.es(d5, "d4") #2 PD 2 PR

comb.data <- rbind(res1, res2, res3, res4) 

names(comb.data)[1] <- "Effect"
names(comb.data)[2] <- "SE"

res <- list()
for (i in 1:length(celltypes)) {
  tmp.data <- comb.data %>% dplyr::filter(celltype == celltypes[i])
  a <- metagen(Effect, SE, studlab = studies, sm = "SMD", data = tmp.data)
  res[[i]] <- tibble(
    SMD = a$TE.random,
    low = a$lower.random,
    up = a$upper.random,
    p = a$pval.random,
    celltype = celltypes[i]
  )
}

#--- plot

wda <- Reduce(rbind, res) %>%
  dplyr::bind_rows(m.res[[1]]) %>%
  dplyr::arrange(desc(SMD)) %>%
  dplyr::mutate(group = 12:1)

wda %>%
  dplyr::mutate(label = ifelse(p < 0.05 & SMD > 0.8, "sig.up", "not")) %>%
  dplyr::mutate(label = ifelse(p < 0.05 & SMD < -0.8, "sig.down", label)) %>%
  ggplot() +
  geom_hline(yintercept = 0, color = "black", lwd = 0.6, linetype = "dashed") +
  geom_segment(aes(x = group, xend = group, y = low, yend = up, color = label), lwd = 0.9) +
  geom_point(aes(group, SMD, color = label), size = 2.6) +
  theme_classic() +
  scale_x_continuous(
    breaks = c(12:1),
    label = wda$celltype
  ) +
  coord_flip() +
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
    axis.text.x = element_text(color="black")
  ) +
  scale_color_manual(values = c("grey70","#00868B", "#CD2626"))

#--- MMP1 fibroblasts


m2 <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/02.MMP1.prop.rds.gz")
m3 <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/10.external.data/04.NC2/03.MMP1.prop.rds.gz")

m1 <- m1 %>% dplyr::mutate(celltype = "MMP1+ fibroblast")
m2 <- m2 %>% dplyr::mutate(celltype = "MMP1+ fibroblast") %>% dplyr::mutate(response = group)
m3 <- m3 %>% dplyr::mutate(celltype = "MMP1+ fibroblast") %>% dplyr::mutate(response = group)


celltypes <- c("MMP1+ fibroblast")

m2.cr <- m2 %>%
  dplyr::filter(data == "CR")

m2.nm <- m2 %>%
  dplyr::filter(data == "NM")


m.res1 <- cal.es(m1, "d1") #4 PD 7 PR
m.res2 <- cal.es(m2.cr, "d2") #6 PD 4 PR
m.res3 <- cal.es(m2.nm, "d3")
m.res5 <- cal.es(m3, "d5") #5 PD 6 PR

comb.res <- rbind(m.res1, m.res2, m.res3, m.res5) 
comb.res


names(comb.res)[1] <- "Effect"
names(comb.res)[2] <- "SE"

m.res <- list()
for (i in 1:length(celltypes)) {
  tmp.data <- comb.res %>% dplyr::filter(celltype == celltypes[i])
  a <- metagen(Effect, SE, studlab = data, sm = "SMD", data = tmp.data)
  m.res[[i]] <- tibble(
    SMD = a$TE.random,
    low = a$lower.random,
    up = a$upper.random,
    p = a$pval.random,
    celltype = celltypes[i]
  )
}
