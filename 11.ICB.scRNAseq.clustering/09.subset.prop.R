
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(ggsci)
library(data.table)
library(reticulate)

#--- functions to calculate cell subset frequency

count.fun <- function(.x){
  as.data.frame(table(paste0(.x$patient, ".",.x$num))) %>%
    as.tibble() %>%
    dplyr::rename(patient = Var1, n = Freq) %>%
    dplyr::mutate(patient = as.character(patient))
}

cal.prop <- function(.x, .y){
  as.data.frame(table(paste0(.x$patient, ".",.x$num), as.character(Idents(.x)))) %>%
    as.tibble() %>%
    dplyr::rename(patient = Var1, cluster = Var2, freq = Freq) %>%
    dplyr::left_join(.y, by = "patient") %>%
    dplyr::mutate(prop = freq/n) %>%
    dplyr::mutate(cluster = as.character(cluster), patient = as.character(patient))
}


add.rp <- function(.x){
  .x %>% dplyr::left_join(rda, by = c("patient" = "sample"))
}
#--- load all data

all.sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/15.ICB.sc/00.all.icb.pro.rds.gz")

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

rda <- rda %>% 
  dplyr::inner_join(rda.pro, by = "sample") %>%
  dplyr::select(sample, timepoint, response, group) %>%
  dplyr::distinct(sample, timepoint, group)

#--- load macrophges and DCs

mac.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/15.ICB.sc/05.macrophage.rds.gz")
dc.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/15.ICB.sc/05.dc.rds.gz")

mye.sce <- all.sce[,Idents(all.sce) %in% c("Myeloid","MAST")]

d2.count <- count.fun(mye.sce)
mac.prop <- cal.prop(mac.sce, d2.count)
dc.prop <- cal.prop(dc.sce, d2.count)

mac.prop <- add.rp(mac.prop)
dc.prop <- add.rp(dc.prop)
endo.prop <- add.rp(endo.prop)
ps.prop <- add.rp(ps.prop)

ps.prop %>%
  ggplot(aes(cluster, prop)) +
  geom_boxplot(aes(color = group))

save(dc.prop, mac.prop, fb.prop, endo.prop, ps.prop, file = "projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/05.mac.dc.prop.rda")

#--- load stromal cell subsets

endo.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/15.ICB.sc/06.endo.rds.gz")
ps.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/15.ICB.sc/07.SMC.pc.rds.gz")


stromal.sce <- all.sce[,Idents(all.sce) %in% c("Endo.","fibroblast","pericyte & SMC")]
d1.count <- count.fun(stromal.sce)
endo.prop <- cal.prop(endo.sce, d1.count)
ps.prop <- cal.prop(ps.sce, d1.count)

fb.prop <- fb.prop %>% dplyr::select(-dataset)

rbind(endo.prop, fb.prop, ps.prop) %>%
  dplyr::mutate(response = ifelse(group %in% c("pre.R","post.R"),"R","NR")) %>%
  dplyr::filter(n > 99) %>%
  dplyr::filter(!(cluster %in% c("NA","other"))) %>%
  ggplot(aes(factor(cluster, levels = c(unique(endo.prop$cluster), unique(fb.prop$cluster), unique(ps.prop$cluster))), prop)) +
  geom_boxplot(aes(color = factor(response, levels = c("R","NR"))), outlier.size = -1) +
  geom_point(aes(color = factor(response, levels = c("R","NR")), group = factor(response, levels = c("R","NR")), shape = `group`), position = position_jitterdodge(jitter.width = 0.23, jitter.height = 0), size = 1.2, stroke = 0.7) +
  scale_shape_manual(values = c(1,1,2)) +
  theme(
    #legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", angle = 45, hjust = 1)
  ) +
  scale_color_manual(values = c("#EEB4B4", "#00868B")) +
  labs(
    x = "",
    y = ""
  )

stromal.prop <- rbind(endo.prop, fb.prop, ps.prop) %>%
  dplyr::mutate(response = ifelse(group %in% c("pre.R","post.R"),"R","NR")) %>%
  dplyr::filter(n > 99) %>%
  dplyr::filter(!(cluster %in% c("NA","other")))

clusters <- unique(stromal.prop$cluster)

for (i in 1:length(clusters)) {
  print(clusters[i])
  print(wilcox.test(stromal.prop$freq[stromal.prop$response == "R" & stromal.prop$cluster == clusters[i]],
         stromal.prop$freq[stromal.prop$response == "NR" & stromal.prop$cluster == clusters[i]]))
}

#--- Dirichlet regression test

counts <- stromal.prop %>%
  dplyr::select(patient, cluster, freq) %>%
  tidyr::spread(key = cluster, value = freq) %>%
  dplyr::select(-patient)

counts[is.na(counts)] <- 0

covariates <- stromal.prop %>%
  dplyr::distinct(patient, response) %>%
  dplyr::rename(sample = patient, condition = response) %>%
  dplyr::mutate(condition = ifelse(condition == "R", "NR", "R"))

counts = as.data.frame(counts)
counts$counts = DR_data(counts)
data = cbind(counts, covariates)
fit = DirichReg(counts ~ condition, data)

# Get p-values
u = summary(fit)
pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
v = names(pvals)
pvals = matrix(pvals, ncol=length(u$varnames))
rownames(pvals) = gsub('condition', '', v[1:nrow(pvals)])
colnames(pvals) = u$varnames
fit$pvals = pvals

pvals

mac.rank <- c("FCN1+","SPP1+","FOLR2+","C1QA+","MARCO+","MT-CO1")
dc.rank <- c("CLEC9A+ DC1","CD1A+ DC2", "CLEC10A+ DC2","LAMP3+ DC", "pDC")

rbind(mac.prop, dc.prop) %>%
  dplyr::mutate(response = ifelse(group %in% c("pre.R","post.R"),"R","NR")) %>%
  dplyr::filter(n > 100) %>%
  dplyr::filter(!(cluster %in% c("NA","other"))) %>%
  ggplot(aes(factor(cluster, levels = c(mac.rank, dc.rank)), prop)) +
  #geom_boxplot(aes(color = response), outlier.size = -1) +
  #geom_point(aes(color = response, group = response, shape = `group`), position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), size = 1.2, stroke = 0.7) +
  geom_boxplot(aes(color = factor(group,levels = c("pre.R","post.R","post.NR"))), outlier.size = -1) +
  geom_point(aes(color = factor(group,levels = c("pre.R","post.R","post.NR")), group = factor(group,levels = c("pre.R","post.R","post.NR"))), position = position_jitterdodge(jitter.width = 0.15, jitter.height = 0), size = 1.5) +
  scale_shape_manual(values = c(1,1,2)) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", angle = 45, hjust = 1)
  ) +
  scale_color_manual(values =  c("#EEB4B4", "#CD2626", "#00868B")) +
  labs(
    x = "",
    y = ""
  )

#--- load T cell subsets

cd8.prop <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/01.cd8.sub.prop.rds.gz")
cd4.prop <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/02.cd4.sub.prop.new.rds.gz")

stromal.prop <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/03.stromal.prop.rds.gz")

cd8.prop %>%
  dplyr::mutate(group = ifelse(response %in% c("pre","PR"),"R","NR")) %>%
  dplyr::mutate(neg = 1-freq) %>%
  dplyr::rename(pos = freq) %>%
  dplyr::select(patient, pos, neg, group, response) %>%
  tidyr::gather(key = cluster, value = prop, -patient, -group, -response) %>%
  dplyr::mutate(response = ifelse(response == "pre","01",response)) %>%
  dplyr::mutate(response = ifelse(response == "PR","02",response)) %>%
  dplyr::mutate(response = ifelse(response == "PD","03",response)) %>%
  ggplot(aes(cluster, prop)) +
  #geom_boxplot(aes(color = group), outlier.size = -1) +
  #geom_point(aes(color = group, group = group, shape = `response`), position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), size = 1.2, stroke = 0.7) +
  geom_boxplot(aes(color = response), outlier.size = -1) +
  geom_point(aes(color = response, group = response), position = position_jitterdodge(jitter.width = 0.15, jitter.height = 0), size = 1.5) +
  scale_shape_manual(values = c(1,1,2)) +
  theme(
    #legend.position = "none",
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

t.test(cd8.prop$freq[cd8.prop$response != "PD"], cd8.prop$freq[cd8.prop$response == "PD"])

cd4.prop %>%
  dplyr::mutate(group = ifelse(response %in% c("pre","PR"),"R","NR")) %>%
  tidyr::gather(key = cluster, value = prop, -patient, -group, -response) %>%
  dplyr::mutate(response = ifelse(response == "pre","01",response)) %>%
  dplyr::mutate(response = ifelse(response == "PR","02",response)) %>%
  dplyr::mutate(response = ifelse(response == "PD","03",response)) %>%
  ggplot(aes(cluster, prop)) +
  geom_boxplot(aes(color = response), outlier.size = -1) +
  geom_point(aes(color = response, group = response), position = position_jitterdodge(jitter.width = 0.15, jitter.height = 0), size = 1.5) +
  scale_shape_manual(values = c(1,1,2)) +
  theme(
    #legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", angle = 45, hjust = 1)
  ) +
  scale_color_manual(values =  c("#EEB4B4", "#CD2626", "#00868B")) +
  labs(
    x = "",
    y = ""
  ) +
  ylim(0,1)
