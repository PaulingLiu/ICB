
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)

#---------- load data -----------#

tcell.sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/15.ICB.sc/03.tcell.rds.gz")
all.tcr <- readr::read_rds("/raid1/pauling/projects/04_lung_immune_therapy/04_TCR_analysis/05.scTCR.processed/merge.tcr.new.rds.gz")
new.tcr <- readr::read_rds("/raid1/pauling/projects/04_lung_immune_therapy/04_TCR_analysis/06.icb.newdata.tcr.rds.gz")

cd8.sce <- tcell.sce[,stringr::str_detect(Idents(tcell.sce),"CD8")]
cd4.sce <- tcell.sce[,stringr::str_detect(Idents(tcell.sce),"CD4")]

all.tcr <- all.tcr %>% dplyr::bind_rows(new.tcr)
cd8.filt <- cd8.sce[,colnames(cd8.sce) %in% all.tcr$cellid]
cd8.tcr <- all.tcr %>% dplyr::filter(cellid %in% colnames(cd8.filt))

#---------- Tex clones ----------#

cd8.tcr <- cd8.tcr %>%
  dplyr::mutate(clone.id = paste(`Identifier(Alpha1)`,`Identifier(Beta1)`, sep = ":")) %>%
  as.tibble()

use.clones <- unique(cd8.tcr$clone.id)

use.clones.cells <- cd8.tcr %>%
  dplyr::select(patient, clone.id, cellid) %>%
  dplyr::group_by(patient, clone.id) %>%
  tidyr::nest() %>%
  dplyr::mutate(num = purrr::map_dbl(data, nrow)) %>%
  dplyr::filter(num >= 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cloneid = paste0("clone.",1:nrow(.)))

mcells <- unique(use.clones.cells$cloneid)

cal.genes <- c("CXCL13","CD8A","CD4")
meta.expr <- Matrix::Matrix(data = 0, nrow = length(mcells), ncol = length(cal.genes))

for (i in 1:length(mcells)) {
  tmp.cell <- use.clones.cells$data[[i]]$cellid
  print(i)
  if (length(tmp.cell) == 1) {
    meta.expr[i,] <- cd8.filt@assays$RNA@scale.data[cal.genes,tmp.cell]
  }else if(length(tmp.cell) > 1){
    meta.expr[i,] <- Matrix::rowMeans(cd8.filt@assays$RNA@scale.data[cal.genes,tmp.cell])
  }
}

pr.matr <- as.matrix(meta.expr)
rownames(pr.matr) <- mcells
colnames(pr.matr) <- cal.genes

tibble(
  CD8A = pr.matr[,"CD8A"],
  CXCL13 = pr.matr[,"CXCL13"],
  CD4 = pr.matr[,"CD4"],
  size = use.clones.cells$num,
  patient = use.clones.cells$patient,
  clone.id = use.clones.cells$clone.id
) -> pda

pda <- as.data.frame(pda)
rownames(pda) <- rownames(pr.matr)

#------------- cutoff determination -------------#

tibble(CD8A = pda$CD8A) %>%
  ggplot(aes(CD8A)) +
  geom_vline(xintercept = 0.6, linetype = "dashed") +
  geom_density(lwd = 1, color = "#EE6363") +
  theme_classic() +
  labs(
    x = "CD8A expression",
    y = "Density"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )

tibble(a = pda$CXCL13) %>%
  ggplot(aes(a)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_density(lwd = 1, color = "#008B8B") +
  theme_classic() +
  labs(
    x = "CXCL13 expression",
    y = "Density"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )

#------------- Tex cell frequency -------------#

tex.clones <- pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CXCL13 > 1, "Tex", "Other")) %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::pull(clone.id)

tex.clones.pro <- pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(clone.id %in% tex.clones, "Tex", "Other"))

response = c('PR',"PD","PR","pre","PR","pre","PD","PR","pre","PR","pre","PR","pre","PR","pre","PD","PD","PD","PD","PD","PD")
used.color <- c("#FA8072", "#6CA6CD", "pink","grey50")

tex.clones.pro %>%
  dplyr::group_by(patient, clone) %>%
  dplyr::summarise(n = sum(size)) %>%
  tidyr::spread(key = clone, value = n) %>%
  dplyr::mutate(freq = Tex/(Tex + Other)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(response = response) %>%
  dplyr::mutate(freq = ifelse(is.na(freq), 0 , freq)) -> pro.pda.cd8

pro.pda.cd8 %>%
  ggplot(aes(factor(response, levels = c("pre","PR","PD","PD.pre")), freq)) +
  geom_boxplot(aes(color = response), outlier.size = -1, lwd = 0.6) +
  geom_jitter(aes(color = response), width = 0.15, size = 2) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "Frequency of CD8+ tumour-reactive\ncells in all T cells"
  ) +
  scale_color_manual(values = c("#009ACD", "#CD2626", "pink","#EE799F"))

t.test(pro.pda.cd8$freq[pro.pda.cd8$response %in% c("PR","pre")], pro.pda.cd8$freq[pro.pda.cd8$response == "PD"])
pro.pda.cd8 %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/01.cd8.sub.prop.rds.gz", compress = "gz")
pda %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/01.cd8.pda.rds.gz", compress = "gz")

pro.pda.cd8 <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/01.cd8.sub.prop.rds.gz")
pda <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/01.cd8.pda.rds.gz")

pda %>%
  dplyr::filter(CXCL13 < 1) %>%
  ggplot(aes(CD8A)) +
  geom_density(fill = "#008B8B", color = "#008B8B") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  xlim(0,6) +
  labs(
    x = "CXCL13 expression"
  )

pda %>%
  dplyr::filter(CXCL13 > 1) %>%
  ggplot(aes(CD8A)) +
  geom_density(fill = "#EE6363", color = "#EE6363", alpha = 0.8) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  xlim(0,6) +
  labs(
    x = "CXCL13 expression"
  )

#----------

cd4.filt <- cd4.sce[,colnames(cd4.sce) %in% all.tcr$cellid]
cd4.tcr <- all.tcr %>% dplyr::filter(cellid %in% colnames(cd4.filt))

#---------- Tex clones ----------#

cd4.tcr <- cd4.tcr %>%
  dplyr::mutate(clone.id = paste(`Identifier(Alpha1)`,`Identifier(Beta1)`, sep = ":")) %>%
  as.tibble()

use.clones <- unique(cd4.tcr$clone.id)

use.clones.cells <- cd4.tcr %>%
  dplyr::select(patient, clone.id, cellid) %>%
  dplyr::group_by(patient, clone.id) %>%
  tidyr::nest() %>%
  dplyr::mutate(num = purrr::map_dbl(data, nrow)) %>%
  dplyr::filter(num >= 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cloneid = paste0("clone.",1:nrow(.)))

mcells <- unique(use.clones.cells$cloneid)

cal.genes <- c("CXCL13","CD8A","CD4","FOXP3","CTLA4")
meta.expr <- Matrix::Matrix(data = 0, nrow = length(mcells), ncol = length(cal.genes))

for (i in 1:length(mcells)) {
  tmp.cell <- use.clones.cells$data[[i]]$cellid
  print(i)
  if (length(tmp.cell) == 1) {
    meta.expr[i,] <- cd4.filt@assays$RNA@scale.data[cal.genes,tmp.cell]
  }else if(length(tmp.cell) > 1){
    meta.expr[i,] <- Matrix::rowMeans(cd4.filt@assays$RNA@scale.data[cal.genes,tmp.cell])
  }
}

pr.matr <- as.matrix(meta.expr)
rownames(pr.matr) <- mcells
colnames(pr.matr) <- cal.genes

tibble(
  CD8A = pr.matr[,"CD8A"],
  CXCL13 = pr.matr[,"CXCL13"],
  CD4 = pr.matr[,"CD4"],
  FOXP3 = pr.matr[,"FOXP3"],
  CTLA4 = pr.matr[,"CTLA4"],
  size = use.clones.cells$num,
  patient = use.clones.cells$patient,
  clone.id = use.clones.cells$clone.id
) -> cd4.pda

cd4.pda <- as.data.frame(cd4.pda)
rownames(cd4.pda) <- rownames(pr.matr)

#------------- cutoff determination -------------#

tibble(CD8A = cd4.pda$CD8A) %>%
  ggplot(aes(CD8A)) +
  geom_vline(xintercept = 0.6, linetype = "dashed") +
  geom_density(lwd = 1, color = "#EE6363") +
  theme_classic() +
  labs(
    x = "CD8A expression",
    y = "Density"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )

tibble(a = cd4.pda$FOXP3) %>%
  ggplot(aes(a)) +
  geom_vline(xintercept = 0.4, linetype = "dashed") +
  geom_density(lwd = 1, color = "#008B8B") +
  theme_classic() +
  labs(
    x = "CXCL13 expression",
    y = "Density"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )

#------------- Tex cell frequency -------------#

tex.clones <- cd4.pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CXCL13 > 0.5, "Tex", "Other")) %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::pull(clone.id)

tex.clones.pro <- cd4.pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(clone.id %in% tex.clones, "Tex", "Other"))

cd4.sda <- cd4.pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CXCL13 > 0.5, "CXCL13+", "CXCL13-")) %>%
  dplyr::mutate(clone = ifelse(FOXP3 > 0.4, "Treg", clone))

used.color <- c("#FA8072", "#6CA6CD", "pink","grey50")

response = c('PR',"PD","PR","pre","PR","pre","PD","PR","pre","PR","pre","PR","pre","PR","pre","PD","PD","PD","PD","PD","PD")

cd4.sda %>%
  dplyr::group_by(patient, clone) %>%
  dplyr::summarise(n = sum(size)) %>%
  dplyr::mutate(f = n/sum(n)) %>%
  dplyr::select(-n) %>%
  tidyr::spread(key = clone, value = f) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(response = response) %>%
  dplyr::mutate(`CXCL13+` = ifelse(is.na(`CXCL13+`), 0, `CXCL13+`)) -> pro.pda.cd4

pro.pda.cd4 %>%
  ggplot(aes(factor(response, levels = c("pre","PR","PD","PD.pre")), `CXCL13+`)) +
  geom_boxplot(aes(color = response), outlier.size = -1, lwd = 0.6) +
  geom_jitter(aes(color = response), width = 0.15, size = 2) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "Frequency of CD8+ tumour-reactive\ncells in all T cells"
  ) +
  scale_color_manual(values = c("#009ACD", "#CD2626", "pink","#EE799F"))

cd4.pda <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/02.cd4.pda.rds.gz")
cd4.pda %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/02.cd4.pda.rds.gz", compress = "gz")
pro.pda.cd4 %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/02.cd4.sub.prop.new.rds.gz", compress = "gz")
pro.pda.cd4 <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/02.cd4.sub.prop.new.rds.gz")

cd4.sda %>%
  dplyr::filter(clone == "CXCL13-") %>%
  ggplot(aes(FOXP3)) +
  geom_density(fill = "#EE6363", color = "#EE6363", alpha = 0.8) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  xlim(0,4.5) +
  labs(
    x = "FOXP3 expression"
  )

cd4.sda %>%
  dplyr::filter(clone == "CXCL13-") %>%
  ggplot(aes(FOXP3)) +
  geom_density(fill = "#008B8B", color = "#008B8B", alpha = 0.8) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  xlim(0,4.5) +
  labs(
    x = "FOXP3 expression"
)

cd4.sda %>%
  dplyr::filter(clone == "Treg") %>%
  ggplot(aes(FOXP3)) +
  geom_density(fill = "#9F79EE", color = "#9F79EE", alpha = 0.8) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  xlim(0,4.5) +
  labs(
    x = "FOXP3 expression"
  )

  

#---

pro.pda.cd4 %>%
  dplyr::inner_join(pro.pda.cd8, by = "patient") %>%
  ggplot(aes(freq.x, freq.y)) +
  geom_point(aes(fill = response.x, color = response.x), size = 2.5, alpha = 1, shape = 21, stroke = 0.5) +
  scale_y_sqrt() +
  scale_x_sqrt() +
  theme_classic() +
  scale_color_manual(values = c("#00868B", "#CD2626","black")) +
  scale_fill_manual(values = c("#00868B","#CD2626","#CD2626")) +
  #geom_hline(yintercept = 0.06, linetype = "dashed") +
  #geom_vline(xintercept = 0.06, linetype = "dashed")  +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  labs(
    x = "",
    y = ""
  )
