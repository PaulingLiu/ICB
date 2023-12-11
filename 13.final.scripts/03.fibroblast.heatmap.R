
#--- library

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)

#--- load data
fb.sce.tn <- readr::read_rds("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/01.clustering/03.fibroblast.new.tn.rds.gz")
col.genes <- rownames(fb.sce.tn)[stringr::str_detect(rownames(fb.sce.tn), "^COL[1-9]")] %>% sort()

fb.sub <- fb.sce.tn[,!(Idents(fb.sce.tn) %in% c("SMC","perycite"))]

#--- collagen genes
col.means <- rowMeans(fb.sub@assays$RNA@scale.data[col.genes,])
use.genes <- names(col.means[col.means > 0.1])


tmp.expr <- as.tibble(t(fb.sub@assays$RNA@scale.data[use.genes,]))
tmp.expr$cluster = as.character(Idents(fb.sub))
tmp.expr %>%
  tidyr::gather(key = "gene", value = "expr", -cluster) %>%
  ggplot(aes(cluster, expr))

fb.sub$cluster = as.character(Idents(fb.sub))

levels(fb.sub) <- c("PI16+","ADH1B","LRRC15+ myCAF","COL18A1+","MMP1+")

VlnPlot(fb.sub, c("COL14A1","COL16A1","COL10A1","COL11A1","COL18A1","COL4A1","COL5A3","COL7A1","COL27A1"), 
        fill.by = "ident", 
        cols = c("#8968CD", "#EEAEEE", "#B9D3EE","#009ACD","#20B2AA"),
        stack = TRUE, 
        sort = F,
        flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Identity on x-axis")
