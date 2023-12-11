
#--- library

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)

# DC heatmap
dc.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/15.ICB.sc/05.dc.rds.gz")
DimPlot(dc.sce)

dc.genes <- c("CLEC9A","RAB7B","C1orf54","CPNE3",
              "CD1A","FCER1A","CD207",
              "CLEC10A","FCGR2B","ANXA1",
              "LAMP3","CCR7","CCL19","CCL22",
              "LILRA4","IRF7","ITM2C","IL3RA","NR3C1")

levels(dc.sce)
levels(dc.sce) <- c("CLEC9A+ DC1","CD1A+ DC2",'CLEC10A+ DC2','LAMP3+ DC', 'pDC')

DoHeatmap(dc.sce, features = dc.genes, size = 3) + 
  NoLegend() + 
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")

# macrophage heatmap
mac.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/15.ICB.sc/05.macrophage.rds.gz")
mac.sce <- ScaleData(mac.sce, do.scale = T, do.center = T, features = rownames(mac.sce))

aa = FindMarkers(mac.sce, ident.1 = "MT-CO1", only.pos = T)

mac.genes <- c("FCN1","S100A8","VCAN","TIMP1","CD55","FPR1",
               "SPP1","FABP5","TPI1",
               "FOLR2","F13A1",
               "APOE","C1QA","C1QB","C1QC",
               "MARCO","MCEMP1","PPARG","FABP4","TREM1","LPL",
               "MT-CO1","MT-ND4L")

levels(mac.sce)
levels(mac.sce) <- c('FCN1+', "SPP1+", "FOLR2+", "C1QA+","MARCO+","MT-CO1")

DoHeatmap(mac.sce, features = mac.genes, size = 3) + 
  scale_fill_gradient2(low = rev(c('white',"white")), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")

# Endothelial cells
endo.sce <- readr::read_rds("projects/05_nsclcpd1.part2/02.data/15.ICB.sc/06.endo.rds.gz")
endo.genes <- c("CCL21","PDPN","TFF3","FABP4",
                "ACKR1","CLU","IL33",
                "INSR","FLT1","IGFBP3","KDR")

levels(endo.sce) 
levels(endo.sce) <- sort(levels(endo.sce))

DoHeatmap(endo.sce, features = endo.genes, size = 3) + 
  NoLegend() + 
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")

