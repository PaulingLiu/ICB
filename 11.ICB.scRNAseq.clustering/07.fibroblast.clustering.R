
#-- library

library(Seurat)
library(tidyverse)
library(ggplot2)
library(reticulate)
library(ggsci)
library(metacell)
library(ks)
library(RColorBrewer)

#-- load data

all.sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/15.ICB.sc/00.all.icb.pro.rds.gz")
DimPlot(all.sce, label = T)

#-- perycite and SMC

ps.sce <- all.sce[,Idents(all.sce) %in% c("pericyte & SMC")]
DimPlot(ps.sce)
FeaturePlot(ps.sce, features = "RGS5")
ps.sce <- FindVariableFeatures(ps.sce, selection.method = "vst", nfeatures = 2000)
ps.sce <- RunPCA(ps.sce, features = VariableFeatures(ps.sce))
ps.sce <- ps.sce[,ps.sce$patient %in% names(table(ps.sce$patient)[table(ps.sce$patient) > 4])]
ps.sce <- bbknn.batch(ps.sce, r = 0.8, n.pc = 30)
DimPlot(ps.sce, label =  T)
FeaturePlot(ps.sce, features = "MYH11")
new.cluster.ids <- ifelse(levels(ps.sce) %in% c(1,4), "SMC","pericyte")
names(new.cluster.ids) <- levels(ps.sce)
ps.sce <- RenameIdents(ps.sce, new.cluster.ids)
ps.sce %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/15.ICB.sc/07.SMC.pc.rds.gz", compress = "gz")

#-- fibroblasts

FeaturePlot(all.sce, features = "PDGFRA")
fb.sce <- all.sce[,Idents(all.sce) == "fibroblast"]
table(fb.sce$patient)

#--- overall cell count

sda <- table(st.sce.tn$patient) %>%
  as.data.frame() %>%
  as.tibble() %>%
  dplyr::rename(patient = Var1, n = Freq)

#--- meta cell calculation

sce = Seurat::as.SingleCellExperiment(fb.sce)

scdb_init("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/02.fibroblast/02.metacell/", force_reinit=T)
scfigs_init("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/02.fibroblast/02.metacell/")

mat1 = scm_import_sce_to_mat(sce)
scdb_add_mat("ut", mat1)

mat <- mat1

nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
ig_genes = c(grep("^IGJ", nms, v=T), 
             grep("^IGH",nms,v=T),
             grep("^IGK", nms, v=T), 
             grep("^IGL", nms, v=T))


bad_genes = unique(c(grep("^MT-", rownames(fb.sce.tn), v=T), grep("^MTMR", rownames(fb.sce.tn), v=T), grep("^MTND", rownames(fb.sce.tn), v=T),c("NEAT1","TMSB4X", "TMSB10"), ig_genes))
mcell_mat_ignore_genes(new_mat_id="ut", mat_id="ut", bad_genes, reverse=F) 

mcell_add_gene_stat(gstat_id="ut", mat_id="ut", force=T)

mcell_gset_filter_varmean(gset_id="ut_feats", gstat_id="ut", T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id = "ut_feats", gstat_id="ut", T_tot=100, T_top3=2)

mcell_add_cgraph_from_mat_bknn(mat_id="ut", 
                               gset_id = "ut_feats", 
                               graph_id="ut_graph",
                               K=10, #30, 80
                               dsamp=F)

mcell_coclust_from_graph_resamp(
  coc_id="ut_coc500", 
  graph_id="ut_graph",
  min_mc_size=5, 
  p_resamp=0.75, n_resamp=500) #500)

mcell_mc_from_coclust_balanced(
  coc_id="ut_coc500", 
  mat_id= "ut",
  mc_id= "ut_mc", 
  K=10, min_mc_size=5, alpha=2)

a <- metacell::scdb_mc("ut_mc")

mda <- tibble(
  metaid = a@mc,
  cellid = names(a@mc)
) %>%
  dplyr::arrange(desc(metaid))


mcells <- unique(mda$metaid)

meta.expr <- list()

for (i in 1:length(mcells)) {
  tmp.cell <- mda$cellid[mda$metaid == mcells[i]]
  meta.expr[[i]] <- Matrix::rowMeans(fb.sce@assays$RNA@scale.data[,tmp.cell])
}

pr.matr <- Reduce(rbind, meta.expr)
pr.matr <- as.matrix(pr.matr)

rownames(pr.matr) <- mcells

tibble(
  PDGFRA = pr.matr[,("PDGFRA")],
  MMP2 = pr.matr[,("MMP2")],
  GRP = pr.matr[,("GRP")],
  MMP1 = pr.matr[,("MMP1")],
  FAP = pr.matr[,("FAP")],
  COL18A1 = pr.matr[,("COL18A1")],
  LRRC15 = pr.matr[,("LRRC15")],
  PI16 = pr.matr[,("PI16")],
  ADH1B = pr.matr[,("ADH1B")],
  CCL19 = pr.matr[,("CCL19")],
  ACTA2 = pr.matr[,("ACTA2")],
  CD34 = pr.matr[,("CD34")],
  COL1A1 = pr.matr[,("COL1A1")],
) -> pda

pda <- as.data.frame(pda)
rownames(pda) <- rownames(pr.matr)

pda %>% readr::write_rds("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/14.fibroblast.metacell/10.ICB.pda.rds.gz", compress = "gz")
mda %>% readr::write_rds("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/14.fibroblast.metacell/10.ICB.mda.rds.gz", compress = "gz")

pda <- readr::read_rds("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/14.fibroblast.metacell/10.ICB.pda.rds.gz")
mda <- readr::read_rds("/raid1/pauling/projects/05_nsclcpd1.part2/02.data/14.fibroblast.metacell/10.ICB.mda.rds.gz")

#--- PI16 and ADH1B

sub.pda <- pda[,c("PI16","CD34","ADH1B","CCL19","ACTA2","MMP1","COL18A1")]
sub.pda$PI16 <- sub.pda$PI16 + sub.pda$CD34
sub.pda$ADH1B <- sub.pda$ADH1B + sub.pda$CCL19

fhat <- kde(x=as.data.frame(sub.pda[,c(1,3)]))
plot(fhat, display="filled.contour2",cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2, col = c("white",pal_material(palette = "teal")(9)))
points(sub.pda[,c(1,3)], pch = 19, cex=0.8, col = rgb(45, 67,121, 100, maxColorValue = 255))

cut1 <- 0.6
cut2 <- 1

abline(h=cut2,lty=2,col="grey50", lwd = 1.5)
abline(v=cut1,lty=2,col="grey50", lwd = 1.5)

meta.cda <- sub.pda[,c(1,3,5,6,7)] %>%
  tibble::rownames_to_column(var = "metaid") %>%
  dplyr::mutate(cluster = ifelse(ACTA2 > 0.5, "03.ACTA2", "other")) %>%
  dplyr::mutate(cluster = ifelse(COL18A1 > 1.5, "04.COL18A1", cluster)) %>%
  dplyr::mutate(cluster = ifelse(MMP1 > 0.6, "05.MMP1", cluster)) %>%
  dplyr::mutate(cluster = ifelse(PI16 >= 0.6 & ADH1B < 1, "01.PI16", cluster)) %>%
  dplyr::mutate(cluster = ifelse(ADH1B >= 1, "02.ADH1B", cluster))

get.cut.frac(sub.pda[,c(1,3)], cut1, cut2) 

#--- MMP1 & COL18A1

sub.pda <- pda[!(rownames(pda) %in% meta.cda$metaid[meta.cda$cluster %in% c('01.PI16',"02.ADH1B")]),c("COL18A1","MMP1")]
fhat <- kde(x=as.data.frame(sub.pda))
plot(fhat, display="filled.contour2",cont=seq(0,90,by=10), cex.lab=1.2, xlim = c(-0.2,3.3), ylim = c(-0.5,5.5), cex.axis=1.2, cex.main=1.2, cex.sub=1.2, col = c("white",pal_material(palette = "teal")(9)))
points(sub.pda, pch = 19, cex=0.8, col = rgb(45, 67,121, 100, maxColorValue = 255))

cut1 <- 1.5
cut2 <- 0.6

abline(h=cut2,lty=2,col="grey50", lwd = 1.5)
abline(v=cut1,lty=2,col="grey50", lwd = 1.5)

get.cut.frac(sub.pda, cut1, cut2)

#--- ACTA2

sub.pda <- pda[!(rownames(pda) %in% meta.cda$metaid[meta.cda$cluster %in% c('01.PI16',"02.ADH1B","05.MMP1","04.COL18A1")]),c("COL1A1","ACTA2")]

sub.pda$COL1A1 <- 1 + rnorm(nrow(sub.pda))/2
fhat <- kde(x=as.data.frame(sub.pda))
plot(fhat, display="filled.contour2",cont=seq(0,90,by=10), cex.lab=1.2, cex.axis=1.2, cex.main=1.2, xlim = c(-2,4), ylim = c(-0.5,4.5), cex.sub=1.2, col = c("white",pal_material(palette = "teal")(9)))
points(sub.pda, pch = 19, cex=0.8, col = rgb(45, 67,121, 100, maxColorValue = 255))

cut1 <- 0.5

abline(h=cut1,lty=2,col="grey50", lwd = 1.5)

get.cut.frac(sub.pda, cut1, cut2)


get.cut.frac <- function(.x, c1, c2){
  colnames(.x) <- c("gene1","gene2")
  n1 <- .x %>% dplyr::filter(gene1 <= c1 & gene2 <= c2) %>% nrow()
  n2 <- .x %>% dplyr::filter(gene1 > c1 & gene2 <= c2) %>% nrow()
  n3 <- .x %>% dplyr::filter(gene1 < c1 & gene2 > c2) %>% nrow()
  n4 <- .x %>% dplyr::filter(gene1 > c1 & gene2 > c2) %>% nrow()
  
  tmp.res1 <- paste0("(", c(n1,n2,n3,n4), " | ", nrow(.x), ")")
  tmp.res2 <- round(100*c(n1,n2,n3,n4)/nrow(.x),1)
  paste0(tmp.res2, "% ", tmp.res1)
}

get.cell.frac <- function(.x, .y, c1, c2){
  colnames(.x) <- c("gene1","gene2")
  tmp.mda <- .y %>% dplyr::mutate(metaid = as.character(metaid)) %>% dplyr::count(metaid)
  tmp.count <- .x %>%
    tibble::rownames_to_column(var = "metaid") %>%
    dplyr::inner_join(tmp.mda, by = "metaid")
  
  n1 <- tmp.count %>% dplyr::filter(gene1 <= c1 & gene2 <= c2) %>% dplyr::pull(n) %>% sum()
  n2 <- tmp.count %>% dplyr::filter(gene1 > c1 & gene2 <= c2) %>% dplyr::pull(n) %>% sum()
  n3 <- tmp.count %>% dplyr::filter(gene1 < c1 & gene2 > c2) %>% dplyr::pull(n) %>% sum()
  n4 <- tmp.count %>% dplyr::filter(gene1 > c1 & gene2 > c2) %>% dplyr::pull(n) %>% sum()
  
  all.n <- sum(tmp.count$n)
  
  tmp.res1 <- paste0("(", c(n1,n2,n3,n4), " | ", all.n, ")")
  tmp.res2 <- round(100*c(n1,n2,n3,n4)/all.n,1)
  paste0(tmp.res2, "% ", tmp.res1)
}

get.cut.frac(sub.pda, cut1, cut2)
get.cell.frac(sub.pda, mda, cut1, cut2)

#--- cell cluster

cda <- tibble(
  cellid = colnames(fb.sce),
  donor = fb.sce$patient,
  num = fb.sce$num
)

cell.cluster <- meta.cda %>%
  dplyr::mutate(metaid = as.numeric(metaid)) %>%
  dplyr::inner_join(mda, by = "metaid")

cell.cluster <- cda %>%
  dplyr::left_join(cell.cluster, by = "cellid") %>%
  as.tibble() %>%
  dplyr::rename(patient = donor) %>%
  dplyr::mutate(cluster = ifelse(is.na(cluster), "NA", cluster))

cell.cluster %>% readr::write_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/14.fibroblast.metacell/08.metacell.cluster/05.icb.rds.gz", compress = "gz")
cell.cluster <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/14.fibroblast.metacell/08.metacell.cluster/05.icb.rds.gz")

#--- heatmap

cell.cluster.filt <- cell.cluster %>% dplyr::filter(cluster != "NA") %>% dplyr::select(cellid, cluster)

stromal.sce <- all.sce[,Idents(all.sce) %in% c("fibroblast","pericyte & SMC")]
ps.sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/15.ICB.sc/07.SMC.pc.rds.gz")
DimPlot(ps.sce)

ps.cell <- tibble(
  cellid = ps.sce$cellid,
  cluster = as.character(Idents(ps.sce))
)


stromal.sce <- stromal.sce[,c(cell.cluster.filt$cellid, ps.cell$cellid)]
Idents(stromal.sce) = c(cell.cluster.filt$cluster, ps.cell$cluster)

fb.genes = c("IGF1",'CD34','WISP2','PI16','CFD','CLU','ADH1B','CCL19','APOD','FAP','ACTA2','TIMP3','COL11A1',
             'COL10A1','THBS2','GJB2','LRRC15','HTRA3',"COL18A1",'F2R','GRP','MMP1','MMP3','CXCL1',
             'CCL11','ADM','RGS5','NDUFA4L2','EGFL6','HIGD1B','ADIRF','MYH11')

stromal.sce <- stromal.sce[,Idents(stromal.sce) != 'other']
levels(stromal.sce)
levels(stromal.sce) <- sort(levels(stromal.sce))


DoHeatmap(stromal.sce, features = fb.genes, size = 3) + 
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")

#--- functions

count.fun <- function(.x){
  as.data.frame(table(paste0(.x$patient, ".",.x$num))) %>%
    as.tibble() %>%
    dplyr::rename(patient = Var1, n = Freq) %>%
    dplyr::mutate(patient = as.character(patient))
}

cal.prop <- function(.x, .y, .z){
  as.data.frame(table(paste0(.x$patient, ".",.x$num), .x$cluster)) %>%
    as.tibble() %>%
    dplyr::rename(patient = Var1, cluster = Var2, freq = Freq) %>%
    dplyr::left_join(.y, by = "patient") %>%
    dplyr::mutate(prop = freq/n) %>%
    dplyr::mutate(dataset = .z) %>%
    dplyr::mutate(cluster = as.character(cluster), patient = as.character(patient))
}

sce <- all.sce[,Idents(all.sce) %in% c("fibroblast","pericyte & SMC")]
d1.count <- count.fun(sce)
d1.prop <- cal.prop(cell.cluster, d1.count, "d1")

d1.prop %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/16.ICB.data.analysis/03.stromal.prop.rds.gz", compress = "gz")

d1.prop <- d1.prop %>%
  dplyr::left_join(rda, by = c("patient" = "sample"))

d1.prop %>%
  dplyr::filter(stringr::str_detect(cluster, "MMP1")) %>%
  as.data.frame() %>%
  dplyr::filter(n > 70) %>%
  dplyr::mutate(response = ifelse(group %in% c("pre.R","post.R"),"R","NR")) %>%
  #dplyr::mutate(response = c("R","R","NR","R","R","R","NR","R","R","R","R","R","R","R","NR","NR","NR","NR","NR")) %>%
  #dplyr::filter(n > 70) %>%
  ggplot(aes(response, prop)) +
  geom_boxplot(outlier.size = -1) +
  geom_jitter(width = 0.2)
