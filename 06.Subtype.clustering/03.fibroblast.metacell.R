
#--- library

library(Seurat)
library(tidyverse)
library(ggplot2)
library(reticulate)
library(ggsci)
library(metacell)
library(ks)
library(RColorBrewer)

sce <- readr::read_rds("/home/pauling/projects/05_nsclcpd1.part2/02.data/10.external.data/02.NM/01.all.sce.rds.gz")
fb.sce.nm <- sce[,Idents(sce) == 14]

fb.sce.nm <- NormalizeData(fb.sce.nm)
fb.sce.nm <- ScaleData(fb.sce.nm, features = rownames(fb.sce.nm), do.center = F, do.scale = F)
fb.sce.nm <- FindVariableFeatures(fb.sce.nm, selection.method = "vst", nfeatures = 2000)
fb.sce.nm <- RunPCA(fb.sce.nm, features = VariableFeatures(fb.sce.nm))
fb.sce.nm <- bbknn.batch(fb.sce.nm, r = 1, n.pc = 30)

DimPlot(fb.sce.nm)
FeaturePlot(fb.sce.pro, features = "MMP2")
FeaturePlot(fb.sce.nm, features = c("COL12A1"))

markers <- FindMarkers(fb.sce.pro, ident.1 = "PI16+", only.pos = T)

#--- metacell analysis

fb.comb <- merge(fb.sce, y = fb.sce.nm)
fb.comb <- NormalizeData(fb.comb)
fb.comb <- ScaleData(fb.comb, do.center = F, do.scale = F, features = rownames(fb.comb))
sce = Seurat::as.SingleCellExperiment(fb.comb)


scdb_init("projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/", force_reinit=T)
scfigs_init("projects/05_nsclcpd1.part2/02.data/10.external.data/01.CR/")

mat1 = scm_import_sce_to_mat(.sce)
scdb_add_mat("ut", mat1)

mat <- mat1

nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
ig_genes = c(grep("^IGJ", nms, v=T), 
             grep("^IGH",nms,v=T),
             grep("^IGK", nms, v=T), 
             grep("^IGL", nms, v=T))


bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes))
mcell_mat_ignore_genes(new_mat_id="ut", mat_id="ut", bad_genes, reverse=F) 

mcell_add_gene_stat(gstat_id="ut", mat_id="ut", force=T)

mcell_gset_filter_varmean(gset_id="ut_feats", gstat_id="ut", T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id = "ut_feats", gstat_id="ut", T_tot=100, T_top3=2)

mcell_add_cgraph_from_mat_bknn(mat_id="ut", 
                               gset_id = "ut_feats", 
                               graph_id="ut_graph",
                               K=.k, #30, 80
                               dsamp=F)

mcell_coclust_from_graph_resamp(
  coc_id="ut_coc500", 
  graph_id="ut_graph",
  min_mc_size=.size, 
  p_resamp=0.75, n_resamp=500) #500)

mcell_mc_from_coclust_balanced(
  coc_id="ut_coc500", 
  mat_id= "ut",
  mc_id= "ut_mc", 
  K=.k, min_mc_size=.size, alpha=2)

a <- metacell::scdb_mc("ut_mc")

mda <- tibble(
  metaid = a@mc,
  cellid = names(a@mc)
) %>%
  dplyr::arrange(desc(metaid))
  

mda <- metacell.fun(sce)

#mda <- readr::read_rds("/raid1/pauling/projects/04_lung_immune_therapy/13.In.silico.FACS/02.meta.info/01.NSCLC.PR.rds.gz")
mcells <- unique(mda$metaid)

meta.expr <- list()

for (i in 1:length(mcells)) {
  tmp.cell <- mda$cellid[mda$metaid == mcells[i]]
  meta.expr[[i]] <- Matrix::rowMeans(fb.comb@assays$RNA@scale.data[,tmp.cell])
}

pr.matr <- Reduce(rbind, meta.expr)
pr.matr <- as.matrix(pr.matr)

rownames(pr.matr) <- mcells

tibble(
  PDGFRA = pr.matr[,"PDGFRA"],
  COL1A1 = pr.matr[,"COL1A1"],
  COL7A1 = pr.matr[,"COL7A1"],
  MMP1 = pr.matr[,"MMP1"],
  MMP2 = pr.matr[,"MMP2"],
  COL18A1 = pr.matr[,"COL18A1"],
  LRRC15 = pr.matr[,"LRRC15"],
  ACTA2 = pr.matr[,"ACTA2"]
) -> pda

pda <- as.data.frame(pda)
rownames(pda) <- rownames(pr.matr)

sub.pda <- pda[,c(1,8)]
fhat <- kde(x=as.data.frame(sub.pda))
plot(fhat, display="filled.contour2",cont=seq(0,90,by=10), xlim = c(-0.2,1.7), ylim = c(-0.2,3), cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2, col = c("white",brewer.pal(n = 9, name = 'GnBu')))
points(sub.pda, pch = 19, cex=0.8, col = rgb(45, 67,121, 100, maxColorValue = 255))

cut1 <- 0.2
cut2 <- 0.5

abline(h=cut2,lty=2,col="grey50", lwd = 1.5)
abline(v=cut1,lty=2,col="grey50", lwd = 1.5)

#--- clustering

sce.comb <- CreateSeuratObject(counts = exp(t(pr.matr))-1)
sce.comb <- NormalizeData(sce.comb)
sce.comb <- ScaleData(sce.comb, features = rownames(sce.comb), do.center = F, do.scale = F)
sce.comb <- FindVariableFeatures(sce.comb, selection.method = "vst", nfeatures = 2000)
sce.comb <- RunPCA(sce.comb, features = c("LRRC15","FAP","MMP2","POSTN","COL11A1","MMP1","RGS5","PDGFRA","MYH11","ACTA2","CLU","PI16","CD34"))
sce.comb <- RunUMAP(sce.comb, dims = 1:3)

DimPlot(sce.comb)
FeaturePlot(sce.comb, features = "MMP2")
