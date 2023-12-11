
# TCR analysis

library(tidyverse)
library(Seurat)
library(ggrepel)
library(ggbeeswarm)
library(reticulate)
library(ggsci)
library(ROGUE)
library(ggpubr)


#------------------------- Functions -----------------------#

process.config <- function(.x){
  .x <- .x %>% 
    dplyr::filter(high_confidence %in% c("TRUE","True") & is_cell == "TRUE" & productive %in% c("TRUE","True")) %>%
    dplyr::filter(full_length == "TRUE" & chain %in% c("TRA","TRB") & umis >= 2) %>%
    .[,c(1,5:16)] %>%
    dplyr::mutate(
      identifier = purrr::pmap_chr(
        .l = list(
          .x = v_gene,
          .y = cdr3_nt,
          .z = j_gene
        ),
        .f = function(.x, .y, .z){
          paste(.x, .y, .z, sep = "_")
        }
      ))
  
  return(.x)
}
process.tcr <- function(cda){
  tibble(cellid = unique(cda$barcode)) %>%
    dplyr::mutate(
      info = purrr::map(
        .x = cellid,
        .f = function(.x){
          TCR_reformat <- TCRReFormat(cda %>% dplyr::filter(barcode == .x))
          return(TCR_reformat)
        }
      )
    ) -> a
  
  aa <- data.frame(do.call(rbind, a$info), check.names = F)
  return(aa)
}

TCRReFormat <- function(data){
  reformat_cellinfo <- as.character(as.vector(data[1,c("barcode")]))
  TRA_data <- data %>% dplyr::filter(chain == "TRA") %>% dplyr::arrange(desc(umis), desc(reads)) %>% head(2)
  TRB_data <- data %>% dplyr::filter(chain == "TRB") %>% dplyr::arrange(desc(umis), desc(reads)) %>% head(2)
  reformat_TRA1 <- TCRextract(TRA_data,1)
  reformat_TRA2 <- TCRextract(TRA_data,2)
  reformat_TRB1 <- TCRextract(TRB_data,1)
  reformat_TRB2 <- TCRextract(TRB_data,2)
  TCR_info <- c(reformat_cellinfo, reformat_TRA1, reformat_TRA2, reformat_TRB1, reformat_TRB2)
  names(TCR_info) <- c("CellName",
                       as.vector(outer(c("Identifier","CDR3","CDR3_nt","V_gene","J_gene","nRead","nUMI","Length","Full_length"), 
                                       c("(Alpha1)","(Alpha2)","(Beta1)","(Beta2)"), 
                                       paste, sep="")
                       )
  )
  return(TCR_info)
}
TCRextract <- function(data, rowid){
  if(rowid > nrow(data)){
    info_extract <- rep(NA,9)
  }else{
    info_extract <- c(data[rowid,"identifier"],
                      as.character(as.vector(data[rowid,c("cdr3","cdr3_nt","v_gene","j_gene","reads","umis","length","full_length")]))
    )
  }
  return(info_extract)
}

#------------------------ load data -----------------------#

tcr.file <- c(
  "/raid1/pauling/projects/01_data/01_lung_immune_therapy/39.P39.post/filtered_contig_annotations.csv",
  "/raid1/pauling/projects/01_data/01_lung_immune_therapy/42.Pxx.WFB.post/filtered_contig_annotations.csv",
  "/raid1/pauling/projects/01_data/01_lung_immune_therapy/43.PXX.ZXZ.post/filtered_contig_annotations.csv"
)

tcr.data <- tibble(
  path = tcr.file,
  patient = c("P39.tr.1","P40.tr.1","P41.tr.1")
) %>%
  dplyr::mutate(
    contig = purrr::map(
      .x = path,
      .f = function(.x){
        readr::read_csv(file = .x)
      }
    )
  ) %>%
  dplyr::select(-path)

#------------------------ Process tcr data -----------------------#

tcr.data %>%
  dplyr::mutate(
    contig = purrr::map2(
      .x = patient,
      .y = contig,
      .f = function(.x, .y){
        tmp1 <- process.config(.y)
        tmp2 <- process.tcr(tmp1)
        tmp2 <- tmp2 %>% 
          dplyr::mutate(patient = .x) %>% 
          dplyr::mutate(cellid = paste0(patient, ".", tmp2$CellName)) %>%
          dplyr::filter(!is.na(`CDR3(Alpha1)`)) %>%
          dplyr::filter(!is.na(`CDR3(Beta1)`))
        
        print(.x)
        return(tmp2)
      }
    )
  ) -> tcr.data

merge.tcr <- Reduce(rbind, tcr.data$contig)
merge.tcr %>% readr::write_rds("/home/pauling/projects/04_lung_immune_therapy/04_TCR_analysis/06.icb.newdata.tcr.rds.gz", compress = "gz")
