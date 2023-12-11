
sce.pro$response = ifelse(sce.pro$patient %in% c("P1","P10","P13","P19","P30","P33","P35"), "R", "NR")
sce.pro$response = ifelse(sce.pro$num == 2, "NR", sce.pro$response)
sce.pro$timepoint = ifelse(sce.pro$num == 0, "pre", "post")
sce.pro$num = ifelse(sce.pro$num == 3, 4, sce.pro$num)
sce.pro$num = ifelse(sce.pro$num == 2, 3, sce.pro$num)
sce.pro$num = ifelse(sce.pro$num == 4, 2, sce.pro$num)
sce.pro$sample = paste(sce.pro$patient, sce.pro$timepoint, sce.pro$num, sep = ".")

sce.pro %>% readr::write_rds("projects/05_nsclcpd1.part2/02.data/15.ICB.sc/09.geo.icb.rds.gz", compress = "gz")
