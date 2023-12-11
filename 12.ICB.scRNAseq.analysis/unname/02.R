
d1.meta %>%
  #dplyr::filter(ct1 %in% c("CD4 T", "CD8 T", "Prolif. T","Myeloid","Plasma B","B cell", "NK", "MAST","pDC")) %>%
  dplyr::count(group, patient, ct1) %>%
  dplyr::group_by(group, patient) %>%
  dplyr::mutate(f = n/sum(n)) %>%
  ggplot(aes(ct1, f)) +
  geom_boxplot(aes(color = factor(group, levels = c("pre.R","post.R","post.NR"))), outlier.size = -1) +
  geom_jitter(aes(color =  factor(group, levels = c("pre.R","post.R","post.NR"))), position = position_dodge(width = 0.75)) +
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

pda.pro <- pda %>% dplyr::filter(gene != "CXCR4")

pda.pro %>%
  ggplot(aes(label, factor(gene, levels = rev(unique(pda.pro$gene))))) +
  geom_tile(aes(fill = zscore)) +
  scale_fill_distiller(palette = "RdBu") +
  theme_classic() +
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
