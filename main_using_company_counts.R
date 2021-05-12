setwd("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/third_party_analyses/bhupinder_bone")

library(tidyverse)
library(tidySummarizedExperiment)
library(tidybulk)
library(tidyHeatmap)
library(patchwork)
library(future)

counts_scaled = 
  read_csv("company_derived_normalised_counts.csv") %>%
  pivot_longer(cols = -Gene_Name, names_to = "sample_id", values_to = "count") %>%
  
  # Add factor of interest
  inner_join(
    tibble(
      sample_id = c("W1", "W2", "W3", "Ko1", "Ko2", "Ko3"),
      type = c(rep("Wild_type", 3), rep("Knock_out", 3))
    ),
    by="sample_id"
  ) %>%
  as_SummarizedExperiment(.sample = sample_id, .transcript = Gene_Name , .abundance = count) %>%
  identify_abundant(factor_of_interest=type) %>%
  scale_abundance() %>%
  reduce_dimensions(method = "PCA")

# Density plot
counts_scaled %>%
  keep_abundant(factor_of_interest=type) %>%
  ggplot(aes(count_scaled+1, group=sample, color=type)) +
  geom_density() +
  scale_x_log10() +
  scale_color_brewer(palette = "Set1") +
  theme_bw()

plot_pca = 
  counts_scaled %>%
  as_tibble() %>%
  nanny::subset(sample) %>%
  ggplot(aes(PC1,PC2, label=sample)) +
  geom_point(aes(color=type)) +
  ggrepel::geom_text_repel() +
  scale_color_brewer(palette = "Set1") +
  theme_bw() 

ggsave(
  "PCA.pdf",
  useDingbats=FALSE,
  units = c("mm"),
  width = 183 ,
  height = 183,
  limitsize = FALSE
)

plot_heatmap = 
  counts_scaled %>%
  keep_variable(.abundance = counts_scaled, top = 500) %>%
  as_tibble() %>%
  heatmap(
    .column = sample,
    .row = feature,
    .value = count_scaled,
    transform = log1p,
    palette_value = circlize::colorRamp2(c(-2, -1, 0, 1, 2), RColorBrewer::brewer.pal(5,"RdBu")) 
  ) %>%
  add_tile(type) 


plot_heatmap %>% save_pdf("heatmap_variable_genes.pdf")

counts_de = 
  counts_scaled %>%
  test_differential_abundance(
    ~0+type, 
    method="edger_robust_likelihood_ratio", 
    .contrasts = "typeKnock_out-typeWild_type", 
    test_above_log2_fold_change = 1,
    omit_contrast_in_colnames=TRUE
  )


library(ppcseq)

counts_de_out = 
  counts_de %>%
  mutate(type = case_when(
    type == "Wild_type" ~ "1_Wild_type",
    type == "Knock_out" ~ "2_Knock_out"
  )) %>%
  as_tibble() %>%
  mutate(check = FDR  < 0.05) %>%
  filter(FDR %>% is.na %>% `!`) %>%
  mutate(count_scaled = as.integer(count_scaled)) %>%
  identify_outliers(
    ~type, 
    sample, 
    feature , 
    count_scaled, 
    .significance =PValue , 
    .do_check = check, 
    approximate_posterior_inference = FALSE
  )

counts_de_out %>% saveRDS("counts_de_out.rds")

counts_de_out %>% 
  filter(`tot deleterious outliers` > 0) %>%
  ppcseq::plot_credible_intervals() %>% 
  left_join(counts_de %>% distinct(feature)) %>%
  mutate(plot = map2(plot, feature, ~ .x + ggtitle(.y))) %>%
  pull(plot) %>%
  wrap_plots()

ggsave(
  "outliers.pdf",
  useDingbats=FALSE,
  units = c("mm"),
  width = 483 ,
  height = 283,
  limitsize = FALSE
)

# DE plotting
counts_de %>%
  filter(FDR < 0.05) %>%
  tidySummarizedExperiment::count(logFC > 0)

# A tibble: 2 x 2
# `logFC > 0`     n
# <lgl>       <int>
# FALSE         522
# TRUE          468

counts_de %>%
  as_tibble() %>%
  nanny::subset(feature) %>%
  arrange(PValue) %>%
  describe_transcript(feature) %>%
  dplyr::select(feature, description, logFC, FDR) %>%
  write_csv("differentially_abundant_gene_features_wild_minus_ko.csv")

# Heatmap
plot_heatmap_de = 
  counts_de %>%
  filter(FDR<0.05)%>% #  & abs(logFC) >=2
  as_tibble() %>%
  heatmap(
    .column = sample,
    .row = feature,
    .value = count_scaled,
    transform = log1p
  ) %>%
  add_tile(type) 

plot_heatmap_de %>% save_pdf("heatmap_differential_genes.pdf", width = 183, height = 183, units = "mm")

topgenes_symbols <-
  counts_de %>%
  as_tibble() %>%
  nanny::subset(feature)  %>%
  filter(FDR<0.05) %>%
  arrange(desc(abs(logFC)))  %>%
  pull(feature)

counts_de %>%
  as_tibble() %>%
  nanny::subset(feature) %>%
  
  # Subset data
  mutate(significant = case_when(
    FDR<0.05 & logFC < 0 ~ "Down-regulated", 
    FDR<0.05 & logFC > 0 ~ "Up-regulated",
    TRUE ~ ""
  )) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  # Plot
  ggplot(aes(x = logFC, y = PValue, label=feature)) +
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color="#8c8c8c") +
  geom_hline(yintercept = mean(c(1.24E-4, 3.74E-4)), linetype="dashed", color="#8c8c8c") +
  geom_point(aes(color = significant,  alpha=significant)) +
  ggrepel::geom_text_repel() +
  
  
  # Custom scales
  xlab("Log 2 fold change") +
  theme_bw() +
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values=c("black", "#1c4966", "#e11f28")) +
  scale_size_discrete(range = c(0, 2))

ggsave(
  "volcano_plot.pdf",
  useDingbats=FALSE,
  units = c("mm"),
  width = 183 ,
  height = 183,
  limitsize = FALSE
)


# Gene enrichment

library(enrichplot)


counts_gene_rank = 
  counts_de %>%
  symbol_to_entrez(feature) %>%
  filter(PValue   %>% is.na %>% `!`) %>%
  test_gene_rank(
    .entrez = entrez,
    .arrange_desc = logFC ,
    species="Homo sapiens",
    gene_collections = c("H", "C2", "C5")
  )


counts_gene_rank %>% saveRDS("counts_gene_rank.rds")

# Top GSEA
counts_gene_rank %>% 
  filter(gs_cat == "C2") %>%
  unnest(test) %>%
  slice(1:30) %>%
  mutate(plot = pmap(
    list(fit, ID, idx_for_plotting, p.adjust), 
    ~ enrichplot::gseaplot2(
      ..1, 
      geneSetID = ..3, 
      title = sprintf("%s \nadj pvalue %s", ..2, round(..4, 2)),
      base_size = 6
    ) 
  )) %>%
  pull(plot)%>%
  wrap_plots()



library(EGSEA)
counts_gene_enrichment = 
  counts_de %>%
  symbol_to_entrez(feature) %>%
  test_gene_enrichment(
    ~0+type, 
    .entrez=entrez, 
    .contrasts = "typeKnock_out-typeWild_type",
    species="human",
    gene_collections = c("H", "C2", "C5")
  )

counts_gene_enrichment %>% saveRDS("counts_gene_enrichment.rds")

my_pathways = 
  counts_gene_enrichment %>%
  #filter(data_base == "c2") %>%
  slice(1:30) %>%
  pull(pathway)

counts_gene_rank %>%
  mutate(test = map(test, ~ .x %>% filter(ID %in% my_pathways))) %>%
  unnest(test) %>%
  mutate(plot = pmap(
    list(fit, ID, idx_for_plotting, p.adjust), 
    ~ enrichplot::gseaplot2(
      ..1, 
      geneSetID = ..3, 
      title = sprintf("%s \nadj pvalue %s", ..2, round(..4, 2)),
      base_size = 6
    ) 
  )) %>%
  pull(plot)%>%
  wrap_plots()

# Plot Vijay list

numeric_list = c(210, 220, 224, 1,                 48,                 69    ,             99,                 48,                 106,                 93,                 7,                 51               ,  46,                 71,                 86,                 367,                 350,                 521,                 549, 179)

p_numeric_list = 
  counts_gene_rank %>%
  filter(gs_cat == "C2") %>%
  unnest(test) %>%
  filter(idx_for_plotting %in% numeric_list) %>%
  mutate(plot = pmap(
    list(fit, ID, idx_for_plotting, p.adjust), 
    ~ enrichplot::gseaplot2(
      ..1, 
      geneSetID = ..3, 
      title = sprintf("%s \nadj pvalue %s", ..2, round(..4, 2)),
      base_size = 6
    ) 
  )) %>%
  pull(plot)%>%
  wrap_plots()

ggsave(
  "p_numeric_list.pdf",
  plot = p_numeric_list,
  useDingbats=FALSE,
  units = c("mm"),
  width = 183 ,
  height = 183,
  limitsize = FALSE
)

# Overall enrichment

counts_gene_rank %>%
  filter(gs_cat == "C2") %>%
  unnest(test) %>%
  mutate(label = if_else(idx_for_plotting %in% numeric_list, ID %>% stringr::str_sub(0, 20), "")) %>%
  filter(p.adjust < 0.05) %>%
  ggplot(aes(forcats::fct_reorder(ID, enrichmentScore), enrichmentScore, label=label )) +
  geom_hline(yintercept = 0, color="#c8c8c8", linetype="dashed") +
  geom_point(aes(size=Count, color=p.adjust)) +
  ggrepel::geom_text_repel(max.overlaps = 100, size = 2) +
  facet_grid(~ gs_cat , scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + # element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_distiller( trans = "log10_reverse", palette = "Spectral") +
  #scale_y_discrete(labels = label_func) +
  ylab("Enrichment score") + 
  xlab(NULL) +
  scale_size(range=c(1, 3)) +
  ggExtra::removeGrid() +
  guides(fill = guide_legend(override.aes = list(size = 1), nrow = 1 ) ) +
  theme(legend.position = 'bottom', text = element_text(size = 8))
  

ggsave(
  "sets_overall.pdf",
  useDingbats=FALSE,
  units = c("mm"),
  width = 183 ,
  height = 100,
  limitsize = FALSE
)
