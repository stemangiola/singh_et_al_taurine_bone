# module load cutadapt; module load trimgalore; ls *fastq.gz | awk '{split($0,a,"_"); print a[3]}' | uniq | parallel --eta -j1 trim_galore --paired SO_9588_{}_R1.fastq.gz SO_9588_{}_R2.fastq.gz --fastqc > trim_galore.log
# STAR --runThreadN 15 --runMode genomeGenerate --genomeDir  ~/third_party_sofware/hg38/karyotypic/ --genomeFastaFiles ~/third_party_sofware/hg38/karyotypic/Homo_sapiens_assembly38.fasta --sjdbGTFfile ~/third_party_sofware/hg38/karyotypic/hg38_ucsc.gtf --sjdbOverhang 100
# module load STAR; for i in $(ls *_val_*fq.gz | rev | cut -c 16- | rev | uniq); do mkdir -p alignment_hg38/$i; cd alignment_hg38/$i; STAR --genomeDir ~/third_party_sofware/hg38/karyotypic/ --readFilesIn ../../$i'_R1_val_1.fq.gz' ../../$i'_R2_val_2.fq.gz' --readFilesCommand zcat --genomeLoad LoadAndKeep --runThreadN 42 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 40000000000 --outReadsUnmapped Fastx; cd ../../; done

setwd("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/third_party_analyses/bhupinder_bone")

library(tidyverse)
library(tidySummarizedExperiment)
library(tidybulk)
library(tidyHeatmap)
library(patchwork)
library(future)
plan(multisession, workers = 15)

# # Counting
# raw_counts =
#   dir(
#     "download.genotypic.co.in/SO_9588_rawdata/alignment_hg38/",
#     recursive = T,
#     full.names = T,
#     pattern = ".bam"
#   ) %>%
#   tidybulk_SAM_BAM(genome = "hg38", isPairedEnd = TRUE)
# 
# raw_counts %>% saveRDS("raw_counts.rds", compress = "gzip")

raw_counts = readRDS("raw_counts.rds")


counts_scaled = 
  raw_counts %>%
  mutate(sample_id = sample) %>%
  extract(sample_id, "sample_id",  regex=".*SO_9588_(.+)/.+") %>%
  
  # Add factor of interest
  inner_join(
    tibble(
      sample_id = c("W1", "W2", "W3", "Ko1", "Ko2", "Ko3"),
      type = c(rep("Wild_type", 3), rep("Knock_out", 3))
    ),
    by="sample_id"
  ) %>%
  identify_abundant(factor_of_interest=type) %>%
  scale_abundance() %>%
  reduce_dimensions(method = "PCA")

counts_scaled %>% saveRDS("counts_scaled.rds")

# Density plot
counts_scaled %>%
  keep_abundant(factor_of_interest=type) %>%
  ggplot(aes(count_scaled+1, group=sample_id, color=type)) +
  geom_density() +
  scale_x_log10() +
  scale_color_brewer(palette = "Set1") +
  theme_bw()

# PCA
# ADD PIVOT SAMPLE 
plot_pca = 
  counts_scaled %>%
  as_tibble() %>%
  nanny::subset(sample) %>%
  ggplot(aes(PC1,PC2, label=sample_id)) +
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

# Heatmap
plot_heatmap = 
  counts_scaled %>%
  keep_variable(.abundance = counts_scaled, top = 500) %>%
  as_tibble() %>%
  heatmap(
    .column = sample_id,
    .row = symbol,
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
  identify_outliers(
    ~type, 
    sample_id, 
    transcript , 
    count, 
    .significance =PValue , 
    .do_check = check, 
    approximate_posterior_inference = FALSE
  )

counts_de_out %>% 
  ppcseq::plot_credible_intervals() %>% 
  left_join(counts_de %>% distinct(symbol, transcript)) %>%
  mutate(plot = map2(plot, symbol, ~ .x + ggtitle(.y))) %>%
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


counts_de %>%
  as_tibble() %>%
  nanny::subset(transcript) %>%
  arrange(PValue) %>%
  describe_transcript(symbol) %>%
  dplyr::select(symbol, description, logFC, FDR) %>%
  write_csv("differentially_abundant_gene_transcripts_wild_minus_ko.csv")

# Heatmap
plot_heatmap_de = 
  counts_de %>%
  filter(FDR<0.05)%>% #  & abs(logFC) >=2
  as_tibble() %>%
  heatmap(
    .column = sample_id,
    .row = symbol,
    .value = count_scaled,
    transform = log1p
  ) %>%
  add_tile(type) 

plot_heatmap_de %>% save_pdf("heatmap_differential_genes.pdf", width = 183, height = 183, units = "mm")

topgenes_symbols <-
  counts_de %>%
  as_tibble() %>%
  nanny::subset(transcript)  %>%
  filter(FDR<0.05) %>%
  arrange(desc(abs(logFC)))  %>%
  pull(symbol)

counts_de %>%
  as_tibble() %>%
  nanny::subset(transcript) %>%
  
  # Subset data
  mutate(significant = FDR<0.05) %>%
  mutate(symbol = ifelse(symbol %in% topgenes_symbols, as.character(symbol), "")) %>%
  
  # Plot
  ggplot(aes(x = logFC, y = PValue, label=symbol)) +
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color="#8c8c8c") +
  geom_hline(yintercept = mean(c(1.24E-4, 3.74E-4)), linetype="dashed", color="#8c8c8c") +
  geom_point(aes(color = significant,  alpha=significant)) +
  ggrepel::geom_text_repel() +

  
  # Custom scales
  xlab("Log 2 fold change") +
  theme_bw() +
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values=c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2))

ggsave(
  "volcano_plot.pdf",
  useDingbats=FALSE,
  units = c("mm"),
  width = 183 ,
  height = 183,
  limitsize = FALSE
)


library(EGSEA)
counts_gene_enrichment = 
  counts_de %>%
  test_gene_enrichment(
    ~0+type, 
    .entrez=transcript, 
    .contrasts = "typeKnock_out-typeWild_type",
    species="human"
  )

counts_gene_enrichment %>% saveRDS("counts_gene_enrichment.rds")

library(enrichplot)


counts_gene_rank = 
  counts_scaled %>%
  test_differential_abundance(
    ~0+type, 
    method="limma_voom", 
    .contrasts = "typeKnock_out-typeWild_type", 
    omit_contrast_in_colnames=TRUE
  ) %>%
  filter(P.Value  %>% is.na %>% `!`) %>%
  mutate(score=logFC) %>% 
  test_gene_rank(
    .entrez = transcript,
    .arrange_desc = logFC ,
    species="Homo sapiens",gene_set = c("H", "C1", "C2", "C3", "C5", "C6", "C7", "C8")
  )

counts_gene_rank %>% saveRDS("counts_gene_rank.rds")

# Vijay sets
my_pathways = 
  counts_gene_enrichment %>%
  slice(c(7, 8, 15, 20, 22, 32, 34, 40, 43, 52, 54, 61, 75, 79, 82, 83, 91)) %>%
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

ggsave(
  "vijay_sets.pdf",
  useDingbats=FALSE,
  units = c("mm"),
  width = 250 ,
  height = 250,
  limitsize = FALSE
)

# Top 16 EGSEA
my_pathways = 
  counts_gene_enrichment %>%
  slice(1:16) %>%
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

ggsave(
  "top_16_EGSEA_sets.pdf",
  useDingbats=FALSE,
  units = c("mm"),
  width = 250 ,
  height = 250,
  limitsize = FALSE
)


# Top 16 GSEA
my_pathways = 
  counts_gene_rank %>%
  filter(gs_cat != "C4") %>%
  select(test) %>%
  unnest(test) %>%
  arrange(p.adjust) %>%
  slice(1:16) %>%
  pull(ID)

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


ggsave(
  "top_16_GSEA_sets.pdf",
  useDingbats=FALSE,
  units = c("mm"),
  width = 250 ,
  height = 250,
  limitsize = FALSE
)

library(enrichplot)



counts_gene_rank %>%
  mutate(test = map(test, ~ .x %>% filter(ID %in% (counts_gene_enrichment %>%
                                            slice(1:100) %>%
                                            pull(pathway))))) %>%
  unnest(test) %>%

  ggplot(aes(forcats::fct_reorder(ID, enrichmentScore), enrichmentScore,size=Count, color=p.adjust )) +
  geom_hline(yintercept = 0, color="#c8c8c8", linetype="dashed") +
  geom_point() +
  facet_grid(~ gs_cat , scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_distiller( trans = "log10_reverse", palette = "Spectral") +
  #scale_y_discrete(labels = label_func) +
  ylab("Enrichment score") + 
  xlab(NULL) +
  scale_size(range=c(3, 8))

ggsave(
  "sets_overall.pdf",
  useDingbats=FALSE,
  units = c("mm"),
  width = 600 ,
  height = 300,
  limitsize = FALSE
)


# Volcano plot non robust

counts_de_limma_voom = 
  counts_scaled %>%
  test_differential_abundance(
    ~0+type, 
    method="limma_voom", 
    .contrasts = "typeKnock_out-typeWild_type", 
    omit_contrast_in_colnames=TRUE
  )

topgenes_symbols <-
  counts_de_limma_voom %>%
  as_tibble() %>%
  nanny::subset(transcript)  %>%
  filter(adj.P.Val<0.05 & abs(logFC) > 1) %>%
  arrange(desc(abs(logFC)))  %>%
  pull(symbol)

counts_de_limma_voom %>%
  as_tibble() %>%
  nanny::subset(transcript) %>%
  
  # Subset data
  mutate(significant = adj.P.Val<0.05 & abs(logFC) > 1) %>%
  mutate(symbol = ifelse(symbol %in% topgenes_symbols, as.character(symbol), "")) %>%
  
  # Plot
  ggplot(aes(x = logFC, y = P.Value, label=symbol)) +
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color="#8c8c8c") +
  geom_hline(yintercept = mean(c(1.24E-4, 3.74E-4)), linetype="dashed", color="#8c8c8c") +
  geom_point(aes(color = significant,  alpha=significant)) +
  ggrepel::geom_text_repel() +
  
  
  # Custom scales
  xlab("Log 2 fold change") +
  theme_bw() +
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values=c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2))

ggsave(
  "volcano_plot_limma_voom.pdf",
  useDingbats=FALSE,
  units = c("mm"),
  width = 183 ,
  height = 183,
  limitsize = FALSE
)
