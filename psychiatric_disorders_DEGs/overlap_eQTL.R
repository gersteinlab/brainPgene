eqtl_df <- read.delim("DEG/eQTL/Brain_Frontal_Cortex_BA9.v7.signif_variant_gene_pairs.txt")
eqtl_df$gene_id <- sub("\\..*", "", eqtl_df$gene_id)

ASD_p <- read.table('DEG/Capstone/pgene/raw/ASD_DE.txt', sep = '\t', header = T)
BD_p <- read.table('DEG/Capstone/pgene/raw/BD_DE.txt', sep = '\t', header = T)
SCZ_p <- read.table('DEG/Capstone/pgene/raw/SCZ_DE.txt', sep = '\t', header = T)

ASD <- read.table('DEG/Capstone/pc/raw/ASD_DE.txt', sep = '\t', header = T)
BD <- read.table('DEG/Capstone/pc/raw/BD_DE.txt', sep = '\t', header = T)
SCZ <- read.table('DEG/Capstone/pc/raw/SCZ_DE.txt', sep = '\t', header = T)

deg_eqtl_overlap_pgene <- function(data, eqtl){
  overlap <- data %>% filter(padj < 0.05) %>%
    inner_join(eqtl_df, by = c("gene" = "gene_id")) %>%
    distinct(gene, disease, type) %>%
    group_by(disease, type) %>%
    summarise(eqtl_overlap_count = n()) %>%
    ungroup()
  total <- data %>% filter(padj < 0.05) %>%
    group_by(disease, type) %>%
    summarise(total_DEG = n()) %>%
    ungroup()
  summary <- left_join(total, overlap, by = c("disease", "type")) %>%
    mutate(eqtl_overlap_count = ifelse(is.na(eqtl_overlap_count), 0, eqtl_overlap_count),
           percent_overlap = round(eqtl_overlap_count / total_DEG * 100, 2))
  
  summary$type <- factor(summary$type, levels = c('upregulated', 'downregulated'))
  print(summary)
  p <- ggplot(summary, aes(x = disease, y = percent_overlap, fill = type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    ylab("Percentage of DEGs overlapping eGenes (%)") + xlab('') +
    scale_fill_manual(values = c("upregulated" = "#1f77b4", "downregulated" = "#ff7f0e")) +
    scale_x_discrete(labels =c('ASD', 'BD', 'SCZ')) + ylim(0,40) +
    theme_bw(base_size = 12) +
    theme(legend.title = element_blank(), axis.text = element_text(color = 'black')) +
    ggtitle("Pseudogene")
  return(p)
}

deg_eqtl_overlap_pc <- function(data, eqtl){
  overlap <- data %>% filter(abs(log2FoldChange) > 0.2, padj < 0.05) %>%
    inner_join(eqtl_df, by = c("gene" = "gene_id")) %>%
    distinct(gene, disease, type) %>%
    group_by(disease, type) %>%
    summarise(eqtl_overlap_count = n()) %>%
    ungroup()
  total <- data %>% filter(abs(log2FoldChange) > 0.2, padj < 0.05) %>%
    group_by(disease, type) %>%
    summarise(total_DEG = n()) %>%
    ungroup()
  summary <- left_join(total, overlap, by = c("disease", "type")) %>%
    mutate(eqtl_overlap_count = ifelse(is.na(eqtl_overlap_count), 0, eqtl_overlap_count),
           percent_overlap = round(eqtl_overlap_count / total_DEG * 100, 2))
  summary$type <- factor(summary$type, levels = c('upregulated', 'downregulated'))
  print(summary)
  p <- ggplot(summary, aes(x = disease, y = percent_overlap, fill = type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    ylab("Percentage of DEGs overlapping eGenes (%)") + xlab('') +
    scale_fill_manual(values = c("upregulated" = "#1f77b4", "downregulated" = "#ff7f0e")) +
    scale_x_discrete(labels =c('ASD', 'BD', 'SCZ')) + ylim(0,40) +
    theme_bw(base_size = 12) +
    theme(legend.title = element_blank(), axis.text = element_text(color = 'black')) +
    ggtitle("Protein-coding gene")
  return(p)
}

p1 <- deg_eqtl_overlap_pgene(rbind(ASD_p, BD_p, SCZ_p), eqtl_df)
p2 <- deg_eqtl_overlap_pc(rbind(ASD, BD, SCZ), eqtl_df)

p <- p1 + p2 + plot_layout(axis_titles = 'collect', guides = 'collect')
ggsave('./DEG/eQTL/perc_DEG_overlap_eQTL.pdf', width = 6, height = 4)
