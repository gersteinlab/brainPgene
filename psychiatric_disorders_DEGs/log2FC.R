setwd('DEG_across_disorders/')
ASD_p <- read.table('DEG/Capstone/pgene/raw/ASD_DE.txt', sep = '\t', header = T)
BD_p <- read.table('DEG/Capstone/pgene/raw/BD_DE.txt', sep = '\t', header = T)
SCZ_p <- read.table('DEG/Capstone/pgene/raw/SCZ_DE.txt', sep = '\t', header = T)
ASD <- read.table('DEG/Capstone/pc/raw/ASD_DE.txt', sep = '\t', header = T)
BD <- read.table('DEG/Capstone/pc/raw/BD_DE.txt', sep = '\t', header = T)
SCZ <- read.table('DEG/Capstone/pc/raw/SCZ_DE.txt', sep = '\t', header = T)

plot_log2FC_comparison <- function(data, padj_threshold = 0.05, log2FC_threshold = 0){
  ASD <- data[1] %>% as.data.frame() %>% filter(padj < padj_threshold, abs(log2FoldChange) > log2FC_threshold) %>% dplyr::mutate(Group = 'Protein-coding Gene')
  BD <- data[2] %>% as.data.frame() %>% filter(padj < padj_threshold, abs(log2FoldChange) > log2FC_threshold)  %>% dplyr::mutate(Group = 'Protein-coding Gene')
  SCZ <- data[3] %>% as.data.frame() %>% filter(padj < padj_threshold, abs(log2FoldChange) > log2FC_threshold) %>% dplyr::mutate(Group = 'Protein-coding Gene')
  ASD_p <-  data[4] %>% as.data.frame() %>% filter(padj < padj_threshold, abs(log2FoldChange) > log2FC_threshold)  %>% dplyr::mutate(Group = 'Pseudogene')
  BD_p <- data[5] %>% as.data.frame() %>% filter(padj < padj_threshold, abs(log2FoldChange) > log2FC_threshold)  %>% dplyr::mutate(Group = 'Pseudogene')
  SCZ_p <- data[6] %>% as.data.frame() %>% filter(padj < padj_threshold, abs(log2FoldChange) > log2FC_threshold)  %>% dplyr::mutate(Group = 'Pseudogene')
  disorder <- rbind(ASD, BD, SCZ, ASD_p, BD_p, SCZ_p)
  disorder$type <- str_to_title(disorder$type)
 
  disorder$type <- factor(disorder$type, levels = c('Upregulated', 'Downregulated'))
  disorder$Group <- factor(disorder$Group, levels = c('Pseudogene', 'Protein-coding Gene'))
#  tbl <- disorder %>% group_by(disease, type) %>% mutate(var = abs(log2FoldChange)) %>% rstatix::wilcox_test(var~Group) %>% rstatix::adjust_pvalue(method = 'fdr')
#  write.table(tbl, './DEG/Capstone/wilcox_test.txt', quote = F, row.names = F, sep = '\t')
  
  p <- disorder %>% ggplot(aes(x = disease, y = abs(log2FoldChange), fill = Group)) + geom_violin(trim = F) + geom_boxplot(width = 0.3, position = position_dodge(width = 0.9)) + 
    facet_grid(~type) + scale_y_log10() + xlab('Brain Disorder') + ylab(TeX('$\\left|\\log_{2}FoldChange\\right|$')) +
    theme_bw(base_size = 12) + theme(legend.position="bottom") + theme(legend.title=element_blank(), axis.text = element_text(colour = 'black')) + 
    scale_fill_manual(values = c('Pseudogene' = '#E64B35CC', 'Protein-coding Gene' = '#3C5488CC'), breaks = c('Pseudogene', 'Protein-coding Gene')) + 
    scale_x_discrete(labels=c('ASD','BD','SCZ')) +
    stat_compare_means(
      aes(group = Group), 
      label = "p.signif",  
      method = "wilcox.test",   
      comparisons = NULL,
      hide.ns = TRUE
    )
  return(p)
}

log2FC <- plot_log2FC_comparison(data = list(ASD, BD, SCZ, ASD_p, BD_p, SCZ_p))
ggsave('./DEG/Capstone/log2FC_comparison.pdf', log2FC, width = 10, height = 8)
