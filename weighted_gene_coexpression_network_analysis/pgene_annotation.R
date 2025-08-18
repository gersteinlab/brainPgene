GOenrichment_hs <- readRDS('WGCNA/human_wgcna_BP_enrichment.rds')
GOenrichment_hs_subset <-GOenrichment_hs %>% filter(grepl('brain|neuron|axon|synaptic|synapse|learning|memory', Description, ignore.case = T)) 
hs_modules <-  readRDS('WGCNA/human_res.rds')  %>% bind_rows()
pgene_expr <- read.table('quantification/BrainSpan/pgene_expression_BrainSpan.tsv', header = T, row.names = 1, sep = '\t')
hs_modules_pgene <- hs_modules %>% mutate(type = ifelse(hs_modules[['gene']] %in% rownames(pgene_expr), 'Pseudogene', 'Protein-coding')) %>% filter(type == 'Pseudogene') %>% merge(GOenrichment_hs_subset, by = 'module') %>% dplyr::select(gene, ID)
GO_terms <- unique(hs_modules_pgene$ID); pgene <- unique(hs_modules_pgene$gene)
df <- data.frame(matrix(ncol = length(pgene), nrow = length(GO_terms)))
colnames(df) <- pgene
rownames(df) <- sapply(pgene,function(x) hs_modules_pgene %>% filter(gene == x) %>% pull(ID)) %>% unlist() %>% unique()
for (i in 1:ncol(df)) {
  for (j in 1:nrow(df)){
    gene_name <- colnames(df)[i]
    GO_name <- rownames(df)[j]
    df[j,i] <- ifelse(hs_modules_pgene %>% filter(ID == GO_name, gene == gene_name) %>% nrow() != 0, 1, 0)
  }
}
saveRDS(df, 'WGCNA/pgene_anno.rds')
library(ComplexHeatmap)
df <- readRDS('WGCNA/pgene_anno.rds')

pdf('WGCNA/pgene_anno3.pdf', width = 6, height = 6)
Heatmap(df, col = c('white', 'black'), column_names_side = 'top', row_names_side = 'right',cluster_columns = T, show_row_names = T, show_column_names = F, show_row_dend = T, cluster_rows = T, row_title = 'GO terms (BP)', column_title = 'Pseudogenes', row_names_gp = grid::gpar(fontsize = 4))
dev.off()

## cluster_1

term_list_1 <- c('GO:0150104', 'GO:0010975', 'GO:0030900', 'GO:0106027', 'GO:0050803', 'GO:0007416', 'GO:0007409', 'GO:0021954', 
                 'GO:0050804','GO:0050808', 'GO:0050807', 'GO:0099560', 'GO:0007612', 'GO:0042551', 'GO:0008038')

count_1 <- c()

for (term in term_list_1) {
  count_1[term] <- df %>% filter(rownames(df) == term) %>% rowSums()
  
}

term_list_2 <- c('GO:0099504', 'GO:0035249', 'GO:0035418', 'GO:0099558', 'GO:0019228')
count_2 <- c()

for (term in term_list_2) {
  count_2[term] <- df %>% filter(rownames(df) == term) %>% rowSums()
}

library(ggplot2)
library(forcats)
p1 <- data.frame(id = rowSums(df) %>% names(), count = rowSums(df)) %>% slice_max(count, n = 25) %>%
  ggplot(aes(x = fct_reorder(id, -count), y = count)) + geom_bar(stat = 'identity', fill = '#4f6d7a') + 
  labs(x = 'Biological Process', y = 'Number of pseudogenes annotated') +
  geom_hline(yintercept = mean(rowSums(df)), linetype = "dashed", color = "#fde4cf", linewidth = 0.5) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(axis.text = element_text(colour = 'black'))

ggsave('WGCNA/num_pgene_annotated.pdf', p1, width = 6, height = 4)


         