plot_triple_venn <- function(data, type, padj_threshold = 0.05, log2FC_threshold = 0){
  ASD <- data[1] %>% as.data.frame() %>% filter(type == !!type, padj < padj_threshold, abs(log2FoldChange) > log2FC_threshold) %>% dplyr::select(gene) %>% pull()
  BD <- data[2] %>% as.data.frame() %>% filter(type == !!type, padj < padj_threshold, abs(log2FoldChange) > log2FC_threshold) %>% dplyr::select(gene) %>% pull()
  SCZ <- data[3] %>% as.data.frame() %>% filter(type == !!type, padj < padj_threshold, abs(log2FoldChange) > log2FC_threshold) %>% dplyr::select(gene) %>% pull()
  model <- euler(c('A' = length(ASD), 'B' = length(BD), 'C' = length(SCZ),
                   'A&B' = Reduce(intersect, list(ASD, BD)) %>% length(), 
                   'A&C' = Reduce(intersect, list(ASD, SCZ)) %>% length(), 
                   'B&C' = Reduce(intersect, list(BD, SCZ)) %>% length(), 
                   'A&B&C' = Reduce(intersect, list(ASD, BD, SCZ)) %>% length()))
  plot(model, quantities = TRUE, labels = c('ASD', 'BD', 'SCZ'), fills = list(fill = pal_npg("nrc")(3), alpha = 0.8), lty = 1, main = str_to_title(type))
}

DE_pc_count <- lapply(seq(0, 1, 0.001), function(x) rbind(ASD, BD, SCZ) %>% dplyr::filter(abs(log2FoldChange) > x, padj < 0.05) %>% summarise(n=n()) %>% pull() )

pdf('./DEG/Capstone/pc/log2FC_threshold.pdf', width = 5, height = 5)
plot(x = seq(0, 1, 0.001), y = DE_pc_count, xlab = 'log2FC', ylab = 'Number of DE genes', main = 'Protein-coding gene', pch = 16)
abline(v = 0.2, col = 'red')
dev.off()
pdf('./DEG/Capstone/pgene/log2FC_threshold.pdf', width = 5, height = 5)
DE_pgene_count <- lapply(seq(0, 1, 0.001), function(x) rbind(ASD_p, BD_p, SCZ_p) %>% dplyr::filter(abs(log2FoldChange) > x, padj < 0.05) %>% summarise(n=n()) %>% pull())
plot(x = seq(0, 1, 0.001), y = DE_pgene_count, xlab = 'log2FC', ylab = 'Number of DE genes', main = 'Pseudogene', pch = 16)
dev.off()

pdf('DEG/Capstone/pgene/pgene_upregulated.pdf', width = 5, height = 5)
plot_triple_venn(list(ASD_p, BD_p, SCZ_p), type =  'upregulated', log2FC_threshold = 0)
dev.off()
pdf('DEG/Capstone/pgene/pgene_downregulated.pdf', width = 5, height = 5)
plot_triple_venn(list(ASD_p, BD_p, SCZ_p), type =  'downregulated', log2FC_threshold = 0)
dev.off()

pdf('DEG/Capstone/pc/pc_upregulated.pdf', width = 5, height = 5)
plot_triple_venn(list(ASD, BD, SCZ), type =  'upregulated', log2FC_threshold = 0.2)
dev.off()
pdf('DEG/Capstone/pc/pc_downregulated.pdf', width = 5, height = 5)
plot_triple_venn(list(ASD, BD, SCZ), type =  'downregulated', log2FC_threshold = 0.2)
dev.off()