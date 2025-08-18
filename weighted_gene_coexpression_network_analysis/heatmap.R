rm(list = ls())
GOenrichment_hs <- readRDS('WGCNA/human_wgcna_BP_enrichment.rds')
GOenrichment_hs_subset <-GOenrichment_hs %>% filter(grepl('brain|neuron|axon|synaptic|synapse|learning|memory', Description, ignore.case = T)) 
metadata <- sapply(unique(GOenrichment_hs_subset$ID), function(x) GOenrichment_hs_subset %>% filter(ID == x) %>% pull(region) %>% unique() %>% length(), USE.NAMES = T) %>% as.list() %>%  data.frame(check.names = F) %>% t() %>% data.frame(check.names = F) %>% add_rownames(var = "ID")
colnames(metadata)[2] <- 'tissue'
library(forcats)
metadata$ID <- fct_reorder(metadata$ID, metadata$tissue, .desc = TRUE)
GOenrichment_hs_subset <- merge(GOenrichment_hs_subset, metadata, by = 'ID')
library(dplyr)
GO_terms <- unique(GOenrichment_hs_subset$ID) ;  regions <- unique(GOenrichment_hs_subset$region)
wide_table <- data.frame(matrix(ncol = length(regions), nrow = length(GO_terms)))


colnames(wide_table) <- c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD')
rownames(wide_table) <- sapply(c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD'), 
                               function(x) GOenrichment_hs_subset %>% filter(region == x) %>% pull(ID)) %>% unlist() %>% unique()

for (i in 1:ncol(wide_table)) {
  for (j in 1:nrow(wide_table)){
    region_name <- colnames(wide_table)[i]
    GO_name <- rownames(wide_table)[j]
    wide_table[j,i] <- ifelse(GOenrichment_hs_subset %>% filter(ID == GO_name, region == region_name) %>% nrow() != 0,
                              -log10(GOenrichment_hs_subset %>% filter(ID == GO_name, region == region_name) %>% pull(p.adjust) %>% min()), NA)
  }
}
library(ComplexHeatmap)
library(circlize)
index <- rowSums(!is.na(wide_table)) == 1
col <- colorRamp2(seq(min(wide_table, na.rm = T), max(wide_table, na.rm = T), length = 3), c("blue", "#EEEEEE", "red"))
pdf('./WGCNA/human_wgcna_BP_enrichment_subset.pdf', width = 8, height = 8)
Heatmap(wide_table[index,], column_names_side = 'top', cluster_columns = F, show_row_names = F, show_row_dend = F, cluster_rows = F, column_names_rot = 0, column_names_centered = T, col = col, na_col = 'gray90', name = '-log10(p.adj)')
Heatmap(wide_table[!index,], column_names_side = 'top', cluster_columns = F, show_row_names = F, show_row_dend = F, cluster_rows = F, column_names_rot = 0, column_names_centered = T, col = col, na_col = 'gray90', name = '-log10(p.adj)')
dev.off()
