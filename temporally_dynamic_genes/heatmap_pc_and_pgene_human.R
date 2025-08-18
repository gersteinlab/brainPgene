pgene <- readRDS('temporal_dynamic_expression/human/pgene_temporallyDynamicGene.rds')
pgene <- sapply(c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD'), 
                function(x) pgene[[x]] %>% rownames(), simplify = F, USE.NAMES = T)
tissues <- names(pgene)
TDGs <- unlist(pgene, use.names = F) %>% unique()
binaryDf <- data.frame(matrix(ncol = length(tissues), nrow = length(TDGs)), row.names = TDGs)
colnames(binaryDf) <- tissues
for (i in 1:ncol(binaryDf)) {binaryDf[,i] <- rownames(binaryDf) %in% pgene[[colnames(binaryDf)[i]]]}
binaryDf_mat <- apply(binaryDf, c(1,2), as.integer) %>% as.matrix()
library(ComplexHeatmap)
pdf('./temporal_dynamic_expression/human_pgene_heatmap.pdf', width = 10, height = 10)
Heatmap(binaryDf_mat[rowSums(binaryDf_mat) == 1,], col = c('black', 'red'), column_names_side = 'top', cluster_columns = F, show_row_names = F, show_row_dend = F, cluster_rows = F, column_names_rot = 0, column_names_centered = T)
Heatmap(binaryDf_mat[rowSums(binaryDf_mat) == 2,], col = c('black', 'red'), column_names_side = 'top', cluster_columns = F, show_row_names = F, show_row_dend = F, cluster_rows = F, column_names_rot = 0, column_names_centered = T)
Heatmap(binaryDf_mat[rowSums(binaryDf_mat) == 3,], col = c('black', 'red'), column_names_side = 'top', cluster_columns = F, show_row_names = F, show_row_dend = F, cluster_rows = F, column_names_rot = 0, column_names_centered = T)
Heatmap(binaryDf_mat[rowSums(binaryDf_mat) >= 4 & rowSums(binaryDf_mat) <=8 ,], col = c('black', 'red'), column_names_side = 'top', cluster_columns = F, show_row_names = F, show_row_dend = F, cluster_rows = F, column_names_rot = 0, column_names_centered = T)
Heatmap(binaryDf_mat[rowSums(binaryDf_mat) >= 9 ,], col = c('black', 'red'), column_names_side = 'top', cluster_columns = F, show_row_names = F, show_row_dend = F, cluster_rows = F, column_names_rot = 0, column_names_centered = T)
dev.off()

pgene_expression <- read.table('quantification/BrainSpan/pgene_expression_BrainSpan.tsv', sep = '\t', header = T, row.names = 1)
pgene_TDG <- list('Non-TD' = base::setdiff(rownames(pgene_expression), rownames(binaryDf_mat)),
                  'TD in one region' = binaryDf_mat[rowSums(binaryDf_mat) == 1,]  %>% rownames(),
                  'TD in multiple regions' = base::setdiff(rownames(binaryDf_mat), binaryDf_mat[rowSums(binaryDf_mat) == 1,]  %>% rownames()))
pgene_TDG_df <- lapply(pgene_TDG, as.data.frame, stringsAsFactors = FALSE) %>% bind_rows(.id = 'Type')
colnames(pgene_TDG_df)[2] <- 'row'
write.table(pgene_TDG_df, './temporal_dynamic_expression/human/pgene_TDG_nonTDG.txt', sep = '\t', quote = F, row.names = F)

pc <- readRDS('temporal_dynamic_expression/human/protein_coding_temporallyDynamicGene.rds')
pc <- sapply(c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD'), 
             function(x) pc[[x]] %>% rownames(), simplify = F, USE.NAMES = T)
tissues <- names(pc)
TDGs <- unlist(pc, use.names = F) %>% unique()
binaryDf <- data.frame(matrix(ncol = length(tissues), nrow = length(TDGs)), row.names = TDGs)
colnames(binaryDf) <- tissues
for (i in 1:ncol(binaryDf)) {binaryDf[,i] <- rownames(binaryDf) %in% pc[[colnames(binaryDf)[i]]]}
binaryDf_mat <- apply(binaryDf, c(1,2), as.integer) %>% as.matrix()
pdf('temporal_dynamic_expression/human_protein_coding_heatmap.pdf', width = 10, height = 10)
Heatmap(binaryDf_mat[rowSums(binaryDf_mat) == 1,], col = c('black', 'red'), column_names_side = 'top', cluster_columns = F, show_row_names = F, show_row_dend = F, cluster_rows = F, column_names_rot = 0, column_names_centered = T, use_raster = F)
Heatmap(binaryDf_mat[rowSums(binaryDf_mat) == 2,], col = c('black', 'red'), column_names_side = 'top', cluster_columns = F, show_row_names = F, show_row_dend = F, cluster_rows = F, column_names_rot = 0, column_names_centered = T, use_raster = F)
Heatmap(binaryDf_mat[rowSums(binaryDf_mat) == 3,], col = c('black', 'red'), column_names_side = 'top', cluster_columns = F, show_row_names = F, show_row_dend = F, cluster_rows = F, column_names_rot = 0, column_names_centered = T, use_raster = F)
Heatmap(binaryDf_mat[rowSums(binaryDf_mat) >= 4 & rowSums(binaryDf_mat) <=8 ,], col = c('black', 'red'), column_names_side = 'top', cluster_columns = F, show_row_names = F, show_row_dend = F, cluster_rows = F, column_names_rot = 0, column_names_centered = T, use_raster = F)
Heatmap(binaryDf_mat[rowSums(binaryDf_mat) >= 9 ,], col = c('black', 'red'), column_names_side = 'top', cluster_columns = F, show_row_names = F, show_row_dend = F, cluster_rows = F, column_names_rot = 0, column_names_centered = T, use_raster = F)
dev.off()
## TDG types
pc_expression <- read.table('quantification/BrainSpan/pc_expression_BrainSpan.tsv', sep = '\t', header = T, row.names = 1)
pc_TDG <- list('Non-TD' = base::setdiff(rownames(pc_expression), rownames(binaryDf_mat)),
               'TD in one region' = binaryDf_mat[rowSums(binaryDf_mat) == 1,]  %>% rownames(),
               'TD in multiple regions' = base::setdiff(rownames(binaryDf_mat), binaryDf_mat[rowSums(binaryDf_mat) == 1,]  %>% rownames()))
pc_TDG_df <- lapply(pc_TDG, as.data.frame, stringsAsFactors = FALSE) %>% bind_rows(.id = 'Type')
colnames(pc_TDG_df)[2] <- 'row'
write.table(pc_TDG_df, './temporal_dynamic_expression/human/pc_TDG_nonTDG.txt', sep = '\t', quote = F, row.names = F)