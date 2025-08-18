hs_modules <- readRDS('WGCNA/human_res.rds')
pgene_expr <- read.table('quantification/BrainSpan/pgene_expression_BrainSpan.tsv', header = T, row.names = 1, sep = '\t')

GOenrichment_module <- function(module_list, pgene_expr, x, orgDB = org.Hs.eg.db ){
  modules <- module_list[[x]] %>% mutate(type = ifelse(module_list[[x]][['gene']] %in% rownames(pgene_expr), 'Pseudogene', 'Protein-coding')) %>%
    group_by(module, type) %>% summarise(count = n()) %>% mutate(proportion= prop.table(count), region = x) %>% 
    filter(type == 'Protein-coding', !grepl("_0", module), proportion != 1, count > 30) %>% as.data.frame() %>% pull(module)
  df <-  module_list[[x]] %>% mutate(type = ifelse(module_list[[x]][['gene']] %in% rownames(pgene_expr), 'Pseudogene', 'Protein-coding')) %>% filter(module %in% modules)
  res <- data.frame()
  for (i in levels(factor(df$module))) {
    print(i)
    tmp <- subset(df, module == i)
    GOres <- enrichGO(tmp$gene, OrgDb = orgDB, keyType = 'ENSEMBL', ont = 'BP', pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 10)
    if(!is.null(GOres)){
      res <- rbind(res, simplify(GOres) %>% as.data.frame() %>% mutate(module = i, region = x))
    }
  }
  return(res)
}
GOenrichment <- lapply(c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD'), 
                       function(x) GOenrichment_module(hs_modules, pgene_expr, x))  %>% bind_rows()
saveRDS(GOenrichment, 'WGCNA/human_wgcna_BP_enrichment.rds')
write.table(GOenrichment, 'WGCNA/human_wgcna_BP_enrichment.txt', sep = '\t', quote = F, row.names = F)

# macaque
mm_modules <- readRDS('WGCNA/macaque_res.rds')
pgene_expr_mm <- read.table('quantification/NHP_macaque/pgene_expression_NHP_macaque.tsv', header = T, row.names = 1, sep = '\t')
GOenrichment_mm <- lapply(c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD'), 
                          function(x) GOenrichment_module(mm_modules, pgene_expr_mm, x))  %>% bind_rows()
saveRDS(GOenrichment_mm, 'WGCNA/macaque_wgcna_BP_enrichment.rds')
write.table(GOenrichment_mm, 'WGCNA/macaque_wgcna_BP_enrichment.txt', sep = '\t', quote = F, row.names = F)
