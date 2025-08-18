# human
library(dplyr)
pgene <- readRDS('temporal_dynamic_expression/human/pgene_temporallyDynamicGene.rds')
pcg  <- readRDS('temporal_dynamic_expression/human/protein_coding_temporallyDynamicGene.rds')

extract_TDgenes <- function(df, region_name) {
  df$gene <- rownames(df)         
  df$region <- region_name        
  df[, c("region", "gene", "type")]
}

pgene_extracted <- do.call(
  rbind,
  lapply(names(pgene), function(region) extract_TDgenes(pgene[[region]], region)))

pcg_extracted <- do.call(
  rbind,
  lapply(names(pcg), function(region) extract_TDgenes(pcg[[region]], region)))

human_TD <- rbind(pgene_extracted, pcg_extracted) 
human_TD <- human_TD %>% arrange(region)

# macaque
pgene <- readRDS('temporal_dynamic_expression/macaque/pgene_temporallyDynamicGene.rds')
pcg  <- readRDS('temporal_dynamic_expression/macaque/protein_coding_temporallyDynamicGene.rds')

extract_TDgenes <- function(df, region_name) {
  df$gene <- rownames(df)         
  df$region <- region_name        
  df[, c("region", "gene", "type")]
}

pgene_extracted <- do.call(
  rbind,
  lapply(names(pgene), function(region) extract_TDgenes(pgene[[region]], region)))

pcg_extracted <- do.call(
  rbind,
  lapply(names(pcg), function(region) extract_TDgenes(pcg[[region]], region)))

macaque_TD <- rbind(pgene_extracted, pcg_extracted)
macaque_TD <- macaque_TD %>% arrange(region)

write.table(human_TD, 'temporal_dynamic_expression/human_TDgene.txt', sep = '\t', quote = F, row.names = F)
write.table(macaque_TD, 'temporal_dynamic_expression/macaque_TDgene.txt', sep = '\t', quote = F, row.names = F)
