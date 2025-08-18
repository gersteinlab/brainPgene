rm(list = ls())
setwd('/Users/yzjiang/Desktop/brain_pgene/')
library(dplyr)
ASD_tbl <- read.table('DEG_across_disorders/DEG/MAGMA/ASD.txt', sep = '\t', header = T) %>% mutate(disease = 'ASD')
BD_tbl <- read.table('DEG_across_disorders/DEG/MAGMA/BD.txt', sep = '\t', header = T) %>% mutate(disease = 'BD')
SCZ_tbl <- read.table('DEG_across_disorders/DEG/MAGMA/SCZ.txt', sep = '\t', header = T) %>% mutate(disease = 'SCZ')
df <- rbind(ASD_tbl, BD_tbl, SCZ_tbl)

pgene_tissue_specific <- read.table('brain_vs._non-brain/tissue_specificity/pgene_tissue_specificity.txt', sep = '\t', header = T) %>%
  mutate(gene_type = 'pseudogene')
pc_tissue_specific <- read.table('brain_vs._non-brain/tissue_specificity/pc_tissue_specificity.txt', sep = '\t', header = T) %>%
  mutate(gene_type = 'protein_coding')
tissue_specific <- rbind(pgene_tissue_specific, pc_tissue_specific)

pgene_region_specific <- read.table('brain_vs._non-brain/region_specificity/pgene_region_specificity.txt', sep = '\t', header = T) %>%
  mutate(gene_type = 'pseudogene')
pc_region_specific <- read.table('brain_vs._non-brain/region_specificity/pc_region_specificity.txt', sep = '\t', header = T) %>%
  mutate(gene_type = 'protein_coding')
region_specific <- rbind(pgene_region_specific, pc_region_specific)

res <- df %>% merge(., y = tissue_specific, by.x = 'Gene', by.y = 'row') %>%
  dplyr::select(Gene, DE, HMAGMA, Type, disease, subtype) %>% dplyr::rename(tissue_specificity = subtype) %>%
  merge(., y = region_specific, by.x = 'Gene', by.y = 'row') %>%
  dplyr::select(Gene, DE, HMAGMA, Type, disease, tissue_specificity, subtype) %>% dplyr::rename(region_specificity = subtype)


df_counts <- res %>% filter(DE == HMAGMA) %>%
  group_by(disease, tissue_specificity, region_specificity) %>%
  summarise(count = n()) %>%
  ungroup()

df_counts_2 <- res %>% filter(DE == HMAGMA) %>%
  group_by(Type, disease, tissue_specificity, region_specificity) %>%
  summarise(count = n()) %>%
  ungroup()

library(ggsci)
p <- ggplot(df_counts_2, aes(x = tissue_specificity, y = region_specificity, color = Type, size = count, shape = disease)) + 
  geom_jitter(width = 0.3, height = 0.3, alpha = 1) +
  scale_size_continuous(range = c(3, 7)) +
  scale_color_igv(alpha = 1) +
  labs(x = 'Tissue specificity', y = 'Region specificity', shape = 'Disease', size ='Count') + coord_fixed() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        legend.position = 'right')

ggsave('DEG_across_disorders/DE_tissue_regin_specificity.pdf', width = 8, height = 8)
