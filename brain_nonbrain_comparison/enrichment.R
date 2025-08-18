library(clusterProfiler)
library(org.Hs.eg.db)

pc_elv_in_brain <- pc_tissueSpecDf %>% dplyr::filter(subtype == 'Elevated in brain') %>% pull(row)
pc_elv_in_nonbrain <- pc_tissueSpecDf %>% dplyr::filter(subtype == 'Elevated in non-brain') %>% pull(row)
pc_elv_in_regions <- pc_regionSpecDf %>% dplyr::filter(subtype == 'Elevated in at least one region') %>% pull(row)
# tissue-specific: brain
ego_elv_in_brain <- enrichGO(gene = pc_elv_in_brain,OrgDb = org.Hs.eg.db, keyType = 'ENSEMBL',ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, readable = TRUE) %>% simplify()
ego_elv_in_brain <- ego_elv_in_brain %>% as.data.frame() %>% slice_min(order_by = p.adjust, n = 5) %>%
  mutate(group = 'Elevated in brain', type = 'Biological Process')
# KEGG
gene.df <- bitr(pc_elv_in_brain,fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
gene.kegg <- bitr_kegg(gene.df$ENTREZID, fromType = 'ncbi-geneid', toType = 'kegg', organism = 'hsa')
ekegg_elv_in_brain <- enrichKEGG(gene.kegg$kegg, organism = "hsa", keyType = "kegg")
ekegg_elv_in_brain <- ekegg_elv_in_brain %>% as.data.frame() %>% slice_min(order_by = p.adjust, n = 5) %>% mutate(group = 'Elevated in brain', type = 'Pathway')

# tissue-specific: non-brain
ego_elv_in_nonbrain <- enrichGO(gene = pc_elv_in_nonbrain,OrgDb = org.Hs.eg.db, keyType = 'ENSEMBL',ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, readable = TRUE) %>% simplify()
ego_elv_in_nonbrain <- ego_elv_in_nonbrain %>% as.data.frame() %>% slice_min(order_by = p.adjust, n = 5) %>%
  mutate(group = 'Elevated in non-brain', type = 'Biological Process')
# KEGG
gene.df <- bitr(pc_elv_in_nonbrain,fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
gene.kegg <- bitr_kegg(gene.df$ENTREZID, fromType = 'ncbi-geneid', toType = 'kegg', organism = 'hsa')
ekegg_elv_in_nonbrain <- enrichKEGG(gene.kegg$kegg, organism = "hsa", keyType = "kegg")
ekegg_elv_in_nonbrain <- ekegg_elv_in_nonbrain %>% as.data.frame() %>% slice_min(order_by = p.adjust, n = 5) %>% mutate(group = 'Elevated in non-brain', type = 'Pathway')

# region-specific
ego_elv_in_regions <- enrichGO(gene = pc_elv_in_regions,OrgDb = org.Hs.eg.db, keyType = 'ENSEMBL',ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, readable = TRUE) %>% simplify()
ego_elv_in_regions<- ego_elv_in_regions %>% as.data.frame() %>% slice_min(order_by = p.adjust, n = 5) %>%
  mutate(group = 'Elevated in at least one region', type = 'Biological Process')

# KEGG
gene.df <- bitr(pc_elv_in_regions,fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
gene.kegg <- bitr_kegg(gene.df$ENTREZID, fromType = 'ncbi-geneid', toType = 'kegg', organism = 'hsa')
ekegg_elv_in_regions <- enrichKEGG(gene.kegg$kegg, organism = "hsa", keyType = "kegg")
ekegg_elv_in_regions <- ekegg_elv_in_regions %>% as.data.frame() %>% slice_min(order_by = p.adjust, n = 5) %>% mutate(group = 'Elevated in at least one region', type = 'Pathway')


# plot
library(stringr)
library(latex2exp)
p <- rbind(ego_elv_in_brain, ego_elv_in_nonbrain, ego_elv_in_regions, ekegg_elv_in_brain, ekegg_elv_in_nonbrain, ekegg_elv_in_regions) %>%
  ggplot(aes(x = group, y = str_to_title(Description), fill = -log10(p.adjust)))+ geom_tile() + facet_grid(type ~. , scales = 'free_y')  + theme_bw() +
  labs(x = '', y = '') + guides(fill=guide_legend(title=TeX(r'($-\log _ {10} FDR$)'))) + theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave('brain_vs._non-brain/enrichment.pdf', p, width = 6, height = 8)