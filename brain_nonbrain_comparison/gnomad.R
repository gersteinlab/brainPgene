cols <- c(
  "Elevated in brain" = "#b5dfcc",
  "Elevated in non-brain" = "#fdcdac",
  "Low tissue-specificity" = "#cbd5e8",  
  "Elevated in at least one region" = "#f2e8db",
  "Low region-specificity" = "#bf8552",  
  "Not detected" = "#5e3023"
)
# pgene
pgene_AF <- read.table('brain_vs._non-brain/gnomad/pgene.AF.prop.txt', header = F, sep = '\t')
colnames(pgene_AF) <- c('row','Proportion')
pgene <- rbind(pgene_tissueSpecDf, pgene_regionSpecDf)
pgene <- pgene_AF %>% merge(., y = pgene, by = 'row')
p1 <- ggplot(pgene, aes(x = subtype, y = Proportion, fill = subtype)) + geom_boxplot() + facet_wrap(~group, scales = 'free', nrow = 1) + scale_fill_manual(values = cols) + 
  labs(x = '', y ='Rare allele proportion', title = 'Pseudogene') + theme_bw(base_size = 14) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=90, hjust=1)) + ylim(c(0,1))

# parent
pc_AF <- read.table('brain_vs._non-brain/1KG/protein_coding.AF.prop.txt', header = F, sep = '\t')
colnames(pc_AF) <- c('row','Proportion')
parent <- rbind(parent_tissueSpecDf, parent_regionSpecDf)
parent <- pc_AF %>% as.data.frame() %>% merge(., y = parent, by = 'row')
p2 <- ggplot(parent, aes(x = subtype, y = Proportion, fill = subtype)) + geom_boxplot() + facet_wrap(~group, scales = 'free', nrow = 1) + scale_fill_manual(values = cols) +
  labs(x = '', y = 'Rare allele proportion', title = 'Parent gene') + theme_bw(base_size = 14) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=90, hjust=1)) +  ylim(c(0,1))

# protein-coding
pc <- rbind(pc_tissueSpecDf, pc_regionSpecDf)
pc <- pc_AF %>% as.data.frame() %>% merge(., y = pc, by = 'row')
p3 <- ggplot(pc, aes(x = subtype, y = Proportion, fill = subtype)) + geom_boxplot() + facet_wrap(~group, scales = 'free', nrow = 1) + scale_fill_manual(values = cols) + 
  labs(x = '', y = 'Rare allele proportion', title = 'Protein-coding gene') + theme_bw(base_size = 14) + theme(legend.position = "none") + theme(axis.text.x=element_text(angle=90, hjust=1)) +  ylim(c(0,1))

p <- p1 / p2 / p3 + plot_layout(axes = 'collect')
ggsave('brain_vs._non-brain/rare_allele_proportion.pdf', p, width = 8, height = 12)


tbl1 <- pgene %>% na.omit %>% group_by(group) %>% rstatix::wilcox_test(Proportion ~ subtype, p.adjust.method = 'fdr')
tbl2 <- parent %>% na.omit %>% group_by(group) %>% rstatix::wilcox_test(Proportion ~ subtype, p.adjust.method = 'fdr')
tbl3 <- pc %>% na.omit %>% group_by(group) %>% rstatix::wilcox_test(Proportion ~ subtype, p.adjust.method = 'fdr')
write.table(tbl1, 'brain_vs._non-brain/tables/rare_prop_pgene_wilcox_test.txt', quote = F, row.names = F, sep = '\t')
write.table(tbl2, 'brain_vs._non-brain/tables/rare_prop_parent_wilcox_test.txt', quote = F, row.names = F, sep = '\t')
write.table(tbl3, 'brain_vs._non-brain/tables/rare_prop_pc_wilcox_test.txt', quote = F, row.names = F, sep = '\t')
