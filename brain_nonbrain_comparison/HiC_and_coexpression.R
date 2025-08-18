# HiC interaction
pgene_parent <- read.table('parent/pgene_dist_HiC.txt', sep = '\t', header = T)
pgene_parent <- pgene_parent %>% filter(Type == 'unprocessed_pseudogene' | Type == 'processed_pseudogene') %>% 
  mutate(Type = if_else(Type == 'unprocessed_pseudogene', 'Duplicated', 'Processed'))
write.table(pgene_parent, 'parent/pgene_dist_HiC_filtered.txt', sep = '\t', row.names = F, col.names = T, quote = F)

library(latex2exp)
options(scipen=200)
p <- pgene_parent %>% filter(Count > 0) %>% mutate(x = log10(Distance/1000), y = log10(Count)) %>%
  ggplot(aes(x = x, y = y, color = Type)) + geom_point() + stat_smooth(method = lm, se = TRUE, fullrange = TRUE) + scale_color_startrek(alpha = 0.8) + 
  labs(x = TeX(r'($\log_{10}(\textrm{Distance (kbp)})$)'), y = TeX(r'($\log_{10}(\textrm{Count})$)')) + facet_wrap(~Type)  +
  theme_bw(base_size = 12) +
  theme(legend.title=element_blank(), legend.position = 'None', axis.text = element_text(colour = 'black'))
ggsave('brain_vs._non-brain/pgene_parent_Hi-C.pdf', p, width = 8, height = 8)

processed_lm <- pgene_parent %>% filter(Type == 'Processed', Count > 0) %>% mutate(x = log10(Distance/1000), y = log10(Count)) %>% lm(y~x, .) %>% summary()
duplicated_lm <- pgene_parent %>% filter(Type == 'Duplicated', Count > 0) %>% mutate(x = log10(Distance/1000), y = log10(Count)) %>% lm(y~x, .) %>% summary()

# co-expression
pgene_expression <- read.table('quantification/Capstone/pgene_expression_Capstone.tsv', header =  T, check.names = F, row.names = 1)
pc_expression <- read.table('quantification/Capstone/pc_expression_Capstone.tsv', header = T, check.names = F, row.names = 1)
# remove rows with all zeros
pgene_expression <- pgene_expression[rowSums(pgene_expression) != 0, ]
pc_expression <- pc_expression[rowSums(pc_expression) != 0, ]

pgene_parent <- read.table('parent/pgene_dist_HiC_filtered.txt', sep = '\t', header = T)
pgene_parent <- pgene_parent[,c(1,2,3,5)] %>%
  mutate(Group = ifelse(is.na(Distance), 'Different Chrom.', ifelse(Distance < 1000000, 'Same Chrom. < 1 Mbp', 'Same Chrom. > 1 Mbp')))

correlation <- as.numeric()
for (i in 1:nrow(pgene_parent)) {
  pgene <- strsplit(pgene_parent[i,1], '\\.')[[1]][1]
  parent <- strsplit(pgene_parent[i,2], '\\.')[[1]][1]
  df1 <- pgene_expression[rownames(pgene_expression) == pgene,]
  df2 <- pc_expression[rownames(pc_expression) == parent,]
  if(nrow(df1) == 0 || nrow(df2) == 0){
    correlation <- c(correlation, NA) 
  }else{
    correlation <- c(correlation, cor(unlist(df1), unlist(df2), method = 'pearson'))
  }
}

pgene_parent <- pgene_parent %>% mutate(cor = correlation) 
my_comparisons <- list(c('Different Chrom.', 'Same Chrom. < 1 Mbp'),  c('Same Chrom. < 1 Mbp', 'Same Chrom. > 1 Mbp'), c('Different Chrom.', 'Same Chrom. > 1 Mbp'))
library(ggpubr)
p <- ggplot(pgene_parent, aes(x = Group, y = cor, fill = Group, color = Group)) + 
  geom_boxplot() + 
  stat_compare_means(method = 'wilcox.test', 
                     comparisons = my_comparisons, p.adjust.method = "fdr", label = 'p.signif', labelhide.ns = TRUE, step.increase = 0.06, tip.length = 0.01, method.args = list(alternative = 'two.sided')) +  
  labs(x = '', y ='Expression correlation') + 
  theme_bw(base_size = 12) +
  theme(legend.position="none", axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), axis.text = element_text(color = 'black')) + 
  facet_grid(~Type) + scale_fill_aaas(alpha = 0.6) + scale_color_aaas(alpha = 0.6)

tbl <- pgene_parent  %>% group_by(Type) %>% rstatix::wilcox_test(cor ~ Group, p.adjust.method = 'fdr') 
write.table(tbl, 'brain_vs._non-brain/expression_correlation_wilcox_test.txt', quote = F, row.names = F, sep = '\t')
ggsave('brain_vs._non-brain/tables/expression_correlation.pdf', p, width = 10, height = 8)
