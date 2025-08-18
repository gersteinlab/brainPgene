# phastcons
grch37_gtf <- import('GTF/gencode.v39lift37.annotation.gtf')
grch37_exon <- as.data.frame(grch37_gtf) %>% dplyr::filter(type == 'exon') %>% dplyr::select(seqnames, start, end, strand, gene_id)
grch37_exon <- grch37_exon %>% dplyr::mutate(row = unlist(lapply(grch37_exon$gene_id, function(x) strsplit(x, '\\.')[[1]][1]))) %>% dplyr::select(-gene_id)

grch37_exon_gr <- makeGRangesFromDataFrame(grch37_exon, starts.in.df.are.0based = F, keep.extra.columns = T)
res <- gscores(phastCons100way.UCSC.hg19, grch37_exon_gr)
res <- as.data.frame(res) %>% na.omit()
res <- res %>% dplyr::group_by(row) %>%
  summarise(weighted_phastcons = sum(default * width) / sum(width))

# 1KG 
pgene_AF <- read.table('brain_vs._non-brain/gnomad/pgene.AF.prop.txt', header = F, sep = '\t')
colnames(pgene_AF) <- c('row','Proportion')
pc_AF <- read.table('brain_vs._non-brain/gnomad/protein_coding.AF.txt', header = F, sep = '\t')
colnames(pc_AF) <- c('row','Proportion')

# TDG and non-TDG
pgene_TDG_df <- read.table('./temporal_dynamic_expression/human/pgene_TDG_nonTDG.txt', header = T, sep = '\t')
pc_TDG_df <- read.table('./temporal_dynamic_expression/human/pc_TDG_nonTDG.txt', header = T, sep = '\t')

# plot
library(ggpubr)
p1 <- res %>% as.data.frame() %>% merge(., y = pgene_TDG_df, by = 'row') %>%
  ggplot(aes(x = Type, y = weighted_phastcons, fill = Type)) + geom_boxplot(notch = T, notchwidth = 0.7) + scale_fill_npg() + 
  labs(x = '', y = 'phastCons 100-way Score', title = 'Pseudogene') + theme_bw(base_size = 12) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(color = 'black')) + 
  stat_compare_means(method = "kruskal.test")

p2 <- pgene_AF %>% merge(y = pgene_TDG_df, by = 'row') %>%
  ggplot(aes(x = Type, y = Proportion, fill = Type)) + geom_boxplot(notch = T, notchwidth = 0.7) + scale_fill_npg() + 
  labs(x = '', y = 'Rare allele proportion', title = 'Pseudogene') + theme_bw(base_size = 12) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(color = 'black'))  +  ylim(c(0,1)) + 
  stat_compare_means(method = "kruskal.test")

p3 <- res %>% as.data.frame() %>% merge(., y = pc_TDG_df, by = 'row') %>%
  ggplot(aes(x = Type, y = weighted_phastcons, fill = Type)) + geom_boxplot(notch = T, notchwidth = 0.7) + scale_fill_npg() + 
  labs(x = '', y = 'phastCons 100-way Score', title = 'Protein-coding gene') + theme_bw(base_size = 12) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(color = 'black')) + 
  stat_compare_means(method = "kruskal.test")

p4 <- pc_AF %>% merge(y = pc_TDG_df, by = 'row') %>%
  ggplot(aes(x = Type, y = Proportion, fill = Type)) + geom_boxplot(notch = T, notchwidth = 0.7) + scale_fill_npg() + 
  labs(x = '', y = 'Rare allele proportion', title = 'Protein-coding gene') + theme_bw(base_size = 12) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(color = 'black'))  +  ylim(c(0,1)) + 
  stat_compare_means(method = "kruskal.test")


p <- p1 + p3 + plot_layout(axis_titles = 'collect')
ggsave('./temporal_dynamic_expression/phastcons.pdf', width = 5, height = 5)

p <- p2 + p4 + plot_layout(axis_titles = 'collect')
ggsave('./temporal_dynamic_expression/gnomad.pdf', width = 6, height = 5)
