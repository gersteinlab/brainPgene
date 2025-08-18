# load the GTF file
grch37_gtf <- import('GTF/gencode.v39lift37.annotation.gtf')
grch37_exon <- as.data.frame(grch37_gtf) %>% dplyr::filter(type == 'exon') %>% dplyr::select(seqnames, start, end, strand, gene_id)
grch37_exon <- grch37_exon %>% dplyr::mutate(row = unlist(lapply(grch37_exon$gene_id, function(x) strsplit(x, '\\.')[[1]][1]))) %>% dplyr::select(-gene_id)
# data
pgene_tissueSpecDf <- read.table('brain_vs._non-brain/tissue_specificity/pgene_tissue_specificity.txt', sep = '\t', header = T)
parent_tissueSpecDf <-read.table('brain_vs._non-brain/tissue_specificity/parent_tissue_specificity.txt', sep = '\t', header = T)
pc_tissueSpecDf <- read.table('brain_vs._non-brain/tissue_specificity/pc_tissue_specificity.txt', sep = '\t', header = T)

pgene_regionSpecDf <- read.table('brain_vs._non-brain/region_specificity/pgene_region_specificity.txt', sep = '\t', header = T)
parent_regionSpecDf <-read.table('brain_vs._non-brain/region_specificity/parent_region_specificity.txt', sep = '\t', header = T)
pc_regionSpecDf <- read.table('brain_vs._non-brain/region_specificity/pc_region_specificity.txt', sep = '\t', header = T)

cols <- c(
  "Elevated in brain" = "#b5dfcc",
  "Elevated in non-brain" = "#fdcdac",
  "Low tissue-specificity" = "#cbd5e8",  
  "Elevated in at least one region" = "#f2e8db",
  "Low region-specificity" = "#bf8552",  
  "Not detected" = "#5e3023"
)

# pgene
library(GenomicScores)
library(phastCons20way.UCSC.hg19)
library(GenomicRanges)
grch37_exon_gr <- makeGRangesFromDataFrame(grch37_exon, starts.in.df.are.0based = F, keep.extra.columns = T)
res <- gscores(phastCons100way.UCSC.hg19, grch37_exon_gr)
res <- as.data.frame(res) %>% na.omit()
res <- res %>% dplyr::group_by(row) %>%
  summarise(weighted_phastcons = sum(default * width) / sum(width))

pgene <- rbind(pgene_tissueSpecDf, pgene_regionSpecDf)
pgene <- pgene %>% as.data.frame() %>% merge(., y = res, by = 'row') %>% distinct()

tbl1 <- pgene %>% group_by(group) %>% wilcox_test(weighted_phastcons~subtype, p.adjust.method = 'fdr')
p1 <- ggplot(pgene, aes(x = subtype, y = weighted_phastcons, fill = subtype)) + 
  geom_boxplot(notch = T, notchwidth = 0.7)  + facet_wrap(~group, scales = 'free', nrow = 1)+ scale_fill_manual(values = cols)+
  labs(x = '', y = 'phastCons 100-way score', title = 'Pseudogene') + theme_bw(base_size = 12) +
  theme(legend.position = "none")  + coord_flip() 

# parent
parent <- rbind(parent_tissueSpecDf, parent_regionSpecDf)
parent <- parent %>% as.data.frame() %>% merge(., y = res, by = 'row') %>% distinct()
tbl2 <- parent %>% group_by(group) %>% wilcox_test(weighted_phastcons~subtype, p.adjust.method = 'fdr')
p2 <- ggplot(parent, aes(x = subtype, y = weighted_phastcons, fill = subtype)) + geom_boxplot(notch = T, notchwidth = 0.7) + facet_wrap(~group, scales = 'free', nrow = 1) + scale_fill_manual(values = cols) + labs(x = '', y = 'phastCons 100-way score', title = 'Parent Gene') + theme_bw(base_size = 14) + theme(legend.position = "none") + coord_flip() 

# protein-coding
pc <- rbind(pc_tissueSpecDf, pc_regionSpecDf)
pc <- pc %>% as.data.frame() %>% merge(., y = res, by = 'row') %>% distinct()
tbl3 <- pc %>% group_by(group) %>% wilcox_test(weighted_phastcons~subtype, p.adjust.method = 'fdr')
p3 <- ggplot(pc, aes(x = subtype, y = weighted_phastcons, fill = subtype)) + geom_boxplot(notch = T, notchwidth = 0.7) + facet_wrap(~group, scales = 'free', nrow = 1) + scale_fill_manual(values = cols) + labs(x = '', y = 'phastCons 100-way score', title = 'Protein-coding Gene') + theme_bw(base_size = 14) + theme(legend.position = "none") + coord_flip() 

p <- p1 + p2 + p3
ggsave('brain_vs._non-brain/phastcons.pdf', p, width = 24, height = 4)

write.table(tbl1, 'brain_vs._non-brain/tables/phastcons100way_pgene_wilcox_test.txt', quote = F, row.names = F, sep = '\t')
write.table(tbl2, 'brain_vs._non-brain/tables/phastcons100way_parent_wilcox_test.txt', quote = F, row.names = F, sep = '\t')
write.table(tbl3, 'brain_vs._non-brain/tables/phastcons100way_pc_wilcox_test.txt', quote = F, row.names = F, sep = '\t')
