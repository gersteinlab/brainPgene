# tissue-specificity
pgene_tissueSpec <- readRDS('brain_vs._non-brain/tissue_specificity/DEG_Pseudogene.rds')
pgene_tissueSpecDf <- data.frame(row = c(pgene_tissueSpec[[1]], pgene_tissueSpec[[2]], pgene_tissueSpec[[3]]),
                                 subtype = c(rep('Elevated in brain', length(pgene_tissueSpec[[1]])),
                                             rep('Elevated in non-brain', length(pgene_tissueSpec[[2]])),
                                             rep('Low tissue-specificity', length(pgene_tissueSpec[[3]]))),
                                 group = 'Tissue-specificity')
write.table(pgene_tissueSpecDf, 'brain_vs._non-brain/tissue_specificity/pgene_tissue_specificity2.txt', quote = F, sep = '\t', row.names = F)



pgene_tissueSpec <- readRDS('brain_vs._non-brain/tissue_specificity/DEG_Pseudogene_nofilter.rds')
pgene_tissueSpecDf <- data.frame(row = c(pgene_tissueSpec[[1]], pgene_tissueSpec[[2]], pgene_tissueSpec[[3]]),
                                 subtype = c(rep('Elevated in brain', length(pgene_tissueSpec[[1]])),
                                             rep('Elevated in non-brain', length(pgene_tissueSpec[[2]])),
                                             rep('Low tissue-specificity', length(pgene_tissueSpec[[3]]))),
                                 group = 'Tissue-specificity')
write.table(pgene_tissueSpecDf, 'brain_vs._non-brain/tissue_specificity/pgene_tissue_specificity_nofilter.txt', quote = F, sep = '\t', row.names = F)

# region specificity
pgene_regionSpec <- readRDS('brain_vs._non-brain/region_specificity/DEG_Pseudogene.rds')
pgene_regionSpecDf <- data.frame(row = c(pgene_regionSpec[[1]], pgene_regionSpec[[2]], pgene_regionSpec[[3]]),
                                 subtype = c(rep('Elevated in at least one region', length(pgene_regionSpec[[1]])),
                                             rep('Low region-specificity', length(pgene_regionSpec[[2]])),
                                             rep('Not detected', length(pgene_regionSpec[[3]]))),
                                 group = 'Region-specificity')
write.table(pgene_regionSpecDf, 'brain_vs._non-brain/region_specificity/pgene_region_specificity.txt', quote = F, sep = '\t', row.names = F)

##
pgene_regionSpec <- readRDS('brain_vs._non-brain/region_specificity/DEG_Pseudogene_nofilter.rds')
pgene_regionSpecDf <- data.frame(row = c(pgene_regionSpec[[1]], pgene_regionSpec[[2]], pgene_regionSpec[[3]]),
                                 subtype = c(rep('Elevated in at least one region', length(pgene_regionSpec[[1]])),
                                             rep('Low region-specificity', length(pgene_regionSpec[[2]])),
                                             rep('Not detected', length(pgene_regionSpec[[3]]))),
                                 group = 'Region-specificity')
write.table(pgene_regionSpecDf, 'brain_vs._non-brain/region_specificity/pgene_region_specificity_nofilter.txt', quote = F, sep = '\t', row.names = F)


# combined with pgene type
pgene_parent <- read.table('parent/pgene_filtered_merged.txt', sep = '\t', header = T) 
pgene_parent <- pgene_parent %>%  mutate(row = unlist(lapply(pgene_parent$Pseudogene, function(x) strsplit(x, '\\.')[[1]][1])))
unitary_pgene <- read.table('parent/fasta/unitary_pgene.bed', sep = '\t', header = F) %>% dplyr::select(V28) 
unitary_pgene <- unitary_pgene %>%  mutate(row = unlist(lapply(unitary_pgene$V28, function(x) strsplit(x, '\\.')[[1]][1])), Type = 'Unitary')
pgene_metadata <- rbind(pgene_parent[,c(ncol(pgene_parent),ncol(pgene_parent)-1)], 
                        unitary_pgene[,c(ncol(unitary_pgene), ncol(unitary_pgene)-1)])

data <- rbind(pgene_tissueSpecDf, pgene_regionSpecDf) %>% merge(., pgene_metadata, by = 'row')
p <- data %>% group_by(subtype, Type, group) %>% summarise(n = n()) %>%
  ggplot(aes(x = subtype, fill = Type, y = n)) +  geom_bar(position='stack', stat='identity') + facet_wrap(~group, scales = 'free', ncol = 1) + scale_fill_igv()  + theme_bw(base_size = 14) + theme(legend.title=element_blank(), axis.text = element_text(colour = 'black')) + labs(y = 'Count', x = '') + coord_flip()
ggsave('brain_vs._non-brain/pgene_subtype_by_specificity.pdf', p, width = 12, height = 4)
