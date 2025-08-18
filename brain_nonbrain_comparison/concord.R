# tissue
pgene_parent <- read.table('parent/pgene_dist_HiC_filtered.txt', sep = '\t', header = T)
pgene_tissueSpecDf <- read.table('brain_vs._non-brain/tissue_specificity/pgene_tissue_specificity.txt', sep = '\t', header = T)
parent_tissueSpecDf <-read.table('brain_vs._non-brain/tissue_specificity/parent_tissue_specificity.txt', sep = '\t', header = T)

pgene_parent <- pgene_parent %>% merge(y = pgene_tissueSpecDf, by.x = 'Pgene', by.y = 'row') %>% mutate(Pseudogene = Pgene, Pseudogene_type = subtype) %>% 
  dplyr::select(-Pgene, -subtype, -group)
pgene_parent <- pgene_parent %>% merge(y = parent_tissueSpecDf, by.x = 'ParentGene', by = 'row') %>% mutate(Parent = ParentGene, Parent_type = subtype) %>% 
  dplyr::select(-ParentGene, -subtype, -group)


pgene_parent <- pgene_parent %>% mutate(name = interaction(pgene_parent$Pseudogene_type, pgene_parent$Parent_type, sep = '/'),
                                        Type = ifelse(Pseudogene_type == Parent_type, 'Consensus', 'Non-consensus')) %>% group_by(name, Type) %>% summarise(n = n())

pgene_parent$name <- factor(pgene_parent$name, levels = c('Elevated in brain/Elevated in brain', 'Elevated in non-brain/Elevated in non-brain', 'Low tissue-specificity/Low tissue-specificity', 'Elevated in brain/Elevated in non-brain', 'Elevated in brain/Low tissue-specificity', 'Elevated in non-brain/Elevated in brain', 'Elevated in non-brain/Low tissue-specificity', 'Low tissue-specificity/Elevated in brain', 'Low tissue-specificity/Elevated in non-brain'))

p <- pgene_parent %>% 
  ggplot(aes(x = 2, y = n, fill = name)) +
  geom_col(aes(alpha = Type), width = 1, size = 0.5, position = position_stack()) +
  coord_polar(theta = 'y') + scale_alpha_manual(values=c(1, 0.5)) +
  xlim(c(0.5, 2 + 0.5)) + scale_fill_brewer(palette = 'Spectral') + ggtitle('Tissue-specificity') +
  geom_text(aes(label = n, group = name), 
            position = position_stack(vjust = 0.5),  # Centers text in each segment
            size = 6, color = "black") +
  theme( panel.background = element_rect(fill = "white"),
         panel.grid = element_blank(),
         axis.title = element_blank(),
         axis.ticks = element_blank(),
         axis.text = element_blank())
ggsave('brain_vs._non-brain/concord_tissue.pdf', p, width = 10, height = 10)
write.table(pgene_parent, 'brain_vs._non-brain/tables/concord_tissue.txt', quote = F, sep = '\t', row.names = F)

# region
pgene_parent <- read.table('parent/pgene_dist_HiC_filtered.txt', sep = '\t', header = T)
pgene_regionSpecDf <- read.table('brain_vs._non-brain/region_specificity/pgene_region_specificity.txt', sep = '\t', header = T)
parent_regionSpecDf <-read.table('brain_vs._non-brain/region_specificity/parent_region_specificity.txt', sep = '\t', header = T)

pgene_parent <- pgene_parent %>% merge(y = pgene_regionSpecDf, by.x = 'Pgene', by.y = 'row') %>% mutate(Pseudogene = Pgene, Pseudogene_type = subtype) %>% 
  dplyr::select(-Pgene, -subtype, -group)
pgene_parent <- pgene_parent %>% merge(y = parent_regionSpecDf, by.x = 'ParentGene', by = 'row') %>% mutate(Parent = ParentGene, Parent_type = subtype) %>% 
  dplyr::select(-ParentGene, -subtype, -group)
pgene_parent <- pgene_parent %>% mutate(name = interaction(pgene_parent$Pseudogene_type, pgene_parent$Parent_type, sep = '/'),
                                        Type = ifelse(Pseudogene_type == Parent_type, 'Consensus', 'Non-consensus')) %>% group_by(name, Type) %>% summarise(n = n())

pgene_parent$name <- factor(pgene_parent$name, levels = c('Elevated in at least one region/Elevated in at least one region', 'Low region-specificity/Low region-specificity', 'Not detected/Not detected', 'Elevated in at least one region/Low region-specificity', 'Elevated in at least one region/Not detected', 'Low region-specificity/Elevated in at least one region', 'Low region-specificity/Not detected', 'Not detected/Elevated in at least one region', 'Not detected/Low region-specificity'))

p2 <- pgene_parent %>% 
  ggplot(aes(x = 2, y = n, fill = name)) +
  geom_col(aes(alpha = Type), width = 1, size = 0.5, position = position_stack()) +
  coord_polar(theta = 'y') + scale_alpha_manual(values=c(1, 0.5)) +
  xlim(c(0.5, 2 + 0.5)) + scale_fill_brewer(palette = 'Spectral') + ggtitle('Region-specificity') + 
  geom_text(aes(label = n, group = name), 
            position = position_stack(vjust = 0.5),  # Centers text in each segment
            size = 6, color = "black") +
  theme( panel.background = element_rect(fill = "white"),
         panel.grid = element_blank(),
         axis.title = element_blank(),
         axis.ticks = element_blank(),
         axis.text = element_blank())
ggsave('brain_vs._non-brain/concord_region.pdf', p2, width = 10, height = 10)
write.table(pgene_parent, 'brain_vs._non-brain/tables/concord_region.txt', quote = F, sep = '\t', row.names = F)
