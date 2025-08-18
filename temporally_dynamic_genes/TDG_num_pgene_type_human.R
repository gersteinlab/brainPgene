unitary_pgene <- read.table('parent/fasta/unitary_pgene.bed', sep = '\t', header = F) %>% dplyr::select(V28) %>% pull()
duplicated_pgene <- read.table('parent/fasta/duplicated_pgene.bed', sep = '\t', header = F) %>% dplyr::select(V28) %>% pull()
processed_pgene <-  read.table('parent/fasta/processed_pgene.bed', sep = '\t', header = F) %>% dplyr::select(V28) %>% pull()
unitary_pgene <- unlist(lapply(unitary_pgene, function(x) strsplit(x, '\\.')[[1]][1])) %>% unique()
duplicated_pgene <- unlist(lapply(duplicated_pgene, function(x) strsplit(x, '\\.')[[1]][1])) %>% unique()
processed_pgene <- unlist(lapply(processed_pgene, function(x) strsplit(x, '\\.')[[1]][1])) %>% unique()
pgene_type <- data.frame(row.names = c(duplicated_pgene, processed_pgene, unitary_pgene), type = c(rep('Duplicated', length(duplicated_pgene)),
                                                                                                   rep('Processed', length(processed_pgene)),
                                                                                                   rep('Unitary', length(unitary_pgene))))
pgene <- readRDS('temporal_dynamic_expression/human/pgene_temporallyDynamicGene.rds')
pgene <- sapply(names(pgene), function(x) merge(pgene[[x]], pgene_type, by = 'row.names'), simplify = F, USE.NAMES = T)
pgene_res <- sapply(names(pgene), function(x) table(factor(pgene[[x]]$type.y, levels = c('Duplicated', 'Processed', 'Unitary'))), simplify = F, USE.NAMES = T)

p <- pgene_res %>% as_tibble() %>% mutate(Group = factor(c('Duplicated', 'Processed', 'Unitary'))) %>% 
  pivot_longer(cols = A1C:VFC, names_to = 'region', values_to = 'count') %>%
  ggplot(., aes(x=region, y=as.numeric(count), fill = Group)) + geom_bar(position="stack", stat="identity") + scale_fill_igv() +
  labs(x = '', y='Number of temporally dynamic pseudogenes') +
  theme_bw(base_size = 12) + theme(legend.title=element_blank(), axis.text = element_text(color = 'black'))
ggsave('temporal_dynamic_expression/human_vs._macaque/subtype/pgene_#human.pdf', p, width = 8, height = 5)


a <- pgene_res %>% as_tibble() %>% mutate(Group = factor(c('Duplicated', 'Processed', 'Unitary'))) %>% 
  pivot_longer(cols = A1C:VFC, names_to = 'region', values_to = 'count')  %>%
  dplyr::filter(Group == "Processed") %>% dplyr::rename(processed_count = count)

b <- pgene_res %>% as_tibble() %>% mutate(Group = factor(c('Duplicated', 'Processed', 'Unitary'))) %>% 
  pivot_longer(cols = A1C:VFC, names_to = 'region', values_to = 'count') %>% 
  group_by(region) %>% 
  summarise(total_count = sum(count, na.rm = TRUE))

a %>%
  left_join(b, by = "region") %>%
  mutate(proportion_processed = processed_count / total_count) %>% summarise(avg = mean(proportion_processed))
