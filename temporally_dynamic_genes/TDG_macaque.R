specimen_meta <- read.table('metadata/NHP_macaque/meta_RNAseq_bam_Macaque.csv', sep = ',', header = T)
edesign_table <- read.table('metadata/NHP_macaque/edesign_RNAseq_bam_Macaque.txt', sep = '\t', header = T)
edesign_table <- edesign_table[edesign_table$specimenID %in% specimen_meta$specimenID,]
# clean the data
edesign_table$specimenID <-  str_replace(edesign_table$specimenID, "(?<=\\d{3}).*(?=_)", '')
edesign_table$specimenID <- gsub('DLPFC', 'DFC', edesign_table$specimenID)
tissue_meta <- read.table('metadata/BrainSpan/tissue_group.csv', sep = ',', header = T)
tissue_meta <- tissue_meta %>% filter(sampleNum > 10)

TDG_acrossTissue <- function(expr_matrix, edesign_table, tissue_metadata, tissueName, degree = 2, rsq = 0.3){
  # subset expression matrix
  expr <- expr_matrix %>% dplyr::select(ends_with(tissueName))
  # subset design matrix
  edesign <- edesign_table[edesign_table[[tissueName]] == 1, ] %>% dplyr::select(specimenID, Time, Replicate, tissueName) %>% mutate(Time = log2(Time))
  rownames(edesign) <- edesign[['specimenID']]
  edesign <- edesign %>% dplyr::select(-specimenID)
  # run maSigPro
  d <- make.design.matrix(edesign, degree = degree)
  fit <- p.vector(expr, d, Q = 0.05, MT.adjust = 'BH', counts = FALSE)
  tstep <- T.fit(fit, step.method = 'two.ways.backward')
  sigs <- get.siggenes(tstep, rsq = rsq, vars = 'all')
  res <- sigs$sig.genes$coefficients
  # assign temporal dynamic gene type based on degree
  if(degree == 2){
    res <- res %>% mutate(axis = ifelse(betaTime2 == 0, NA, -betaTime/(2*betaTime2)))
    min <- min(edesign$Time); max <- max(edesign$Time)
    df <- assignTemporalGeneType(res, min, max)
  }else if(degree == 1){
    df <- res %>% mutate(type = ifelse(betaTime > 0, 'Rising', 'Falling'))
  }
  return(df)
}

pgene_expression <- read.table('quantification/NHP_macaque/pgene_expression_NHP_macaque.tsv', sep = '\t', header = T, row.names = 1)
# clean the data
colNames <- gsub('mRNA_ployA_.*_','',colnames(pgene_expression), perl = T)
colNames <- sapply(colNames, function(x) strsplit(x, '_')[[1]][5:4] %>% paste(collapse = '_'), simplify = T, USE.NAMES = F)
colNames <- gsub('DLPFC', 'DFC', colNames)
colNames <- str_replace(colNames, "(?<=\\d{3}).*(?=_)", '')
colnames(pgene_expression) <- colNames
# pgene
pgene <- sapply(tissue_meta$tissue, function(x) TDG_acrossTissue(pgene_expression, edesign_table, tissue_meta, x), simplify = FALSE, USE.NAMES = TRUE)
pgene_res <- sapply(names(pgene), function(x) 
  table(factor(pgene[[x]]$type, levels = c('Falling', 'Rising', 'U-shaped', 'Inverted U-shaped'))), simplify = F, USE.NAMES = T)
saveRDS(pgene, 'temporal_dynamic_expression/macaque/pgene_temporallyDynamicGene.rds')


# protein coding
pc_expression <- read.table('quantification/NHP_macaque/pc_expression_NHP_macaque.tsv', sep = '\t', header = T, row.names = 1)
colNames <- gsub('mRNA_ployA_.*_','',colnames(pc_expression), perl = T)
colNames <- sapply(colNames, function(x) strsplit(x, '_')[[1]][5:4] %>% paste(collapse = '_'), simplify = T, USE.NAMES = F)
colNames <- gsub('DLPFC', 'DFC', colNames)
colNames <- str_replace(colNames, "(?<=\\d{3}).*(?=_)", '')
colnames(pc_expression) <- colNames
protein_coding <- sapply(tissue_meta$tissue, function(x) TDG_acrossTissue(pc_expression, edesign_table, tissue_meta, x), simplify = FALSE, USE.NAMES = TRUE)
protein_coding_res <- sapply(names(protein_coding), function(x) 
  table(factor(protein_coding[[x]]$type, levels = c('Falling', 'Rising', 'U-shaped', 'Inverted U-shaped'))), simplify = F, USE.NAMES = T)
saveRDS(protein_coding, 'temporal_dynamic_expression/macaque/protein_coding_temporallyDynamicGene.rds')



total_counts <- pgene_res %>%
  as_tibble() %>%
  pivot_longer(cols = A1C:VFC, names_to = 'region', values_to = 'count') %>%
  group_by(region) %>%
  summarise(total_count = sum(as.numeric(count), na.rm = TRUE))

pgene_plot <- pgene_res %>% as_tibble() %>% mutate(Group = factor(c('Falling', 'Rising', 'U-shaped', 'Inverted U-shaped'),
                                                                 levels = c('Falling', 'Rising', 'U-shaped', 'Inverted U-shaped'))) %>% 
  pivot_longer(cols = A1C:VFC, names_to = 'region', values_to = 'count') %>%
  ggplot(., aes(x=region, y=as.numeric(count), fill = Group)) + geom_bar(position="fill", stat="identity") + 
  geom_text(data = total_counts, aes(x = region, y = 1.01, label = total_count), inherit.aes = FALSE, 
            vjust = -0.5, size = 6) +
  scale_fill_simpsons() + labs(x = '', y='', title = 'Pseudogene') +
  theme_bw(base_size = 12) + theme(legend.title=element_blank(), axis.text = element_text(color = 'black'))

##
total_counts <- protein_coding_res %>%
  as_tibble() %>%
  pivot_longer(cols = A1C:VFC, names_to = 'region', values_to = 'count') %>%
  group_by(region) %>%
  summarise(total_count = sum(as.numeric(count), na.rm = TRUE))

pc_plot <- protein_coding_res %>% as_tibble() %>% mutate(Group = factor(c('Falling', 'Rising', 'U-shaped', 'Inverted U-shaped'),
                                                                       levels = c('Falling', 'Rising', 'U-shaped', 'Inverted U-shaped'))) %>% 
  pivot_longer(cols = A1C:VFC, names_to = 'region', values_to = 'count') %>%
  ggplot(., aes(x=region, y=as.numeric(count), fill = Group)) + geom_bar(position="fill", stat="identity") + 
  geom_text(data = total_counts, aes(x = region, y = 1.01, label = total_count), inherit.aes = FALSE, 
            vjust = -0.5, size = 4) +
  scale_fill_simpsons() + labs(x = '', y='', title = 'Protein-coding gene') + 
  theme_bw(base_size = 12) + theme(legend.title=element_blank(), axis.text = element_text(color = 'black'))


p <- pgene_plot + pc_plot + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave('temporal_dynamic_expression/human_vs._macaque/macaque.pdf', p, width = 18, height = 8)
