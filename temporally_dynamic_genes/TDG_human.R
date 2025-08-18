rm(list = ls())
# metadata
specimen_meta <- read.table('metadata/BrainSpan/meta_RNAseq_bam_BrainSpan.csv', sep = ',', header = T)
edesign_table <- read.table('metadata/BrainSpan/edesign_RNAseq_bam_BrainSpan.txt', sep = '\t', header = T)
edesign_table <- edesign_table[edesign_table$specimenID %in% specimen_meta$specimenID,]
tissue_meta <- read.table('metadata/BrainSpan/tissue_group.csv', sep = ',', header = T)
tissue_meta <- tissue_meta %>% filter(sampleNum > 10)

# temporal dynamic genes across tissues
assignTemporalGeneType <- function(df, min, max){
  df <- df %>% mutate(axis = ifelse(betaTime2 == 0, NA, -betaTime/(2*betaTime2)))
  types <- c()
  for (i in 1:nrow(df)) {
    if(is.na(df[i,4])){
      if(df[i,2] > 0){type <- 'Rising'}
      if(df[i,2] < 0){type <- 'Falling'}
    }
    else{
      if(df[i,3] > 0){
        if(df[i,4] <= min){type <- 'Rising'}
        if(df[i,4] >= max){type <- 'Falling'}
        if(df[i,4] > min & df[i,4] < max){type <- 'U-shaped'}
      }
      if(df[i,3] < 0){
        if(df[i,4] <= min){type <- 'Falling'}
        if(df[i,4] >= max){type <- 'Rising'}
        if(df[i,4] > min & df[i,4] < max){type <- 'Inverted U-shaped'}
      }
    }
    types <- c(types,type)
  }
  df <- df %>% mutate(type = types)
  return(df)
}

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
  # assign temporal dynamic gene type based on degree (rising, falling, rising then falling, falling then rising)
  if(degree == 2){
    res <- res %>% mutate(axis = ifelse(betaTime2 == 0, NA, -betaTime/(2*betaTime2)))
    min <- min(edesign$Time); max <- max(edesign$Time)
    df <- assignTemporalGeneType(res, min, max)
  }else if(degree == 1){
    df <- res %>% mutate(type = ifelse(betaTime > 0, 'Rising', 'Falling'))
  }
  return(df)
}

# pseudogene
pgene_expression <- read.table('quantification/BrainSpan/pgene_expression_BrainSpan.tsv', sep = '\t', header = T, row.names = 1)
pgene <- sapply(tissue_meta$tissue, function(x) TDG_acrossTissue(pgene_expression, edesign_table, tissue_meta, x), simplify = FALSE, USE.NAMES = TRUE)
pgene_res <- sapply(names(pgene), function(x) table(factor(pgene[[x]]$type, levels = c('Falling', 'Rising', 'U-shaped', 'Inverted U-shaped'))), simplify = F, USE.NAMES = T)
saveRDS(pgene, 'temporal_dynamic_expression/human/pgene_temporallyDynamicGene.rds')

# protein coding
pc_expression <- read.table('quantification/BrainSpan/pc_expression_BrainSpan.tsv', sep = '\t', header = T, row.names = 1)
protein_coding <- sapply(tissue_meta$tissue, function(x) TDG_acrossTissue(pc_expression, edesign_table, tissue_meta, x), simplify = FALSE, USE.NAMES = TRUE)
protein_coding_res <- sapply(names(protein_coding), function(x) table(factor(protein_coding[[x]]$type, levels = c('Falling', 'Rising', 'U-shaped', 'Inverted U-shaped'))), simplify = F, USE.NAMES = T)
saveRDS(protein_coding, 'temporal_dynamic_expression/human/protein_coding_temporallyDynamicGene.rds')

# visualization
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

#
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
ggsave('temporal_dynamic_expression/human_vs._macaque/human.pdf', p, width = 18, height = 8)
