rm(list = ls())
# human
hs_modules <- readRDS('WGCNA/human_res.rds')
hs_modules %>% bind_rows() %>% group_by(region) %>% summarise(num_modules = n_distinct(module)) %>% summarise(avg = mean(num_modules))
colors_16 <- paletteer_d("ggthemes::Tableau_20")[1:16]
pgene_expr <- read.table('quantification/BrainSpan/pgene_expression_BrainSpan.tsv', header = T, row.names = 1, sep = '\t')

## number
lapply(c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD'), 
       function(x) hs_modules[[x]] %>% mutate(type = ifelse(hs_modules[[x]][['gene']] %in% rownames(pgene_expr), 'Pseudogene', 'Protein-coding')) %>%
         group_by(module, type) %>% filter(type == 'Pseudogene', !grepl("_0", module))) %>% bind_rows() %>% pull(gene) %>% unique()

p1 <- lapply(c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD'), 
             function(x) hs_modules[[x]] %>% mutate(type = ifelse(hs_modules[[x]][['gene']] %in% rownames(pgene_expr), 'Pseudogene', 'Protein-coding')) %>%
               group_by(module, type) %>% summarise(count = n()) %>% mutate(proportion= prop.table(count), region = x) %>% 
               filter(type == 'Pseudogene', !grepl("_0", module)) %>% as.data.frame())  %>% bind_rows() %>%
  ggplot(aes(x = factor(region, levels = c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD')), y = proportion, 
             fill = factor(region, levels = c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD')))) + 
  geom_boxplot(alpha = 1, outlier.alpha = 1)  + labs(x = '', y = 'Proportion of pseudogenes in modules', title = 'Human') + scale_y_continuous(trans = "sqrt", breaks = scales::trans_breaks("sqrt", function(x) x^2, n = 5)) +
  scale_fill_manual(values = colors_16) +
  theme_bw(base_size = 12) + theme(legend.position='none', axis.text = element_text(color = 'black'))
# macaque
mm_modules <- readRDS('WGCNA/macaque_res.rds')
mm_modules %>% bind_rows() %>% group_by(region) %>% summarise(num_modules = n_distinct(module)) %>% summarise(avg = mean(num_modules))
pgene_expr <- read.table('quantification/NHP_macaque/pgene_expression_NHP_macaque.tsv', header = T, row.names = 1, sep = '\t')

#  num
lapply(c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD'), 
       function(x) mm_modules[[x]] %>% mutate(type = ifelse(mm_modules[[x]][['gene']] %in% rownames(pgene_expr), 'Pseudogene', 'Protein-coding')) %>%
         group_by(module, type) %>% filter(type == 'Pseudogene', !grepl("_0", module))) %>% bind_rows() %>% pull(gene) %>% unique()


p2 <- lapply(c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD'), 
             function(x) mm_modules[[x]] %>% mutate(type = ifelse(mm_modules[[x]][['gene']] %in% rownames(pgene_expr), 'Pseudogene', 'Protein-coding')) %>%
               group_by(module, type) %>% summarise(count = n()) %>% mutate(proportion= prop.table(count), region = x) %>% 
               filter(type == 'Pseudogene', !grepl("_0", module)) %>% as.data.frame())  %>% bind_rows() %>%
  ggplot(aes(x = factor(region, levels = c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD')), 
             y = proportion, fill = factor(region, levels = c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD')))) + 
  geom_boxplot(alpha = 1, outlier.alpha = 1)  + labs(x = '', y = 'Proportion of pseudogenes in modules', title = 'Macaque') + scale_y_continuous(trans = "sqrt", breaks = scales::trans_breaks("sqrt", function(x) x^2, n = 5)) +
  scale_fill_manual(values = colors_16) +
  theme_bw(base_size = 12) + theme(legend.position='none', axis.text = element_text(color = 'black'))

library(patchwork)
ggsave('WGCNA/pgene_proportion_in_modules_human.pdf', p1, width = 10, height = 4)
ggsave('WGCNA/pgene_proportion_in_modules_macaque.pdf', p2, width = 10, height = 4)
