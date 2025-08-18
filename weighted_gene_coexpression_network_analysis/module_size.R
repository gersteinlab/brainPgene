hs_modules <- readRDS('WGCNA/human_res.rds') %>% bind_rows()
mm_modules <- readRDS('WGCNA/macaque_res.rds')   %>% bind_rows()
colors_16 <- paletteer_d("ggthemes::Tableau_20")[1:16]

p1 <- hs_modules %>% group_by(region,module) %>% summarise(n = n()) %>% filter(!grepl("_0", module)) %>%
  ggplot(aes(x = n, 
             y = factor(region, levels = c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD')), 
             fill = factor(region, levels = c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD')) )) + 
  geom_density_ridges2(stat = "binline", bins = 30, scale = 0.95, draw_baseline = TRUE) + labs(x = 'Module size', y = '', title = 'Human') + scale_fill_manual(values = colors_16) + 
  scale_x_continuous(trans = "sqrt", breaks = scales::trans_breaks("sqrt", function(x) x^2, n = 5)) +   coord_cartesian(clip = "on") +
  theme_bw(base_size = 12) +theme(legend.position='none', axis.text = element_text(color = 'black'))

p2 <- mm_modules %>% group_by(region,module) %>% summarise(n = n()) %>% filter(!grepl("_0", module)) %>%
  ggplot(aes(x = n, 
             y = factor(region, levels = c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD')), 
             fill = factor(region, levels = c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD')) )) + 
  geom_density_ridges2(stat = "binline", bins = 30, scale = 0.95, draw_baseline = TRUE) + labs(x = 'Module size', y = '', title = 'Macaque') + scale_fill_manual(values = colors_16) + 
  scale_x_continuous(trans = "sqrt", breaks = scales::trans_breaks("sqrt", function(x) x^2, n = 5)) + coord_cartesian(clip = "on") +
  theme_bw(base_size = 12) +theme(legend.position='none', axis.text = element_text(color = 'black'))

ggsave('WGCNA/module_size_human.pdf', p1, width = 10, height = 6)
ggsave('WGCNA/module_size_macaque.pdf', p2, width = 10, height = 6)

