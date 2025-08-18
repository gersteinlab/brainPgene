library(circlize)
library(RColorBrewer)
# define the color
grid.col = c('Elevated in at least one region' = '#f2e8db', 
             'Low region-specificity' = '#bf8552',
             'Not detected' = '#5e3023',
             'Elevated in brain' = '#b5dfcc',
             'Elevated in non-brain' = '#fdcead',
             'Low tissue-specificity' = '#ccd5e8')
# get the legend
leg <- data.frame(x = c('Elevated in at least one region', 'Low region-specificity','Not detected', 'Elevated in brain', 'Elevated in non-brain','Low tissue-specificity'), y = 1) %>% ggplot(aes(as.factor(x), y, color = x)) + 
  geom_point(size = 5, pch = 15, alpha = 0.8) + scale_color_brewer(palette = 'Set3') + theme_minimal() +  theme(legend.title=element_blank())
leg <- get_legend(leg)
leg <- as_ggplot(leg)
ggsave('brain_vs._non-brain/chrod_diagram/legend.pdf', leg, height = 5, width = 5)

## draw the chrod diagram
plot_ChrodDiagram <- function(tissueSpecDf, regionSpecDf){
  combination <- expand.grid(levels(as.factor(tissueSpecDf$subtype)), levels(as.factor(regionSpecDf$subtype)))
  count <- c()
  for (i in 1:nrow(combination)) {
    gene_1 <- tissueSpecDf %>% dplyr::filter(subtype == combination[i,1]) %>% pull(row)
    gene_2 <- regionSpecDf %>% dplyr::filter(subtype == combination[i,2]) %>% pull(row)
    count <- c(count, length(intersect(gene_1, gene_2)))
  }
  df <- data.frame(from = combination[,1], to = combination[,2], value = count)
  chordDiagram(df, grid.col = grid.col, transparency = 0.2, big.gap = 10)
  circos.clear()
}

# data
pgene_tissueSpecDf <- read.table('brain_vs._non-brain/tissue_specificity/pgene_tissue_specificity.txt', sep = '\t', header = T)
parent_tissueSpecDf <-read.table('brain_vs._non-brain/tissue_specificity/parent_tissue_specificity.txt', sep = '\t', header = T)
pc_tissueSpecDf <- read.table('brain_vs._non-brain/tissue_specificity/pc_tissue_specificity.txt', sep = '\t', header = T)

pgene_regionSpecDf <- read.table('brain_vs._non-brain/region_specificity/pgene_region_specificity.txt', sep = '\t', header = T)
parent_regionSpecDf <-read.table('brain_vs._non-brain/region_specificity/parent_region_specificity.txt', sep = '\t', header = T)
pc_regionSpecDf <- read.table('brain_vs._non-brain/region_specificity/pc_region_specificity.txt', sep = '\t', header = T)

pdf('brain_vs._non-brain/chrod_diagram/pgene.pdf', width = 8, height = 8)
plot_ChrodDiagram(pgene_tissueSpecDf, pgene_regionSpecDf) ## dataframes from pgene_type_by_specificity.R
dev.off()

pdf('brain_vs._non-brain/chrod_diagram/parent.pdf', width = 8, height = 8)
plot_ChrodDiagram(parent_tissueSpecDf, parent_regionSpecDf) ## dataframes from pc_and_parent_by_specificity.R
dev.off()

pdf('brain_vs._non-brain/chrod_diagram/protein_coding.pdf', width = 8, height = 8)
plot_ChrodDiagram(pc_tissueSpecDf, pc_regionSpecDf)  ## dataframes from pc_and_parent_by_specificity.R
dev.off()