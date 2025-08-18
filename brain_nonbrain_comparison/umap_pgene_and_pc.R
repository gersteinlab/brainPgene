setwd('Desktop/brain_pgene/')
library(dplyr)
t_exprMatrix <- function(file, attribute){
  data <- read.table(file, sep = '\t', header = T)
  temp <- t(data[,-1]) %>% as.data.frame()
  colnames(temp) <- data[,1]
  temp <- temp %>% mutate(Group = attribute)
  return(temp)
}

pc_entex <- t_exprMatrix('brain_vs._non-brain/data/pc_expression_ENTEx.tsv', 'Non-brain')
pgene_entex <- t_exprMatrix('brain_vs._non-brain/data/pgene_expression_ENTEx.tsv', 'Non-brain')
pc_brainspan <- t_exprMatrix('brain_vs._non-brain/data/pc_expression_BrainSpan.tsv', 'Brain')
pgene_brainspan <- t_exprMatrix('brain_vs._non-brain/data/pgene_expression_BrainSpan.tsv', 'Brain')
pc <- rbind(pc_entex, pc_brainspan)
pgene <- rbind(pgene_entex, pgene_brainspan)

set.seed(612)
pgene.umap <- umap(pgene[,-ncol(pgene)])
pgene.umap.res <- data.frame(pgene.umap$layout, Group = pgene$Group)
p1 <- ggplot(pgene.umap.res, aes(X1, X2, color = Group)) + geom_point(size = 2) +  theme_bw(base_size = 12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle('Pseudogene') + labs(x = "UMAP1", y = "UMAP2", color = "Tissue") + scale_color_lancet()

set.seed(612)
pc.umap <- umap(pc[,-ncol(pc)])
pc.umap.res <- data.frame(pc.umap$layout, Group = pc$Group)
p2 <- ggplot(pc.umap.res, aes(X1, X2, color = Group)) + geom_point(size = 2) + theme_bw(base_size = 12) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle('Protein-coding gene') + labs(x = "UMAP1", y = "UMAP2", color = "Tissue") + scale_color_lancet()

library(patchwork)
p <- p2 + p1 + plot_layout(guides = "collect", ax = 'collect')
ggsave('brain_vs._non-brain/umap_pgene_and_pc.pdf', p, width = 8, height = 4)
