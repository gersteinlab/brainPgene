### ORA Enrichment
ORA_GO_analysis <- function(data, type, disease, simplify, padj_threshold = 0.05, log2FC_threshold = 0.2){
  data <- data %>% filter(padj < padj_threshold, abs(log2FoldChange) > log2FC_threshold)
  ego <- enrichGO(gene = data %>% filter(type == !!type) %>% dplyr::select(gene) %>% pull(),
                  OrgDb=org.Hs.eg.db, ont = 'BP', pAdjustMethod = 'BH',
                  minGSSize = 3, pvalueCutoff = 0.1, keyType='ENSEMBL', readable = T)
  if(simplify){
    ego <- clusterProfiler::simplify(ego) %>% as.data.frame()
  }
  else{
    ego <- ego %>% as.data.frame()
  }
  if(nrow(ego) != 0){
    ego[,'Disease'] <- disease
    ego[,'Type'] <- str_to_sentence(type)
  }
  return(ego)
}

ASD <- read.table('DEG/Capstone/pc/raw/ASD_DE.txt', sep = '\t', header = T)
BD <- read.table('DEG/Capstone/pc/raw/BD_DE.txt', sep = '\t', header = T)
SCZ <- read.table('DEG/Capstone/pc/raw/SCZ_DE.txt', sep = '\t', header = T)
GO_ASD_up <- ORA_GO_analysis(ASD, simplify = T, 'upregulated', 'ASD')
GO_ASD_down <- ORA_GO_analysis(ASD, simplify = T, 'downregulated', 'ASD')
GO_BD_up <- ORA_GO_analysis(BD, simplify = T, 'upregulated', 'BD')
GO_BD_down <- ORA_GO_analysis(BD, simplify = T, 'downregulated', 'BD')
GO_SCZ_up <- ORA_GO_analysis(SCZ, simplify = T, 'upregulated', 'SCZ')
GO_SCZ_down <- ORA_GO_analysis(SCZ, simplify = T, 'downregulated', 'SCZ')

D <- rbind(GO_ASD_up, GO_ASD_down, GO_BD_up, GO_BD_down, GO_SCZ_up, GO_SCZ_down)
D$Disease <- factor(D$Disease, levels = c('ASD', 'BD', 'SCZ'))
D_sub <-  D %>% group_by(Disease, Type) %>% slice_min(order_by = p.adjust, n = 7)
D_sub$Type <- factor(D_sub$Type, levels = c('Upregulated', 'Downregulated'))

p1 <- ggplot(D_sub, aes(x = -log10(p.adjust), y = Description)) + geom_point(aes(size = Count, color = Disease)) + geom_vline(xintercept = 1, linetype = 'dotted') + 
  xlab(TeX('$-\\log_{10}FDR$')) + ylab(NULL) + theme_bw(base_size = 12) + scale_color_npg(alpha = 0.8) + facet_grid( ~ Type) + scale_x_continuous(limits = c(0,8)) +
  theme(axis.text = element_text(color = 'black'))

ggsave('./DEG/Capstone/pc/BP_enrichment.pdf', p1, width = 12, height = 10)
