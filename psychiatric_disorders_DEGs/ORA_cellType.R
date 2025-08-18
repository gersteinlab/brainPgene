library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(ggsci)
library(readxl)
ORA_cellType_analysis <- function(data, type, disease, padj_threshold = 0.05, log2FC_threshold = 0.2){
  data <- data %>% dplyr::filter(type == !!type, padj < padj_threshold, abs(log2FoldChange) > log2FC_threshold)
  mappings <- bitr(data$gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  entrezID <- mappings %>% na.omit() %>% dplyr::select(ENTREZID) %>% dplyr::pull()
  cell_marker_data <- read_excel("Cell_marker_Human.xlsx")
  cells <- cell_marker_data %>% dplyr::select(cell_name, GeneID) %>% dplyr::mutate(geneID = strsplit(GeneID, ', ')) %>% tidyr::unnest(c(geneID))
  res <- enricher(entrezID, TERM2GENE = cells, pvalueCutoff = 0.1)
  res <- res@result %>% na.omit() %>% dplyr::filter(p.adjust < 0.1)
  res[,'Disease'] <- disease
  res[,'Type'] <- str_to_sentence(type)
  return(res)
}

cellType_ASD_up <- ORA_cellType_analysis(ASD, 'upregulated', 'ASD')
cellType_ASD_down <- ORA_cellType_analysis(ASD, 'downregulated', 'ASD')
cellType_BD_up <- ORA_cellType_analysis(BD, 'upregulated', 'BP')
cellType_BD_down <- ORA_cellType_analysis(BD, 'downregulated', 'BP')
cellType_SCZ_up <- ORA_cellType_analysis(SCZ, 'upregulated', 'SCZ')
cellType_SCZ_down <- ORA_cellType_analysis(SCZ, 'downregulated', 'SCZ')

D <- rbind(cellType_ASD_up, cellType_ASD_down, cellType_BD_up, cellType_BD_down, cellType_SCZ_up, cellType_SCZ_down)
D$Type <- factor(D$Type, levels = c('Upregulated', 'Downregulated'))
D_sub <-  D %>% group_by(Disease, Type) %>% slice_min(order_by = p.adjust, n = 7)

p2 <- ggplot(D_sub, aes(y = Disease, x = Description, fill = -log10(p.adjust))) + 
  geom_tile(lwd = 1.5,linetype = 1) + facet_wrap( ~ Type, nrow =2) +scale_fill_gsea(alpha = 1) + 
  theme_bw(base_size = 14) + coord_fixed() + xlab(NULL) + ylab('Brain disorder') +
  guides(fill = guide_colourbar(title = TeX('$-\\log_{10}FDR$'), barwidth = 1, barheight = 10)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.text = element_text(color = 'black'))
ggsave('./DEG/Capstone/pc/cellType.pdf', p2, width = 10, height = 5)
