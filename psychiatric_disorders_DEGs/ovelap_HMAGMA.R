pgene_gtf <- import('../GTF/pgene/human/pexon_human.gtf') %>% as.data.frame()
pc_gtf <- import('../GTF/protein_coding/human/protein_coding_human.gtf') %>% as.data.frame()
# ASD
ASD_HMAGMA <- read.table('DEG/MAGMA/ASD_H-MAGMA.txt',header = T) %>% dplyr::mutate(padj = p.adjust(P, method = 'BH'))  %>%dplyr::filter(padj < 0.05)
# BD
BD_HMAGMA <- read.table('DEG/MAGMA/BD_H-MAGMA.txt', header = T) %>% dplyr::mutate(padj = p.adjust(P, method = 'BH'))  %>% dplyr::filter(padj < 0.05)
# SCZ
SCZ_HMAGMA <- read.table('DEG/MAGMA/SCZ_H-MAGMA.txt', header = T) %>% dplyr::mutate(padj = p.adjust(P, method = 'BH'))  %>% dplyr::filter(padj < 0.05)

plot_UpSet_DEG_GWAS <- function(DEG_pc, DEG_pgene, HMAGMA_df, protein_coding_gtf, pseudogene_gtf, title, padj_threshold = 0.05, log2FC_threshold = 0.2){
  pseudogene <- pseudogene_gtf %>% dplyr::select(gene_id) %>% dplyr::mutate(gene = gsub('\\..*', '', gene_id)) %>%
    distinct(gene) %>% pull()
  protein_coding <- protein_coding_gtf %>% dplyr::select(gene_id) %>%dplyr::mutate(gene = gsub('\\..*', '', gene_id)) %>%
    distinct(gene) %>% pull()
  
  DEG_1 <- DEG_pc %>% dplyr::filter(padj < padj_threshold, abs(log2FoldChange) > log2FC_threshold) %>% dplyr::select(gene) %>% pull()
  DEG_2 <- DEG_pgene %>% dplyr::filter(padj < padj_threshold) %>% dplyr::select(gene) %>% pull()
  DEG <- c(DEG_1, DEG_2)
  HMAGMA_gene <- HMAGMA_df %>% dplyr::filter(GENE %in% c(pseudogene, protein_coding)) %>% dplyr::select(GENE) %>% pull()
  
  geneList <- unique(c(HMAGMA_gene, DEG))
  
  DE <- geneList %in% DEG
  HMAGMA <- geneList %in% HMAGMA_gene
  Type <- ifelse(geneList %in% pseudogene, 'Pseudogene', 'Protein-coding gene')
  
  df <- data.frame(Gene = geneList, DE = DE, HMAGMA = HMAGMA, Type = Type)
  group <- colnames(df)[2:3]
  
  write.table(df, paste0('DEG/MAGMA/', title, '.txt'), quote = F, row.names = F, sep = '\t')
  
  p <- upset(df, group, name = '', width_ratio = 0.3,
             base_annotations = list(),
             annotations = list(
               'Count'=(ggplot(mapping=aes(fill=Type)) + 
                          geom_bar(stat='count', position='dodge') + coord_trans(y = "sqrt") +
                          geom_text(stat='count', aes(label=..count..), color='black', size = 3, position = position_dodge(width = 0.9), vjust = -1) +
                          ylab('Count') + scale_fill_igv(alpha = 1) + theme(legend.position = 'top', legend.title = element_blank())
               )
             ), 
             set_sizes = (upset_set_size(geom = geom_bar(aes(fill = Type, x = group), width = 0.2), position='left') + 
                            scale_fill_igv(alpha = 1) +  theme(axis.text.x=element_text(angle=90)) + theme(legend.position="none")), wrap = TRUE) + ggtitle(title)
  return(p)
}

p1 <- plot_UpSet_DEG_GWAS(ASD, ASD_p, ASD_HMAGMA, pc_gtf, pgene_gtf, 'ASD')
p2 <- plot_UpSet_DEG_GWAS(BD, BD_p, BD_HMAGMA, pc_gtf, pgene_gtf, 'BD')
p3 <- plot_UpSet_DEG_GWAS(SCZ, SCZ_p, SCZ_HMAGMA, pc_gtf, pgene_gtf, 'SCZ')
ggsave('DEG/MAGMA/ASD.pdf', p1, width = 8, height = 8)
ggsave('DEG/MAGMA/BD.pdf', p2, width = 8, height = 8)
ggsave('DEG/MAGMA/SCZ.pdf', p3, width = 8, height = 8)



