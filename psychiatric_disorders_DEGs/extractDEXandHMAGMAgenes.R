ASD_pgene <- read.table('DEG_across_disorders/DEG/Capstone/pgene/raw/ASD_DE.txt', sep = '\t', header = T) %>% filter(padj < 0.05, abs(log2FoldChange) > 0)
ASD_pcg <- read.table('DEG_across_disorders/DEG/Capstone/pc/raw/ASD_DE.txt', sep = '\t', header = T) %>% filter(padj < 0.05, abs(log2FoldChange) > 0.2)
BD_pgene <- read.table('DEG_across_disorders/DEG/Capstone/pgene/raw/BD_DE.txt', sep = '\t', header = T) %>% filter(padj < 0.05, abs(log2FoldChange) > 0)
BD_pcg <- read.table('DEG_across_disorders/DEG/Capstone/pc/raw/BD_DE.txt', sep = '\t', header = T) %>% filter(padj < 0.05, abs(log2FoldChange) > 0.2)
SCZ_pgene <- read.table('DEG_across_disorders/DEG/Capstone/pgene/raw/SCZ_DE.txt', sep = '\t', header = T) %>% filter(padj < 0.05, abs(log2FoldChange) > 0)
SCZ_pcg <- read.table('DEG_across_disorders/DEG/Capstone/pc/raw/SCZ_DE.txt', sep = '\t', header = T) %>% filter(padj < 0.05, abs(log2FoldChange) > 0.2)

library(dplyr)
df <- rbind(ASD_pgene, ASD_pcg, BD_pgene, BD_pcg, SCZ_pgene, SCZ_pcg) %>% select(gene, log2FoldChange, padj, disease) %>% arrange(disease)
write.table(df, 'DEG_across_disorders/DEXgenes_across_disorders.txt', sep = '\t', quote = F, row.names = F)



ASD <- read.table('DEG_across_disorders/DEG/MAGMA/ASD_H-MAGMA.txt', sep = '\t', header = T) %>% mutate(DISEASE = 'ASD')
BD <- read.table('DEG_across_disorders/DEG/MAGMA/BD_H-MAGMA.txt', sep = '\t', header = T) %>% mutate(DISEASE = 'BD')
SCZ <- read.table('DEG_across_disorders/DEG/MAGMA/SCZ_H-MAGMA.txt', sep = '\t', header = T) %>% mutate(DISEASE = 'SCZ')

HMAGMA <- rbind(ASD, BD, SCZ)
write.table(df, 'DEG_across_disorders/HMAGMA.txt', sep = '\t', quote = F, row.names = F)
