# Protein-coding
## tissue-specificity
pc_tissueSpec <- readRDS('brain_vs._non-brain/tissue_specificity/DEG_Protein-coding_Gene.rds')
pc_tissueSpecDf <- data.frame(row = c(pc_tissueSpec[[1]], pc_tissueSpec[[2]], pc_tissueSpec[[3]]),
                              subtype = c(rep('Elevated in brain', length(pc_tissueSpec[[1]])),
                                          rep('Elevated in non-brain', length(pc_tissueSpec[[2]])),
                                          rep('Low tissue-specificity', length(pc_tissueSpec[[3]]))),
                              group = 'Tissue-specificity')
write.table(pc_tissueSpecDf, 'brain_vs._non-brain/tissue_specificity/pc_tissue_specificity.txt', quote = F, sep = '\t', row.names = F)

## region-specificity
pc_regionSpec <- readRDS('brain_vs._non-brain/region_specificity/DEG_Protein-coding_Gene.rds')
pc_regionSpecDf <- data.frame(row = c(pc_regionSpec[[1]], pc_regionSpec[[2]], pc_regionSpec[[3]]),
                              subtype = c(rep('Elevated in at least one region', length(pc_regionSpec[[1]])),
                                          rep('Low region-specificity', length(pc_regionSpec[[2]])),
                                          rep('Not detected', length(pc_regionSpec[[3]]))),
                              group = 'Region-specificity')
write.table(pc_regionSpecDf, 'brain_vs._non-brain/region_specificity/pc_region_specificity.txt', quote = F, sep = '\t', row.names = F)


# Parent 
## tissue-specificity
parent_tissueSpec <- readRDS('brain_vs._non-brain/tissue_specificity/DEG_Parent_Gene.rds')
parent_tissueSpecDf <- data.frame(row = c(parent_tissueSpec[[1]], parent_tissueSpec[[2]], parent_tissueSpec[[3]]),
                                  subtype = c(rep('Elevated in brain', length(parent_tissueSpec[[1]])),
                                              rep('Elevated in non-brain', length(parent_tissueSpec[[2]])),
                                              rep('Low tissue-specificity', length(parent_tissueSpec[[3]]))),
                                  group = 'Tissue-specificity')
write.table(parent_tissueSpecDf, 'brain_vs._non-brain/tissue_specificity/parent_tissue_specificity.txt', quote = F, sep = '\t', row.names = F)

## region-specificity
parent_regionSpec <- readRDS('brain_vs._non-brain/region_specificity/DEG_Parent_Gene.rds')
parent_regionSpecDf <- data.frame(row = c(parent_regionSpec[[1]], parent_regionSpec[[2]], parent_regionSpec[[3]]),
                                  subtype = c(rep('Elevated in at least one region', length(parent_regionSpec[[1]])),
                                              rep('Low region-specificity', length(parent_regionSpec[[2]])),
                                              rep('Not detected', length(parent_regionSpec[[3]]))),
                                  group = 'Region-specificity')
write.table(parent_regionSpecDf, 'brain_vs._non-brain/region_specificity/parent_region_specificity.txt', quote = F, sep = '\t', row.names = F)
