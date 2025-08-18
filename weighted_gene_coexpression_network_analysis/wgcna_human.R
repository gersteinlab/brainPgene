library(dplyr)
library(WGCNA)
library(stringr)

setwd('~/Desktop/brain_pgene/')
pgene_expr <- read.table('quantification/BrainSpan/pgene_expression_BrainSpan.tsv', header = T, row.names = 1, sep = '\t')
pc_expr <- read.table('quantification/BrainSpan/pc_expression_BrainSpan.tsv', header = T, row.names = 1, sep = '\t')
pgene_expr <- pgene_expr[rowSums(pgene_expr > 0.5) >= 200,] # threshold 0.5
pc_expr <- pc_expr[rowSums(pc_expr > 1) >= 200, ] # threshold 1
expr <- rbind(pgene_expr, pc_expr)

regions <- c('DFC', 'MFC', 'OFC', 'VFC', 'M1C', 'IPC', 'S1C', 'A1C', 'ITC', 'STC', 'V1C', 'AMY', 'HIP', 'STR', 'CBC', 'MD')

WGCNA_acrossRegions <- function(expression, region, species = 'HS'){
  print(region)
  expr_df <- expression %>% select(ends_with(region))
  dataExpr <- as.data.frame(t(expr_df)) # reformat
  # check genes and samples
  gsg <- goodSamplesGenes(dataExpr, verbose=3)
  if (!gsg$allOK){
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
    if(sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:",paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
    dataExpr <- dataExpr[gsg$goodSamples, gsg$goodGenes]
  }
  nGenes <- ncol(dataExpr); nSamples <- nrow(dataExpr)
  dataExpr <- log1p(dataExpr) # transfrorm to log
  sampleTree <- hclust(dist(dataExpr), method = "average")
  # Outlier dectection
  pdf(paste0('WGCNA/human/tree/', region, '.pdf'), width = 10, height = 10)
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
  dev.off()
  
  powers <- c(1:20)
  allowWGCNAThreads(nThreads = 16)
  sft <- pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5)
  net <- blockwiseModules(dataExpr, power = sft$powerEstimate, maxBlockSize = nGenes,
                          TOMType = 'signed', minModuleSize = 30, networkType="signed", 
                          mergeCutHeight = 0.25, detectCutHeight = 0.995, deepSplit = 2, 
                          minCoreKME = 0.7, minCoreKMESize=3, minKMEtoStay=0.7,
                          numericLabels = TRUE, pamRespectsDendro = TRUE,
                          saveTOMs=TRUE,
                          saveTOMFileBase = paste0('WGCNA/human/tom/', region),
                          verbose = 3, nThreads = 16)
  # WGCNA plot
  moduleLabels <- net$colors
  moduleColors <- labels2colors(moduleLabels)
  pdf(paste0('WGCNA/human/plot/', region, '.pdf'), width = 20, height = 20)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  modulesDf <- net$colors %>% as.list() %>% data.frame() %>% t() %>% as.data.frame() %>% tibble::rownames_to_column('gene')  %>%   mutate(region = region, module = paste(species, region, V1,sep = '_')) %>% select(gene, region, module)
  return(modulesDf)
}

res <- sapply(regions, function(x) WGCNA_acrossRegions(expr, x), USE.NAMES = T, simplify = F)
saveRDS(res, 'WGCNA/human_res.rds')