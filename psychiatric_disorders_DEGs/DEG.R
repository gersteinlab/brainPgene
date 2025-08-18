filter_exprData <- function(meta, expr){
  sampleID <- rownames(meta)
  df <- read.table(expr, sep  = '\t', header = TRUE, check.names = FALSE)
  df_filtered <- df[,c(TRUE,colnames(df[,-1]) %in% sampleID)]
  rownames(df_filtered) <- df_filtered$Gene
  df_filtered <- df_filtered[,-1]
  df_filtered <- df_filtered[ ,order(colnames(df_filtered))]
  return(df_filtered)
}

DEG_pairwiseAnalysis <- function(meta, count, fpkm, case, pgene, fpkm_threshold = 0.1, ctrl = 'Control'){
  # select samples
  meta.samples <- meta[meta$diagnosis %in% c(case,ctrl), ]
  # convert to factor
  meta.samples$diagnosis <- gsub(' ', '', meta.samples$diagnosis)
  meta.samples$diagnosis <- as.factor(meta.samples$diagnosis)
  # relevel
  meta.samples$diagnosis <- relevel(meta.samples$diagnosis, ref = ctrl)
  count.samples <- count[ ,colnames(count) %in% rownames(meta.samples)]
  fpkm.samples <- fpkm[ ,colnames(fpkm) %in% rownames(meta.samples)]
  if(pgene){
    expressed_genes <- fpkm.samples[apply(fpkm.samples, 1, function(x) quantile(x, 0.9) > fpkm_threshold),] %>% rownames()
    print(paste0(length(expressed_genes), ' genes remained.'))
    count.samples <- count.samples[rownames(count.samples) %in% expressed_genes, ]
  }
  count.samples <- as.matrix(count.samples)
  # run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = count.samples, colData = meta.samples, design = ~ diagnosis)
  dds <- estimateSizeFactors(dds)
  dat <- counts(dds, normalized=TRUE)
  idx <- rowMeans(dat) > 1
  dat <- dat[idx,]
  # SVA
  mod <- model.matrix(~ diagnosis, colData(dds))
  mod0 <- model.matrix(~ 1, colData(dds))
  svseq <- svaseq(dat, mod, mod0, n.sv = 5)
  ddssva <- dds
  ddssva$SV1 <- svseq$sv[,1]
  ddssva$SV2 <- svseq$sv[,2]
  ddssva$SV3 <- svseq$sv[,3]
  ddssva$SV4 <- svseq$sv[,4]
  ddssva$SV5 <- svseq$sv[,5]
  design(ddssva) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + diagnosis
  ddssva <- DESeq(ddssva, parallel = T)
  ddssva <- ddssva[which(mcols(ddssva)$betaConv),] # exclude genes that are not converge
  res <- results(ddssva, parallel = T, tidy = T) %>% na.omit()
  # compute ratio of case to ctrl
  case_samples <- meta.samples[meta.samples$diagnosis == gsub(' ', '', case), ] %>% rownames()
  ctrl_samples <- meta.samples[meta.samples$diagnosis == ctrl, ] %>% rownames()
  fpkm.samples.case <- fpkm.samples[ , colnames(fpkm.samples) %in% case_samples]
  fpkm.samples.ctrl <- fpkm.samples[ , colnames(fpkm.samples) %in% ctrl_samples]
  ratio <- rowMeans(fpkm.samples.case) / rowMeans(fpkm.samples.ctrl) %>% as.matrix()
  ratio <- ratio[rownames(ratio) %in% res$row, ] %>% as.data.frame()
  # reformat output
  DEG_res <- res[ ,c(1,3,6,7)]
  DEG_res <- cbind(DEG_res, ratio)
  colnames(DEG_res) <- c('gene', 'log2FoldChange', 'pvalue', 'padj', 'ratio')
  DEG_res$disease <- replicate(nrow(DEG_res), case)
  DEG_res$type <- ifelse(DEG_res$ratio > 1, 'upregulated', 'downregulated')
  return(DEG_res)
}

# pgene
meta <- read.table('DEG/Capstone/Capstone_metadata.tsv', header = TRUE, sep = '\t', row.names = 1, check.names = F)
count <- read.table('DEG/Capstone/pgene/Capstone_pgene_count.tsv', sep = '\t', header = T, row.names = 1, check.names = F)
fpkm <- read.table('DEG/Capstone/pgene/Capstone_pgene_fpkm.tsv', sep = '\t', header = T, row.names = 1, check.names = F)
ASD_p <- DEG_pairwiseAnalysis(meta = meta, count = count, fpkm = fpkm, pgene = T, case = 'Autism Spectrum Disorder')
BD_p <- DEG_pairwiseAnalysis(meta = meta, count = count, fpkm = fpkm, pgene = T, case = 'Bipolar Disorder')
SCZ_p <- DEG_pairwiseAnalysis(meta = meta, count = count, fpkm = fpkm, pgene = T, case = 'Schizophrenia')
write.table(ASD_p, 'DEG/Capstone/pgene/raw/ASD_DE.txt', sep = '\t', row.names = F, quote = F)
write.table(BD_p, 'DEG/Capstone/pgene/raw/BD_DE.txt', sep = '\t', row.names = F, quote = F)
write.table(SCZ_p, 'DEG/Capstone/pgene/raw/SCZ_DE.txt', sep = '\t', row.names = F, quote = F)

# protein-coding
meta <- read.table('DEG/Capstone/Capstone_metadata.tsv', header = TRUE, sep = '\t', row.names = 1, check.names = FALSE)
count <- read.table('DEG/Capstone/pc/Capstone_pc_count.tsv', sep = '\t', header = T, row.names = 1, check.names = FALSE)
fpkm <- read.table('DEG/Capstone/pc/Capstone_pc_fpkm.tsv', sep = '\t', header = T, row.names = 1, check.names = FALSE)
ASD <- DEG_pairwiseAnalysis(meta = meta, count = count, fpkm = fpkm, pgene = F, fpkm_threshold = 1, case = 'Autism Spectrum Disorder')
BP <- DEG_pairwiseAnalysis(meta = meta, count = count, fpkm = fpkm, pgene = F,fpkm_threshold = 1, case = 'Bipolar Disorder')
SCZ <- DEG_pairwiseAnalysis(meta = meta, count = count, fpkm = fpkm, pgene = F,fpkm_threshold = 1, case = 'Schizophrenia')
write.table(ASD, 'DEG/Capstone/pc/raw/ASD_DE.txt', sep = '\t', row.names = F, quote = F)
write.table(BP, 'DEG/Capstone/pc/raw/BD_DE.txt', sep = '\t', row.names = F, quote = F)
write.table(SCZ, 'DEG/Capstone/pc/raw/SCZ_DE.txt', sep = '\t', row.names = F, quote = F)
