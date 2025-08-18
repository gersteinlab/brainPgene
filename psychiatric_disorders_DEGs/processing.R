metatable <- read.csv('metadata/Capstone/Capstone_RNAseq_meta.csv', na.strings=c('','NA'))
metatable$ageDeath <- gsub('\\+','', metatable$ageDeath)
metatable$ageDeath <- str_replace_all(metatable$ageDeath, '(\\d+)\\s+PCW', '-(1-\\1/40)')
metatable$ageDeath <- unlist(lapply(metatable$ageDeath, function(x) unlist(eval(parse(text = x)))))
metatable$ethnicity <- toupper(metatable$ethnicity)
metatable$ethnicity <- as.factor(metatable$ethnicity)
metatable$study <- as.factor(metatable$study)
metatable$diagnosis <- as.factor(metatable$diagnosis)
metatable$libraryPrep <- as.factor(metatable$libraryPrep)
metatable$sex <- as.factor(metatable$sex)
metatable$platform <- as.factor(metatable$platform)
metatable$runType <- as.factor(metatable$runType)
# select covariates
metadata <- select(metatable, c('fileName','sex', 'ageDeath', 'diagnosis', 'study', 'PMI', 'RIN', 'libraryPrep','isStranded', 'platform'))
# check the missing data
md.pattern(metadata, plot = FALSE)
metadata <- na.omit(metadata)
metadata <- mutate(metadata, ageDeathSquared = ageDeath^2, RINSquared = RIN^2)
rownames(metadata) <- unlist(lapply(metadata$fileName, function(x) unlist(strsplit(x, '.', fixed =TRUE))[1]))
metadata <- metadata[,-1]
metadata <- metadata[order(rownames(metadata)),]
write.table(metadata, 'DEG/Capstone/Capstone_metadata.tsv', sep = '\t', row.names = TRUE, quote = FALSE)