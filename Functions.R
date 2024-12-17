## function for DESeq analysis
DEG <- function(countdata,
                meta,
                include.batch = FALSE,
                ref = ref) {
  dseq_res <- data.frame()
  All_res <- data.frame()
  
  if (include.batch) {
    cat("Including batch as covariate\n")
    design_formula <- ~ Batch + Genotype
  }
  else{
    design_formula <- ~ Genotype
  }
  
  dat2 <-
    as.matrix(countdata[, colnames(countdata) %in% rownames(meta)])
  ddsHTSeq <-
    DESeqDataSetFromMatrix(countData = dat2,
                           colData = meta,
                           design = design_formula)
  ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq)) >= 10,]
  ddsHTSeq$Genotype <- relevel(ddsHTSeq$Genotype, ref = ref)
  dds <- DESeq(ddsHTSeq, parallel = TRUE)
  
  res <- results(dds, alpha = 0.05)
  #summary(res)
  res_sub <- subset(res[order(res$padj),], padj < 0.05)
  
  res_sub$symbol <- map_function.df(res_sub, "ENSEMBL", "SYMBOL")
  res$symbol <- map_function.df(res, "ENSEMBL", "SYMBOL")
  
  res_sub$EntrezGene <-
    map_function.df(res_sub, "ENSEMBL", "ENTREZID")
  res$EntrezGene <- map_function.df(res, "ENSEMBL", "ENTREZID")
  
  dseq_res <<- as.data.frame(res_sub[, c(7:8, 1:6)])
  
  All_res <<- as.data.frame(res[, c(7:8, 1:6)])
  
}
