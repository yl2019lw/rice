#!/usr/bin/env Rscript

library("stringr")

setwd("/home/longwei/Project/first/out")

lfcs <- c("HP15vsHP21.txt","HP15vsHP21_0d.txt", "HP15vsHP21_6d.txt", "HP15vsHP21_12d.txt",
          "HP16vsHP21.txt", "HP16vsHP21_0d.txt", "HP16vsHP21_6d.txt", "HP16vsHP21_12d.txt")
                                     
ncbi_dir <- file.path(getwd(), "ncbi/deseq2/hisat")
ensemble_dir <- file.path(getwd(), "ensemble/deseq2/salmon")
msu_dir <- file.path(getwd(), "msu/deseq2/salmon")

for (lfc in lfcs) {
  nf <- file.path(ncbi_dir, lfc)
  ef <- file.path(ensemble_dir, lfc)
  uf <- file.path(msu_dir, lfc)
  
  ndf <- read.delim(nf, header = TRUE)
  edf <- read.delim(ef, header = TRUE)
  udf <- read.delim(uf, header = TRUE)
  
  ndf1 <- subset(ndf, padj < 0.01, select = c("NCBI_ID", "ENSEMBLE_ID"))
  edf1 <- subset(edf, padj < 0.01, select = c("ENSEMBLE_ID", "padj", "log2FoldChange"))
  #edf1$ENSEMBLE_ID <- sapply(edf1$ENSEMBLE_ID, function(x) {
  #  str_replace(x, "OS", "Os")
  #})
  #edf1$ENSEMBLE_ID <- sapply(edf1$ENSEMBLE_ID, function(x) {
  #  str_replace(x, "G", "g")
  #})

  #add base mean in 2018-1-22
  udf1 <- subset(udf, padj < 0.01, select = c("ENSEMBLE_ID", "MSU_ID", "MSU_NAME", "baseMean"))
  
  mdf <- merge(ndf1, edf1, by = "ENSEMBLE_ID")
  mdf <- merge(udf1, mdf, by = "ENSEMBLE_ID")
  
  outf <- file.path(getwd(), "mergeResult", lfc)
  write.table(mdf[order(mdf$padj), ], file = outf, quote = F, sep = "\t", row.names = FALSE)
  
  pdfpath <- file.path(getwd(), "mergeResult", str_replace(lfc, "txt", "pdf"))
  pdf(pdfpath)
  plot(mdf$padj, 0 - mdf$log2FoldChange, xlab = "adjusted p value", ylab = "- log2FoldChange")
  dev.off()
}