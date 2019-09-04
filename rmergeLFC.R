#!/usr/bin/env Rscript

library("stringr")

setwd("/home/longwei/Project/first/out")
ncbi2ensemble_f <- file.path(getwd(), "idconvert", "ncbiID2ensemble.csv")
ncbi2ensemble <- read.table(ncbi2ensemble_f, header = FALSE, sep = " ")
colnames(ncbi2ensemble) <- c("NCBI_ID", "ENSEMBLE_ID")

msu2ensemble_f <- file.path(getwd(), "idconvert", "msu2ensemble_pcsd.csv")
msu2ensemble <- read.table(msu2ensemble_f, header = TRUE, sep = ",")
colnames(msu2ensemble) <- c("MSU_ID", "ENSEMBLE_ID")

msu2name_f <- file.path(getwd(), "idconvert", "msu2name.csv")
msu2name <- read.delim(msu2name_f, header = TRUE)

ncbi_dir <- file.path(getwd(), "ncbi/deseq2/hisat")
ensemble_dir <- file.path(getwd(), "ensemble/deseq2/salmon")
msu_dir <- file.path(getwd(), "msu/deseq2/salmon")


ncbi_files <- file.path(ncbi_dir, grep("*.txt", list.files(ncbi_dir), value = TRUE))
ensemble_files <- file.path(ensemble_dir, grep("*.txt", list.files(ensemble_dir), value = TRUE))
msu_files <- file.path(msu_dir, grep("*.txt", list.files(msu_dir), value = TRUE))

for (nf in ncbi_files) {
  ndf <- read.table(nf, header = TRUE)
  mdf <- merge(ndf, ncbi2ensemble, by = "NCBI_ID", all.x = TRUE)
  #mdf <- merge(mdf, msu2ensemble, by = "ENSEMBLE_ID", all.x = TRUE)
  #mdf <- merge(mdf, msu2name, by = "MSU_ID", all.x = TRUE) 
  write.table(mdf[order(mdf$padj), ], file = nf, quote = F, sep = "\t", row.names = FALSE)
}

for (ef in ensemble_files) {
  edf <- read.table(ef, header = TRUE)
  edf$ENSEMBLE_ID <- sapply(edf$ENSEMBLE_ID, function(x) {
    str_replace(x, "OS", "Os")
  })
  edf$ENSEMBLE_ID <- sapply(edf$ENSEMBLE_ID, function(x) {
    str_replace(x, "G", "g")
  })
  
  mdf <- merge(edf, msu2ensemble, by = "ENSEMBLE_ID", all.x = TRUE)
  mdf <- merge(mdf, msu2name, by = "MSU_ID", all.x = TRUE)
  write.table(mdf[order(mdf$padj), ], file = ef, quote = F, sep = "\t", row.names = FALSE) 
}

for (uf in msu_files) {
  udf <- read.table(uf, header = TRUE)
  mdf <- merge(udf, msu2ensemble, by = "MSU_ID", all.x = TRUE)
  mdf <- merge(mdf, msu2name, by = "MSU_ID", all.x = TRUE)
  write.table(mdf[order(mdf$padj), ], file = uf, quote = F, sep = "\t", row.names = FALSE)
}