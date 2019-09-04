#!/usr/bin/env Rscript
# Author:long wei

library("tximport")
library("stringr")
library("DESeq2")
library("vsn")
library("pheatmap")
library("RColorBrewer")
library("BiocParallel")

register(MulticoreParam(8))

#init directory
home <- Sys.getenv("HOME")
project <- file.path(home, "Project/first")
ds <- "msu" #for rice annotation data source, change ds value.
countdir <- file.path(project, "out", ds, "count") #location of my htseq counts file
outdir <- file.path(project, "out", ds, "deseq2/hisat") #location of output of result

#construct sample info
sampleFiles <- grep("*.count", list.files(countdir), value = TRUE)
sampleNames <- str_replace(sampleFiles, ".count", "") #(HP15_0d_rep1)
sampleGroup <- str_replace(sampleNames, "(HP\\d+)_(\\d+d)_(rep\\d+)", "\\1")  #(HP15,HP16,HP21)
sampleTime <- str_replace(sampleNames, "(HP\\d+)_(\\d+d)_(rep\\d+)", "\\2")   #(0d,6d,12d)
sampleRep <- str_replace(sampleNames, "(HP\\d+)_(\\d+d)_(rep\\d+)", "\\3")    #(rep1, rep2)

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          sampleGroup = sampleGroup,
                          sampleTime = sampleTime,
                          sampleRep = sampleRep)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = countdir,
                                       design = ~ sampleGroup)

dds <- DESeq(ddsHTSeq, parallel = TRUE)

rld <- rlogTransformation(dds, blind = TRUE)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
ntd <- normTransform(dds)

#plot effects of transformations on the variance
pdf(file.path(outdir, "meanSd.pdf"))
meanSdPlot(assay(ntd), ylab = "sd(ntd)")
meanSdPlot(assay(vsd), ylab = "sd(vsd)")
meanSdPlot(assay(rld), ylab = "sd(rld)")
dev.off()
#end of transformation variance

#heatmap of sample distance 
pdf(file.path(outdir, "distance.pdf"))
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- sampleNames
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, main = "sample to sample distance heatmap")
dev.off()

#heatmap of sample count matrix
pdf(file.path(outdir, "countheat.pdf"))
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("sampleGroup","sampleTime")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main = "ntd count heatmap of top 20")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main = "vsd count heatmap of top 20")
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main = "rld count heatmap of top 20")
dev.off()
#end of sample count matrix heatmap

#plot PCA
pdf(file.path(outdir, "pca.pdf"))
plotPCA(vsd, intgroup=c("sampleGroup"))
plotPCA(vsd, intgroup=c("sampleGroup", "sampleTime"))
plotPCA(vsd, intgroup=c("sampleGroup", "sampleTime", "sampleRep"))
dev.off()
#end of plot PCA

#plot outliers
pdf(file.path(outdir, "outliers.pdf"))
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()
#end of plot outliers

#plot dispersion
pdf(file.path(outdir, "dispersion.pdf"))
plotDispEsts(dds)
dev.off()
#end of plot dispersion

#write diff expression result
writeResult <- function(result, lfcresult, subject) {
  lfcresult <- lfcresult[order(lfcresult$padj), ]
  write.table(as.data.frame(lfcresult), 
              file = file.path(outdir, paste0(subject, ".txt")), 
              quote = F, sep = "\t")
  pdf(file.path(outdir, paste0(subject, ".pdf")))
  par(mfrow=c(2,2))
  plotMA(result, ylim=c(-6, 6), main=subject)
  plotMA(lfcresult, ylim=c(-6, 6), main=paste0(subject, "_shrink"))
  par(mfrow=c(1,1))
  dev.off()
}
#end of writeResult

#contrast by sample group
res_15vs21 <- results(dds, contrast = c("sampleGroup", "HP15", "HP21"))
reslfc_15vs21 <- lfcShrink(dds,contrast = c("sampleGroup", "HP15", "HP21"), parallel = TRUE)
#reslfc_15vs21 <- reslfc_15vs21[order(reslfc_15vs21$padj),]

res_16vs21 <- results(dds, contrast = c("sampleGroup", "HP16", "HP21"))
reslfc_16vs21 <- lfcShrink(dds, contrast = c("sampleGroup", "HP16", "HP21"), parallel = TRUE)

#contrast by sample group & time(multi factor)
ddsMF <- dds
ddsMF$groupTime <- factor(paste0(ddsMF$sampleGroup, "_", ddsMF$sampleTime))
design(ddsMF) <- ~ groupTime
gtdds <- DESeq(ddsMF, parallel = TRUE)

res_15vs21_0d <- results(gtdds, contrast = c("groupTime", "HP15_0d", "HP21_0d"), parallel = TRUE)
reslfc_15vs21_0d <- lfcShrink(gtdds, contrast = c("groupTime", "HP15_0d", "HP21_0d"), parallel = TRUE)

res_15vs21_6d <- results(gtdds, contrast = c("groupTime", "HP15_6d", "HP21_6d"), parallel = TRUE)
reslfc_15vs21_6d <- lfcShrink(gtdds, contrast = c("groupTime", "HP15_6d", "HP21_6d"), parallel = TRUE)

res_15vs21_12d <- results(gtdds, contrast = c("groupTime", "HP15_12d", "HP21_12d"), parallel = TRUE)
reslfc_15vs21_12d <- lfcShrink(gtdds, contrast = c("groupTime", "HP15_12d", "HP21_12d"), parallel = TRUE)

res_16vs21_0d <- results(gtdds, contrast = c("groupTime", "HP16_0d", "HP21_0d"), parallel = TRUE)
reslfc_16vs21_0d <- lfcShrink(gtdds, contrast = c("groupTime", "HP16_0d", "HP21_0d"), parallel = TRUE)

res_16vs21_6d <- results(gtdds, contrast = c("groupTime", "HP16_6d", "HP21_6d"), parallel = TRUE)
reslfc_16vs21_6d <- lfcShrink(gtdds, contrast = c("groupTime", "HP16_6d", "HP21_6d"), parallel = TRUE)

res_16vs21_12d <- results(gtdds, contrast = c("groupTime", "HP16_12d", "HP21_12d"), parallel = TRUE)
reslfc_16vs21_12d <- lfcShrink(gtdds, contrast = c("groupTime", "HP16_12d", "HP21_12d"), parallel = TRUE)

writeResult(res_15vs21, reslfc_15vs21, "HP15vsHP21")
writeResult(res_16vs21, reslfc_16vs21, "HP16vsHP21")
writeResult(res_15vs21_0d, reslfc_15vs21_0d, "HP15vsHP21_0d")
writeResult(res_15vs21_6d, reslfc_15vs21_6d, "HP15vsHP21_6d")
writeResult(res_15vs21_12d, reslfc_15vs21_12d, "HP15vsHP21_12d")
writeResult(res_16vs21_0d, reslfc_16vs21_0d, "HP16vsHP21_0d")
writeResult(res_16vs21_6d, reslfc_16vs21_6d, "HP16vsHP21_6d")
writeResult(res_16vs21_12d, reslfc_16vs21_12d, "HP16vsHP21_12d")
