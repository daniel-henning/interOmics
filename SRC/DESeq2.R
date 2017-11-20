## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#setwd()
## ----loadDESeq2, echo=FALSE----------------------------------------------
library(DESeq2)
library(ggplot2)
library(dplyr)

htseq.count<- read.table("GSE94655_htseq_count_table.txt",header=TRUE,row.names=1,check.names = FALSE)

exprSet <- as.matrix(htseq.count)

group <- factor(c('Kc','Fb','Mc','Mc','Kc','skin','Mc','skin','Mc','Fb','Kc','Kc'))
colData <- data.frame(row.names=colnames(exprSet), condition=group)
#df=data.frame(Treatment = group)
dds <- DESeqDataSetFromMatrix(countData = exprSet,colData = colData,design = ~ condition)
#pre-filter low count genes before running the DESeq2
dds <- dds[rowSums(counts(dds)) > 4,]
#to define which comparison to make using the contrast argument
dds$condition <- relevel(dds$condition, ref="skin")
nrow(dds)

ds2 <- DESeq(dds) resultsNames(dds2)
#es <-  results(dds2, contrast=c("Tcondition","Mc,"cskin))
resOrdered <- res[order(res$pvalue),]
resOrdered =  as.data.frame(resOrdered)
write.csv(resOrdered,'Mc-skin.csv')


############################
prSet_new=assay(rld)
par(cex = 0.7)
n.sample=ncol(exprSet)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
boxplot(exprSet, col = cols,main="expression value",las=2)
boxplot(exprSet_new, col = cols,main="expression value",las=2)
hist(exprSet)
hist(exprSet_new)
## ----parallel, eval=FALSE------------------------------------------------
library("BiocParallel")
register(SnowParam(6))
## ----resOrder------------------------------------------------------------
resOrdered <- res[order(res$padj),]
## ----sumRes--------------------------------------------------------------
summary(res)
## ----sumRes01------------------------------------------------------------
sum(res$padj < 0.1,na.rm=TRUE)
## ----resAlpha05----------------------------------------------------------
res05 <- results(dds2, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
## ----IHW-----------------------------------------------------------------
library("IW")resIHW <- results(dds2, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult
## ----MA, fig.width=4.5, fig.height=4.5-----------------------------------
plotMA(res, main="DESeq2", ylim=c(-2,2))
## ----resMLE--------------------------------------------------------------
resMLE <- results(dds2, addMLE=TRUE)
head(resMLE, 4)
## ----MANoPrior, fig.width=4.5, fig.height=4.5----------------------------
plotMA(resMLE, MLE=TRUE, main="unshrunken LFC", ylim=c(-2,2))

## ----plotCounts, dev="pdf", fig.width=4.5, fig.height=5------------------
plotCounts(dds2, gene=which.min(res$padj), intgroup="grgroup_list")

## ----plotCountsAdv, dev="pdf", fig.width=3.5, fig.height=3.5-------------
d <- plotCounts(dds2, gene=which.min(res$padj), intgroup="group", 
                returnData=TRUE)
library("ggplot2")
ggplot(dds, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

## ----metadata------------------------------------------------------------
mcols(res)$description

## ----export, eval=FALSE--------------------------------------------------
## write.csv(as.data.frame(resOrdered),
##           file="condition_treated_results.csv")
## ----subset--------------------------------------------------------------
resSig <- subset(resOrdered, padj < 0.1)

resSig

## ----multifactor---------------------------------------------------------
colData(dds)

## ----copyMultifactor-----------------------------------------------------
ddsMF <- dds

## ----replaceDesign-------------------------------------------------------
design(ddsMF) <- formula(~ type + group)
ddsMF <- DESeq(ddsMF)

## ----multiResults--------------------------------------------------------
resMF <- results(ddsMF)
head(resMF)

## ----multiTypeResults----------------------------------------------------
resMFType <- results(ddsMF,contrast=c("group"))
head(resMFType)

## ----rlogAndVST----------------------------------------------------------
rld <- rlog(dds2, blind=FALSE)


vsd <- varianceStabilizingTransformation(dds2, blind=FALSE)
vsd.fast <- vst(dds2, blind=FALSE)
head(assay(rld), 3)

## ----vsd1, echo=FALSE, fig.width=4.5, fig.height=4.5, fig.show="asis", fig.small=TRUE, fig.pos="!bt", fig.cap="VST and log2. Graphs of the variance stabilizing transformation for sample 1, in blue, and of the transformation $f(n) = \\log_2(n/s_1)$, in black. $n$ are the counts and $s_1$ is the size factor for the first sample.\\label{figure/vsd1-1}"----
px     <- counts(dds2)[,1] / sizeFactors(dds2)[1]
ord    <- order(px)
ord    <- ord[px[ord] < 150]
ord    <- ord[seq(1, length(ord), length=50)]
last   <- ord[length(ord)]
vstcol <- c("blue", "black")
matplot(px[ord],
        cbind(assay(vsd)[, 1], log2(px))[ord, ],
        type="l", lty=1, col=vstcol, xlab="n", ylab="f(n)")
legend("bottomright",
       legend = c(
         expression("variance stabilizing transformation"),
         expression(log[2](n/s[1]))),
       fill=vstcol)

## ----meansd, fig.width=4, fig.height=3, fig.show="asis", fig.wide=TRUE, fig.pos="tb", out.width=".32\\linewidth", fig.cap="Per-gene standard deviation (taken across samples), against the rank of the mean. {\\bfhelvet(a)} for the shifted logarithm $\\log_2(n+1)$, the regularized log transformation {\\bfhelvet(b)} and the variance stabilizing transformation {\\bfhelvet(c)}.\\label{fig:meansd}", fig.subcap=""----
library("vsn")
notAllZero <- (rowSums(counts(dds2))>0)
meanSdPlot(log2(counts(dds2,normalized=TRUE)[notAllZero,] + 1))
meanSdPlot(assay(rld[notAllZero,]))
meanSdPlot(assay(vsd[notAllZero,]))

## ----heatmap, dev="pdf", fig.width=5, fig.height=7-----------------------
library("pheatmap")
select <- order(rowMeans(counts(dds2,normalized=TRUE)),
                decreasing=TRUE)[1:2000]

nt <- normTransform(dds2) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds2)[,c("group_list")])
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE)

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE)

## ----Sample distances---------------------------------------------------------
sampleDists <- dist(t(assay(rld)))

## ----figHeatmapSamples, dev="pdf", fig.width=7, fig.height=7, fig.show="asis", fig.small=TRUE, fig.pos="tb", fig.cap="Sample-to-sample distances.  Heatmap showing the Euclidean distances between the samples as calculated from the regularized log transformation.\\label{figure/figHeatmapSamples-1}"----
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$group_list)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

## ----figPCA, dev="pdf", fig.width=5, fig.height=3------------------------
plotPCA(rld, intgroup=c("group_list"))

## ----figPCA2, dev="pdf", fig.width=5, fig.height=3-----------------------
data <- plotPCA(rld, intgroup=c("group_list"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=group_list)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
