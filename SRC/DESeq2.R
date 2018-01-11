## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
## ----loadDESeq2, echo=FALSE----------------------------------------------
library(DESeq2)
library(ggplot2)
library(dplyr)
#
library("BioParallel")
register(MulticoreParam(4)) ##申请
#rawdata<- read.csv("zebrafish_gene_count_matrix.csv",header=TRUE,row.names=1,check.names = FALSE)
rawdata<- read.csv("./data/TCGA_SKCM/SKCM_mRNA_htseq_count.txt",header=TRUE,sep = '\t',row.names = 1,check.names = FALSE)
rawdata <- trimming(rawdata)

clust <- read.table('./data/clust_diff/mrna_tsne_ahc_clust.txt',header = T,sep = '\t')
##
exprSet <- as.matrix(rawdata)
##
group_list <- factor(clust$sub)
#patient.id <- factor(c(clinical$colnames.mrna.f.))
colData <- data.frame(row.names = colnames(exprSet), group_list = group_list)
#
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)
dds2 <- DESeq(dds)  ##第二步，直接用DESeq函数即可

resultsNames(dds)
#
res <-  results(dds, contrast=c("groupList","1","2"))
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)

write.table(resOrdered,'mrna_tsne-aghc_cluster1-cluster2_DE_genes.txt',quote = F,sep = '\t')

#
res <-  results(dds, contrast=c("groupList","1","3"))
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)

write.table(resOrdered,'./cluster1-cluster3_DE_genes.txt',quote = F,sep = '\t')
#
res <-  results(dds, contrast=c("groupList","1","4"))
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)

write.table(resOrdered,'./cluster1-cluster4_DE_genes.txt',quote = F,sep = '\t')
#
res <-  results(dds, contrast=c("groupList","2","3"))
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)

write.table(resOrdered,'./cluster2-cluster3_DE_genes.txt',quote = F,sep = '\t')

#
res <-  results(dds, contrast=c("groupList","2","4"))
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)

write.table(resOrdered,'./cluster2-cluster4_DE_genes.txt',quote = F,sep = '\t')

#
res <-  results(dds, contrast=c("groupList","3","4"))
## 提取你想要的差异分析结果，我们这里是treated组对untreated组进行比较
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
write.table(resOrdered,'./cluster3-cluster4_DE_genes.txt',quote = F,sep = '\t')

##############################################################################
##############################################################################
DE1 <- read.table('cluster1-cluster2_DE_genes.txt',header = T,sep = '\t')
DE1<- DE1[which(DE1$log2FoldChange > 1 | DE1$log2FoldChange < -1 & DE1$padj < 0.005),]
genes1=DE1$GeneSymbol
eg1 = bitr(genes1, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
gene1 <- eg1[,2]
#DEplot(gene1,dir = 'cluster1-cluster2')

####################
DE2 <- read.table('cluster1-cluster3_DE_genes.txt',header = T,sep = '\t')
DE2<- DE2[which(DE2$log2FoldChange > 1 | DE2$log2FoldChange < -1 & DE2$padj < 0.005),]
genes2=DE2$GeneSymbol
eg2 = bitr(genes2, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
gene2 <- eg2[,2]
#DEplot(gene2,dir = 'cluster1-cluster3')

####################
DE3 <- read.table('cluster1-cluster4_DE_genes.txt',header = T,sep = '\t')
DE3<- DE3[which(DE3$log2FoldChange > 1 | DE3$log2FoldChange < -1 & DE3$padj < 0.005),]
genes3=DE3$GeneSymbol
eg3 = bitr(genes3, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
gene3 <- eg3[,2]
#DEplot(gene3,dir = 'cluster1-cluster4')

####################
DE4 <- read.table('cluster2-cluster3_DE_genes.txt',header = T,sep = '\t')
DE4<- DE4[which(DE4$log2FoldChange > 1 | DE4$log2FoldChange < -1 & DE4$padj < 0.005),]
genes4=DE$GeneSymbol
eg4 = bitr(genes4, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
gene4 <- eg4[,2]
#DEplot(gene4,dir = 'cluster2-cluster3')

################################
DE5 <- read.table('cluster2-cluster4_DE_genes.txt',header = T,sep = '\t')
DE5<- DE5[which(DE5$log2FoldChange > 2 | DE5$log2FoldChange < -2 & DE5$padj < 0.005),]
genes5=DE5$GeneSymbol
####convert gene symbol to entrezID####
eg5 = bitr(genes5, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
gene5 <- eg5[,2]
#DEplot(DE5,dir = 'cluster2-cluster4')

################################
DE6 <- read.table('cluster3-cluster4_DE_genes.txt',header = T,sep = '\t')
DE6<- DE6[which(DE6$log2FoldChange > 1 | DE6$log2FoldChange < -1 & DE6$padj < 0.005),]
genes6=DE6$GeneSymbol
####convert gene symbol to entrezID####
eg6 = bitr(genes6, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
gene6 <- eg6[,2]
#DEplot(DE6,dir = 'cluster3-cluster4')


gene <- gene5
dir = 'cluster2-cluster4'

geneList <- as.data.frame(t(subset(as.data.frame(t(genes)), select = as.factor(eg$SYMBOL))))
geneList <- as.data.frame(subset(geneList,select = c(GeneSymbol,log2FoldChange)))
#DEplot <- function(gene,dir=''){
  ##########################################enrichGO###############################################
  mfgo <- enrichGO(gene = gene, OrgDb="org.Hs.eg.db",ont = "MF", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  bpgo <- enrichGO(gene = gene, OrgDb="org.Hs.eg.db",ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  ccgo <- enrichGO(gene = gene, OrgDb="org.Hs.eg.db",ont = "CC", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  mf <- as.data.frame(mfgo)
  m <- c('MF')
  mf_t <- rep(m,times=length(mfgo$ID))
  mf$GO_term <- mf_t
  bp <- as.data.frame(bpgo)
  b <- c('BP')
  bp_t <- rep(b,times=length(bpgo$ID))
  bp$GO_term <- bp_t
  cc <- as.data.frame(ccgo)
  c <- c('CC')
  cc_t <- rep(c, times=length(ccgo$ID))
  cc$GO_term <- cc_t
  GO <- rbind(mf,bp,cc)
  #
  write.csv(GO,paste(dir,'mrna_tsne-aghc_GO_up.csv',sep = '/'),row.names=F)
############################
  df <- data.frame(x = GO$Description,y = GO$pvalue,z = GO$GO_term)
  df$y <- -log2(with(df, y))
  df <- df[with(df,order(y,decreasing = TRUE)),]
  df <- df[1:50,]
  png(paste(dir,'GO_plot.png',sep = '/'),height = 1000,width = 1000)
  ggplot(df,aes(x = interaction(x,z),y = y,fill = z))+ geom_bar(stat = "identity")+coord_flip()+xlab("GO Term Description (top100)")+ylab("P-value (-log2)")+labs(fill="GO_term")+coord_flip()+theme(axis.text.x=element_text(size=10,colour="black",face="bold"),axis.text.y=element_text(colour = "black",face="bold",vjust=0.8))
  dev.off()
##########################
  png(paste(dir,'MF_GO_plot.png',sep = '/'),height = 800,width = 600)
  dotplot(mfgo, showCategory=10,colorBy="pvalue")
  dev.off()
  png(paste(dir,'BP_GO_plot.png',sep = '/'),height = 800,width = 600)
  dotplot(bpgo, showCategory=10,colorBy="pvalue")
  dev.off()
  png(paste(dir,'CC_GO_plot.png',sep = '/'),height = 800,width = 600)
  dotplot(ccgo, showCategory=10,colorBy="pvalue")
  dev.off()
  
  ####label
  pdf(paste(dir,'mfGOenrichment.pdf',sep = '/'))
  plotGOgraph(mfgo)
  dev.off()
  pdf(paste(dir,'bpGOenrichment.pdf',sep = '/'))
  plotGOgraph(bpgo)
  dev.off()
  pdf(paste(dir,'ccGOenrichment.pdf',sep = '/'))
  plotGOgraph(ccgo)
  dev.off()
  #########################################KEGG analysis############################################
  ekk <- enrichKEGG(gene,organism="hsa",keyType='kegg', pvalueCutoff = 0.05,pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,qvalueCutoff = 0.05)
  write.csv(as.data.frame(ekk),paste(dir,"KEGG-enrich.csv",sep = '/'),row.names =F)
#}

  immuno.gene <- read.table('clipboard',header = F)
  immuno.gene.mrna <- as.data.frame(subset(mrna.t,select = as.array(immuno.gene[,1])))
  immuno.gene.mrna$label <- clust$ck3
  immuno.gene.mrna <- immuno.gene.mrna[order(immuno.gene.mrna$label),]
  immuno.gene.mrna <- subset(immuno.gene.mrna,select = -c(label))
  clust <- read.table('./mrna_tsne-aghc_clust.txt',header = T)
  sub <- sort(clust$ck3)
  ann_col = data.frame(sub)
  rownames(ann_col) = rownames(immuno.gene.mrna)
  ann_col$sub=paste('C',ann_col$sub,sep='')
  ann_col$sub=as.factor(ann_col$sub)
  #mrna.de <- subset(mrna.de, select = -c(sub))
  
  pheatmap(log10(t(immuno.gene.mrna)+0.1),col=colorRampPalette(c('blue','white','red'))(50),fontsize=8,cluster_col=F,cluster_row=T,annotation = ann_col)
  
  
