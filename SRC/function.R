library(ggfortify)

###################################################
### Define a zero-center function to preprocess mRNA,miRNA,methylation data et al.
###################################################
zero.cent <- function(data_input){
  if(dim(data_input)[1]>dim(data_input)[2]){
    #Because mRNA and methylation data have a p>n feature,
    # so this step is to check if the data_input is defined as nxp or pxn.
    data.new <- apply(data_input,2,function(x)x-mean(x))
  }else{
    data.new <- apply(data_input,1,function(x)x-mean(x))
  }
  return(data.new)
}


###################################################
###
###################################################
normalize <- function(x) {return ((x - min(x)) / (max(x) - min(x))) }

###################################################
### Define a function to transfer fpkm to tpm
###################################################

fpkm2tpm <- function(rawdat){
  if(dim(rawdat)[1]<dim(rawdat)[2]){
    dat <- apply(rawdat,1,function(x)x/sum(x))*1000000
  }else{
    dat <- apply(rawdat,2,function(x)x/sum(x))*1000000
  }
  return(t(dat))
}

###################################################
### Define a fuction to trimming expression data with number of fpkm > 1().
###################################################
trimming <- function(dat,cutoff=1,sub_group=4){
  if(nrow(dat)>ncol(dat)){
    for(i in 1:nrow(dat)){
      if(length(dat[i,][which(as.numeric(dat[i,])>cutoff)]) < sub_group){
        dat[i,1] <- NA
      }
    }
  }else{
    dat <- t(dat)
    for(i in 1:nrow(dat)){
      if(length(dat[i,][which(as.numeric(dat[i,])>cutoff)]) < sub_group){
        dat[i,1] <- NA
      }
    }
  }
  dat <- na.omit(dat)
  return(dat)
}

###################################################
### define a function to caculate coefficient of variation
###################################################
c.v <- function(data){
  if(dim(data)[1]<dim(data)[2]){
    c.v  <- apply(data,2,function(x)sd(x)/mean(x))
    mean <- apply(data,2,function(x)mean(x))
    sd   <- apply(data,2,function(x)sd(x))
  }else{
    c.v  <- apply(data,1,function(x)sd(x)/mean(x))
    mean <- apply(data,1,function(x)mean(x))
    sd   <- apply(data,1,function(x)sd(x))
  }
  return(list(c.v=c.v,sd=sd,mean=mean))
}
################################################################################
#####loading data matrix
tcga.mrna <- read.table('./data/TCGA_SKCM/TCGA_SKCM_mrna_fpkm.txt',header = T,row.names = 1,check.names = F)
#mrna.t <- as.data.frame(t(mrna))
#mrna.t.sort <- mrna.t[order(mrna.t$V1),]
#write.table(t(mrna.t.sort),'TCGA_SKCM_mrna_fpkm.txt',row.names = F,quote = F,sep = '\t')
tcga.mrna.nromal <- subset(tcga.mrna,select = c(TCGA_GN_A4U8_11))
tcga.mrna.tumor <- subset(tcga.mrna,select = -c(TCGA_GN_A4U8_11))
tcga.mrna.tumor.t <- as.data.frame(t(trimming(tcga.mrna.tumor,cutoff = 1,sub_group = 20)))
tcga.mrna.tumor.t.scaled <- scale(tcga.mrna.tumor.t)


##############Plotting PCA (Principal Component Analysis)
rawdat.log <- t(log2(mRNA.trimmed+0.01))
pc <- pca(rawdat.log)
autoplot(pc,data = rawdat.log, label = TRUE, label.size = 3)


#get top 20% most-variable genes.
fpkm.top <- sort(mRNA.cv$c.v,decreasing = T)[1:floor(length(mRNA.cv$c.v)*0.2)]
mRNA.top <- mRNA.trimmed[names(fpkm.top),]


#gene_filter
gene.ft <- function(data){
  f1 <- kOverA(5, 1)
  ffun <- filterfun(f1)
  wh1 <- genefilter(exprs(sample.ExpressionSet), ffun)
  sum(wh1)
}

#<<===============================================================================
###################################clustering method##############################
#===============================================================================>>

###################################################
### 期望最大化聚类(Expectation Maximization Algorithm)
###################################################
library(mclust)
mc <- Mclust(as.matrix(tpm.trimmed), G=1:20)
m.best <- dim(mc$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")
plot(d_clust)


###################################################
### k-means (k-均值聚类)
###################################################
library(stats)
kmeans.result<-kmeans(tpm.trimmed,4)
table(group_info,kmeans.result$cluster)

plot(iris[,-5],col=kmeans.result$cluster)
# plot cluster centers
points(kmeans.result$centers[,c("Sepal.Length","Sepal.Width")],col=1:3,pch=8,cex=2)
#define a function to see the results with different k-means
wss <- (nrow(tpm.trimmed)-1)*sum(apply(tpm.trimmed,2,var))
for(i in 2:30){wss[i] <- sum(kmeans(tpm.trimmed,centers=i)$withinss)}
###这里的wss(within-cluster sum of squares)是组内平方和
plot(1:20, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")


###################################################
### K-Medoids (k-中心点聚类)
###################################################
library(fpc)
pamk.result <- pamk(tpm.trimmed)
pamk.result$nc
plot(pamk.result$pamobject)
#======#
library(cluster)
pam.result <- pam(tpm.trimmed,3)
table(pam.result$clustering)
plot(pam.result)



###################################################
### DBSCAN (密度聚类)
###################################################
library(cluster)#做聚类的包  
library(fpc)#有dbscan  
library("dbscan")
data("iris")
x <- as.matrix(iris[, 1:4])
#Run DBSCAN
db <- dbscan(x, eps = .4, minPts = 10)
db
#Visualize results (noise is shown in black)
pairs(x, col = db$cluster + 1L)
#Calculate LOF (local outlier factor) and visualize (larger bubbles in the visualization have a larger LOF)
lof <- lof(x, k = 4)
pairs(x, cex = lof)
#Run OPTICS
opt <- optics(x, eps = 1, minPts = 10)
opt
#Extract DBSCAN-like clustering from OPTICS and create a reachability plot (extracted DBSCAN clusters at eps_cl=.4 are colored)
opt <- extractDBSCAN(opt, eps_cl = .4)
plot(opt)
#Extract a hierarchical clustering using the Xi method (captures clusters of varying density)
opt <- extractXi(opt, xi = .05)
opt
plot(opt)
#Run HDBSCAN (captures stable clusters)
hdb <- hdbscan(x, minPts = 4)
hdb
#Visualize the results as a simplified tree
plot(hdb, show_flat = T)
#See how well each point corresponds to the clusters found by the model used
colors <- mapply(function(col, i) adjustcolor(col, alpha.f = hdb$membership_prob[i]), 
                 palette()[hdb$cluster+1], seq_along(hdb$cluster))
plot(x, col=colors, pch=20)


###################################################
### h-clust (层次聚类)
###################################################
#d为待处理数据集样本间的距离矩阵,可用dist()函数计算得到;
#method参数用于选择聚类的具体算法,可供选择的有ward,single及complete等7种,默认选择complete方法;
#参数members用于指出每个待聚类样本点/簇是由几个单样本构成,
#该参数默认值为NULL,表示每个样本点本身即为单样本.
#cutree()函数则可以对hclust()函数的聚类结果进行剪枝
library(stats)
hh<-hclust(dist(scale(input_data)),"complete")
plot(hh,labels=w[,1],cex=0.6)
#自动分成5类
rect.hclust(hh,k=5)


############################################
###MCLUST
###########################################
#install.packages("mclust")
library(mclust)
#You can then perform model-based clustering on the iris dataset using Mclust:
mb = Mclust(iris[,-5])
#or specify number of clusters
mb3 = Mclust(iris[,-5], 3)
# optimal selected model
mb$modelName
# optimal number of cluster
mb$G
# probality for an observation to be in a given cluster
head(mb$z)
# get probabilities, means, variances
summary(mb, parameters = TRUE)
#Compare amount of the data within each cluster
table(iris$Species, mb$classification)
# vs
table(iris$Species, mb3$classification)
#After the data is fit into the model, we plot the model based on clustering results.
plot(mb, what=c("classification"))






#<<===============================================================================
##############################Dimension Reduction method##########################
#===============================================================================>>


###################################################
### 
###################################################





###################################################
### 2nd: t-SNE
###################################################
mrna.t <- as.data.frame(t(mrna.f))

##########################################
# load the Rtsne package
library(Rtsne)
# run Rtsne with default parameters
rtsne_out <- Rtsne(as.data.frame(t(mrna.f)), dims = 3, perplexity=30)

# plot the output of Rtsne into d:\\barneshutplot.jpg file of 2400x1800 dimension
jpeg("barneshutplot.jpg", width=1000, height=1000)
plot(rtsne_out$Y, t='n', main="BarnesHutSNE")
text(rtsne_out$Y, labels=rownames(t(mirna.f)))
dev.off()
#############################################

#############################################
wss <- (nrow(rtsne_out$Y)-1)*sum(apply(rtsne_out$Y,2,var))
for (i in 2:15)
  wss[i] <- sum(kmeans(rtsne_out$Y,centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

#############################################
# mydata is the matrix loaded into R previously. Run the following command if it isn't yet loaded into R. 
# Cluster Plot against 1st 2 principal components
library(cluster)
# K-Means Clustering with 8 clusters
fit <- kmeans(rtsne_out$Y, 2)
#
clusplot(rtsne_out$Y, fit$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)
## plot the output of Rtsne
colors = rainbow(length(unique(fit$cluster)))
names(colors) = unique(fit$cluster)
plot(rtsne_out$Y, t='n', main="tSNE-kmeans")
text(rtsne_out$Y, labels=fit$cluster,col = colors[fit$cluster])

############################################
group_m <- as.data.frame(colnames(mrna.f))
group_m$mrna_clust <- fit$cluster
#
group_mi <- as.data.frame(colnames(mirna.f))
group_mi$mirna_clust <- fit$cluster


#

library(pcaMethods)
library(pcaReduce)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
set.seed(1234567)

iris <- plotTSNE(iris[,-4], rand_seed = 1, return_SCE = TRUE)



wss <- (nrow(rtsne_out$Y)-1)*sum(apply(rtsne_out$Y,2,var))
for(i in 2:15){wss[i] <- sum(kmeans(rtsne_out$Y,centers=i)$withinss)}
###这里的wss(within-cluster sum of squares)是组内平方和
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

mirna.clust <-  kmeans(rtsne_out$Y, centers = 2)$clust

plot(rtsne_out$Y, t='n', main="tsne")
text(rtsne_out$Y, mirna.clust, col=colors[mirna.clust])


library(fpc)

min.max.norm <- function(x){
  (x-min(x))/(max(x)-min(x))
}
raw.data <- iris[,1:4]
norm.data <- data.frame(sl = min.max.norm(raw.data[,1]),
                        sw = min.max.norm(raw.data[,2]),
                        pl = min.max.norm(raw.data[,3]),
                        pw = min.max.norm(raw.data[,4]))


K <- 2:10
round <- 30 # 每次迭代30次，避免局部最优
rst <- sapply(K, function(i){
  print(paste("K=",i))
  mean(sapply(1:round,function(r){
    print(paste("Round",r))
    result <- kmeans(rtsne_out$Y, i)
    stats <- cluster.stats(dist(rtsne_out$Y), result$cluster)
    stats$avg.silwidth
  }))
})

plot(K,rst,type='l',main='轮廓系数与K的关系', ylab='轮廓系数')


old.par <- par(mfrow = c(1,2))
k = 2 # 根据上面的评估 k=2最优
clu <- kmeans(mirna.f,k)
mds = cmdscale(dist(mirna.f,method="euclidean"))
plot(mds, col=clu$cluster, main='kmeans聚类 k=2', pch = 19)
plot(mds, col=iris$Species, main='原始聚类', pch = 19)
par(old.par)


