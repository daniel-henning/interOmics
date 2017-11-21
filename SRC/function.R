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
trimming <- function(dat,sub_group=20){
  if(nrow(dat)>ncol(dat)){
    for(i in 1:nrow(dat)){
      if(length(dat[i,][which(dat[i,]>1)]) < sub_group){
        dat[i,1] <- NA
      }
    }
  }else{
    dat <- t(dat)
    for(i in 1:nrow(dat)){
      if(length(dat[i,][which(dat[i,]>1)]) < sub_group){
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

######################loading data
rawdat <- read.table('data//mRNA_htseq_FPKM.txt',header = T,row.names = 1,check.names = F,sep = '\t')

#mRNA.tpm <- fpkm2tpm(rawdat)
#mRNA.tpm.trimmed <- trimming(mRNA.tpm,20,1)
mRNA.trimmed <- trimming(rawdat)
mRNA.cv <- c.v(mRNA.trimmed)
write.csv(mRNA.cv,'c.v.csv')


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
#city <- read.csv("中国城市坐标.csv")  
x <- city[,c(3,2)]#行，列  
#ds <- dbscan(x, 2)#2是半径，最小点数默认为5  
ds <- dbscan(tpm.trimmed, 2,6)
#ds <- dbscan(x,1,3)  
ds#可以看border，seed数  
str(ds)#可以看列数  
par(bg="grey")  
plot(ds, tpm.trimmed)


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




#<<===============================================================================
##############################Dimension Reduction method##########################
#===============================================================================>>


###################################################
### 
###################################################





###################################################
### 2nd: t-SNE
###################################################
# load the tsne package
library(tsne)

# initialize counter to 0
x <- 0
epc <- function(x) {
  x <<- x + 1
  filename <- paste(".//plot//", x, "jpg", sep=".")
  cat("> Plotting TSNE to ", filename, " ")
  
  # plot to d:\\plot.x.jpg file of 2400x1800 dimension
  jpeg(filename, width=2400, height=1800)
  
  plot(x, t='n', main="T-SNE")
  text(x, labels=rownames(tpm.trimmed))
  dev.off()
}
# run tsne (maximum iterations:500, callback every 100 epochs, target dimension k=5)
tsne_data <- tsne(tpm.trimmed, k=4, epoch_callback=epc, max_iter=500, epoch=100)
##################################
# load the Rtsne package
library(Rtsne)
# run Rtsne with default parameters
rtsne_out <- Rtsne(as.matrix(t(mRNA.trimmed)))
# plot the output of Rtsne into d:\\barneshutplot.jpg file of 2400x1800 dimension
jpeg("barneshutplot.jpg", width=2400, height=1800)
plot(rtsne_out$Y, t='n', main="BarnesHutSNE")
text(rtsne_out$Y, labels=rownames(t(mRNA.trimmed)))
##################################
# mydata is the matrix loaded into R previously. Run the following command if it isn't yet loaded into R. 
mydata <- read.table("d:\\samplewordembedding.csv", header=TRUE, sep=",")
# K-Means Clustering with 20 clusters
fit <- kmeans(mydata, 20)
# Cluster Plot against 1st 2 principal components
library(cluster)
clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)




#############################
dat <- loadDataSet("3D S Curve")

## use the S4 Class directly:
fastica <- FastICA()
emb <- fastica@fun(tpm.trimmed, pars = list(ndim = 2))

## simpler, use embed():
emb2 <- embed(dat, "FastICA", ndim = 2)


plot(emb@data@data)

###############################
# dat <- loadDataSet("3D S Curve")
# 
# ## Use the S4 Class directly:
mds <- MDS()
emb <- mds@fun(tpm.trimmed, mds@stdpars)
# 
# ## use embed():
emb2 <- embed(scale(tpm.trimmed), "MDS", d = function(x) exp(stats::dist(x)))
# 
# 
plot(emb, type = "2vars")
plot(emb2, type = "2vars")
