#Define a zero-center function to preprocess mRNA,miRNA,methylation data et al.
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

#Define a function to transfer fpkm to tpm
fpkm2tpm <- function(rawdat){
  if(dim(data_input)[1]<dim(data_input)[2]){
    dat <- apply(rawdat,1,function(x)x/sum(x))*1000000
  }else{
    dat <- apply(rawdat,2,function(x)x/sum(x))*100000
  }
  return(dat)
}

#Define a fuction to trimming expression data with number of fpkm > 1,
#the default number was set as n>6 with 470 samples.
trimming <- function(rawdat,n){
  for(i in 1:nrow(rawdat)){
    if(length(rawdat[i,][rawdat[i,]>1]) <= n){
      rawdat[i,] <- NA
    }
  }
  new.dat <- na.omit(rawdat)
  return(t(new.dat))
}

#<<===============================================================================
###################################clustering method##############################
#===============================================================================>>
########
#1st####
#期望最大化聚类(Expectation Maximization Algorithm)
library(mclust)
mc <- Mclust(input_data)

########
#2nd####
#k-means (k-均值聚类)
library(stats)
kmeans.result<-kmeans(input_data,4)
table(group_info,kmeans.result$cluster)

plot(iris2[c language="("][/c],col=kmeans.result$cluster)
# plot cluster centers
points(kmeans.result$centers[,c(“Sepal.Length”,”Sepal.Width“)],col=1:3,pch=8,cex=2)
########
#3rd####
#K-Medoids (k-中心点聚类)


########
#4th####
#DBSCAN (密度聚类)
library(cluster)#做聚类的包  
library(fpc)#有dbscan  
city <- read.csv("中国城市坐标.csv")  
x <- city[,c(3,2)]#行，列  
#ds <- dbscan(x, 2)#2是半径，最小点数默认为5  
ds <- dbscan(x, 2,6)  
#ds <- dbscan(x,1,3)  
ds#可以看border，seed数  
str(ds)#可以看列数  
par(bg="grey")  
plot(ds, x)  

########
#5th####
#h-clust (层次聚类)
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



#define a function to caculate coefficient of variation
##################################
c.v <- function(data){
  c.v  <- apply(data,2,function(x)sd(x)/mean(x))
  mean <- apply(data,2,function(x)mean(x))
  sd   <- apply(data,2,function(x)sd(x))
  return(list(c.v=c.v,sd=sd,mean=mean))
}



######################
setwd("C:\\Users\\Ning\\Desktop\\melanoma_data/")
rawdat <- read.table('mRNA.fpkm.sorted.txt',header = T,row.names = 1,check.names = F,sep = ',')

#
library(fpc)
pamk.result <- pamk(tpm)
pamk.result$nc




