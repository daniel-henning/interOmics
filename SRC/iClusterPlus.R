library(iClusterPlus)
library(GenomicRanges)
library(gplots)
library(lattice)
data(gbm)
###################################3
#The somatic mutation data should be stored in binary matrix (1: mutation, 0: no mutation)
#with the rows and columns corresponding to the samples and genes, respectively.
mut.rate=apply(gbm.mut,2,mean)
gbm.mut2 = gbm.mut[,which(mut.rate>0.02)]
gbm.mut2[1:10,1:8]

#For gene expression data, we recommend using the top variable genes for integrative clustering analysis, 
#which can be obtained by variance filtering. For example, we use the top 1740 genes for our iCluster analysis.
dim(gbm.exp)
# the rows and columns corresponding to the samples and genes respectively
gbm.exp[1:3,1:8]

#We reduce the GBM copy number regions to 5K by removing the redundant regions using function CNregions.
dim(gbm.seg)
gbm.seg[1:3,] #gbm.cn is the segmentation results produced by DNAcopy

data(variation.hg18.v10.nov.2010)
gbm.cn=CNregions(seg=gbm.seg,epsilon=0,adaptive=FALSE,rmCNV=TRUE,
                   cnv=variation.hg18.v10.nov.2010[,3:5],
                   frac.overlap=0.5, rmSmallseg=TRUE,nProbes=5)
dim(gbm.cn)
gbm.cn[1:3,1:5]
#Sort gbm.cn to make sure all the samples are in the same order.
gbm.cn=gbm.cn[order(rownames(gbm.cn)),]
# check if all the samples are in the same order for the three data sets
all(rownames(gbm.cn)==rownames(gbm.exp))
all(rownames(gbm.cn)==rownames(gbm.mut2))

#### #3 Integrative clustering analysis
fit.single=iClusterPlus(dt1=gbm.mut2,dt2=gbm.cn,dt3=gbm.exp,
                        type=c("binomial","gaussian","gaussian"),
                        lambda=c(0.04,0.61,0.90),K=2,maxiter=10)
#### #4 Model tuning using tune.iClusterPlus
set.seed(123)
date()

for(k in 1:5){
  cv.fit = tune.iClusterPlus(cpus=12, dt1=gbm.mut2, dt2=gbm.cn, dt3=gbm.exp, 
                             type=c("binomial","gaussian","gaussian"),K=k, 
                             n.lambda=185, scale.lambda=c(1,1,1), maxiter=20)
  save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
}
date()

#### #5 Model selection
output=alist()
files=grep("cv.fit",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)

#Now we get the ID for the lambda vector at which the BIC is minimum. Then we obtain
#the deviance ratio of the lambda vector at which the BIC is minimum.
minBICid = apply(BIC,2,which.min)
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}

#The optimal k (number of latent variables) is where the curve of %Explained variation levels off.
clusters=getClusters(output)
rownames(clusters)=rownames(gbm.exp)
colnames(clusters)=paste("K=",2:(length(output)+1),sep="")
#write.table(clusters, file="clusterMembership.txt",sep='\t',quote=F)
k=2
best.cluster=clusters[,k]
best.fit=output[[k]]$fit[[which.min(BIC[,k])]]

#### #6 Generate heatmap
plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
       ylab="%Explained Variation")

chr=unlist(strsplit(colnames(gbm.cn),"\\."))
chr=chr[seq(1,length(chr),by=2)]
chr=gsub("chr","",chr)
chr=as.numeric(chr)
#truncate the values for a better image plot
cn.image=gbm.cn
cn.image[cn.image>1.5]=1.5
cn.image[cn.image< -1.5]= -1.5
exp.image=gbm.exp
exp.image[exp.image>2.5]=2.5
exp.image[exp.image< -2.5]= -2.5

#### #7 Feature selection
#Select the top features based on lasso coecient estimates for the 3-cluster solution.
features = alist()
features[[1]] = colnames(gbm.mut2)
features[[2]] = colnames(gbm.cn)
features[[3]] = colnames(gbm.exp)
sigfeatures=alist()
for(i in 1:3){
  rowsum=apply(abs(best.fit$beta[[i]]),1, sum)
  upper=quantile(rowsum,prob=0.75)
  sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
}
names(sigfeatures)=c("mutation","copy number","expression")
#print a few examples of selected features
head(sigfeatures[[1]])
head(sigfeatures[[2]])
head(sigfeatures[[3]])

#We notice that no gene is selected from the mutation data, which indicates that the selected
#lambda value is too large and it is not in the same scale as those for the copy number and
#gene expression data. To solve this problem, we need to set the scale.lambda (an argument
#of tune.iClusterPlus) to a value between 0 and 1.

bw.col = colorpanel(2,low="white",high="black")
col.scheme = alist()
col.scheme[[1]] = bw.col
col.scheme[[2]] = bluered(256)
col.scheme[[3]] = bluered(256)
plotHeatmap(fit=best.fit,datasets=list(gbm.mut2,cn.image,exp.image),
            type=c("binomial","gaussian","gaussian"), col.scheme = col.scheme,
            row.order=c(F,F,T),chr=chr,plot.chr=c(F,T,F),sparse=c(T,F,T),cap=c(F,T,F))

#
set.seed(123)
date()
for(k in 1:5){
  cv2.fit = tune.iClusterPlus(cpus=12,dt1=gbm.mut2,dt2=gbm.cn,dt3=gbm.exp,
                              type=c("binomial","gaussian","gaussian"),K=k,n.lambda=185,
                              scale.lambda=c(0.05,1,1),maxiter=20)
  save(cv2.fit, file=paste("cv2.fit.k",k,".Rdata",sep=""))
}

date()

output2=alist()
files=grep("cv2.fit",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output2[[i]]=cv2.fit
}
nLambda = nrow(output2[[1]]$lambda)
nK = length(output2)
BIC = getBIC(output2)
devR = getDevR(output2)
minBICid = apply(BIC,2,which.min)
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}

#The optimal k (number of latent variables) is where the curve of %Explained variation levels off.
clusters=getClusters(output2)
rownames(clusters)=rownames(gbm.exp)
colnames(clusters)=paste("K=",2:(length(output2)+1),sep="")
#write.table(clusters, file="clusterMembership.txt",sep='\t',quote=F)
k=2
best.cluster=clusters[,k]
best.fit=output2[[k]]$fit[[which.min(BIC[,k])]]
#Select the top features based on lasso coecient estimates for the 3-cluster solution.
features = alist()
features[[1]] = colnames(gbm.mut2)
features[[2]] = colnames(gbm.cn)
features[[3]] = colnames(gbm.exp)
sigfeatures=alist()
for(i in 1:3){
  rowsum=apply(abs(best.fit$beta[[i]]),1, sum)
  upper=quantile(rowsum,prob=0.75)
  sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
}
names(sigfeatures)=c("mutation","copy number","expression")
#print a few examples of selected features
head(sigfeatures[[1]])
head(sigfeatures[[2]])
head(sigfeatures[[3]])
plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
     ylab="%Explained Variation")

plotHeatmap(fit=best.fit,datasets=list(gbm.mut2,cn.image,exp.image),
            type=c("binomial","gaussian","gaussian"), col.scheme = col.scheme,
            row.order=c(F,F,T),chr=chr,plot.chr=c(F,T,F),sparse=c(T,F,T),cap=c(F,T,F))