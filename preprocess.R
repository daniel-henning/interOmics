#define a function to caculate coefficient of variation
##################################
c.v <- function(data){
  c.v  <- apply(data,2,function(x)sd(x)/mean(x))
  mean <- apply(data,2,function(x)mean(x))
  sd   <- apply(data,2,function(x)sd(x))
  return(list(c.v=c.v,sd=sd,mean=mean))
}


#define a fuction to trimming expression data with number of fpkm > 1,
#the default number was set as n>6 with 470 samples.
trimming <- function(dat,n){
  for(i in 1:nrow(dat)){
    if(length(dat[i,][dat[i,]>1]) <= n){
      dat[i,] <- NA
    }
  }
  new.dat <- na.omit(dat)
  return(t(new.dat))
}

#define a function to transfer fpkm to tpm
fpkm2tpm <- function(rawdat){
  dat <- apply(rawdat,1,function(x)x/sum(x))*1000000
  return(dat)
}


######################
setwd("C:\\Users\\Ning\\Desktop\\melanoma_data/")
rawdat <- read.table('mRNA.fpkm.sorted.txt',header = T,row.names = 1,check.names = F,sep = ',')

#
library(fpc)
pamk.result <- pamk(tpm)
pamk.result$nc


hc <- hclust()
