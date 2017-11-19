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


