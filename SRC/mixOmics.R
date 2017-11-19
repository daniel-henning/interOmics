################################################case study#######################################
#######PCA
library(mixOmics)
data(multidrug)
X <- multidrug$ABC.trans 
dim(X) # check dimension of data

trans.pca <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)
trans.pca

plot(trans.pca)

trans.pca2 <- pca(X, ncomp = 48, center = TRUE, scale = TRUE) 
# some warnings may appear as we are asking for many comp and the algo may not converge
plot(trans.pca2)

#samples plots
plotIndiv(trans.pca, comp = c(1, 2), ind.names = TRUE, 
          group = multidrug$cell.line$Class, 
          legend = TRUE, title = 'Multidrug transporter, PCA comp 1 - 2')

#variable plots
plotVar(trans.pca, comp = c(1, 2), var.names = TRUE, 
        title = 'Multidrug transporter, PCA comp 1 - 2')
#Biplots allow to both samples and variables to be graphically displayed simultaneously.
biplot(trans.pca, cex = 0.7,
       xlabs = paste(multidrug$cell.line$Class, 1:nrow(X)))



###############ipca
#Case study Independent PCA analysis (IPCA) using Liver Toxicity data set
data(liver.toxicity)
X <- liver.toxicity$gene
#Preliminary analysis with PCA
liver.pca<- pca(X, ncomp = 3, scale = FALSE)
plotIndiv(liver.pca, ind.names = liver.toxicity$treatment[, 3], group= liver.toxicity$treatment[,4], legend = TRUE, title = 'Liver, PCA')

#IPCA combines the advantages of both PCA and Independent Component Analysis (ICA) to reveal insightful patterns in the data
liver.ipca <- ipca(X, ncomp = 3, mode="deflation", scale = FALSE)
#sample plots
plotIndiv(liver.ipca, ind.names = liver.toxicity$treatment[, 3], group= liver.toxicity$treatment[,4], legend = TRUE, title = 'Liver, IPCA')

#Variable Plots
head(selectVar(liver.ipca, comp = 1)$value)

head(selectVar(liver.pca, comp = 1)$value)

plotVar(liver.ipca, pch = 20)


#The kurtosis measure is used to order the loading vectors to order the Independent Principal Components
liver.ipca$kurtosis

#Sparse Independent PCA analysis (sIPCA)
liver.sipca <- sipca(X, ncomp = 3, mode = "deflation",
                     scale = FALSE, keepX = c(50,50,50))
#sample plot
plotIndiv(liver.sipca, ind.names = liver.toxicity$treatment[, 3], group= liver.toxicity$treatment[,4], legend = TRUE, title = 'Liver, IPCA')

#Variable Plots

head(selectVar(liver.sipca, comp = 1)$value)


plotVar(liver.sipca, pch = 20)


#Case study with PCA and sPLS on Liver Toxicity data set




