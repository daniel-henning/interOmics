#################################################################################
library(ggfortify)
setwd('C://Users/Ning/Desktop/melanoma_data/')
##############Plotting PCA (Principal Component Analysis)
data <- read.delim('mRNA_htseq_FPKM.txt',check.names = FALSE,header = T,row.names = 1)
data <- read.csv('clipboard',check.names = FALSE,header = T,row.names = 1,sep = '\t')

df <- t(data)
df <- na.omit(data)
autoplot(prcomp(df[,-25269]))
autoplot(prcomp(df),data = t(Group), colour = 'group')
autoplot(prcomp(df), data = t(Group), colour = 'group', label = TRUE, label.size = 3)
autoplot(prcomp(df), data = t(Group), colour = 'group', shape = FALSE, label.size = 3)
autoplot(prcomp(df), data = t(Group), colour = 'group', loadings = TRUE)
autoplot(prcomp(df), data = t(Group), colour = 'group',label = TRUE,loadings = TRUE, 
         loadings.colour = 'blue',loadings.label = TRUE, loadings.label.size = 2)
#
df.log <- log(df)

mu = colMeans(df.log)

X = df.log
Xpca = prcomp(X)
nComp = 5
Xhat = Xpca$x[,1:nComp] %*% t(Xpca$rotation[,1:nComp])
Xhat = scale(Xhat, center = mu, scale = FALSE)

gene <- sort(Xhat[1,],decreasing = TRUE)

par(mar = c(4,4,4,8))
plot(Xpca$x[,c(1,2)],pch=16,cex=1,main = "PCA map")
text(Xpca$x[,c(1,2)],row.names(Xpca$x),col="black",pos=3,cex=0.8)
legend("right",legend=Group,ncol = 1,xpd=T,inset = -0.15,
       pch=16,cex=1,col=rainbow(length(Group)),bty="n")

dev.off()
#<<==================================================================
#Plotting Factor Analysis
#==================================================================>>
#{ggfortify} supports stats::factanal object as the same manner as PCAs. Available opitons are the same as PCAs.
#Important You must specify scores option when calling factanal to calcurate sores (default scores = NULL). 
#Otherwise, plotting will fail.

d.factanal <- factanal(state.x77, factors = 3, scores = 'regression')
autoplot(d.factanal, data = state.x77, colour = 'Income')

autoplot(d.factanal, label = TRUE, label.size = 3,
         loadings = TRUE, loadings.label = TRUE, loadings.label.size  = 3)
#<<==================================================================
#Plotting K-means
#==================================================================>>
#{ggfortify} supports stats::kmeans class. 
#You must explicitly pass original data to autoplot function via data keyword !!!
set.seed(1)
autoplot(kmeans(USArrests, 3), data = USArrests)

autoplot(kmeans(USArrests, 3), data = USArrests, label = TRUE, label.size = 3)
#<<==================================================================
#Plotting cluster package
#==================================================================>>
library(cluster)
autoplot(clara(iris[-5], 3))
#Specifying frame = TRUE in autoplot for stats::kmeans and cluster::* draws convex for each cluster.
autoplot(fanny(iris[-5], 3), frame = TRUE)
#If you want probability ellipse, {ggplot2} 1.0.0 or later is required. 
#Specify whatever supported in ggplot2::stat_ellipse's type keyword via frame.type option.
autoplot(pam(iris[-5], 3), frame = TRUE, frame.type = 'norm')
#<<==================================================================
#Plotting Local Fisher Discriminant Analysis with {lfda} package
#==================================================================>>
#{lfda} package supports a set of Local Fisher Discriminant Analysis methods. 
#You can use autoplot to plot the analysis result as the same manner as PCA.
library(lfda)
# Local Fisher Discriminant Analysis (LFDA)
model <- lfda(iris[-5], iris[, 5], r = 3, metric="plain")
autoplot(model, data = iris, frame = TRUE, frame.colour = 'Species')
# Semi-supervised Local Fisher Discriminant Analysis (SELF)
model <- self(iris[-5], iris[, 5], beta = 0.1, r = 3, metric="plain")
autoplot(model, data = iris, frame = TRUE, frame.colour = 'Species')


#<<==================================================================
#Plotting Multidimensional Scaling
#==================================================================>>
#Before Plotting
#Even though MDS functions returns matrix or list (not specific class), 
#{ggfortify} can infer background class from list attribute and perform autoplot.
# NOTE: Inference from matrix is not supported.
# NOTE: {ggfortify} can plot stats::dist instance as heatmap.

autoplot(eurodist)
#<<==================================================================
#Plotting Classical (Metric) Multidimensional Scaling
#stats::cmdscale performs Classical MDS and returns point coodinates as matrix, 
#thus you can not use autoplot in this case. However, either eig = TRUE, 
#add = True or x.ret = True is specified, stats::cmdscale return list instead of matrix. 
#In these cases, {ggfortify} can infer how to plot it via autoplot. 
#Refer to help(cmdscale) to check what these options are.
autoplot(cmdscale(eurodist, eig = TRUE))
#Specify label = TRUE to plot labels.

autoplot(cmdscale(eurodist, eig = TRUE), label = TRUE, label.size = 3)
#====================================================================
#Plotting Non-metric Multidimensional Scaling
#MASS::isoMDS and MASS::sammon perform Non-metric MDS and return list which contains point coordinates.
#Thus, autoplot can be used.
# NOTE: On background, autoplot.matrix is called to plot MDS. See help(autoplot.matrix) to check available options.
library(MASS)
autoplot(isoMDS(eurodist), colour = 'orange', size = 4, shape = 3)
#Passing shape = FALSE makes plot without points. In this case, label is turned on unless otherwise specified.
autoplot(sammon(eurodist), shape = FALSE, label.colour = 'blue', label.size = 3)
#==================================================================>>


#**************************************************************************************
#--------------------------------------------------------------------------------------
#**************************************************************************************


#<<==================================================================
#Plotting with survival package
#==================================================================>>

#{ggfortify} let {ggplot2} know how to draw survival curves. After loading {ggfortify}, 
#you can use ggplot2::autoplot function for survfit objects.
library(ggfortify)
library(survival)
fit <- survfit(Surv(time, status) ~ sex, data = lung)
autoplot(fit)
#There are some options to change survival curve output. 
#Use help(autoplot.survfit) (or help(autoplot.*) for any other objects) to check available options.

autoplot(fit, surv.linetype = 'dashed', conf.int = FALSE,
         censor.shape = '*', censor.size = 5, facets = TRUE, ncol = 2)

autoplot(survfit(Surv(time, status) ~ 1, data = lung), surv.colour = 'orange', censor.colour = 'red')

autoplot(survfit(Surv(time, status) ~ sex, data = lung), fun = 'event')

d.coxph <- survfit(coxph(Surv(time, status) ~ sex, data = lung))
autoplot(d.coxph, surv.linetype = 'dashed', surv.colour = 'blue',
         conf.int.fill = 'dodgerblue3', conf.int.alpha = 0.5, censor = FALSE)

#Also, you can use autoplot for survival::aareg instance.
autoplot(aareg(Surv(time, status) ~ age + sex + ph.ecog, data = lung))


#<<==================================================================
#Plotting {glmnet}
#==================================================================>>

#{ggfortify} supports {glmnet} package which supports Regularized Generalized Linear Models (Ridge, Lasso and Elastic-net).

library(glmnet)
data(QuickStartExample)
fit <- glmnet::glmnet(x, y)
autoplot(fit)

fit <- glmnet::cv.glmnet(x, y)
autoplot(fit, colour = 'blue')


#Plotting Diagnostics for Linear Models

#{ggfortify} let {ggplot2} know how to interpret lm objects. 
#After loading {ggfortify}, you can use ggplot2::autoplot function for lm objects.

library(ggfortify)
autoplot(lm(Petal.Width ~ Petal.Length, data = iris), label.size = 3)

#You can select desired plot by which option as the same manner as standard plot. 
#Also, ncol and nrow allows you to specify the number of subplot columns and rows.

par(mfrow = c(1, 2))
m <- lm(Petal.Width ~ Petal.Length, data = iris)

autoplot(m, which = 1:6, ncol = 3, label.size = 3)

#Plotting Diagnostics for Generalized Linear Models
#t also suppotgs glm instance.

m <- glm(Murder ~ Assault + UrbanPop + Rape,
         family = gaussian, data = USArrests)

autoplot(m, which = 1:6, label.size = 3)

#Decorating Plots

#Because {ggplot2} itself cannot handle different kinds of plots in a single instance, 
#{ggfortify} handle them using its original class named ggmultiplot. You can use + operator to decorate ggmultiplot.

class(autoplot(m))
## [1] "ggmultiplot"
## attr(,"package")
## [1] "ggfortify"
autoplot(m, label.size = 3) + theme_bw()
#Specifing Plotting Options
#Some properties can be changed by passing corresponding keywords. 
#For example, colour keyword is for data points, smooth.colour is for smoothing lines and 
#ad.colour is for additional auxiliary lies. Also, ncol and nrow control facet layout. 
#Use help(autoplot.lm) (or help(autoplot.*) for any other objects) to check available options.

autoplot(m, which = 1:6, colour = 'dodgerblue3',
         smooth.colour = 'black', smooth.linetype = 'dashed',
         ad.colour = 'blue',
         label.size = 3, label.n = 5, label.colour = 'blue',
         ncol = 3)

#Also, you can use column names for these properties. 
#Note that lm and glm instances doesn't retain original data, 
#you should pass original data via data keyword to use column names not included in the model.

autoplot(lm(Petal.Width ~ Petal.Length, data = iris), data = iris,
colour = 'Species', label.size = 3)
