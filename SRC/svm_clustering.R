#SVM example with Iris Data in R

#Use library e1071, you can install it using install.packages(??e1071??). Load library

library("e1071")
#Using Iris data
data(iris)
head(iris,5)
rawdat <- read.table('mRNA.fpkm.sorted.txt',header = T,row.names = 1,check.names = FALSE,sep = ',')
#Attach the Data
attach(rawdat)
#Divide Iris data to x (containt the all features) and y only the classes

x <- subset(iris, select=-Species)
y <- Species

#Create SVM Model and show summary
svm_model <- svm(Species ~ ., data=iris)
summary(svm_model)


#Or you can use command like this
#Create SVM Model and show summary
x <- rawdat
y <- read.csv('cluster.csv')
svm_model1 <- svm(x,y)
summary(svm_model1)


##  setosa versicolor virginica
#Run Prediction and you can measuring the execution time in R
pred <- predict(svm_model1,x)
system.time(pred <- predict(svm_model1,x))

#See the confusion matrix result of prediction, using command table to compare the result of SVM prediction and the class data in y variable.
table(pred,y)


svm_tune <- tune(svm, train.x=x, train.y=y, 
              kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))

print(svm_tune)

#After you find the best cost and gamma, you can create svm model again and try to run again
svm_model_after_tune <- svm(Species ~ ., data=iris, kernel="radial", cost=1, gamma=0.5)
summary(svm_model_after_tune)
## 


#Run Prediction again with new model
pred <- predict(svm_model_after_tune,x)
system.time(predict(svm_model_after_tune,x))


#See the confusion matrix result of prediction, using command table to compare the result of SVM prediction and the class data in y variable.
table(pred,y)
