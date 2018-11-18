library(tidyverse)
library(reshape2)
library(caret)
library(ROCR)
library(pROC)
library(FFTrees)
library(boot)

data <- read.csv("data/final_data.csv")

set.seed(40) # Reproducibility

data$Group <- as.factor(data$Group)

# Cross-validation https://www.youtube.com/watch?v=p5rDg4OVBuA. Need to create validation set 
indices <- createDataPartition(data$Group, p = 0.7, list = FALSE) #70% of data will be for training
train <- data[indices,]
test <- data[-indices,]

# Shuffling data:

train <- train[sample(nrow(train)),]
# After cross validation, testing on training data. Then get best model --> predict for test data

controlParameters <- trainControl(
  method = "LOOCV", #Leave one out cross validation
  savePrediction = TRUE,
  classProbs = TRUE
)

parameterGrid <- data.frame(mtry = c(2,3,4)) # Parameter for optimization: https://stats.stackexchange.com/questions/102867/random-forest-mtry-question

accuracyMetric <- function (conMatrix) {
  tp <- conMatrix[1,1]
  fp <- conMatrix[1,2]
  fn <- conMatrix[2,1]
  tn <- conMatrix[2,2]
  
  metric <- (tp + tn)/(tp + fp + fn + tn)
  return(metric)
}

mccMetric <- function (conMatrix) {
  tp <- conMatrix[1,1]
  fp <- conMatrix[1,2]
  fn <- conMatrix[2,1]
  tn <- conMatrix[2,2]
  
  metric <- ((tp * tn) - (fp * fn))/(sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
  return(metric)
}

# Training random forest model
rfmodel <- train(Group ~.,
  data = train,
  method = "rf",
  trControl = controlParameters,
  tuneGrid = parameterGrid
)


plot(varImp(rfmodel), main = "Variable Importance - Random Forest")


rfPredictions <- predict(rfmodel, test, type = "prob")
confMatrixrf <- table(predictions = rfPredictions, actual = test$Group)
rfAccuracy <- accuracyMetric(confMatrixrf)
rfMCC <- mccMetric(confMatrixrf)

#Model: loadRDS("random_forest.rds") OR loadRDS("loocv_rf.rds")

# TUNING PARAMETERS: https://www.analyticsvidhya.com/blog/2016/08/practicing-machine-learning-techniques-in-r-with-mlr-package/

# Training logistic regression model. No tuning parameters?

logmodel <- train(Group ~ .,
  data = train,
  method = "glm",
  family = "binomial",
  trControl = controlParameters)

logmodel <- readRDS("log_model.rds")

plot(varImp(logmodel), main = "Variable Importance - Logistic Regression")

logPredictions <- predict(logmodel, test)
confMatrixlog <- table(predictions = logPredictions, actual = test$Group)
logAccuracy <- accuracyMetric(confMatrixlog)
logMCC <- mccMetric(confMatrixlog)
# MAKE ROC curves

# Training CART model
CARTmodel <- train(Group ~.,
  data = train,
  method = "rpart",
  trControl = controlParameters
)

CARTmodel <- readRDS("CARTmodel.rds")

rpart.plot::rpart.plot(CARTmodel$finalModel)

cartPredictions <- predict(CARTmodel, test, type = "prob")
confMatrixCART <- table(predictions = cartPredictions, actual = test$Group)
cartAccuracy <- accuracyMetric(confMatrixCART)
cartMCC <- mccMetric(confMatrixCART)

# Training FFT
fft.train <- train
fft.test <- test
fft.train$target <- NULL
fft.test$target <- NULL

fft.train$target[fft.train$Group == "AD"] <- 1
fft.test$target[fft.test$Group == "AD"] <- 1
fft.train$target[fft.train$Group == "CONTROL"] <- 0
fft.test$target[fft.test$Group == "CONTROL"] <- 0

fft.train <- fft.train[,2:6]
fft.test <- fft.test[,2:6]

fftModel <- FFTrees(formula = target ~.,
                    data = fft.train,
                    data.test = fft.test)

fftModel <- readRDS("fft_model.rds")

plot(fftModel, main = "Fast and Frugal Tree", decision.labels = c("AD", "Control"))

# Training KNN model
parameterGrid <- expand.grid(k = c(2,3,4,5,6,7,8,9,10))

knnModel <- train(Group ~ .,
  data = train,
  method = "knn",
  trControl = controlParameters,
  tuneGrid = parameterGrid)
knnModel <- readRDS("knn_model.RDS")

knnPredictions <- predict(knnModel, test, type = "prob")
confMatrixKNN <- table(predictions = knnPredictions, actual = test$Group)
knnAccuracy <- accuracyMetric(confMatrixKNN)
knnMCC <- mccMetric(confMatrixKNN)

# Training SVM model (radial kernel??)
parameterGrid <- expand.grid(C = c(0.1,0.2,0.3,0.4,0.5,0.6,0.6,0.7,0.8,0.9,1), sigma = 0.5)

svmModel <- train (Group ~ .,
  data = train,
  method = "svmRadial",
  trControl = controlParameters,
  tuneGrid = parameterGrid)
svmModel <- readRDS("svm_model.RDS")

svmPredictions <- predict(svmModel, test, type = "prob")
confMatrixSVM <- table(predictions = svmPredictions, actual = test$Group)
svmAccuracy <- accuracyMetric(confMatrixSVM)
svmMCC <- mccMetric(confMatrixSVM)

#Stacking models: https://www.youtube.com/watch?v=k7sTiTWWCXM&t=228s

## Correlation between models; fix this
#model_corr <- list(mod1 = CARTmodel, mod2 = fftModel, mod3 = knnModel, mod4 = logmodel, mod5 = rfmodel, mod6 = svmModel)
#modelCor(model_corr)

cartPredictions <- predict(CARTmodel, train)
knnPredictions <- predict(knnModel, train)
logPredictions <- predict(logmodel, train)
rfPredictions <- predict(rfmodel, train)
svmPredictions <- predict(svmModel, train)

test$cartPredictions <- predict(CARTmodel, test)
test$knnPredictions <- predict(knnModel, test)
test$logPredictions <- predict(logmodel, test)
test$rfPredictions <- predict(rfmodel, test)
test$svmPredictions <- predict(svmModel, test)

predDF <- data.frame(cartPredictions, knnPredictions, logPredictions, rfPredictions, svmPredictions, class = train$Group)

combineModel <- train(as.factor(class) ~., method = "rf", data = predDF, trControl = controlParameters)

ensemblePredct <- predict(combineModel, test, type = "prob")
confMatrixEnsemble <- table(predictions = ensemblePredct, actual = test$Group)
ensembleAccuracy <- accuracyMetric(confMatrixEnsemble)
ensembleMCC <- mccMetric(confMatrixEnsemble)

# Confidence intervals for random forest classifier
controlParameters <- trainControl(
  method = "boot",
  number = 1000
)

parameterGrid <- expand.grid(k = c(2,3,4,5,6,7,8,9,10))

knnModel <- train(Group ~.,
                  data = train,
                  method = "knn",
                  trControl = controlParameters,
                  tuneGrid = parameterGrid
)

predictions <- predict(knnModel, test)
cmatrix <- table(predictions = predictions, actual = test$Group)
accMetric <- accuracyMetric(cmatrix)
mccMet <- mccMetric (cmatrix)

# ROC 
#rocKNN <- roc(test[,1], knnPredictions[,2])
#rocCART <- roc(test[,1], cartPredictions[,2])
#rocSVM <- roc(test[,1], svmPredictions[,2])
#rocRF <- roc(test[,1], rfPredictions[,2])
#rocLog <- roc(test[,1], logPredictions[,2])
#rocEnsemble <- roc(test[,1], ensemblePredct[,2])
#ggroc(list(KNN = rocKNN, CART = rocCART, SVM = rocSVM, RF = rocRF, Logistic = rocLog, Ensemble = rocEnsemble)) +
  #ggtitle("ROC Curves of All Models")

logPred <- prediction(logPredictions[,2], test$Group)
knnPred <- prediction(knnPredictions[,2], test$Group)
cartPred <- prediction(cartPredictions[,2], test$Group)
rfPred <- prediction(rfPredictions[,2], test$Group)
svmPred <- prediction(svmPredictions[,2], test$Group)
ensemblePred <- prediction(ensemblePredct[,2], test$Group)

perfLog <- performance(logPred, "tpr", "fpr")
perfKnn <- performance(knnPred, "tpr", "fpr")
perfCart <- performance(cartPred, "tpr", "fpr")
perfrf <- performance(rfPred, "tpr", "fpr")
perfsvm <- performance(svmPred, "tpr", "fpr")
perfensemble <- performance(ensemblePred, "tpr", "fpr")

plot(perfLog, col = "red")
plot(perfKnn, col = "orange", add = TRUE)
plot(perfCart, col = "yellow", add = TRUE)
plot(perfrf, col = "green", add = TRUE)
plot(perfsvm, col = "blue", add = TRUE)
plot(perfensemble, col = "purple", add = TRUE)
title(main = "Receiver operator charateristic curves")

legend(x = "bottomright", legend=c("Logistic Regression", "K-Nearest Neighbors", "CART", "Random Forest", "Support Vector Machine", "Ensemble"),
       col=c("red", "orange", "yellow", "green", "blue", "purple"), pch=15, cex = 0.50)
#AUROC
aucLog <- performance(logPred, "auc") #0.846
aucKNN <- performance(knnPred, "auc") #0.942
aucCART <- performance(cartPred, "auc") #0.977
aucRF <- performance(rfPred, "auc") #1
aucSVM <- performance(svmPred, "auc")#1
aucEnsemble <- performance(ensemblePred, "auc") #1
