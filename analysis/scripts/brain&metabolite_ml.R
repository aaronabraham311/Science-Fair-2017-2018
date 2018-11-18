# Including libraries
library(tidyverse)
library(reshape2)
library(caret)

# Loading in dataset
data <- read.csv("data/brain_volume_mlData.csv")

# Dropping subject id column
drops <- c("Subject.ID")
data <- data[,!(names(data) %in% drops)]

#Setting seed for reproducibility 
set.seed(123)

#Factoring classes and splitting data into testing and training
data$Name <- as.factor(data$Name)

indices <- createDataPartition(data$Name, p = 0.7, list = FALSE) #70% of data will be for training
train <- data[indices,]
test <- data[-indices,]

train <- train[sample(nrow(train)),] # Mixing up rows to remove any patterns

# Control parameters
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

# Random forest
rfmodel <- train(Name ~.,
                 data = train,
                 method = "rf",
                 trControl = controlParameters,
                 tuneGrid = parameterGrid
)

rfPredictions <- predict(rfmodel, test)
confMatrixrf <- table(predictions = rfPredictions, actual = test$Name)
rfAccuracy <- accuracyMetric(confMatrixrf)
rfMCC <- mccMetric(confMatrixrf)

plot(varImp(rfmodel), main = "Variable Importance - Random Forest with Metabolites and Brain Volume")

# Logistic regression
logmodel <- train(Name ~ .,
                  data = train,
                  method = "glm",
                  family = "binomial",
                  trControl = controlParameters)
logPredictions <- predict(logmodel, test)
confMatrixlog <- table(predictions = logPredictions, actual = test$Name)
logAccuracy <- accuracyMetric(confMatrixlog)
logMCC <- mccMetric(confMatrixlog)

#CART Model
CARTmodel <- train(Name ~.,
                   data = train,
                   method = "rpart",
                   trControl = controlParameters
)

rpart.plot::rpart.plot(CARTmodel$finalModel)

cartPredictions <- predict(CARTmodel, test)
confMatrixCART <- table(predictions = cartPredictions, actual = test$Name)
cartAccuracy <- accuracyMetric(confMatrixCART)
cartMCC <- mccMetric(confMatrixCART)

#KNN Model
parameterGrid <- expand.grid(k = c(2,3,4,5,6,7,8,9,10))

knnModel <- train(Name ~ .,
                  data = train,
                  method = "knn",
                  trControl = controlParameters,
                  tuneGrid = parameterGrid)

knnPredictions <- predict(knnModel, test)
confMatrixKNN <- table(predictions = knnPredictions, actual = test$Name)
knnAccuracy <- accuracyMetric(confMatrixKNN)
knnMCC <- mccMetric(confMatrixKNN)

#SVM Model
parameterGrid <- expand.grid(C = c(0.1,0.2,0.3,0.4,0.5,0.6,0.6,0.7,0.8,0.9,1), sigma = 0.5)

svmModel <- train (Name ~ .,
                   data = train,
                   method = "svmRadial",
                   trControl = controlParameters,
                   tuneGrid = parameterGrid)
svmModel <- readRDS("svm_model.RDS")

svmPredictions <- predict(svmModel, test)
confMatrixSVM <- table(predictions = svmPredictions, actual = test$Name)
svmAccuracy <- accuracyMetric(confMatrixSVM)
svmMCC <- mccMetric(confMatrixSVM)

#Ensemble
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

predDF <- data.frame(cartPredictions, knnPredictions, logPredictions, rfPredictions, svmPredictions, class = train$Name)

combineModel <- train(as.factor(class) ~., method = "rf", data = predDF, trControl = controlParameters)

ensemblePredct <- predict(combineModel, test)
confMatrixEnsemble <- table(predictions = ensemblePredct, actual = test$Name)
ensembleAccuracy <- accuracyMetric(confMatrixEnsemble)
ensembleMCC <- mccMetric(confMatrixEnsemble)
