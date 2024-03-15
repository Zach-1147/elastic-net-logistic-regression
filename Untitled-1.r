

### --------------- SET PARAMETERS ------------------

Input_Dir <- "C:/Users/zachr/Desktop/6970_Asn_3"

Patient_filename <- "Immunologic profiles of patients with COVID-19.xlsx"

### --------------- LIBRARIES ------------------------
library(tidyverse)
library(readxl)
library(ggplot2)
library(DAAG)
library(plotly)
library(glmnet)
library(caret)

## -------------- DATA CLEANING ----------------------------

#Read in data
setwd(Input_Dir) 
data <- as.data.frame(read_xlsx(Patient_filename))

#Convert Sex to numeric binary outcome
data$SEX <- ifelse(data$SEX == 'F', 1, 0)

#Convert response variable to numeric binary outcome, severe = 1.
data$Severirty <- ifelse(data$Severirty == 'Severe', 1, 0)

predictors <- data %>%
dplyr::select(-`Patient Number`, -Severirty)

response <- data$Severirty

#Scale the predictors
predictors <- as.data.frame(scale(predictors))

view(predictors)


## -------------- DATA SPLITTING ----------------------------

set.seed(1717)  #Keep it reproducibile

#Sample indices for 75% of data to uss in training

test_size <- round(dim(predictors)[1] * 0.25)
test_index <- sample.int(n = dim(predictors)[1], size = test_size, replace = FALSE)

#Subset training from selected indices
X_train <- as.matrix(predictors[-test_index, ])
Y_train <- response[-test_index]

#Subset testing data as everything not selected in sampled indices for training
X_test <- as.matrix(predictors[test_index, ])
Y_test <- response[test_index]

view(X_test)

## -------------- MODEL TRAINING ----------------------------

alpha_values <- seq(0, 1, by=0.1) 

#Initialize vectors to store AUC and threshold values for models
train_auc <- NA  
test_auc <- NA

thresholds_train <- NA
thresholds_test <- NA

#For loop for each alpha value
for (i in seq_along(alpha_values)) {
  alpha <- alpha_values[i]
  
  #Fit model using cross-validation
  cvx <- cv.glmnet(X_train, Y_train, nfolds = 10, family="binomial", alpha=alpha, type.measure = "auc")
  
  #Make predictions on both training and test sets
  prds.train <- predict(cvx, newx = X_train, type = "response", s="lambda.min")
  prds.test <- predict(cvx, newx = X_test, type = "response", s="lambda.min")
  
  #Compute and store AUC for training and test sets
  roc_train <- roc(Y_train, prds.train[,1])
  roc_test <- roc(Y_test, prds.test[,1])
  
  #Apply min-max approach for training set to select best threshold
  snsp_train <- cbind(roc_train$sensitivities, roc_train$specificities)
  index_train_thresh <- which.max(apply(snsp_train, 1, min))
  thresholds_train[i] <- roc_train$thresholds[index_train_thresh]
  
  #Apply min-max approach on testing set to select best threshold
  snsp_test <- cbind(roc_test$sensitivities, roc_test$specificities)
  index_test_thresh <- which.max(apply(snsp_test, 1, min))
  thresholds_test[i] <- roc_test$thresholds[index_test_thresh]

  #Extract AUC's
  train_auc[i] <- auc(roc_train)
  test_auc[i] <- auc(roc_test)
}


