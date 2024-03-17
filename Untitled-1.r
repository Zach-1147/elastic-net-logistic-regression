

### --------------- SET PARAMETERS ------------------

Input_Dir <- "C:/Users/zachr/Desktop/6970_Asn_3"

Patient_filename <- "Immunologic profiles of patients with COVID-19.xlsx"

### --------------- LIBRARIES ------------------------
library(tidyverse)
library(readxl)
library(ggplot2)
library(DAAG)
library(pROC)
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

#GET SOME FOLD IDs

n <- nrow(X_train)  # Number of observations
nfolds <- 10  # Number of folds
fold_ids <- sample(rep(seq_len(nfolds), length.out = n))

#CONSIDER USING DEVIANCE

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
  
  #Make predictions on both training and test sets using best lambda value for each alpha
  prds.train <- predict(cvx, newx = X_train, type = "response", s=cvx$lambda.min)
  prds.test <- predict(cvx, newx = X_test, type = "response", s=cvx$lambda.min)
  
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

## OUR FINAL MODEL IS SELECTED AS THE BEST AUC OF ALL THE ALPHA VALUES WITH THEIR RESPECTIVE LAMBDA MINS. FROM HERE, OUR THRESHOLD SELECTION INDICATES THE BEST TRESHOLD FOR THAT SPECIFIC MODEL. THIS IS THE BEST THRESHOLD FOR THE BEST MODEL, AND WE WILL LOOK AT THE COEFFICIENTS TO DISCUSS. 

## WE WILL EXAMINE RESULTS FOR BOTH TRAINING AND TESTING, BUT THE TESING IS OUR FINAL MODEL OF INTEREST FOR DISCUSSION. WE CAN COMPARE TRAINING AND TESTING CURVES FOR INTEREST AND DICSCUSSION (TESTING BETTER PERF THAN TRAINING? --> SMALL DATASET SYMPTOMS.)

#Find the index of the best model based on the test and train AUCs
best_model_index <- which.max(test_auc)
best_train_index <- which.max(train_auc)

#Look at coefficients for best testing model
best_coefficients <- coefficients_list[[best_model_index]]
best_coefficients <- as.data.frame(best_coefficients)

view(best_coefficients)

#Check stats for best test model
best_alpha <- alpha_values[best_model_index]
best_test_auc <- test_auc[best_model_index]
best_alpha
best_test_auc

#And for best in training
best_coefficients_train <- coefficients_list[[best_train_index]]
best_coefficients <- as.data.frame(best_coefficients_train)
view(best_coefficients_train)

best_alpha_train <- alpha_values[best_train_index]
best_train_auc <- test_auc[best_train_index]

best_alpha_train
best_train_auc




### NOTE: WE NEED TO DO LITERATURE RESEARCH ON THE BIOLOGY OF THE CYTOKINES FOR SECOND PART OF FIRST QUESTION



## -------------- PLOTTING MODELS ----------------------------

#Extract the AUC and sensitivity-specificity pairs for the best model
best_auc_train <- train_auc[best_model_index]
best_auc_test <- test_auc[best_model_index]
best_snsp_train <- cbind(roc_train$sensitivities, roc_train$specificities)
best_snsp_test <- cbind(roc_test$sensitivities, roc_test$specificities)
best_indx_train <- which.min(apply(best_snsp_train, 1, function(x) abs(x[1] - x[2])))
best_indx_test <- which.min(apply(best_snsp_test, 1, function(x) abs(x[1] - x[2])))

#Plot ROC curves and selected thresholds for the best model
par(mfrow=c(1,2))
plot(roc_train, main = "ROC Curve - Training", col = "blue")
abline(h = best_snsp_train[best_indx_train, 1], v = best_snsp_train[best_indx_train, 2], col = "red", lty = 2)
plot(roc_test, main = "ROC Curve - Test", col = "blue")
abline(h = best_snsp_test[best_indx_test, 1], v = best_snsp_test[best_indx_test, 2], col = "red", lty = 2)
par(mfrow=c(1,1))




