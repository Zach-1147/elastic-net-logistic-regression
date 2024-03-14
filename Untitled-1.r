

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
select(-Patient Number, -Severirty)

response <- data$Severirty

#Scale the predictors
predictors <- as.data.frame(scale(predictors))

view(predictors)

#Reconnect response variable to scaled predictors df
data$Severirty <- response

## -------------- DATA SPLITTING ----------------------------

set.seed(1717)  #Keep it reproducibile

#Sample indices for 75% of data to uss in training

test_size <- round(dim(predictors)[1] * 0.25)
test_index <- sample.int(n = dim(predictors)[1], size = test_size, replace = FALSE)

#Subset training from selected indices
X_train <- predictors[-test_index, ]
Y_train <- response[-test_index]

#Subset testing data as everything not selected in sampled indices for training
X_test <- predictors[test_index, ]
Y_test <- response[test_index]

view(X_test)

## -------------- MODEL TRAINING ----------------------------


