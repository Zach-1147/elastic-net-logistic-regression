
rm(list=ls())
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

packageVersion("pROC")

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


## -------------- DATA SPLITTING ----------------------------

set.seed(1717)  #Keep it reproducible

#Sample indices for 25% of data for use as testing

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

#Dry run cross validation to obtain fold IDs to stroe for use throughout analysis.
fold_id_run <- cv.glmnet(X_train, Y_train, nfolds = 10, family="binomial", type.measure = "auc", keep = TRUE)

#Save the fold IDs to a variable to set later
std_foldid <- fold_id_run$foldid

alpha_values <- seq(0, 1, by=0.01) 

#Initialize vectors to store AUC, threshold values, and coefficients for models
train_auc <- list() 
test_auc <- list()

thresholds_train <- list()
thresholds_test <- list()

coefficients_list <- list() #to access coefficients by index
dev_plots <- list() #store pots as well

#For loop for each alpha value
for (i in seq_along(alpha_values)) {
  alpha <- alpha_values[i]
  
  #Fit model using cross-validation
  cvx <- cv.glmnet(X_train, Y_train, nfolds = 10, family="binomial", alpha=alpha, type.measure = "auc", foldid = std_foldid)
  
  dev_plots[[i]] <- plot(cvx) #store current model's plot

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

  #Extract coefficients
  coefficients_for_model <- coef(cvx, s = cvx$lambda.min)
  coefficients_list[[i]] <- coefficients_for_model

}


#Find the index of the best model based on the test and train AUCs
best_model_index <- which.max(test_auc)
best_train_index <- which.max(train_auc)

best_alpha_plot <- dev_plots[[best_model_index]]
best_alpha_plot

#Look at coefficients for best testing model
best_coefficients <- coefficients_list[[best_model_index]]
best_coefficients <- as.data.frame(as.matrix(best_coefficients))
view(best_coefficients)

#Check stats for best test model
best_alpha <- alpha_values[best_model_index]
best_test_auc <- test_auc[best_model_index]
best_lambda <- 
best_alpha
best_test_auc

#And for best in training
best_coefficients_train <- coefficients_list[[best_train_index]]
best_coefficients_train <- as.data.frame(as.matrix(best_coefficients_train))
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


par(mfrow = c(1, 2), cex.main = 1.1, cex.lab = .95, cex.axis = .95)
plot(roc_train, main = "ROC Curve - Training", col = "violet")
abline(h = best_snsp_train[best_indx_train, 1], v = best_snsp_train[best_indx_train, 2], col = "darkorange", lty = 2)
plot(roc_test, main = "ROC Curve - Test", col = "violet")
abline(h = best_snsp_test[best_indx_test, 1], v = best_snsp_test[best_indx_test, 2], col = "darkorange", lty = 2)
par(mfrow=c(1,1))

##### Problem 2 ----

library(GGally)

p2_data <- load("geneexpression2.rda")
gene_exp_data <- get(p2_data)
str(gene_exp_data)
view(gene_exp_data)

# Eigenvalues Using prcomp() and covariance matrix
mat <- var(gene_exp_data)
eigen(mat)$values
prcomp(gene_exp_data)$sdev^2

# Using eigen()

eg.vals <- eigen(mat)$values
eg.vcs <- eigen(mat)$vectors

result_matrix <- matrix(NA, ncol=156, nrow = 158)

result_matrix[157,] <- eg.vals

result_matrix[158,] <- cumsum(eg.vals/sum(eg.vals)*100)

rownames(result_matrix) <- c(colnames(gene_exp_data), "lambda.hat", "cumulative %")

result_matrix <- as.data.frame(result_matrix)

for (i in 1:ncol(result_matrix)) {
  colnames(result_matrix)[i] <- paste("PC", i, sep = "")
}

result_matrix[1:156,1:156] <- eg.vcs

round(result_matrix, 3)

plot(as.ts(eg.vals), ylab="lambda-hat", xlab="PC", main = "The Scree Plot")

result_matrix[,1:2]

pc.scores <- matrix(NA, ncol=156, nrow = 30)

score.fn <- function (i,eigen.vector) {	
  as.numeric(t(eigen.vector)%*%i)}

for (j in 1:156) {
  pc.scores[,j] <- apply(gene_exp_data, 1, score.fn, eigen.vector = eg.vcs[,j])}

var(pc.scores)
diag(var(pc.scores))
eg.vals

plot (pc.scores[,1] , pc.scores[,2], main = " PC2  Vs.  PC1" )
text(pc.scores[, 1], pc.scores[, 2], labels = rownames(gene_exp_data), pos = 3)

# Plotting PCAs using ggfortify (based on ggplot2)
pca_result <- prcomp(gene_exp_data)

pc_scores <- as.data.frame(pca_result$x)

ggplot(pc_scores, aes(x = PC1, y = PC2, label = rownames(pc_scores))) +
  geom_point() +
  geom_text(hjust = 0, vjust = 0) +  # Adjust label position
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()

# Eigenvalues Using prcomp() and covariance matrix
mat <- var(gene_exp_data)
eigen(mat)$values
prcomp(gene_exp_data)$sdev^2

# Using eigen()

eg.vals <- eigen(mat)$values
eg.vcs <- eigen(mat)$vectors

result_matrix <- matrix(NA, ncol=156, nrow = 158)

result_matrix[157,] <- eg.vals

result_matrix[158,] <- cumsum(eg.vals/sum(eg.vals)*100)

rownames(result_matrix) <- c(colnames(gene_exp_data), "lambda.hat", "cumulative %")

result_matrix <- as.data.frame(result_matrix)

for (i in 1:ncol(result_matrix)) {
  colnames(result_matrix)[i] <- paste("PC", i, sep = "")
}

result_matrix[1:156,1:156] <- eg.vcs

round(result_matrix, 3)

plot(as.ts(eg.vals), ylab="lambda-hat", xlab="PC", main = "The Scree Plot")

result_matrix[,1:2]

pc.scores <- matrix(NA, ncol=156, nrow = 30)

score.fn <- function (i,eigen.vector) {	
  as.numeric(t(eigen.vector)%*%i)}

for (j in 1:156) {
  pc.scores[,j] <- apply(gene_exp_data, 1, score.fn, eigen.vector = eg.vcs[,j])}

var(pc.scores)
diag(var(pc.scores))
eg.vals

plot (pc.scores[,1] , pc.scores[,2], main = " PC2  Vs.  PC1" )
text(pc.scores[, 1], pc.scores[, 2], labels = rownames(gene_exp_data), pos = 3)

# Plotting PCAs using ggfortify (based on ggplot2)
pca_result <- prcomp(gene_exp_data)

pc_scores <- as.data.frame(pca_result$x)

ggplot(pc_scores, aes(x = PC1, y = PC2, label = rownames(pc_scores))) +
  geom_point() +
  geom_text(hjust = 0, vjust = 0) +  # Adjust label position
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()

