
rm(list=ls())


### ----------------- ###
### --- PROBLEM 1---- ###
###------------------ ###

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


## -------------- MODEL TRAINING 10-FOLD CV ----------------------------

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

lambda_mins <- list() #store lambda values

#For loop for each alpha value
for (i in seq_along(alpha_values)) {
  
  suppressMessages({
  alpha <- alpha_values[i]
  
  #Fit model using cross-validation
  cvx <- cv.glmnet(X_train, Y_train, nfolds = 10, family="binomial", alpha=alpha, type.measure = "auc", foldid = std_foldid)
  
  plot(cvx, main="Binomial Deviance vs. Lambda: 10-Fold Cross-Validation with Elastic Net Regularization", line = 2.5)
  dev_plots[[i]] <- recordPlot()

  
  
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
  
  lambda_mins[[i]] <- cvx$lambda.min #Append lambda min for each iteration

  #Extract coefficients
  coefficients_for_model <- coef(cvx, s = cvx$lambda.min)
  coefficients_list[[i]] <- coefficients_for_model

  }
)}


#Find the index of the best model based on the test and train AUCs
best_model_index <- which.max(test_auc)
best_train_index <- which.max(train_auc)

best_alpha_plot <- dev_plots[[best_model_index]]

#Look at coefficients for best testing model
best_coefficients <- coefficients_list[[best_model_index]]
best_coefficients <- as.data.frame(as.matrix(best_coefficients))
view(best_coefficients)

#Check stats for best test model
best_alpha <- alpha_values[best_model_index]
best_test_auc <- test_auc[best_model_index]
best_lambda <- lambda_mins[[best_model_index]]
best_lambda
best_alpha
best_test_auc

#And for best in training
best_coefficients_train <- coefficients_list[[best_train_index]]
best_coefficients_train <- as.data.frame(as.matrix(best_coefficients_train))

best_alpha_train <- alpha_values[best_train_index]
best_train_auc <- test_auc[best_train_index]

best_alpha_train
best_train_auc


## -------------- PLOTTING MODELS FOR 10-FOLD CV ----------------------------


replayPlot(best_alpha_plot)


#Extract the AUC and sensitivity-specificity pairs for the best model
best_auc_train <- train_auc[best_model_index]
best_auc_test <- test_auc[best_model_index]
best_snsp_train <- cbind(roc_train$sensitivities, roc_train$specificities)
best_snsp_test <- cbind(roc_test$sensitivities, roc_test$specificities)
best_indx_train <- which.min(apply(best_snsp_train, 1, function(x) abs(x[1] - x[2])))
best_indx_test <- which.min(apply(best_snsp_test, 1, function(x) abs(x[1] - x[2])))

par(mfrow = c(1, 2), cex.main = 1.1, cex.lab = .95, cex.axis = .95)
plot(roc_train, main = "ROC Curve - Training - 10-fold CV", col = "violet")
abline(h = best_snsp_train[best_indx_train, 1], v = best_snsp_train[best_indx_train, 2], col = "darkorange", lty = 2)
plot(roc_test, main = "ROC Curve - Test - 10-fold CV", col = "violet")
abline(h = best_snsp_test[best_indx_test, 1], v = best_snsp_test[best_indx_test, 2], col = "darkorange", lty = 2)
par(mfrow=c(1,1))

## -------------- MODEL TRAINING 20-FOLD CV ----------------------------

#Dry run cross validation to obtain fold IDs to store for use throughout analysis, for 20 fold now.
fold_id_run <- cv.glmnet(X_train, Y_train, nfolds = 20, family="binomial", type.measure = "auc", keep = TRUE)

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

lambda_mins <- list() #store lambda values

#For loop for each alpha value
for (i in seq_along(alpha_values)) {
  suppressMessages({
  alpha <- alpha_values[i]
  
  #Fit model using cross-validation
  cvx <- cv.glmnet(X_train, Y_train, nfolds = 20, family="binomial", alpha=alpha, type.measure = "auc", foldid = std_foldid)
  
  plot(cvx, main="Binomial Deviance vs. Lambda: 20-Fold Cross-Validation with Elastic Net Regularization", line = 2.5)
  dev_plots[[i]] <- recordPlot()
  
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
  
  lambda_mins[[i]] <- cvx$lambda.min #Append lambda min for each iteration
  
  #Extract coefficients
  coefficients_for_model <- coef(cvx, s = cvx$lambda.min)
  coefficients_list[[i]] <- coefficients_for_model
  
  }
)}

#Find the index of the best model based on the test and train AUCs
best_model_index_twf <- which.max(test_auc)
best_train_index_twf <- which.max(train_auc)

best_alpha_plot_twf <- dev_plots[[best_model_index_twf]]

#Look at coefficients for best testing model
best_coefficients_twf <- coefficients_list[[best_model_index_twf]]
best_coefficients_twf <- as.data.frame(as.matrix(best_coefficients_twf))
view(best_coefficients_twf)

#Check stats for best test model
best_alpha_twf <- alpha_values[best_model_index_twf]
best_test_auc_twf <- test_auc[best_model_index_twf]
best_lambda_twf <- lambda_mins[[best_model_index_twf]]
best_lambda_twf
best_alpha_twf
best_test_auc_twf

#And for best in training
best_coefficients_train <- coefficients_list[[best_train_index]]
best_coefficients_train <- as.data.frame(as.matrix(best_coefficients_train))

best_alpha_train <- alpha_values[best_train_index]
best_train_auc <- test_auc[best_train_index]

best_alpha_train
best_train_auc


## -------------- PLOTTING MODELS FOR 20-FOLD CV ----------------------------

replayPlot(best_alpha_plot)

#Extract the AUC and sensitivity-specificity pairs for the best model
best_auc_train <- train_auc[best_model_index]
best_auc_test <- test_auc[best_model_index]
best_snsp_train <- cbind(roc_train$sensitivities, roc_train$specificities)
best_snsp_test <- cbind(roc_test$sensitivities, roc_test$specificities)
best_indx_train <- which.min(apply(best_snsp_train, 1, function(x) abs(x[1] - x[2])))
best_indx_test <- which.min(apply(best_snsp_test, 1, function(x) abs(x[1] - x[2])))

par(mfrow = c(1, 2), cex.main = 1.1, cex.lab = .95, cex.axis = .95)
plot(roc_train, main = "ROC Curve - Training - 20-fold CV", col = "violet")
abline(h = best_snsp_train[best_indx_train, 1], v = best_snsp_train[best_indx_train, 2], col = "darkorange", lty = 2)
plot(roc_test, main = "ROC Curve - Test - 20-fold CV", col = "violet")
abline(h = best_snsp_test[best_indx_test, 1], v = best_snsp_test[best_indx_test, 2], col = "darkorange", lty = 2)
par(mfrow=c(1,1))


### ----------------- ###
### --- PROBLEM 2---- ###
###------------------ ###

library(GGally)

# Load the RDA data file and get the data into an R object.
p2_data <- load("geneexpression2.rda")
gene_exp_data <- get(p2_data)

# Analyze the structure of the data take a look at the data by using the view() function.
str(gene_exp_data)
view(gene_exp_data)

# We will construct PCAs using ggfortify (based on ggplot2).
# Use the prcomp() function to run PCA on the gene_exp_data.
pca_result <- prcomp(gene_exp_data)

# Extract the PC scores as a data frame.
pc_scores <- as.data.frame(pca_result$x)

# Extract the eigenvalues for the scree plot.
eigenvalues <- pca_result$sdev^2

# Construct the scree plot.
plot(1:length(eigenvalues), eigenvalues, 
     type = "b", 
     xlab = "Principal Component", 
     ylab = "Eigenvalue",
     main = "Scree Plot")

# Extract the subject health and cell types information from the names of the observations and add them to new columns in the data frame.
pc_scores$health <- substr(rownames(gene_exp_data), 1, 3)
pc_scores$cell_type <- sapply(strsplit(rownames(gene_exp_data), "_"), function(x) x[2])

# Calculate the total variances explained by PC1 and PC2.
variance_pc1 <- round(100 * (pca_result$sdev[1]^2 / sum(pca_result$sdev^2)), 2)
variance_pc2 <- round(100 * (pca_result$sdev[2]^2 / sum(pca_result$sdev^2)), 2)

# Construct the PCA plot, plotting PC1 on the x-axis and PC2 on the y-axis. Differentiate the healthy and melanoma patients by shape. Differentiate the cell type by color.
ggplot(pc_scores, aes(x = PC1, y = PC2, shape = health, color = cell_type)) +
  geom_point(size = 2) +
  labs(x = paste("PC1 (", variance_pc1, "% of total variance)", sep = ""),
       y = paste("PC2 (", variance_pc2, "% of total variance)", sep = ""),
       title = "PC1 vs PC2") +
  scale_shape_manual(values = c(16, 17)) +
  scale_color_manual(values = c("red", "blue", "green")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# We will construct another PCA using the eigen() function directly. We want to confirm that the PCA plots are similar to verify that the PCA was performed properly.
# Construct the covariance matrix that we will work with.
mat <- var(gene_exp_data)

# Using the eigen() function, we will extract the eigenvalues and eigenvectors.
eg.vals <- eigen(mat)$values
eg.vcs <- eigen(mat)$vectors

# Initialize the matrix to contain our results.
result_matrix <- matrix(NA, ncol=156, nrow = 158)

# Add the eigenvalues into the 157th row.
result_matrix[157,] <- eg.vals

# Add the cumulative sum of PC variances into the last row.
result_matrix[158,] <- cumsum(eg.vals/sum(eg.vals)*100)

# Update the row names for the result_matrix. The names of the differentially expressed genes are used to name the first 156 rows.
rownames(result_matrix) <- c(colnames(gene_exp_data), "lambda.hat", "cumulative %")

# Convert to a data frame.
result_matrix <- as.data.frame(result_matrix)

# This for loop renames the columns of the results_matrix to "PC" then the number of the PC.
for (i in 1:ncol(result_matrix)) {
  colnames(result_matrix)[i] <- paste("PC", i, sep = "")
}

# Input the eigenvectors into the result_matrix.
result_matrix[1:156,1:156] <- eg.vcs

# Round the values in the matrix to the thousandth decimal place.
round(result_matrix, 3)

# Construct the scree plot to analyze the variance explained by the PC's.
plot(as.ts(eg.vals), ylab="lambda-hat", xlab="PC", main = "The Scree Plot")

# Initialize the matrix to contain the PC scores for the observations.
pc.scores <- matrix(NA, ncol=156, nrow = 30)

# Write a function that calculate the PC scores of the observations.
score.fn <- function (i,eigen.vector) {	
  as.numeric(t(eigen.vector)%*%i)}

# This for loop uses the function created above to calculate the PC scores.
for (j in 1:156) {
  pc.scores[,j] <- apply(gene_exp_data, 1, score.fn, eigen.vector = eg.vcs[,j])}

# Check the computation by comparing the eigenvalues.
var(pc.scores)
diag(var(pc.scores))
eg.vals

#Construct the PCA plot and compare to the previous PCA plot. Analyze the point labels to ensure that the PCA plots are similar.
plot (-pc.scores[,1] , -pc.scores[,2], main = " PC2  Vs.  PC1")
text(-pc.scores[, 1], -pc.scores[, 2], labels = rownames(gene_exp_data), pos = 3)
