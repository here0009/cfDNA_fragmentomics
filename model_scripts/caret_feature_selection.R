#!/usr/bin/env Rscript
# ref: https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
# ref: https://www.machinelearningplus.com/machine-learning/feature-selection/
# ref: https://www.datacareer.ch/blog/ridge-and-lasso-in-r/
library(GenomicRanges)
library(glmnet)
library(caret)
library(pROC)
library(doSNOW)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)


input_file <- args[1]
output_dir <- args[2]
sample_info_file <- args[3]
threads <- as.integer(args[4])

lines <- c()
repeated_num <- 10

get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

output_file <- paste0(output_dir, "/feature_selection.csv")
threads <- 10


cl <- makeCluster(threads, type = "SOCK")
# Register cluster so that caret will know to train in parallel.
registerDoSNOW(cl)
lines <- cbind(lines, paste(sep=",", collapse=",",c('method', 'feature_number', 'features')))
input_table <- read_csv(input_file)
input_table <- as.data.frame(input_table)
rownames(input_table) <- input_table$sample_id
input_table$type <- ifelse(input_table$type == 0, "Healthy", "Cancer")
input_table$type <- relevel(as.factor(input_table$type), "Healthy", levels= c("Healthy", "Cancer"))
input_table <- input_table[, -1] # remove column sample_id
sample_info_table <- read_csv(sample_info_file)



# recursive feature elimination
# define the control using a random forest selection function
subsets <- c(1:10) * 25
control <- rfeControl(functions=rfFuncs, method="repeatedcv", number=10, saveDetails=TRUE, repeats=repeated_num, verbose=TRUE, rerank=TRUE)
# function options : linear regression (in the object lmFuncs), random forests (rfFuncs), naive Bayes (nbFuncs), bagged trees (treebagFuncs) and functions that can be used with caret’s train function (caretFuncs)
# run the RFE algorithm
 
results <- rfe(input_table[,-dim(input_table)[2]], input_table$type, rfeControl=control, preProcess=c("center", "scale","corr", "nzv"), sizes=subsets)
print(results)
fifty_predictors <- unique(results$variables[results$variables$Variables==50,]$var)
lines <- cbind(lines, paste(sep=',', collapse=',', c('rfe_fifty_predictors', as.character(length(fifty_predictors)), paste(sep=';',  collapse=';',fifty_predictors))))
lines <- cbind(lines, paste(sep=',', collapse=',', c('rfe_best', as.character(length(predictors(results))), paste(sep=';', collapse=';',predictors(results)))))


### glmnet feature selection
tunegrid=expand.grid(
              .alpha=seq(0,1,0.2),
              .lambda=seq(0.0001, 1, length = 100))
ctrl <- trainControl(verboseIter = TRUE, classProbs = TRUE, 
                     summaryFunction = twoClassSummary, method = "repeatedcv",
                     repeats = repeated_num)
results <- train(input_table[,-dim(input_table)[2]], input_table$type, method = "glmnet", metric = "ROC", trControl = ctrl,preProcess=c("center", "scale","corr", "nzv"), tuneGrid=tunegrid)
print(results)
lines <- cbind(lines, paste(sep=',', collapse=',', c('glmnet', as.character(length(predictors(results))), paste(sep=';', collapse=';',predictors(results)))))


# Univariate Filters
ctrl <- sbfControl(functions = rfSBF, method = "repeatedcv", repeats = repeated_num, saveDetails=TRUE)
results <- sbf(input_table[,-dim(input_table)[2]], input_table$type, sbfControl=ctrl)
print(results)


lines <- cbind(lines, paste(sep=',', collapse=',', c('univariate-filter', as.character(length(predictors(results))), paste(sep=';', collapse=';',predictors(results)))))
writeLines(lines, output_file, sep='\n')