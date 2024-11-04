#!/usr/bin/env Rscript
# source : https://rpubs.com/phamdinhkhanh/389752
library(tidyverse)
library(GenomicRanges)
library(randomForest)
library(xgboost)
library(kernlab)
library(gbm)
library(pROC)
library(caret)
library(doSNOW)
library(doParallel)
args <- commandArgs(trailingOnly = TRUE)


input_file <- args[1]
output_dir <- args[2]
sample_info_file <- args[3]
threads <- as.integer(args[4])

input_table <- read_csv(input_file)
input_table <- as.data.frame(input_table)
rownames(input_table) <- input_table$sample_id
input_table$type <- ifelse(input_table$type == 0, "Healthy", "Cancer")
input_table$type <- relevel(as.factor(input_table$type), "Healthy", levels= c("Healthy", "Cancer"))
input_table <- input_table[, -1] # remove column sample_id
sample_info_table <- read_csv(sample_info_file)
set.seed(123)

train_idx <- createDataPartition(input_table$type, p=0.8,list=FALSE)
trainning <- input_table[train_idx,]

ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 5,
                     verboseIter = FALSE,
                     savePredictions=TRUE,
                     classProbs=TRUE,
                     summaryFunction = twoClassSummary)

# rfGrid <- data.frame(n.trees=200,interaction.depth=1,shrinkage=0.2,n.minobsinnode=5)
# rfGrid <- data.frame(.mtry = 50)
# rfGrid <- expand.grid(nodesize = c(1,2,3,4), .mtry = c(1:8)*10, ntree = c(1:7)*100) # we can only tune .mtry
cl <- makeCluster(threads, type = "SOCK")
registerDoSNOW(cl)
cl <- makePSOCKcluster(threads)
registerDoParallel(cl)
# Register cluster so that caret will know to train in parallel.
# tunegrid=expand.grid(
#               .alpha=seq(0,1,0.2),
#               .lambda=seq(0.0001, 1, length = 100))
# model_glmnet <- caret::train(type ~ .,
#                         data = trainning,
#                         method='glmnet', 
#                         metric = "ROC",
#                         trControl = ctrl,
#                         preProcess = c("corr", "nzv", "center", "scale"),
#                         tuneGrid = tunegrid
#                         )
tunegrid=expand.grid(
              .nprune=seq(1, 20),
              .degree=seq(1, 10))
model_earth <- caret::train(type ~ .,
                        data = trainning, 
                        method='earth', 
                        tuneLength=3,
                        metric = "ROC",
                        trControl = ctrl,
                        preProcess = c("corr", "nzv", "center", "scale"),
                        tuneGrid = tunegrid
                        )
model_rf <- caret::train(type ~ .,
                          data = trainning,
                          method = 'rf',
                          preProcess = c("corr", "nzv", "center", "scale"),
                          metric = "ROC",
                          tuneLength  = 5,
                          trControl = ctrl)
model_svmRadial <- caret::train(type ~ .,
                          data = trainning,
                          method = 'svmRadial',
                          preProcess = c("corr", "nzv", "center", "scale"),
                          metric = "ROC",
                          tuneLength  = 5,
                          nthread = 1,
                          trControl = ctrl)
model_svmLinear <- caret::train(type ~ .,
                          data = trainning,
                          method = 'svmLinear',
                          preProcess = c("corr", "nzv", "center", "scale"),
                          metric = "ROC",
                          tuneLength  = 5,
                          nthread = 1,
                          trControl = ctrl)
model_gbm <- caret::train(type ~ .,
                          data = trainning,
                          method = 'gbm',
                          preProcess = c("corr", "nzv", "center", "scale"),
                          metric = "ROC",
                          tuneLength  = 5,
                          trControl = ctrl)
stopCluster(cl)
models_compare <- resamples(list(MARS=model_earth, RF=model_rf, SVM_radial=model_svmRadial, SVM_linear=model_svmLinear, GBM=model_gbm))
modle_list <- list('MARS'=model_earth ,'RF'=model_rf, 'SVM_radial'=model_svmRadial, 'SVM_linear'=model_svmLinear, 'GBM'=model_gbm)
summary(models_compare)
saveRDS(modle_list, file = paste0(output_dir, "/model_list.rds"))

scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(models_compare, scales=scales)

library(caretEnsemble)
stackControl <- trainControl(method="repeatedcv", 
                             number=10, 
                             repeats=3,
                             savePredictions=TRUE, 
                             classProbs=TRUE)

# Ensemble the predictions of `models` to form a new combined prediction based on glm
stack.glm <- caretStack(modle_list, method="glm", metric="Accuracy", trControl=stackControl)
results <- resamples(models)
summary(results)
print(stack.glm)

test <- input_table[-train_idx,]
# test <- predict(preProcValues, raw_test)
test.tbl <- select(test, c("type"))
stack_predicteds <- predict(stack.glm, newdata=test)
head(stack_predicteds)