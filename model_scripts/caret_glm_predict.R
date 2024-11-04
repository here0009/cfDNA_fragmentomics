#!/usr/bin/env Rscript
# source : https://rpubs.com/phamdinhkhanh/389752
# source: https://www.kaggle.com/code/pelkoja/visual-xgboost-tuning-with-caret/report
library(tidyverse)
library(GenomicRanges)
library(xgboost) 
library(pROC)
library(kernlab)
# library(doParallel)
library(doSNOW)
library(caret)
args <- commandArgs(trailingOnly = TRUE)


input_file <- args[1]
output_dir <- args[2]
sample_info_file <- args[3]
threads <- as.integer(args[4])

get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}


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
cl <- makeCluster(threads, type = "SOCK")
# cl <- makePSOCKcluster(threads)
# Register cluster so that caret will know to train in parallel.
registerDoSNOW(cl)
# registerDoParallel(cl)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 3,
                     verboseIter = FALSE,
                     savePredictions=TRUE,
                     classProbs=TRUE,
                     summaryFunction = twoClassSummary,
                     allowParallel = FALSE)

# svmGrid <- expand.grid(sigma = seq(0, 0.3, 0.05), C = seq(0.5, 5, 0.25))
svmGrid <- expand.grid(C = seq(0.2, 3, length=30)) # svmlinear, C larger , less prone to overfitting
# svmGrid <- expand.grid(sigma = seq(0, 0.3, 0.1), C = seq(0.5, 5.5, 1))


model_xgbTree <- caret::train(type ~ .,
                          data = trainning,
                          method = 'adaboost',
                          # method = 'xgbLinear',
                          preProcess = c("corr", "nzv", "center", "scale"),
                          metric = "ROC",
                          tuneLength = 10,
                          nthread = 3,
                          trControl = ctrl)
grid_default <- expand.grid(
  nrounds = 100,
  max_depth = 6,
  eta = 0.3,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

train_control <- caret::trainControl(
  method = "none",
  verboseIter = FALSE, # no training log
  allowParallel = TRUE # FALSE for reproducible results 
)


# model_xgbTree <- caret::train(type ~ .,
#                           data = trainning,
#                           method = 'svmPoly',
#                           preProcess = c("corr", "nzv", "center", "scale"),
#                           tuneLength = 5,
#                           tuneGrid = svmGrid,
#                           metric = "ROC",
#                           # nthread = 1,
#                           trControl = ctrl)
stopCluster(cl)
print(model_xgbTree)
# plot summary of model
png(file.path(output_dir, paste0("model_summary.png")))
plot(model_xgbTree)
dev.off()
# save model
saveRDS(model_xgbTree, file.path(output_dir, "model.rds"))
# model_xgbTree <- readRDS('model.rds')
# model_xgbTree <- readRDS(file.path(output_dir, "end_6_models_list.rds"))
# check the performance of model of trainning data
# model_xgbTree$pred %>% filter(model)
# pred.tbl <- model_xgbTree$pred %>% filter(C==model_xgbTree$best$C, sigma==model_xgbTree$best$sigma) %>%
  # group_by(rowIndex) %>% dplyr::summarize(obs=obs[1], Cancer=mean(Cancer))
pred.tbl <- model_xgbTree$pred %>% filter(C==model_xgbTree$best$C) %>% group_by(rowIndex) %>% dplyr::summarize(obs=obs[1], Cancer=mean(Cancer))
pred.tbl$sample_id <- rownames(trainning)
pred.tbl <- inner_join(pred.tbl, sample_info_table, by=c("sample_id"="Sample_ID"))

options(digits=3)
print(model_xgbTree$finalModel)
varImp(model_xgbTree)

print("confusion matrix of trainning data")
confusionMatrix(model_xgbTree)
print(get_best_result(model_xgbTree))
write.csv(pred.tbl,file.path(output_dir, paste0("train_predictions.csv")),row.names=FALSE, quote=FALSE)

# check the performance of model of test data
test <- input_table[-train_idx,]
# test <- predict(preProcValues, raw_test)
test.tbl <- select(test, c("type"))
test_pred <- predict(model_xgbTree, test[,-dim(test)[2]]) # remove last column of test
# model_xgbTree$results %>% filter(n.trees==model_xgbTree$best$n.trees, interaction.depth==model_xgbTree$best$interaction.depth, shrinkage==model_xgbTree$best$shrinkage, n.minobsinnode==model_xgbTree$best$n.minobsinnode)
print("confusion matrix of test data")
confusionMatrix(reference = test$type, data = test_pred, mode='everything', positive='Cancer')
# confusion_table <- data.frame(confusionMatrix(test_pred, test$type)$table)
# print(confusion_table)
test.tbl$predicted <- test_pred
test.tbl$sample_id <-rownames(test.tbl)
test.tbl <- inner_join(test.tbl, sample_info_table, by=c("sample_id"="Sample_ID"))
write.csv(test.tbl,file.path(output_dir, paste0("test_predictions.csv")),row.names=FALSE, quote=FALSE)
