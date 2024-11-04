#!/usr/bin/env Rscript
# source : https://rpubs.com/phamdinhkhanh/389752
# source: https://stackoverflow.com/questions/57939453/building-a-randomforest-with-caret
library(tidyverse)
library(GenomicRanges)
library(kernlab) # SVM
library(pROC)
# library(doSNOW)
library(doParallel)
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
# cl <- makeCluster(threads, type = "SOCK")
cl <- makePSOCKcluster(threads)
# Register cluster so that caret will know to train in parallel.
# registerDoSNOW(cl)
registerDoParallel(cl)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 5,
                     verboseIter = FALSE,
                     savePredictions=TRUE,
                     classProbs=TRUE,
                     summaryFunction = twoClassSummary
                     )

# svmGrid <- expand.grid(sigma = seq(0, 0.3, 0.05), C = seq(0.5, 5, 0.25))
svmGrid <- expand.grid(C = seq(0.2, 3, length=30)) # svmlinear, C larger , less prone to overfitting
# svmGrid <- expand.grid(sigma = seq(0, 0.3, 0.1), C = seq(0.5, 5.5, 1))


# medhod: svmLinear, svmPoly, svmRadial
model_svm <- caret::train(type ~ .,
                          data = trainning,
                          method = 'svmLinear',
                          preProcess = c("corr", "nzv", "center", "scale"),
                          metric = "ROC",
                          # tuneLength = 10,
                          tuneGrid = svmGrid,
                          # nthread = 1,
                          trControl = ctrl)

# model_svm <- caret::train(type ~ .,
#                           data = trainning,
#                           method = 'svmPoly',
#                           preProcess = c("corr", "nzv", "center", "scale"),
#                           tuneLength = 5,
#                           tuneGrid = svmGrid,
#                           metric = "ROC",
#                           # nthread = 1,
#                           trControl = ctrl)
stopCluster(cl)
print(model_svm)
# plot summary of model
png(file.path(output_dir, paste0("model_summary.png")))
plot(model_svm)
dev.off()
# save model
saveRDS(model_svm, file.path(output_dir, "model.rds"))
# model_svm <- readRDS('model.rds')
# model_svm <- readRDS(file.path(output_dir, "end_6_models_list.rds"))
# check the performance of model of trainning data
# model_svm$pred %>% filter(model)
# pred.tbl <- model_svm$pred %>% filter(C==model_svm$best$C, sigma==model_svm$best$sigma) %>%
  # group_by(rowIndex) %>% dplyr::summarize(obs=obs[1], Cancer=mean(Cancer))
pred.tbl <- model_svm$pred %>% filter(C==model_svm$best$C) %>% group_by(rowIndex) %>% dplyr::summarize(obs=obs[1], Cancer=mean(Cancer))
pred.tbl$sample_id <- rownames(trainning)
pred.tbl <- inner_join(pred.tbl, sample_info_table, by=c("sample_id"="Sample_ID"))

options(digits=3)
print(model_svm$finalModel)
varImp(model_svm)

print("confusion matrix of trainning data")
confusionMatrix(model_svm)
print(get_best_result(model_svm))
write.csv(pred.tbl,file.path(output_dir, paste0("train_predictions.csv")),row.names=FALSE, quote=FALSE)

# check the performance of model of test data
test <- input_table[-train_idx,]
# test <- predict(preProcValues, raw_test)
test.tbl <- select(test, c("type"))
test_pred <- predict(model_svm, test[,-dim(test)[2]]) # remove last column of test
# model_svm$results %>% filter(n.trees==model_svm$best$n.trees, interaction.depth==model_svm$best$interaction.depth, shrinkage==model_svm$best$shrinkage, n.minobsinnode==model_svm$best$n.minobsinnode)
print("confusion matrix of test data")
confusionMatrix(reference = test$type, data = test_pred, mode='everything', positive='Cancer')
# confusion_table <- data.frame(confusionMatrix(test_pred, test$type)$table)
# print(confusion_table)
test.tbl$predicted <- test_pred
test.tbl$sample_id <-rownames(test.tbl)
test.tbl <- inner_join(test.tbl, sample_info_table, by=c("sample_id"="Sample_ID"))
write.csv(test.tbl,file.path(output_dir, paste0("test_predictions.csv")),row.names=FALSE, quote=FALSE)
