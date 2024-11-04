#!/usr/bin/env Rscript
# source : https://rpubs.com/phamdinhkhanh/389752
# source: https://stackoverflow.com/questions/57939453/building-a-randomforest-with-caret
library(tidyverse)
library(GenomicRanges)
library(randomForest)
library(pROC)
library(caret)
library(doSNOW)
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
# input_file <- 'results/end_motif_2/end_motif_6.csv'
# output_dir <- 'results/models_220803/end_motif_6_model/rf'
# sample_info_file <- "sample_info/sample_info_csv/220708_total_sample_info.csv"
# threads <- 10

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
rfGrid <- expand.grid(.mtry = c(1:10)*6)
cl <- makeCluster(threads, type = "SOCK")
# Register cluster so that caret will know to train in parallel.
registerDoSNOW(cl)
# grid search
# model_rf <- caret::train(type ~ .,
#                           data = trainning,
#                           method = 'rf',
#                           preProcess = c("corr", "nzv", "center", "scale"),
#                           tuneGrid = rfGrid,
#                           metric = "ROC",
#                           trControl = ctrl)
# random search
model_rf <- caret::train(type ~ .,
                          data = trainning,
                          method = 'rf',
                          preProcess = c("corr", "nzv", "center", "scale"),
                          metric = "ROC",
                          tuneGrid = rfGrid,
                          # tuneLength  = 5,
                          trControl = ctrl)
stopCluster(cl)
print(model_rf)
# plot summary of model
png(file.path(output_dir, paste0("model_summary.png")))
plot(model_rf)
dev.off()
# save model
saveRDS(model_rf, file.path(output_dir, "model.rds"))
# model_rf <- readRDS('model.rds')
# model_rf <- readRDS(file.path(output_dir, "end_6_models_list.rds"))
# check the performance of model of trainning data
# model_rf$pred %>% filter(model)
pred.tbl <- model_rf$pred %>% filter(mtry==model_rf$best$mtry) %>%
  group_by(rowIndex) %>% dplyr::summarize(obs=obs[1], Cancer=mean(Cancer))

pred.tbl$sample_id <- rownames(trainning)
pred.tbl <- inner_join(pred.tbl, sample_info_table, by=c("sample_id"="Sample_ID"))

options(digits=3)
print(model_rf$finalModel)
varImp(model_rf)

print("confusion matrix of trainning data")
confusionMatrix(model_rf)
print(get_best_result(model_rf))
write.csv(pred.tbl,file.path(output_dir, paste0("train_predictions.csv")),row.names=FALSE, quote=FALSE)

# check the performance of model of test data
test <- input_table[-train_idx,]
# test <- predict(preProcValues, raw_test)
test.tbl <- select(test, c("type"))
test_pred <- predict(model_rf, test[,-dim(test)[2]]) # remove last column of test
# model_rf$results %>% filter(n.trees==model_rf$best$n.trees, interaction.depth==model_rf$best$interaction.depth, shrinkage==model_rf$best$shrinkage, n.minobsinnode==model_rf$best$n.minobsinnode)
print("confusion matrix of test data")
confusionMatrix(reference = test$type, data = test_pred, mode='everything', positive='Cancer')
# confusion_table <- data.frame(confusionMatrix(test_pred, test$type)$table)
# print(confusion_table)
test.tbl$predicted <- test_pred
test.tbl$sample_id <-rownames(test.tbl)
test.tbl <- inner_join(test.tbl, sample_info_table, by=c("sample_id"="Sample_ID"))
write.csv(test.tbl,file.path(output_dir, paste0("test_predictions.csv")),row.names=FALSE, quote=FALSE)
