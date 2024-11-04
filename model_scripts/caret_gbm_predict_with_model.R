#!/usr/bin/env Rscript
library(tidyverse)
library(GenomicRanges)
library(caret)
library(pROC)
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
# output_dir <- 'results/models_220803/end_motif_6_model'
# sample_info_file <- "sample_info/sample_info_csv/220708_total_sample_info.csv"
# threads <- 10

input_file <- 'results/end_motif_220805/tcr/end_motif_6.csv'
sample_info_file <- 'sample_info/sample_info_csv/220805_tcr_sample_info.csv'
model_file <- 'results/selected_features/glmnet/breakpoint_motif/model.rds'

input_table <- read_csv(input_file)
input_table <- as.data.frame(input_table)
rownames(input_table) <- input_table$sample_id
input_table$type <- ifelse(input_table$type == 0, "Healthy", "Cancer")
input_table$type <- relevel(as.factor(input_table$type), "Healthy", levels= c("Healthy", "Cancer"))
input_table <- input_table[, -1] # remove column sample_id
sample_info_table <- read_csv(sample_info_file)
set.seed(101)




cl <- makeCluster(threads, type = "SOCK")
# Register cluster so that caret will know to train in parallel.
registerDoSNOW(cl)


stopCluster(cl)

# save model

model_gbm <- readRDS(model_file)

library(gbm)
options(digits=3)

# check the performance of model of test data
# raw_test <- input_table[-train_idx,]
# test <- predict(preProcValues, raw_test)
test.tbl <- select(input_table, c("type"))
test_pred <- predict(model_gbm, input_table[,-dim(input_table)[2]]) # remove last column of test

confusion_table <- data.frame(confusionMatrix(test_pred, input_table$type)$table)
# print(confusion_table)
test.tbl$predicted <- test_pred
test.tbl$sample_id <-rownames(test.tbl)
test.tbl <- inner_join(test.tbl, sample_info_table, by=c("sample_id"="Sample_ID"))
write.csv(test.tbl,file.path(output_dir, paste0("test_predictions.csv")),row.names=FALSE, quote=FALSE)
