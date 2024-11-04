#!/usr/bin/env Rscript
# remove the quotes in the input file
library(tidyverse)
library(GenomicRanges)
library(caret)
library(pROC)
library(doSNOW)
args <- commandArgs(trailingOnly = TRUE)
# input_file <- args[1]
# sample_info_file <- args[2]
# summary_file <- args[3]
# outdir <- args[4]

get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

# input_file <- "results/delfi/total/total_results.tsv"
# sample_info_file <- "sample_info/sample_info_csv/220708_total_sample_info.csv"
# summay_file <- "results/delfi/total/summary_correlations.txt"
# outdir <- "results/delfi/total/model2"

input_file <- "results/delfi_220809/all/total_results.tsv"
sample_info_file <- "sample_info/sample_info_csv/220809_total_sample_info.csv"
summay_file <- "results/delfi_220809/all/summary_correlations.txt"
outdir <- "results/delfi_220809/tcr_control/model"

threads <- 10
dir.create(outdir, showWarnings = FALSE)

input_table <- read_tsv(input_file)
sample_info_table <- read_csv(sample_info_file)
summary_table <- read_tsv(summay_file)
selected_ids <- sample_info_table[sample_info_table$tag == 'Exp',] 
sample_info_table %>% filter(tag == 'Exp' )  %>% filter(type == 'control') 
selected_ids <- (sample_info_table %>% filter(!(tag == 'Exp' & type == 'cancer')))$Sample_ID  # only reserve control and tcr samples
input_table <- input_table %>% filter(sample_id %in% selected_ids)
sample_info_table <- sample_info_table %>% filter(Sample_ID %in% selected_ids)

df.fr3 <- input_table %>% filter(binsize==50)
df.fr3 <- df.fr3 %>% group_by(sample_id) %>% mutate(bin = 1:length(sample_id))
df.fr3 <- inner_join(df.fr3, sample_info_table, by=c("sample_id"="Sample_ID"))

summary.df <- inner_join(summary_table, sample_info_table, by=c("sample_id"="Sample_ID"))
summary.df$`type` = relevel(as.factor(summary.df$`type`), "control")


features.cov <- df.fr3  %>% ungroup() %>%
  select(nfrags_corrected, sample_id, bin) %>%
  spread(sample_id, nfrags_corrected) %>%
  select(-bin) %>% 
  na.omit() %>%
  scale() %>%
  t() %>%
  as.data.frame()

features.short <- df.fr3  %>% ungroup() %>%
  select(short_corrected, sample_id, bin) %>%
  spread(sample_id, short_corrected) %>%
  select(-bin) %>% 
  na.omit() %>%
  scale() %>%
  t() %>%
  as.data.frame()

features.sl <- cbind(features.cov, features.short)
colnames(features.sl) <- c(paste0("total", 1:dim(features.cov)[2]), paste0("short", 1:dim(features.short)[2]))
features.sl$type <- ifelse(summary.df$type == "control", "Healthy", "Cancer")

features <- cbind(features.sl,
                  as.matrix(summary.df %>% ungroup() %>%
                              select(contains("Z Score"))))
set.seed(77)
train_idx <- createDataPartition(features$type, p=0.8,list=FALSE)
trainning <- features[train_idx,]
test <- features[-train_idx,]
test_summary <- summary.df[-train_idx,]

ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     verboseIter = FALSE,
                     savePredictions=TRUE,
                     classProbs=TRUE,
                     summaryFunction = twoClassSummary)

# gbmGrid <- expand.grid(interaction.depth = c(1,3,5), n.trees=(1:4)*50, shrinkage=c(0.1,0.2), n.minobsinnode=c(5,10))
gbmGrid <- expand.grid(interaction.depth = c(1,3,5), n.trees=(1:5)*50, shrinkage=c(0.1,0.2,0.3), n.minobsinnode=c(5,10,15))

cl <- makeCluster(threads, type = "SOCK")
registerDoSNOW(cl)
# model_gbm <- caret::train(type ~ .,
#                           data = trainning,
#                           method = 'gbm',
#                           tuneGrid=data.frame(n.trees=150,
#                                               interaction.depth=3,
#                                               shrinkage=0.1,
#                                               n.minobsinnode=5),
#                           preProcess = c("corr", "nzv"),
#                           trControl = ctrl)

model_gbm <- caret::train(type ~ .,
                          data = trainning,
                          method = 'gbm',
                          preProcess = c("corr", "nzv"),
                          tuneGrid = gbmGrid,
                          trControl = ctrl)
# The final values used for the model were n.trees = 150, interaction.depth =
#  3, shrinkage = 0.2 and n.minobsinnode = 5.
stopCluster(cl)
print(model_gbm)
# trellis.par.set(caretTheme())
png("model_gbm_summary.png")
plot(model_gbm)
dev.off()
# png("model_gbm_kappa.png")
# plot(model_gbm, metric = "ROC", plotType = "level",scales = list(x = list(rot = 90)))
# dev.off()
####### Only short/total coverage
# set.seed(1234)
# model_sl <- caret::train(type ~ .,
#                          data = features.sl,
#                          method = 'gbm',
#                          tuneGrid=data.frame(n.trees=150, 
#                                              interaction.depth=3,
#                                              shrinkage=0.1,
#                                              n.minobsinnode=5),
#                          preProcess = c("corr", "nzv"),
#                          trControl = ctrl)

###### Only z-scores
# features.z <- summary.df %>% ungroup() %>% select(contains("Z Score"))
# summary.df$`type` = relevel(as.factor(summary.df$`type`), "control")
# features.z$type <- ifelse(summary.df$type == "control", "Healthy", "Cancer")
# features.z$type <- relevel(as.factor(summary.df$`type`), "control")
# set.seed(1234)
# model_z <- caret::train(type ~ .,
#                         data = features.z,
#                         method = 'gbm',
#                         tuneGrid=data.frame(n.trees=150,
#                                             interaction.depth=3,
#                                             shrinkage=0.1,
#                                             n.minobsinnode=5),
#                         preProcess = c("corr", "nzv"),
#                         trControl = ctrl)

#### Save
# models.list <- list("all"=model_gbm, "SL"=model_sl, "z"=model_z)
# models.list <- list("all"=model_gbm, "SL"=model_sl)
models.list <- list("all"=model_gbm)
saveRDS(models.list, file.path(outdir, "models_list.rds"))
# model_gbm <- readRDS("models_list.rds")$all
test_type <- relevel(as.factor(test$type), "Healthy", levels= c("Cancer", "Healthy"))
test.tbl <- select(test, c("type"))
test_pred <- predict(model_gbm, test[, -dim(test)[2]])

confusionMatrix(reference = test_type, data = test_pred, mode='everything', positive='Cancer')
confusion_table <- data.frame(confusionMatrix(test_pred, test_type)$table)

test.tbl$predicted <- test_pred
test.tbl$sample_id <-rownames(test.tbl)
# test.tbl <- inner_join(test.tbl, summary.df)



pred.tbl <- model_gbm$pred %>% filter(n.trees==model_gbm$best$n.trees, interaction.depth==model_gbm$best$interaction.depth, shrinkage==model_gbm$best$shrinkage, n.minobsinnode==model_gbm$best$n.minobsinnode) %>%
  group_by(rowIndex) %>% dplyr::summarize(obs=obs[1], Cancer=mean(Cancer))
pred.tbl$sample_id <- rownames(trainning)
pred.tbl <- inner_join(pred.tbl, summary.df)
library(gbm)
options(digits=3)
print(model_gbm$finalModel)
varImp(model_gbm)

print("confusion matrix of trainning data")
confusionMatrix(model_gbm)
print(get_best_result(model_gbm))
## 95% specificity
# cutoff <- (pred.tbl %>% filter(type=="Healthy") %>%
#              arrange(desc(Cancer)))$Cancer[11]
# cutoff98 <- (pred.tbl %>% filter(type=="Healthy") %>%
#                arrange(desc(Cancer)))$Cancer[5]
# ## 90% specificity cutoff to be used in tissue prediction.
# cutoff90 <- (pred.tbl %>% filter(type=="Healthy") %>%
#                arrange(desc(Cancer)))$Cancer[21]

# pred.tbl <- pred.tbl %>%
#   mutate(detected95 = ifelse(Cancer > cutoff, "Detected", "Not detected"),
#          detected98 = ifelse(Cancer > cutoff98, "Detected", "Not detected"),
#          detected90 = ifelse(Cancer > cutoff90, "Detected", "Not detected"),
#          stage = gsub("A|B|C", "", `Stage at Diagnosis`))
# pred.tbl <- pred.tbl %>%
#   mutate(detected95 = ifelse(Cancer > cutoff, "Detected", "Not detected"),
#          detected98 = ifelse(Cancer > cutoff98, "Detected", "Not detected"),
#          detected90 = ifelse(Cancer > cutoff90, "Detected", "Not detected"),
#          )

# write.csv(inner_join(summary.df %>% select(-contains("Z Score")), pred.tbl %>%
#                        select(rowIndex, sample_id, stage, Cancer, detected95, detected98),
#                      by=c("sample_id"="sample_id")),file.path(outdir, paste0("predictions_gbm.csv"),
#                      row.names=FALSE))
write.csv(inner_join(summary.df %>% select(-contains("Z Score")), pred.tbl %>%
                       select(rowIndex, sample_id, Cancer),
                     by=c("sample_id"="sample_id")),file.path(outdir, paste0("tranning_predictions.csv")),
                     row.names=FALSE, quote=FALSE)

write.csv(inner_join(summary.df %>% select(-contains("Z Score")), test.tbl %>%
                       select(sample_id, predicted),
                     by=c("sample_id"="sample_id")),file.path(outdir, paste0("test_predictions.csv")),
                     row.names=FALSE, quote=FALSE)