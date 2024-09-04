####------------------ fit MRFs to CODY data -----------------------------####

# load data
load("Data_preparation/Daten_cody_gold_sub_wide.RData")
head(Daten_wide)

load("Feature_engineering/pca_scores.RData")

# load fit lists
fit_list_short <- readRDS("Data_analysis/fit_list_short.rds")
fit_list <- readRDS("Data_analysis/fit_list.rds")

library(randomForestSRC)
# library(mrfs)
devtools::load_all()
# library(doParallel)
library(doMPI)

####--------------------- define variables for model types ----------------------####

# predictor variables
G <- 6 #number of gold coin tests

# list with entries for each status test
predictors_basic <- predictors_training <- vector("list", length = G)
for (g in 1:G) {
  predictors_basic[[g]] <- c("Score_cody_1", "Grade", "Grade_cody_1")
  predictors_training[[g]] <- c(predictors_basic[[g]],
                                colnames(pca_scores[[g]])[-c(1)]) # without Id
}
predictors_training3 <- vector("list", length = G)
for (g in 1:G) predictors_training3[[g]] <- predictors_training[[3]]
predictors <- list("basic" = predictors_basic, "training" = predictors_training, "training3" = predictors_training3)

# outcome variables
gold_uni <- lapply(1:G, function(g) paste0("Score_gold_sum_", g))
gold_uni_long <- lapply(1:G, function(g) paste0("Score_gold_sum_", 1:G))
#list, all others vector
outcomes <- list("uni" = gold_uni, "multi_long" = gold_uni_long)

####-------------------------- split intro training and test data with cross-validation -----####

outer_folds <- 10
data_splitted <- vector("list", outer_folds)

set.seed(3103)
for (f in 1:nrow(fit_list_short)) {
  outcomes_f <- outcomes[[as.character(fit_list_short[[f, "outcome"]])]][[fit_list_short[f, "timepoint"]]]
  predictors_f <- predictors[[as.character(fit_list_short[[f, "predictors"]])]][[fit_list_short[f, "timepoint"]]]

  Daten_wide_pca <- merge(Daten_wide, pca_scores[[fit_list_short[f, "timepoint"]]], by = "Id", all = TRUE)

  data_complete <- na.omit(Daten_wide_pca[, c("Id", outcomes_f, predictors_f)])
  # z-standardize y
  data_complete[, outcomes_f] <- scale(data_complete[, outcomes_f])

  data_splitted[[f]] <- split_crossval(data_complete, outer_folds)
}


####-------------------------- real fitting with hyperparameter tuning ----------------------####

get_f <- function(outcome, predictors, timepoint, fl = fit_list_short) {
  which(fl$outcome == outcome &
          fl$predictors == predictors &
          fl$timepoint == timepoint)
}


# files <- list.files("results_reg")
# fs <- gsub("results_reg_f|.rds", "", files)
# missing_rows <- setdiff(1:nrow(fit_list), fs)

# cl <- makeCluster(neval)
# registerDoParallel(cl)

cl <- startMPIcluster()
registerDoMPI(cl)

sinkWorkerOutput(paste0("worker_iter_glmnet.out"))

# define chunkSize so that each cluster worker gets a single task chunk
chunkSize <- ceiling(nrow(fit_list)/getDoParWorkers())
mpiopts <- list(chunkSize=chunkSize)


res <-  foreach (f=1:nrow(fit_list), .verbose=T, .packages=c("randomForestSRC", "glmnet", "mrfs"),
                 .inorder=TRUE, .errorhandling="pass", .options.mpi=mpiopts) %dopar% {

                   devtools::load_all()
                   set.seed(1401 + f)

                   outcomes_f <- outcomes[[as.character(fit_list[f, "outcome"])]][[fit_list[f, "timepoint"]]]
                   predictors_f <- predictors[[as.character(fit_list[f, "predictors"])]][[fit_list[f, "timepoint"]]]

                   data_complete <- data_splitted[[get_f(fit_list[f, "outcome"],
                                                         fit_list[f, "predictors"],
                                                         fit_list[f, "timepoint"])]][[fit_list[f, "fold"]]]

                   if(is.null(predictors_f)) { # predictors status for timepoint 1
                     fit <- NULL
                   } else {
                     fit <- fit_1_glmnet(data = data_complete, y_names = outcomes_f, x_names = predictors_f,
                                         multi = length(outcomes_f) > 1)
                   }

                   saveRDS(fit, file = paste0("results_reg/results_reg_f",f,".rds"))
                   fit
                 }

saveRDS(res, file="results_reg.rds")

# stopCluster(cl)
closeCluster(cl)
mpi.quit()

# mse_test_mean <- do.call(cbind, lapply(fits, function(fit) fit$mse_test_mean))
# colnames(mse_test_mean) <- names(y_names)
# mse_test_mean <- cbind(mse_test_mean, gold1_avg = rowMeans(mse_test_mean[, c("gold1_add", "gold1_sub", "gold1_row")]))
#
# mse_long <- data.frame("MSE" = c(mse_test_mean), "outcome" = rep(colnames(mse_test_mean), each = nrow(mse_test_mean)))
# mse_long$outcome <- factor(mse_long$outcome, levels = colnames(mse_test_mean)[c(1,2,6,3:5)])
#
#
# ####--------------------------------- plot MSEs --------------------------------####
# # violin plot
# library(ggplot2)
# library(colorspace)
#
# ggplot(data=mse_long, aes(y=MSE, x=outcome, fill=outcome)) +
#   geom_violin(show.legend=FALSE) +
#   labs(y="MSE", x="Outcome") +
#   scale_x_discrete(labels = c('gold1_add' = "Addition",
#                               'gold1_sub'   = "Subtraction",
#                               "gold1_row" = "Number Rows",
#                               "gold1_multi" = "Multivariate",
#                               "gold1_uni" = "Univariate",
#                               "gold1_avg" = "Uni Average")) +
#   scale_fill_manual(values=qualitative_hcl(ncol(mse_test_mean))) +
#   theme(axis.text=element_text(size=11),
#         axis.title=element_text(size=11))

