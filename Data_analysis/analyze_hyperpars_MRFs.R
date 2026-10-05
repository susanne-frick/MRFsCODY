####------------- analyze hyperparameter selections ----------------####

library(psych)

setwd("~/Dokumente/FAIR/Reha/")

hyperpars <- readRDS("hyperpars_MRFs_neval100_k5.rds")
hyperpars_mse <- readRDS("hyperpars_mse_MRFs_neval100_k5.rds")

# check
lapply(hyperpars, dim)
do.call(c, lapply(hyperpars, is.null))

# re-build design list
G <- 6
fit_list <- expand.grid("timepoint" = 1:G,
                        "outcome" = c("multi", "uni", "row", "sub", "add"),
                        "predictors" = c("basic", "training", "status"))

####-------------------- histograms of hyperparameter distributions -----------------------------####
plot_hist_hyper <- function(hyper, selected = FALSE) {
  par(mfrow = c(ncol(hyper[[1]]) - 1, length(hyper)))
  for(p in colnames(hyper[[1]])){
    if(p == "mse_mean") next
    for(o in names(hyper)) {
        print(hist(hyper[[o]][, p], xlab = p, ylab = "Frequency", main = o))
    }
  }
  par(mfrow = c(1,1))
}

# basic (across timepoints, because the predictor set stay constant)
# average across timepoints
hyperpars_basic <- lapply(unique(fit_list$outcome), function(o, fb, hy) {
  do.call(rbind, hy[which((fb$outcome == o) & (fb$predictors == "basic"))])
},fb = fit_list, hy = hyperpars)
names(hyperpars_basic) <- unique(fit_list$outcome)

plot_hist_hyper(hyperpars_basic)

# training, timepoint 6
hyperpars_training_T6 <- hyperpars[which((fit_list$timepoint == 6) & (fit_list$predictors == "training"))]
names(hyperpars_training_T6) <- fit_list[which((fit_list$timepoint == 6) & (fit_list$predictors == "training")), "outcome"]
plot_hist_hyper(hyperpars_training_T6)

# training, timepoint 1
hyperpars_training_T1 <- hyperpars[which((fit_list$timepoint == 1) & (fit_list$predictors == "training"))]
names(hyperpars_training_T1) <- fit_list[which((fit_list$timepoint == 1) & (fit_list$predictors == "training")), "outcome"]
plot_hist_hyper(hyperpars_training_T1)

# status, timepoint 6
hyperpars_status_T6 <- hyperpars[which((fit_list$timepoint == 6) & (fit_list$predictors == "status"))]
names(hyperpars_status_T6) <- fit_list[which((fit_list$timepoint == 6) & (fit_list$predictors == "status")), "outcome"]
plot_hist_hyper(hyperpars_status_T6)

####------------------------------ bivariate plots: MSE on hyperparamters --------------------####
plot_mse_hyper <- function(hyper_mse, hyper = NULL, pars = c("ntree", "mtry", "nodesize")) {
  par(mfrow = c(length(pars), length(hyper_mse)))
  for(p in pars){
    if(p == "mse_mean") next
    for(o in names(hyper_mse)) {
      print(plot(hyper_mse[[o]][, p], hyper_mse[[o]][, "mse_mean"], xlab = p, ylab = "MSE", main = o))
      if (is.null(hyper) == FALSE) {
        print(points(hyper[[o]][, p], hyper[[o]][, "mse_mean"], col = "red"))
      }
    }
  }
  par(mfrow = c(1,1))
}

plot_mse_hyper(hyperpars_basic)
plot_mse_hyper(hyperpars_status_T6)

# with all mses
hyperpars_mse_basic <- lapply(unique(fit_list$outcome), function(o, fb, hy) {
  do.call(rbind, hy[which((fb$outcome == o) & (fb$predictors == "basic"))])
},fb = fit_list, hy = hyperpars_mse)
names(hyperpars_mse_basic) <- unique(fit_list$outcome)

# status, timepoint 6
hyperpars_mse_status_T6 <- hyperpars_mse[which((fit_list$timepoint == 6) & (fit_list$predictors == "status"))]
names(hyperpars_mse_status_T6) <- fit_list[which((fit_list$timepoint == 6) & (fit_list$predictors == "status")), "outcome"]

plot_mse_hyper(hyperpars_mse_basic)
plot_mse_hyper(hyperpars_mse_status_T6)

plot_mse_hyper(hyperpars_mse_basic, hyperpars_basic)
plot_mse_hyper(hyperpars_mse_status_T6, hyperpars_status_T6)
