####------------ summarize variable importance measures --------------------####
# for each condition in fit_list_short
# mean importance across folds
# for permutation importance forest, trees, count variable used in splits

# load fit lists
fit_list_short <- readRDS("Data_analysis/fit_list_short.rds")
fit_list <- readRDS("Data_analysis/fit_list.rds")


get_f <- function(outcome, predictors, timepoint, fl = fit_list) {
  which(fl$outcome == outcome &
          fl$predictors == predictors &
          fl$timepoint == timepoint)
}

get_var_used <- function(fl) {
  fit_f <- readRDS(fl)
  var_used <- predict(fit_f$fit, var.used = "all.trees")$var.used
  res <- data.frame("var_used" = 100 * var_used / sum(var_used))
  res
}

read_results <- function(files_condition, readFUN) {
  res <- lapply(files_condition, function(fl, fu) tryCatch({
    as.matrix(fu(fl))
  }, error = function(e) NULL)
  , fu = readFUN)
  res
}


avg_results <- function(res, mean_res) {
  res_all <- array(do.call(c, res), dim = c(dim(res[[1]]), length(res)))
  mean_res[[f]] <- data.frame(rowMeans(res_all, dims = 2, na.rm = TRUE))
  dimnames(mean_res[[f]]) <- dimnames(res[[1]])
  mean_res
}


setwd("results_MRFs")

files <- list.files()

mean_imp_forest <- mean_imp_trees <- mean_var_used <- vector("list", length = nrow(fit_list_short))

for(f in 1:nrow(fit_list_short)) {

  fs <- get_f(fit_list_short[f, "outcome"],
              fit_list_short[f, "predictors"],
              fit_list_short[f, "timepoint"])

  pats <- paste0("results_MRF_f", fs, ".rds", collapse  = "|")
  files_condition <- grep(pats, files, value = TRUE)
  if(length(files_condition) > 0) {

    imp_forest <- read_results(files_condition, function(fl) readRDS(fl)$importance_forest)
    imp_trees <- read_results(files_condition, function(fl) readRDS(fl)$importance_trees)
    var_used <- read_results(files_condition, get_var_used)

    mean_imp_forest <- avg_results(imp_forest, mean_imp_forest)
    mean_imp_trees <- avg_results(imp_trees, mean_imp_trees)
    mean_var_used <- avg_results(var_used, mean_var_used)

    saveRDS(list("imp_forest" = mean_imp_forest,
                 "imp_trees" = mean_imp_trees,
                 "var_used" = mean_var_used),
            file = "../mean_importance.rds")

  }
}
