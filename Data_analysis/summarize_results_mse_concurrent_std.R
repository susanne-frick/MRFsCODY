####------------ summarize results --------------------####

# load fit lists
fit_list_short <- readRDS("Data_analysis/fit_list_short.rds")
fit_list <- readRDS("Data_analysis/fit_list.rds")


get_f <- function(outcome, predictors, timepoint, fl = fit_list) {
  which(fl$outcome == outcome &
          fl$predictors == predictors &
          fl$timepoint == timepoint)
}

setwd("results_MRFs_con_std")

files <- list.files()

mean_mse <- vector("list", length = nrow(fit_list_short))

for(f in 1:nrow(fit_list_short)) {

  fs <- get_f(fit_list_short[f, "outcome"],
              fit_list_short[f, "predictors"],
              fit_list_short[f, "timepoint"])

  pats <- paste0("results_MRF_con_std_f", fs, ".rds", collapse  = "|")
  files_condition <- grep(pats, files, value = TRUE)
  if(length(files_condition) > 0) {
    ids <- readRDS(files_condition[1])$mse_test[, "Id"]
    res <- lapply(files_condition, function(fl) tryCatch({
      as.matrix(readRDS(fl)$mse_test[, -c(1)])
    }, error = function(e) NULL)
    )
    res_all <- array(do.call(c, res), dim = c(dim(res[[1]]), length(res)))
    print(fit_list_short[f,])
    print(length(res))
    mean_mse[[f]] <- data.frame(ids, rowMeans(res_all, dims = 2, na.rm = TRUE))
    saveRDS(mean_mse, file = "../mean_mse_MRFs_con_std.rds")

    fits <- lapply(files_condition, function(fl) tryCatch({
      readRDS(fl)$fit
    }, error = function(e) NULL)
    )
    saveRDS(fits, file = paste0("../results_MRFs_fit_con_std/fit_MRF_con_std_f", f, ".rds"))
  }
}
