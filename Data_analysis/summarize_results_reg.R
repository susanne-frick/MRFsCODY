####------------ summarize results --------------------####

res <- readRDS("Data_analysis/results_reg.rds")

G <- 6
# load fit lists
fit_list_short <- readRDS("Data_analysis/fit_list_short.rds")
fit_list <- readRDS("Data_analysis/fit_list.rds")


get_f <- function(outcome, predictors, timepoint, fl = fit_list) {
  which(fl$outcome == outcome &
          fl$predictors == predictors &
          fl$timepoint == timepoint)
}

mean_mse <- vector("list", length = nrow(fit_list_short))

for(f in 1:nrow(fit_list_short)) {

  fs <- get_f(fit_list_short[f, "outcome"],
              fit_list_short[f, "predictors"],
              fit_list_short[f, "timepoint"])

  ids <- res[[fs[1]]]$mse_test[, 1]
  resf <- res[fs]

  resf <-lapply(resf, function(rf)
    tryCatch({
      as.matrix(rf$mse_test[, -c(1)])
    }, error = function(e) NULL)
  )

  res_all <- array(do.call(c, resf), dim = c(dim(res[fs][[1]]$mse_test[, -c(1)]), length(resf)))
  print(fit_list_short[f,])
  print(length(do.call(c, lapply(resf, class))))
  mean_mse[[f]] <- data.frame(ids, rowMeans(res_all, dims = 2, na.rm = TRUE))
  saveRDS(mean_mse, file = "Data_analysis/mean_mse_reg.rds")
}
