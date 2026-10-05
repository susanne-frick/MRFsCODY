###------------------------- partial dependence plots for univariate forests ------------------------####

library(ggplot2)
library(gridExtra)
devtools::load_all()
devtools::load_all("~/Dokumente/packages/DataAnalysisSimulation/")

dir_manuscript <- "~/Dokumente/FAIR/Reha/paper/mrfs_cody_man/"

G <- 6
# load fit lists
fit_list_short <- readRDS("Data_analysis/fit_list_short.rds")
fit_list <- readRDS("Data_analysis/fit_list.rds")

get_f <- function(outcome, predictors, timepoint, fl = fit_list_short) {
  which(fl$outcome == outcome &
          fl$predictors == predictors &
          fl$timepoint == timepoint)
}

# partial dependence plots, only for those variables with high importance values
# get the whole results object

####------------------------ variable importance to select relevant predictors -------------------------#####

imp <- readRDS("Data_analysis/mean_importance_con_std.rds")

# check
length(imp)
length(imp$imp_forest)

# prepare
imp <- lapply(1:length(imp$imp_forest), function(i, ip) {
  res <- data.frame(ip$imp_trees[[i]], ip$imp_forest[[i]], ip$var_used[[i]])
  if(nrow(res) > 0) {
    colnames(res) <- c(paste0("imp_trees_", colnames(ip$imp_trees[[i]])),
                       paste0("imp_forest_", colnames(ip$imp_trees[[i]])),
                       "var_used")
  }
  res
}, ip = imp)

rm_0 <- function(imp, cut = 0) round(imp[rowSums(round(imp[, -c(ncol(imp))], 2) <= cut) < (ncol(imp) - 1),], 2)

prep_imp <- function(tp, imp, cut, outcome = "multi_long", predictors = "training") {
  imp_t <- rm_0(imp[[get_f(outcome, predictors, tp)]], cut = cut)
  rownames(imp_t) <- gsub("_", " ", rownames(imp_t))
  if(outcome == "multilong") {
    colnames(imp_t) <-  c(paste0("T", rep(1:6, 2)), "\\% ")
  } else {
    colnames(imp_t) <-  c(paste0("T", rep(tp, 2)), "\\% ")
  }
  imp_t
}

imp_uni_training3 <- lapply(4:6, prep_imp, imp = imp, cut = 0,
                              outcome = "uni", predictors = "training3") #cut = .05
imp_uni_training3

imp_uni_basic <- lapply(1:6, prep_imp, imp = imp, cut = 0,
                            outcome = "uni", predictors = "basic") #cut = .05
imp_uni_basic


name_variables <- function(vec) {
  vec <- gsub("up", "Up", vec)
  vec <- gsub("leveldif", "Level Diff.", vec)
  vec <- gsub("score", "Level", vec)
  vec <- gsub("time", "Time", vec)
  vec <- gsub("discrep", "Discrepancy", vec)
  vec <- gsub("Grade cody 1", "Pretest x Grade", vec)
  return(vec)
}
#
# imp_uni_training3 <- lapply(imp_uni_training3, function(imp) {
#   rownames(imp) <- name_variables(rownames(imp))
#   imp
# })
#
# imp_uni_basic <- lapply(imp_uni_basic, function(imp) {
#   rownames(imp) <- name_variables(rownames(imp))
#   imp
# })


####-------------------------------- training until T3 ----------------------------------------####
# pd_plot delivers mean prediction across folds if a list of fits is given
# fit_f_single is needed to supply the dataset

plots_training3 <- vector("list", 3)

for (tp in c(4,5,6)) {
  fit_f <- readRDS(paste0("results_MRFs_con_std/fit_MRF_con_std_f", get_f("uni", "training3", tp), ".rds"))
  fit_f_single <- readRDS(paste0("results_MRFs_con_std/results_MRF_con_std_f", get_f("uni", "training3", tp, fl = fit_list)[1], ".rds"))

  plot_list <- pd_plot(fit = fit_f,
                       data_plot = rbind(fit_f_single$data$train, fit_f_single$data$test),
                       x_plot = c("Grade_cody_1", "score_PC1", "leveldif_PC2"),
                       legend_limits = c(-1.5, 1.5), legend_name = "Status Test")


  x_y_names <- plot_list$x_y_names
  x_y_names[, "Predictor1"] <- recode.df(x_y_names[, "Predictor1"],
                                         c("Grade", "Grade_cody_1",
                                           "leveldif_PC1", "leveldif_PC2", "leveldif_PC3",
                                           "score_PC1", "score_PC2"),
                                         c("Grade", "Grade x Pre-Test",
                                           "Level Difference PC1", "Level Difference PC2", "Level Difference PC3",
                                           "Level PC1", "Level PC2"))
  x_y_names[, "Predictor2"] <- recode.df(x_y_names[, "Predictor2"],
                                         c("Grade", "Grade_cody_1",
                                           "leveldif_PC1", "leveldif_PC2", "leveldif_PC3",
                                           "score_PC1", "score_PC2"),
                                         c("Grade", "Grade x Pre-Test",
                                           "Level Difference PC1", "Level Difference PC2", "Level Difference PC3",
                                           "Level PC1", "Level PC2"))

  # legend only on last plot, add title
  for(p in 1:length(plot_list$plot_list)) {
    plot_list$plot_list[[p]] <- plot_list$plot_list[[p]] + ggtitle(paste0("O(T", tp, ")")) +
      xlab(x_y_names[p, "Predictor1"]) + ylab(x_y_names[p, "Predictor2"])
    if(! (tp %in% 6)) plot_list$plot_list[[p]] <- plot_list$plot_list[[p]] + theme(legend.position = "none")
  }

  plots_training3[[tp - 3]] <- plot_list
}

plots_training3_unlist <- list(plots_training3[[1]]$plot_list[[1]],
                               plots_training3[[1]]$plot_list[[2]],
                               plots_training3[[1]]$plot_list[[3]],
                               plots_training3[[2]]$plot_list[[1]],
                               plots_training3[[2]]$plot_list[[2]],
                               plots_training3[[2]]$plot_list[[3]],
                               plots_training3[[3]]$plot_list[[1]],
                               plots_training3[[3]]$plot_list[[2]],
                               plots_training3[[3]]$plot_list[[3]])

ggsave(marrangeGrob(plots_training3_unlist,
                    layout_matrix = matrix(1:9, ncol = 3, nrow = 3, byrow = FALSE),
                    top = NULL, widths = c(rep(1,2), 1.4)),
       filename = paste0("plots/PD_plot_uni_training3.jpeg"),
       width = 23, height = 20, units = "cm")

file.copy(from = "plots/PD_plot_uni_training3.jpeg",
          to = paste0(dir_manuscript, "figures"),
          overwrite = TRUE)

####-------------------------------- basic predictors ----------------------------------------####
# pd_plot delivers mean prediction across folds if a list of fits is given
# fit_f_single is needed to supply the dataset

plots_basic <- vector("list", 6)

for (tp in c(1:6)) {
  fit_f <- readRDS(paste0("results_MRFs_con_std/fit_MRF_con_std_f", get_f("uni", "basic", tp), ".rds"))
  fit_f_single <- readRDS(paste0("results_MRFs_con_std/results_MRF_con_std_f", get_f("uni", "basic", tp, fl = fit_list)[1], ".rds"))

  plot_list <- pd_plot(fit = fit_f,
                       data_plot = rbind(fit_f_single$data$train, fit_f_single$data$test),
                       x_plot = c("Score_cody_1", "Grade", "Grade_cody_1"),
                       legend_limits = c(-1.5, 1.5), legend_name = "Status Test")


  x_y_names <- plot_list$x_y_names
  x_y_names[, "Predictor1"] <- recode.df(x_y_names[, "Predictor1"],
                                         c("Grade", "Grade_cody_1", "Score_cody_1",
                                           "leveldif_PC1", "leveldif_PC2", "leveldif_PC3",
                                           "score_PC1", "score_PC2"),
                                         c("Grade", "Grade x Pre-Test", "Pre-Test",
                                           "Level Difference PC1", "Level Difference PC2", "Level Difference PC3",
                                           "Level PC1", "Level PC2"))
  x_y_names[, "Predictor2"] <- recode.df(x_y_names[, "Predictor2"],
                                         c("Grade", "Grade_cody_1", "Score_cody_1",
                                           "leveldif_PC1", "leveldif_PC2", "leveldif_PC3",
                                           "score_PC1", "score_PC2"),
                                         c("Grade", "Grade x Pre-Test", "Pre-Test",
                                           "Level Difference PC1", "Level Difference PC2", "Level Difference PC3",
                                           "Level PC1", "Level PC2"))

  # legend only on last plot, add title
  for(p in 1:length(plot_list$plot_list)) {
    plot_list$plot_list[[p]] <- plot_list$plot_list[[p]] + ggtitle(paste0("O(T", tp, ")")) +
      xlab(x_y_names[p, "Predictor1"]) + ylab(x_y_names[p, "Predictor2"])
    if(! (tp %in% 6)) plot_list$plot_list[[p]] <- plot_list$plot_list[[p]] + theme(legend.position = "none")
  }

  plots_basic[[tp]] <- plot_list
}

plots_basic_unlist <- list(plots_basic[[1]]$plot_list[[1]],
                           plots_basic[[1]]$plot_list[[2]],
                           plots_basic[[1]]$plot_list[[3]],
                           plots_basic[[2]]$plot_list[[1]],
                           plots_basic[[2]]$plot_list[[2]],
                           plots_basic[[2]]$plot_list[[3]],
                           plots_basic[[3]]$plot_list[[1]],
                           plots_basic[[3]]$plot_list[[2]],
                           plots_basic[[3]]$plot_list[[3]],
                           plots_basic[[4]]$plot_list[[1]],
                           plots_basic[[4]]$plot_list[[2]],
                           plots_basic[[4]]$plot_list[[3]],
                           plots_basic[[5]]$plot_list[[1]],
                           plots_basic[[5]]$plot_list[[2]],
                           plots_basic[[5]]$plot_list[[3]],
                           plots_basic[[6]]$plot_list[[1]],
                           plots_basic[[6]]$plot_list[[2]],
                           plots_basic[[6]]$plot_list[[3]])

ggsave(marrangeGrob(plots_basic_unlist,
                    layout_matrix = matrix(1:18, ncol = 6, nrow = 3, byrow = FALSE),
                    top = NULL, widths = c(rep(1,5), 1.4)),
       filename = paste0("plots/PD_plot_uni_basic.jpeg"),
       width = 45, height = 20, units = "cm")

file.copy(from = "plots/PD_plot_uni_basic.jpeg",
          to = paste0(dir_manuscript, "figures"),
          overwrite = TRUE)
