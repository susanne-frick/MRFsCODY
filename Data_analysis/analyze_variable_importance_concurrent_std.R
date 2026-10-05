#####-------------------------------- analyze variable importance --------------------------------------####

library(ggplot2)
library(gridExtra)
devtools::load_all()
devtools::load_all("~/Dokumente/packages/DataAnalysisSimulation/")

# dir_presentation <- "~/Dokumente/FAIR/Reha/presentations/Kolloquium_Bamberg_June2023/"
dir_manuscript <- "~/Dokumente/FAIR/Reha/paper/mrfs_cody_man/"

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


G <- 6
# load fit lists
fit_list_short <- readRDS("Data_analysis/fit_list_short.rds")
fit_list <- readRDS("Data_analysis/fit_list.rds")


# which ones did converge?
rows_converged <- which(do.call(c, lapply(imp, function(i) nrow(i) > 0)))
fit_list_short[rows_converged, ]
nrow(fit_list_short[rows_converged, ])
# all converged

get_f <- function(outcome, predictors, timepoint, fl = fit_list_short) {
  which(fl$outcome == outcome &
          fl$predictors == predictors &
          fl$timepoint == timepoint)
}

rm_0 <- function(imp, cut = 0) round(imp[rowSums(round(imp[, -c(ncol(imp))], 2) <= cut) < (ncol(imp) - 1),], 2)

prep_imp <- function(tp, imp, cut, outcome = "multi_long", predictors = "training") {
  imp_t <- rm_0(imp[[get_f(outcome, predictors, tp)]], cut = cut)
  rownames(imp_t) <- gsub("_", " ", rownames(imp_t))
  colnames(imp_t) <-  c(paste0("T", rep(1:(ncol(imp_t)/2), 2)), "\\% ")
  imp_t
}

####---------------------- multivariate longitudinal ----------------------------####

# training predictors

lapply(1:6, prep_imp, imp = imp, cut = .1)
# > .1

# Grade (T1-T2)
# Grade x CODY (T1-T4)
# Score_PC1 (T1, T4-T5)

# leveldif_PC1 (T1)
# Score_PC2 (T3)
# leveldif_PC3 (T4)
# leveldif_PC2 (T5, T6)

imp_long_timepoints <- lapply(1:6, prep_imp, imp = imp, cut = 0) #cut = .05
imp_long_timepoints
# > .05

# overall: grade x cody1, score PC1,
# score PC2, leveldif PC1, leveldifPC2, timePC1

# Grade (T1-T4)
# Grade x CODY (T1-T6)
# Score PC1 (T1-T6)
# Score PC2 (T2-T5)
# leveldif PC1 (T1-T2, T4-T5)
# leveldif PC2 (T5-T6)
# leveldif PC3 (T4)
# time PC1 (T1-T3, T6)
# time PC3 (T3)
# time PC2 (T4)

name_variables <- function(vec) {
  vec <- gsub("up", "Up", vec)
  vec <- gsub("leveldif", "Level Diff.", vec)
  vec <- gsub("score", "Level", vec)
  vec <- gsub("time", "Time", vec)
  vec <- gsub("discrep", "Discrepancy", vec)
  vec <- gsub("Grade cody 1", "Pretest x Grade", vec)
  return(vec)
}

imp_long_timepoints <- lapply(imp_long_timepoints, function(imp) {
  rownames(imp) <- name_variables(rownames(imp))
  imp
})

pdf("plots/variable_importances.pdf", width = 13, height = 7)
par(mfrow = c(2,3), mar = c(5,4,4,12), xpd = TRUE)
for (tp in 1:G) {
  sub_imp <- imp_long_timepoints[[tp]][, 1:G]
  sub_imp
  sub_imp <- sub_imp[order(rowMeans(sub_imp), decreasing = TRUE),]

  cols <- rev(RColorBrewer::brewer.pal(min(nrow(sub_imp), 9), name = "YlGnBu")[-c(1:2)]) # too few colors
  if (nrow(sub_imp) > 7) cols <- c(cols, rep(cols[7], nrow(sub_imp) - 7))

  ltypes <- rep(1:6, times = ceiling(nrow(sub_imp)/6))[1:nrow(sub_imp)]

  plot(1:G, sub_imp[1,], ylim = c(0, max(sub_imp)), type = "l",
       xlab = "Status Test", ylab = "Importance Tree",
       col = cols[1], main = paste0("P(T", tp, ")"), xaxt = "none")
  axis(side = 1, at = 1:G, labels = paste0("O(T", 1:G, ")"))
  for(v in 2:nrow(sub_imp)) lines(1:G, sub_imp[v, ], col = cols[v], lty = ltypes[v])
  legend("topright", inset = c(-0.65, 0), legend = rownames(sub_imp), col = cols, lty = ltypes) #box.col = "transparent"
  # or: ggplot grouped bars
}
dev.off()
file.copy(from = "plots/variable_importances.pdf", to = paste0(dir_manuscript, "/figures"), overwrite = TRUE)

####----------------------------- unidimensional --------------------------------------####

lapply(1:6, prep_imp, imp = imp, cut = .05, outcome = "uni", predictors = "basic")
# similar to multilong

# with training predictors up to T3
lapply(4:6, prep_imp, imp = imp, cut = .05, outcome = "uni", predictors = "training3")
# also similar

####---------------------------------------- partial dependence plots -------------------------------------####
# partial dependence plots, only for those variables with high importance values
# get the whole results object

####-------------------------------- longitudinal ----------------------------------------####
# pd_plot delivers mean prediction across folds if a list of fits is given
# fit_f_single is needed to supply the dataset
imp_long_timepoints <- lapply(1:6, prep_imp, imp = imp, cut = 0)
imp_long_timepoints


for (tp in 3:6) {
  fit_f <- readRDS(paste0("results_MRFs_con_std/fit_MRF_con_std_f", get_f("multi_long", "training", tp), ".rds"))
  fit_f_single <- readRDS(paste0("results_MRFs_con_std/results_MRF_con_std_f", get_f("multi_long", "training", tp, fl = fit_list)[1], ".rds"))

  imp_tp <- imp_long_timepoints[[tp]]
  imp_tp <- imp_tp[order(rowSums(imp_tp[, 1:G]), decreasing = TRUE)[1:3],]
  xp <- gsub(" ", "_", rownames(imp_tp))
  xp <- xp[order(xp)]

  plot_list <- pd_plot(fit = fit_f,
                       data_plot = rbind(fit_f_single$data$train, fit_f_single$data$test),
                       x_plot = xp,
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
    plot_list$plot_list[[p]] <- plot_list$plot_list[[p]] + ggtitle(paste0("O(T", rep(1:6, each = 3)[p], ")")) +
      xlab(x_y_names[p, "Predictor1"]) + ylab(x_y_names[p, "Predictor2"])
    if(! (p %in% 16:18)) plot_list$plot_list[[p]] <- plot_list$plot_list[[p]] + theme(legend.position = "none")
  }

  ggsave(marrangeGrob(plot_list$plot_list, layout_matrix = matrix(1:18, ncol = 6, nrow = 3, byrow = FALSE),
                      top = NULL, widths = c(rep(1,5), 1.6)),
         filename = paste0("plots/PD_plot_multi-long_training_T", tp, ".jpeg"),
         width = 45, height = 20, units = "cm")
}

file.copy(from = paste0("plots/", list.files("plots", pattern = "PD_plot_multi-long_training_T")),
          to = paste0(dir_manuscript, "figures"),
          overwrite = TRUE)

# how are the predictions distributed?
preds <- predict(fit_f$fit, newdata = fit_f$data$test)
summary(preds$predicted)
hist(preds$predicted)
plot(preds$predicted, fit_f$data$test[, fit_f$fit$yvar.names], xlab = "Predicted", ylab = "True")
summary(fit_f$mse_test[, "mean"])



####----------------------------- example data ---------------------------####
# to explain permutation importance in presentations

load("Aufbereitung/Daten_cody_gold_sub_wide.RData")
sub <- head(Daten_wide)[, 105:109]
colnames(sub) <- paste0("$x_", 1:ncol(sub), "$")
sub
set.seed(132)
sub_permuted <- data.frame(sub[, 1:2], sub[, 3][sample(1:nrow(sub))], sub[, 4:5])
colnames(sub_permuted) <- colnames(sub)

sub_permuted[, 3] <- paste("\\textcolor{red}{", sub_permuted[, 3], "}")
colnames(sub_permuted)[3] <- paste("\\textcolor{red}{", colnames(sub_permuted)[3], "}")

print(xtable::xtable(sub, digits=0), include.colnames = T, include.rownames=F, hline.after=c(-1, 0),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){gsub("_", " ", x)},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp",
      caption.placement = "top", latex.environments = NULL,
      file="presentations/Update_Meeting_March2023/tables/textable_example_data.tex")

print(xtable::xtable(sub_permuted, digits=0, align = "rrrrrr"), include.colnames = T, include.rownames=F, hline.after=c(-1, 0),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){gsub("_", " ", x)},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp",
      caption.placement = "top", latex.environments = NULL,
      file="presentations/Update_Meeting_March2023/tables/textable_example_data_permuted.tex")

