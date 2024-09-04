#####-------------------------------- analyze variable importance --------------------------------------####

library(ggplot2)
library(gridExtra)
devtools::load_all()
devtools::load_all("~/Dokumente/packages/DataAnalysisSimulation/")

# dir_presentation <- "~/Dokumente/FAIR/Reha/presentations/Kolloquium_Bamberg_June2023/"
dir_manuscript <- "~/Dokumente/FAIR/Reha/paper/MRFs CODY/"

imp <- readRDS("Data_analysis/mean_importance.rds")

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
  imp_t <- rm_0(imp[[get_f("multi_long", predictors, tp)]], cut = cut)
  rownames(imp_t) <- gsub("_", " ", rownames(imp_t))
  colnames(imp_t) <-  c(paste0("T", rep(1:6, 2)), "\\% ")
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

imp_long_timepoints <- lapply(1:6, prep_imp, imp = imp, cut = .05)
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


# combine into one table
tab_imp_long_timepoints <- do.call(rbind, lapply(1:length(imp_long_timepoints), function(ti, imp) {
  tb <- cbind(rownames(imp[[ti]]), imp[[ti]])
  rownames(tb) <- NULL
  tb <- rbind(c(paste0("up to T", ti), rep("", ncol(tb) - 1)), tb)
  tb
}, imp = imp_long_timepoints)
)
colnames(tab_imp_long_timepoints)[1] <- "Predictor"
tab_imp_long_timepoints

tab_imp_long_timepoints$Predictor <- gsub("up", "Up", tab_imp_long_timepoints$Predictor)
tab_imp_long_timepoints$Predictor <- gsub("leveldif", "Level Diff.", tab_imp_long_timepoints$Predictor)
tab_imp_long_timepoints$Predictor <- gsub("score", "Level", tab_imp_long_timepoints$Predictor)
tab_imp_long_timepoints$Predictor <- gsub("time", "Time", tab_imp_long_timepoints$Predictor)
tab_imp_long_timepoints$Predictor <- gsub("Grade cody 1", "Pretest x Grade", tab_imp_long_timepoints$Predictor)
tab_imp_long_timepoints

header <- list()
header$pos <- list(-1)
header$command <- c("\\hline \n & \\multicolumn{6}{c}{Importance Trees} & \\multicolumn{6}{c}{Importance Forest} &  Splits\\\\")

print(xtable::xtable(tab_imp_long_timepoints,
                     digits = c(0, 0, rep(2, 12), 0),
                     caption = "Variable Importance Measures for the Training Predictors in the Multivariate Random Forests. Only predictors with importances > .05 are shown.",
                     label = "tb:importance_timepoints"),
      include.colnames = T, include.rownames=F,
      hline.after=c(-1, 0, nrow(tab_imp_long_timepoints), grep("Up to", tab_imp_long_timepoints$Predictor) - 1),
      sanitize.rownames.function=function(x){x},
      sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file = paste0(dir_manuscript, "tables/textable_importance_MRF_long-timepoints-05.tex"))


# basic predictors
imp_long_basic_T1 <- prep_imp(1, imp = imp, cut = .05, predictors = "basic")
imp_long_basic_T1
rownames(imp_long_basic_T1) <- c("Pre-test", "Grade", "Pretest $\\times$ Grade")

header <- list()
header$pos <- list(-1)
header$command <- c("\\hline \n Predictor & \\multicolumn{6}{c}{Importance Trees} & \\multicolumn{6}{c}{Importance Forest} &  Splits\\\\")

print(xtable::xtable(imp_long_basic_T1,
                     digits = c(0, rep(2, 12), 0),
                     caption = "Variable Importance Measures for the Basic Predictors in the Multivariate Random Forests.",
                     label = "tb:importance_basic"),
      include.colnames = T, include.rownames=T,
      hline.after=c(-1, 0, nrow(imp_long_basic_T1)),
      sanitize.rownames.function=function(x){x},
      sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file = paste0(dir_manuscript, "tables/textable_importance_MRF_long-basic-T1.tex"))

#-#-#- T4
imp_long_training_T4 <- prep_imp(4, imp = imp, cut = .05)
rownames(imp_long_training_T4) <- gsub("Grade cody 1", "Pre-test x Grade", rownames(imp_long_training_T4))
rownames(imp_long_training_T4) <- gsub("score", "Level", rownames(imp_long_training_T4))
rownames(imp_long_training_T4) <- gsub("leveldif", "Level difference", rownames(imp_long_training_T4))
rownames(imp_long_training_T4) <- gsub("time", "Response time", rownames(imp_long_training_T4))
imp_long_training_T4

header <- list()
header$pos <- list(-1)
header$command <- c("\\hline \n Predictor & \\multicolumn{6}{c}{Importance Trees} & \\multicolumn{6}{c}{Importance Forest} &  Splits\\\\")

print(xtable::xtable(imp_long_training_T4,
                     digits = c(0, rep(2, 12), 0),
                     caption = "Variable Importance Measures for the Training Predictors up to T4. Only predictors with importances > .05 are shown.",
                     label = "tb:importance_T4"),
      include.colnames = T, include.rownames=T,
      hline.after=c(-1, 0, nrow(imp_long_training_T4)),
      sanitize.rownames.function=function(x){x},
      sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file = paste0(dir_manuscript, "tables/textable_importance_MRF_long-training-T4.tex"))

####---------------------------- basic and long T4 - shortened for presentations ------------------------------####

# shortened for presentation (Tübingen)
imp_long_training_T4_short <- imp_long_training_T4[, c(1:6, 13)]
names_rw <- rownames(imp_long_training_T4_short)
rw_large <- apply(imp_long_training_T4_short[, 1:6], 1, function(rw) any(rw > .1))
imp_long_training_T4_short <- apply(imp_long_training_T4_short, 2, as.character)
imp_long_training_T4_short[rw_large, ] <- paste0("\\textcolor{green}{", imp_long_training_T4_short[rw_large, ], "}")
rownames(imp_long_training_T4_short) <- ifelse(rw_large, paste0("\\textcolor{green}{", names_rw, "}"), names_rw)
imp_long_training_T4_short

header <- list()
header$pos <- list(-1)
header$command <- c("\\hline \n Predictor & \\multicolumn{6}{c}{Importance Trees} & Splits\\\\")

print(xtable::xtable(imp_long_training_T4_short,
                     digits = c(0, rep(2, 6), 0),
                     caption = NULL,
                     label = "tb:importance_T4"),
      include.colnames = T, include.rownames=T,
      hline.after=c(-1, 0, nrow(imp_long_training_T4)),
      sanitize.rownames.function=function(x){x},
      sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file = paste0(dir_manuscript, "tables/textable_importance_MRF_long-training-T4_shortened.tex"))

# basic
imp_long_basic_T1
imp_long_basic_T1_short <- imp_long_basic_T1[, c(1:6, 13)]
names_rw <- rownames(imp_long_basic_T1_short)
rw_large <- apply(imp_long_basic_T1_short[, 1:6], 1, function(rw) any(rw > .1))
imp_long_basic_T1_short <- apply(imp_long_basic_T1_short, 2, as.character)
imp_long_basic_T1_short[rw_large, ] <- paste0("\\textcolor{green}{", imp_long_basic_T1_short[rw_large, ], "}")
rownames(imp_long_basic_T1_short) <- ifelse(rw_large, paste0("\\textcolor{green}{", names_rw, "}"), names_rw)
imp_long_basic_T1_short

header <- list()
header$pos <- list(-1)
header$command <- c("\\hline \n Predictor & \\multicolumn{6}{c}{Importance Trees} & Splits\\\\")

print(xtable::xtable(imp_long_basic_T1_short,
                     digits = c(0, rep(2, 6), 0),
                     caption = NULL,
                     label = "tb:importance_basic"),
      include.colnames = T, include.rownames=T,
      hline.after=c(-1, 0, nrow(imp_long_basic_T1)),
      sanitize.rownames.function=function(x){x},
      sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file = paste0(dir_manuscript, "tables/textable_importance_MRF_long-basic-T1_shortened.tex"))


####----------------------------- unidimensional --------------------------------------####

lapply(1:6, prep_imp, imp = imp, cut = .05, outcome = "uni")
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

for (tp in c(2,3,4,6)) {
  fit_f <- readRDS(paste0("results_MRFs/fit_MRF_f", get_f("multi_long", "training", tp), ".rds"))
  fit_f_single <- readRDS(paste0("results_MRFs/results_MRF_f", get_f("multi_long", "training", tp, fl = fit_list)[1], ".rds"))

  imp <- imp_long_timepoints[[tp]]
  imp <- imp[order(rowSums(imp[, 1:G]), decreasing = TRUE)[1:3],]
  xp <- gsub(" ", "_", rownames(imp))
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
    plot_list$plot_list[[p]] <- plot_list$plot_list[[p]] + ggtitle(paste0("T", rep(1:6, each = 3)[p])) +
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

