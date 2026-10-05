####---------- examine MRFs (preliminary) -----------------####

library(ggplot2)
library(colorspace)
library(gridExtra)

mean_mse <- readRDS("Data_analysis/mean_mse_MRFs.rds")

# check
lapply(mean_mse, dim)
do.call(c, lapply(mean_mse, is.null))
table(do.call(c, lapply(mean_mse, is.null)))

# re-build design list
G <- 6
# without predictors (on x-axis in plot)
fit_list_short <- rbind(expand.grid("timepoint" = 1:G,
                                    "outcome" = c("multi", "uni", "avg", "row", "sub", "add")),
                        expand.grid("timepoint" = 1:G,
                                    "outcome" = "multi_long") # multi_long without status predictors
)

fit_list_sum <- rbind(expand.grid("timepoint" = 1:G,
                                    "outcome" = c("multi", "uni", "row", "sub", "add"),
                                    "predictors" = c("basic", "training", "status")),
                        expand.grid("timepoint" = 1:G,
                                    "outcome" = "multi_long",
                                    "predictors" = c("basic", "training")) # multi_long without status predictors
)

# descriptives on n (drop-out)
ns <- data.frame(fit_list_sum, "n" = do.call(c, lapply(mean_mse, function(ms) ifelse(is.null(ms), 0, nrow(ms)))))
ns
ns[ns$outcome %in% c("uni", "multi_long"),]

mean_mse_long <- do.call(rbind, lapply(1:nrow(fit_list_sum), function(f, fl, mm) {
  if (is.null(mm[[f]])) {
    m <- cbind(0,0)
  }  else {
    m <- mm[[f]][, c(1, ncol(mm[[f]]))]
  }
  df <- data.frame(fl[f,], m)
  colnames(df) <- c(colnames(fl[f,]), c("ids", "mse"))
  return(df)
} , mm = mean_mse, fl = fit_list_sum))

mean_mse_long$log_mse <- log(mean_mse_long$mse)
mean_mse_long$log_mse[mean_mse_long$mse == 0] <- 0
nlme_mse <- nlme::lme(log_mse ~ timepoint * outcome * predictors, random = ~1 | ids,
                      data = mean_mse_long[mean_mse_long$outcome %in% c("multi", "uni"),])
summary(nlme_mse)

# multi_long vs. uni
# subset without status predictors
mean_mse_long_sub <- mean_mse_long[mean_mse_long$predictors != "status", ]
nlme_mse <- nlme::lme(log_mse ~ timepoint * outcome * predictors, random = ~1 | ids,
                      data = mean_mse_long_sub[mean_mse_long_sub$outcome %in% c("multi_long", "uni"),])
summary(nlme_mse)
# according to p-value: multilong has a smaller MSE than uni

# select only column with mean (if multi)
mean_mse <- lapply(mean_mse, function(mm) mm[, ncol(mm)])

####---------------------------- average across add, sub, row ------------------------####

fit_list_avg <- expand.grid("timepoint" = 1:G,
                            "predictors" = c("basic", "training", "status"))
mean_mse_avg <- vector("list", length = nrow(fit_list_avg))

for (d in 1:nrow(fit_list_avg)) {
  if((fit_list_avg[d, "timepoint"] == 1) & (fit_list_avg[d, "predictors"] == "status")) {
    next
  } else {
  fs <- which((fit_list_sum$timepoint == fit_list_avg[d, "timepoint"]) & (fit_list_sum$predictors == fit_list_avg[d, "predictors"]) &
                (fit_list_sum$outcome %in% c("row", "sub", "add")))
  fs <- fs[fs <= length(mean_mse)]
  # average only if the same amount of observations
  if(length(unique(do.call(c, lapply(mean_mse[fs], length)))) == 1) {
    mean_mse_avg[[d]] <- rowMeans(do.call(cbind, mean_mse[fs]))
  }
  }
}

####------------------------- plots all conditions ---------------------------------####
plot_list <- vector("list", length = nrow(fit_list_short))


for (d in 1:nrow(fit_list_short)) {

  if((fit_list_sum[d, "timepoint"] == 1) & (fit_list_sum[d, "predictors"] == "status")) {
    next
  } else {

  if(fit_list_short[d, "outcome"] != "avg") {
    fs <- which((fit_list_sum$timepoint %in% fit_list_short[d, "timepoint"]) & (fit_list_sum$outcome %in% fit_list_short[d, "outcome"]))
    mm <- mean_mse
    fl <- fit_list_sum
  } else {
    fs <- which((fit_list_avg$timepoint %in% fit_list_short[d, "timepoint"]))
    mm <- mean_mse_avg
    fl <- fit_list_avg
  }

  if(all(do.call(c, lapply(mm[fs], is.null)))) next

  preds <- do.call(c, lapply(fs, function(f, fl, mm) rep(fl[f, "predictors"], length(mm[[f]])),
                             fl = fl, mm = mm))
  mse_long <- data.frame(MSE = do.call(c, mm[fs]), predictors = preds)

  mse_long$predictors <- factor(mse_long$predictors, levels = c("basic","training", "status"))
  # remove any empty rows
  mse_long <- mse_long[is.na(mse_long$MSE) == FALSE, ]

  title <- paste0("T", fit_list_short[d, "timepoint"], " ", fit_list_short[d, "outcome"])

  # violin plot

  plot_list[[d]] <- ggplot(data=mse_long, aes(y=MSE, x=predictors, fill=predictors)) +
    geom_violin(show.legend=FALSE) +
    labs(y="MSE", x="Predictors") +
    scale_fill_manual(values=qualitative_hcl(3)) +
#    scale_y_continuous(limits=c(0, 8)) +
    theme(axis.text=element_text(size=11),
          axis.title=element_text(size=11),
          title=element_text(size = 11)) +
    ggtitle(label = title)
  }
}

# save plots
ggsave(marrangeGrob(plot_list, layout_matrix = matrix(1:42, 7, 6, byrow = TRUE), top = NULL),
       file="plots/plots_mse.pdf",
       width=50, height=50, units="cm")
saveRDS(plot_list, file = "Data_analysis/plot_list_MRFs.rds")

# only multi_long and uni
plot_list <- readRDS("Data_analysis/plot_list_MRFs.rds")
ds <- which(fit_list_short$outcome %in% c("multi_long", "uni"))
ggsave(marrangeGrob(lapply(plot_list[ds], function(pl) pl + scale_y_continuous(limits = c(0, 10))),
                    layout_matrix = matrix(1:12, 2, 6, byrow = TRUE),
                    top = NULL),
       file="plots/plots_mse_longitudinal-uni.pdf",
       width=50, height=17, units="cm")


####--------------------  plots for presentation FAIR Update Seminar Mary 2023 ------------------------------####

# multi, uni, avg x T2, T6
# f_outcomes <- which((fit_list_short$timepoint %in% c(2, 6)) & (fit_list_short$outcome %in% c("multi", "uni", "avg")))
# f_outcomes <- f_outcomes[c(1,3,5,2,4,6)]
# fit_list_short[f_outcomes, ]
# plots_outcomes <- plot_list[f_outcomes]
# plots_outcomes[1:3] <- lapply(plots_outcomes[1:3], function(pl) pl + scale_y_continuous(limits=c(0, 8)))
# plots_outcomes[4:6] <- lapply(plots_outcomes[4:6], function(pl) pl + scale_y_continuous(limits=c(0, 8)))
# plots_outcomes <- lapply(1:length(plots_outcomes), function(p, pp, out) pp[[p]] + ggtitle(out[p]),
#                          pp = plots_outcomes, out = c("Multivariate", "Univariate (sum)", "Average over univariate",
#                                                       "", "", ""))
# ggsave(marrangeGrob(marrangeGrob(plots_outcomes, nrow = 1, ncol = 3, left = "T2", top = NULL),
#                     marrangeGrob(plots_outcomes, nrow = 1, ncol = 3, left = "T6", top = NULL),
#                     nrow = 2, ncol = 1, top = NULL),
#        file="~/Dokumente/FAIR/Reha/presentations/Update_Seminar_May2023/figures/plots_mse_outcomes.pdf",
#        width=25, height=15, units="cm")

####-------------------- plots for Kolloquium Bamberg June 2023 ------------------------------------------- ###
# similar to Advisory Board Meeting January 2023

# timpoints: uni # multi
f_timepoints <- which(fit_list_short$outcome == "uni")
plots_timepoints <- plot_list[f_timepoints]
plots_timepoints <- lapply(plots_timepoints, function(pl) pl + scale_y_continuous(limits=c(0, 8)))
plots_timepoints <- lapply(1:length(plots_timepoints), function(p, pp, tp) pp[[p]] + ggtitle(paste0("Timepoint ", tp[p])),
                           pp = plots_timepoints, tp = 1:6)
ggsave(marrangeGrob(plots_timepoints, layout_matrix = matrix(1:6, 2, 3, byrow = TRUE), top = NULL),
       file="plots/plots_mse_timepoints.pdf",
       width=25, height=15, units="cm")
file.copy(from = "plots/plots_mse_timepoints.pdf",
          to = "~/Dokumente/FAIR/Reha/presentations/Kolloquium_Bamberg_June2023/figures/plots_mse_timepoints.pdf",
          overwrite = TRUE)

# predictors, multi, T2 and T4
f_predictors <- which((fit_list_short$timepoint %in% c(2,4)) & (fit_list_short$outcome == "multi"))
plots_predictors <- plot_list[f_predictors]
plots_predictors <- lapply(plots_predictors, function(pl) pl + scale_y_continuous(limits=c(0, 6)))
plots_predictors <- lapply(1:length(plots_predictors), function(p, pp, tp) pp[[p]] + ggtitle(paste0("Timepoint ", tp[p])),
                           pp = plots_predictors, tp = c(2,4))
ggsave(marrangeGrob(plots_predictors, nrow = 1, ncol = 2, top = NULL),
       file="plots/plots_mse_predictors.pdf",
       width=20, height=10, units="cm")
file.copy(from = "plots/plots_mse_predictors.pdf",
          to = "~/Dokumente/FAIR/Reha/presentations/Kolloquium_Bamberg_June2023/figures/plots_mse_predictors.pdf",
          overwrite = TRUE)


# outcomes: T3 (without row, sub, add)
f_outcomes <- which(fit_list_short$timepoint == 3 & fit_list_short$outcome %in% c("multi", "uni", "avg"))
plots_outcomes <- plot_list[f_outcomes]
plots_outcomes <- lapply(plots_outcomes, function(pl) pl + scale_y_continuous(limits=c(0, 6)))
plots_outcomes <- lapply(1:length(plots_outcomes), function(p, pp, out) pp[[p]] + ggtitle(out[p]),
                           pp = plots_outcomes, out = c("Multivariate", "Univariate (sum)", "Average over univariate"))
ggsave(marrangeGrob(plots_outcomes, layout_matrix = matrix(1:3, 1, 3, byrow = TRUE), top = NULL),
       file="plots/plots_mse_outcomes.pdf",
       width=25, height=8, units="cm")
file.copy(from = "plots/plots_mse_outcomes.pdf",
          to = "~/Dokumente/FAIR/Reha/presentations/Kolloquium_Bamberg_June2023/figures/plots_mse_outcomes.pdf",
          overwrite = TRUE)

# MRFs versus regession, T2 multi
plot_list_reg <- readRDS("Data_analysis/plot_list_reg.rds")
f_models <- which((fit_list_short$timepoint %in% c(2)) & (fit_list_short$outcome == "multi"))
plots_models <- list(plot_list[f_models][[1]], plot_list_reg[f_models][[1]])
plots_models <- lapply(plots_models, function(pl) pl + scale_y_continuous(limits=c(0, 6)))
plots_models <- lapply(1:length(plots_models), function(p, pp, mod) pp[[p]] + ggtitle(mod[p]),
                         pp = plots_models, mod = c("Forest", "Lasso Regression"))

ggsave(marrangeGrob(plots_models, nrow = 1, ncol = 2, top = NULL),
       file="plots/plots_mse_models.pdf",
       width=20, height=10, units="cm")
file.copy(from = "plots/plots_mse_models.pdf",
          to = "~/Dokumente/FAIR/Reha/presentations/Kolloquium_Bamberg_June2023/figures/plots_mse_models.pdf",
          overwrite = TRUE)
