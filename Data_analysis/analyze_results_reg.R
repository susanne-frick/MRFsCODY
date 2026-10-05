####---------- examine MRFs (preliminary) -----------------####

library(ggplot2)
library(colorspace)
library(gridExtra)

mean_mse <- readRDS("Data_analysis/mean_mse_reg.rds")

# check
lapply(mean_mse, dim)
do.call(c, lapply(mean_mse, is.null))

G <- 6
# load fit lists
fit_list_short <- readRDS("Data_analysis/fit_list_short.rds")
fit_list <- readRDS("Data_analysis/fit_list.rds")


# only non-null results
# mean_mse <- mean_mse[! do.call(c, lapply(mean_mse, is.null))]

# select only column with mean (if multi)
mean_mse <- lapply(mean_mse, function(mm) mm[, ncol(mm)])

####------------------------- plots all conditions ---------------------------------####
plot_list <- vector("list", length = nrow(fit_list_short))

for (fs in 1:nrow(fit_list_short)) {
  # fs <- which((fit_list_short$timepoint %in% fit_list_short[d, "timepoint"]) & (fit_list_short$outcome %in% fit_list_short[d, "outcome"]))
  # mm <- mean_mse
  # fl <- fit_list_short
  #
  # fs <- fs[fs <= length(mm)]
  preds <- do.call(c, lapply(fs, function(f, fl, mm) rep(fl[f, "predictors"], length(mm[[f]])),
                             fl = fit_list_short, mm = mean_mse))
  mse_long <- data.frame(MSE = do.call(c, mean_mse[fs]), predictors = preds)

  mse_long$predictors <- factor(mse_long$predictors, levels = c("basic","training","training3"))
  # remove any empty rows
  mse_long <- mse_long[is.na(mse_long$MSE) == FALSE, ]

  title <- paste0("T", fit_list_short[fs, "timepoint"], " ", fit_list_short[fs, "outcome"])

  # violin plot

  plot_list[[fs]] <- ggplot(data=mse_long, aes(y=MSE, x=predictors, fill=predictors)) +
    geom_violin(show.legend=FALSE) +
    labs(y="MSE", x="Predictors") +
    scale_fill_manual(values=qualitative_hcl(3)) +
    #    scale_y_continuous(limits=c(0, 8)) +
    theme(axis.text=element_text(size=11),
          axis.title=element_text(size=11),
          title=element_text(size = 11)) +
    ggtitle(label = title)
}

# save plots
ggsave(marrangeGrob(plot_list, layout_matrix = matrix(1:16, 2, 8, byrow = TRUE), top = NULL),
       file="plots/plots_mse_reg.pdf",
       width=50, height=50, units="cm")
saveRDS(plot_list, file = "Data_analysis/plot_list_reg.rds")
