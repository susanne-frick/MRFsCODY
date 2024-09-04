####-------------------- analyse results MRFs longitudinal outcome ------------------####

library(ggplot2)
library(colorspace)
library(gridExtra)
library(colorspace)

mean_mse <- readRDS("Data_analysis/mean_mse_MRFs.rds")

dir_manuscript <- "~/Dokumente/FAIR/Reha/paper/MRFs CODY/"

# re-build design list
G <- 6
fit_list_short <- readRDS("Data_analysis/fit_list_short.rds")

print_tab <- function(tab, name_rownames = NULL, file, addr = NULL, ...) {
  # change rownames to first column if rownames title is given
  if(isFALSE(is.null(name_rownames))) {
    cn <- colnames(tab)
    tab <- data.frame(rownames(tab), tab)
    colnames(tab) <- c(name_rownames, cn)
    hl <- c(0, nrow(tab))
  } else {
    hl <- c(-1, 0, nrow(tab))
  }

  print(xtable::xtable(tab, digits=3, ...),
        include.colnames = T, include.rownames = is.null(name_rownames),
        hline.after = hl,
        add.to.row = addr,
        sanitize.rownames.function=function(x){x},
        sanitize.colnames.function = function(x){x},
        sanitize.text.function = function(x){x},
        NA.string = "", table.placement = "htp",
        caption.placement = "top", latex.environments = NULL,
        file=file)
}


####------------------- basic predictors -----------------------####
f <- which(fit_list_short$outcome == "multi_long" & fit_list_short$predictors == "basic")
head(mean_mse[[f]])
mse_long_basic <- data.frame("id" = rep(mean_mse[[f]]$ids, G),
                             "timepoint" = rep(1:G, each = nrow(mean_mse[[f]])),
                             "mse" = do.call(c, mean_mse[[f]][ , 2:(G + 1)]))
mse_long_basic$timepoint <- factor(mse_long_basic$timepoint)
head(mse_long_basic)


plot_basic_long <- ggplot(data = mse_long_basic, aes(y = mse, x=timepoint)) +
  geom_violin(show.legend=FALSE, fill = qualitative_hcl(3)[1]) +
  labs(y="MSE", x="Timepoint") +
  # scale_fill_manual(values=qualitative_hcl(3)[1])
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11),
        title=element_text(size = 11)) +
  ggtitle(label = "Multivariate, Pre-test & Grade")


####---------------------- training predictors ----------------------------####

fs <- which(fit_list_short$outcome == "multi_long" & fit_list_short$predictors == "training")
fit_list_short[fs, ]

plot_list <- vector("list", length = G)

for (f in fs) {
  mse_long <- data.frame("id" = rep(mean_mse[[f]]$ids, G),
                               "timepoint" = rep(1:G, each = nrow(mean_mse[[f]])),
                               "mse" = do.call(c, mean_mse[[f]][ , 2:(G + 1)]))
  mse_long$timepoint <- factor(mse_long$timepoint)
  head(mse_long)


  plot_list[[f - fs[1] + 1]] <- ggplot(data = mse_long, aes(y = mse, x=timepoint)) +
    geom_violin(show.legend=FALSE, fill = qualitative_hcl(3)[2]) +
    labs(y="MSE", x="Timepoint") +
    theme(axis.text=element_text(size=11),
          axis.title=element_text(size=11),
          title=element_text(size = 11)) +
    scale_y_continuous(limits = c(0, 6.5)) +
    ggtitle(label = paste0("Training T", f - fs[1] + 1))

}

ggsave(marrangeGrob(plot_list, layout_matrix = matrix(1:G, 2, 3, byrow = TRUE), top = NULL),
       file="plots/plots_mse_longitudinal.pdf",
       width=30, height=15, units="cm")
file.copy(from = "plots/plots_mse_longitudinal.pdf", to = paste0(dir_manuscript, "/figures"), overwrite = TRUE)

ggsave(grid.arrange(plot_basic_long + scale_y_continuous(limits = c(0,8)),
                    plot_list[[2]]  + scale_y_continuous(limits = c(0,8)),
                    plot_list[[6]] + scale_y_continuous(limits = c(0,8)),
                    nrow = 1, ncol = 3),
       file = "plots/plot_longitudinal_basic-training-T6.pdf",
       width = 30, height = 9, units = "cm")
file.copy(from = "plots/plot_longitudinal_basic-training-T6.pdf", to = paste0(dir_manuscript, "/figures"), overwrite = TRUE)

# Kolmogorov-Smirnov test for T6, data until T6 vs. T5 to T1
lapply(mean_mse[fs], nrow)


# Kolmogorov-Smirnov test and Wasserstein distance for predictors up to T1 to T5 vs. T6
# Wasserstein with p = 2 (squared distance)
# for each outcome (T1-T6) separately
table_w_timepoints <- table_p_timepoints <- matrix(NA, 5, 6)
rownames(table_p_timepoints) <- rownames(table_w_timepoints) <- paste0("T", 1:(G-1))
colnames(table_p_timepoints) <- colnames(table_w_timepoints) <- paste0("T", 1:G)

for(o in 1:6) {
  for(d in 1:5) {
    ks <- ks.test(mean_mse[[fs[6]]][, o + 1], mean_mse[[fs[d]]][, o + 1])
    table_w_timepoints[d, o] <- transport::wasserstein1d(mean_mse[[fs[6]]][, o + 1], mean_mse[[fs[d]]][, o + 1], p = 2)
    table_p_timepoints[d, o] <- ks$p.value
  }
}

table_timepoints <- cbind(round(table_p_timepoints, 3),
                          round(table_w_timepoints, 3)
)

# exemplarily for up to T1 versus up to T6
desc_timepoints <- cbind(apply(mean_mse[[fs[1]]][, -c(1,8)], 2, summary),
                         apply(mean_mse[[fs[6]]][, -c(1,8)], 2, summary)
)
desc_timepoints <- round(desc_timepoints, 3)
colnames(desc_timepoints) <- rep(paste0("T", 1:G), 2)
rownames(desc_timepoints) <- c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum")
desc_timepoints


print_tab(table_timepoints,
          name_rownames = "Up to",
          file = paste0(dir_manuscript, "/tables/table_timepoints.tex"),
          addr = list("pos" = list(-1),
          "command" = c("\\hline \n & \\multicolumn{6}{c}{Kolmogorov-Smirnov p-value} & \\multicolumn{6}{c}{Wasserstein Distance} \\\\ \\cmidrule{2-7} \\cmidrule{8-13}")
          ),
          caption = "Tests for Differences in the Distributions of Mean Squared Error for the Training Predictors up to T6 versus T1 to T5 for the Multivariate Random Forests.",
          label = "tb:test_mse_timepoints")

print_tab(desc_timepoints,
          name_rownames = "Statistic",
          file = paste0(dir_manuscript, "/tables/descriptives_timepoints.tex"),
          addr = list("pos" = list(-1),
                      "command" = c("\\hline \n & \\multicolumn{6}{c}{Training Predictors up to T1} & \\multicolumn{6}{c}{Training Predictors up to T6} \\\\ \\cmidrule{2-7} \\cmidrule{8-13}")
          ),
          caption = "Descriptives on the Distributions of Mean Squared Error for the Training Predictors up to T1 versus T6 for the Multivariate Random Forests.",
          label = "tb:desc_mse_timepoints")

####-------------- basic versus training predictors ----------------------------------####

f_basic <- which(fit_list_short$outcome == "multi_long" & fit_list_short$predictors == "basic")
fs_training <- which(fit_list_short$outcome == "multi_long" & fit_list_short$predictors == "training")

# Kolmogorov-Smirnov test and Wasserstein distance
table_w_basic <- table_p_basic <- matrix(NA, 6, 6)
rownames(table_p_basic) <- rownames(table_w_basic) <- paste0("T", 1:G)
colnames(table_p_basic) <- colnames(table_w_basic) <- paste0("T", 1:G)

for(o in 1:6) {
  for(d in 1:6) {
    ks <- ks.test(mean_mse[[f_basic]][, o + 1], mean_mse[[fs_training[d]]][, o + 1])
    table_w_basic[d, o] <- transport::wasserstein1d(mean_mse[[f_basic]][, o + 1], mean_mse[[fs_training[d]]][, o + 1], p = 2)
    table_p_basic[d, o] <- ks$p.value
  }
}
# ties in the training data

table_basic <- cbind(round(table_p_basic, 3),
                          round(table_w_basic, 3)
)

# exemplarily for up to T1 versus up to T6
desc_basic <- apply(mean_mse[[f_basic]][, -c(1,8)], 2, summary)
desc_basic <- round(desc_basic, 3)
colnames(desc_basic) <- paste0("T", 1:G)
rownames(desc_basic) <- c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum")
desc_basic

print_tab(table_basic,
          name_rownames = "Up to",
          file = paste0(dir_manuscript, "/tables/table_basic.tex"),
          addr = list("pos" = list(-1),
                      "command" = c("\\hline \n & \\multicolumn{6}{c}{Kolmogorov-Smirnov p-value} & \\multicolumn{6}{c}{Wasserstein Distance} \\\\ \\cmidrule{2-7} \\cmidrule{8-13}")
          ),
          caption = "Tests for Differences in the Distributions of Mean Squared Error for the Basic versus the Training Predictors for the Multivariate Random Forests.",
          label = "tb:test_mse_basic")

print_tab(desc_basic,
          name_rownames = "Statistic",
          file = paste0(dir_manuscript, "/tables/descriptives_basic.tex"),
          addr = list("pos" = list(-1),
                      "command" = c("\\hline \n")
          ),
          caption = "Descriptives on the Distributions of Mean Squared Error for the Basic Predictors for the Multivariate Random Forests.",
          label = "tb:desc_mse_basic")

####-------------- comparison with univariate sum score, separately for the 6 timepoints --------####

# basic predictors
fs <- which(fit_list_short$outcome == "uni" & fit_list_short$predictors == "basic")
lapply(mean_mse[fs], nrow)

mse_uni_basic <- data.frame("id" = do.call(c, lapply(mean_mse[fs], function(m) m$ids)),
                             "timepoint" = do.call(c, lapply(1:G, function(g, m) rep(g, nrow(m[[g]])), m = mean_mse[fs])),
                             "mse" = do.call(c, lapply(mean_mse[fs], function(m) m[,2])))
mse_uni_basic$timepoint <- factor(mse_uni_basic$timepoint)
head(mse_uni_basic)


plot_basic_uni <- ggplot(data = mse_uni_basic, aes(y = mse, x=timepoint)) +
  geom_violin(show.legend=FALSE, fill = qualitative_hcl(3)[1]) +
  labs(y="MSE", x="Timepoint") +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11),
        title=element_text(size = 11)) +
  ggtitle(label = "Univariate, Pre-test & Grade")

# training predictors
# ! up to T3 for T4-T6
fs <- which(fit_list_short$outcome == "uni" & fit_list_short$predictors == "training3")
lapply(mean_mse[fs], nrow)

mse_uni_training <- data.frame("id" = do.call(c, lapply(mean_mse[fs], function(m) m$ids)),
                            "timepoint" = do.call(c, lapply(1:3, function(g, m) rep(g, nrow(m[[g]])), m = mean_mse[fs])),
                            "mse" = do.call(c, lapply(mean_mse[fs], function(m) m[,2])))
mse_uni_training$timepoint <- mse_uni_training$timepoint + 3
mse_uni_training$timepoint <- factor(mse_uni_training$timepoint)
head(mse_uni_training)


plot_training_uni <- ggplot(data = mse_uni_training, aes(y = mse, x=timepoint)) +
  geom_violin(show.legend=FALSE, fill = qualitative_hcl(3)[2]) +
  labs(y="MSE", x="Timepoint") +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11),
        title=element_text(size = 11)) +
  ggtitle(label = "Univariate, Training up to T3")

f <- which(fit_list_short$outcome == "multi_long" & fit_list_short$predictors == "training" & fit_list_short$timepoint == 3)
fit_list_short[f, ]

mse_long_T3 <- data.frame("id" = rep(mean_mse[[f]]$ids, G),
                          "timepoint" = rep(1:G, each = nrow(mean_mse[[f]])),
                          "mse" = do.call(c, mean_mse[[f]][ , 2:(G + 1)]))
mse_long_T3 <- mse_long_T3[mse_long_T3$timepoint > 3,]
mse_long_T3$timepoint <- factor(mse_long_T3$timepoint)
head(mse_long_T3)

plot_long_T3 <- ggplot(data = mse_long_T3, aes(y = mse, x=timepoint)) +
  geom_violin(show.legend=FALSE, fill = qualitative_hcl(3)[2]) +
  labs(y="MSE", x="Timepoint") +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11),
        title=element_text(size = 11)) +
  ggtitle(label = "Multivariate, Training up to T3")

ggsave(grid.arrange(plot_basic_long + scale_y_continuous(limits = c(0,8)),
                    plot_basic_uni  + scale_y_continuous(limits = c(0,8)),
                    plot_long_T3 + scale_y_continuous(limits = c(0, 6.5)),
                    plot_training_uni + scale_y_continuous(limits = c(0, 6.5)),
                    nrow = 2, ncol = 2),
       file = "plots/plot_long-uni_basic-training.pdf",
       width = 20, height = 15, units = "cm")
file.copy(from = "plots/plot_long-uni_basic-training.pdf", to = paste0(dir_manuscript, "/figures"), overwrite = TRUE)

# Kolmogorov-Smirnov test and Wasserstein distance for predictors up to T1 to T5 vs. T6
# Wasserstein with p = 2 (squared distance)
# for each outcome (T1-T6) separately
table_w_uni <- table_p_uni <- matrix(NA, 2, 6)
rownames(table_p_uni) <- rownames(table_w_uni) <- c("Basic", "Training T3")
colnames(table_p_uni) <- colnames(table_w_uni) <- paste0("T", 1:G)

for(o in 1:6) {
  ks <- ks.test(mse_uni_basic[mse_uni_basic$timepoint == o, "mse"], mse_long_basic[mse_long_basic$timepoint == o, "mse"])
  # there is a few ties
  table_p_uni[1, o] <- ks$p.value
  table_w_uni[1, o] <- transport::wasserstein1d(mse_uni_basic[mse_uni_basic$timepoint == o, "mse"], mse_long_basic[mse_long_basic$timepoint == o, "mse"], p = 2)
}

for(o in 4:6) {
  ks <- ks.test(mse_uni_training[mse_uni_training$timepoint == o, "mse"], mse_long_T3[mse_long_T3$timepoint == o, "mse"])
  table_p_uni[2, o] <- ks$p.value
  table_w_uni[2, o] <- transport::wasserstein1d(mse_uni_training[mse_uni_training$timepoint == o, "mse"], mse_long_T3[mse_long_T3$timepoint == o, "mse"], p = 2)
}

table_uni <- cbind(round(table_p_uni, 3),
                   round(table_w_uni, 3)
)

# basic
desc_uni <- rbind(
  cbind(round(do.call(cbind, tapply(mse_uni_basic$mse, mse_uni_basic$timepoint, summary)), 3),
        round(do.call(cbind, tapply(mse_long_basic$mse, mse_long_basic$timepoint, summary)), 3)
  ),
  #training
  cbind(matrix(NA, 6, 3),
        round(do.call(cbind, tapply(mse_uni_training$mse, mse_uni_training$timepoint, summary)), 3),
        matrix(NA, 6, 3),
        round(do.call(cbind, tapply(mse_long_T3$mse, mse_long_T3$timepoint, summary)), 3)
  )
)
desc_uni <- round(desc_uni, 3)
colnames(desc_uni) <- rep(paste0("T", 1:G), 2)
rownames(desc_uni) <- rep(c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum"), 2)
desc_uni

print_tab(table_uni,
          name_rownames = "Predictors",
          file = paste0(dir_manuscript, "/tables/table_univariate.tex"),
          addr = list("pos" = list(-1),
                      "command" = c("\\hline \n & \\multicolumn{6}{c}{Kolmogorov-Smirnov p-value} & \\multicolumn{6}{c}{Wasserstein Distance} \\\\ \\cmidrule{2-7} \\cmidrule{8-13}")
          ),
          caption = "Tests for Differences in the Distributions of Mean Squared Error for the Multivariate versus the Univariate Random Forests.",
          label = "tb:test_mse_univariate")

print_tab(desc_uni,
          name_rownames = "",
          file = paste0(dir_manuscript, "/tables/descriptives_univariate.tex"),
          addr = list("pos" = list(-1, 0, 6),
                      "command" = c("\\hline \n & \\multicolumn{6}{c}{Univariate} & \\multicolumn{6}{c}{Multivariate} \\\\ \\cmidrule{2-7} \\cmidrule{8-13}",
                                   "\\hline \n Basic &&&&&&&&&&&& \\\\",
                                   "\\hline \n Training T3 &&&&&&&&&&&& \\\\ \\hline \n ")
          ),
          caption = "Descriptives on the Distributions of Mean Squared Error for the Univariate versus Multivariate Random Forests.",
          label = "tb:desc_mse_univariate")

####-------------------- MRFs versus regression -------------------####

mean_mse_reg <- readRDS("Data_analysis/mean_mse_reg.rds")

fs <- which(fit_list_short$outcome == "multi_long" & fit_list_short$predictors == "training")
fit_list_short[fs, ]

plot_list_reg <- vector("list", length = G)

for (f in fs) {
  mse_long <- data.frame("id" = rep(mean_mse_reg[[f]]$ids, G),
                         "timepoint" = rep(1:G, each = nrow(mean_mse_reg[[f]])),
                         "mse" = do.call(c, mean_mse_reg[[f]][ , 2:(G + 1)]))
  mse_long$timepoint <- factor(mse_long$timepoint)
  head(mse_long)


  plot_list_reg[[f - fs[1] + 1]] <- ggplot(data = mse_long, aes(y = mse, x=timepoint)) +
    geom_violin(show.legend=FALSE, fill = qualitative_hcl(3)[2]) +
    labs(y="MSE", x="Timepoint") +
    theme(axis.text=element_text(size=11),
          axis.title=element_text(size=11),
          title=element_text(size = 11)) +
    scale_y_continuous(limits = c(0, 6.5)) +
    ggtitle(label = paste0("Lasso Regression, Training T", f - fs[1] + 1))

}

ggsave(marrangeGrob(plot_list, layout_matrix = matrix(1:G, 2, 3, byrow = TRUE), top = NULL),
       file="plots/plots_mse_longitudinal_reg.pdf",
       width=30, height=15, units="cm")
file.copy(from = "plots/plots_mse_longitudinal.pdf", to = paste0(dir_manuscript, "/figures"), overwrite = TRUE)

ggsave(grid.arrange(plot_list[[3]] + scale_y_continuous(limits = c(0, 6.5)),
                    plot_list_reg[[3]] + scale_y_continuous(limits = c(0, 6.5)),
                    plot_list[[6]] + scale_y_continuous(limits = c(0, 6.5)),
                    plot_list_reg[[6]] + scale_y_continuous(limits = c(0, 6.5)),
                    nrow = 2, ncol = 2),
       file = "plots/plot_long_MRF-reg.pdf",
       width = 20, height = 15, units = "cm")
file.copy(from = "plots/plot_long_MRF-reg.pdf", to = paste0(dir_manuscript, "/figures"), overwrite = TRUE)


# Kolmogorov-Smirnov test and Wasserstein distance
# Wasserstein with p = 2 (squared distance)
# for each outcome (T1-T6) separately
table_w_reg <- table_p_reg <- matrix(NA, 6, 6)
rownames(table_p_reg) <- rownames(table_w_reg) <- paste0("T", 1:G)
colnames(table_p_reg) <- colnames(table_w_reg) <- paste0("T", 1:G)

for(o in 1:6) {
  for(d in 1:6) {
    ks <- ks.test(mean_mse[[fs[d]]][, o + 1], mean_mse_reg[[fs[d]]][, o + 1])
    table_w_reg[d, o] <- transport::wasserstein1d(mean_mse[[fs[d]]][, o + 1], mean_mse_reg[[fs[d]]][, o + 1], p = 2)
    table_p_reg[d, o] <- ks$p.value
  }
}

table_reg <- cbind(round(table_p_reg, 3),
                   round(table_w_reg, 3)
)

# exemplarily for up to T1 versus up to T6
desc_reg <- cbind(apply(mean_mse_reg[[fs[1]]][, -c(1,8)], 2, summary),
                         apply(mean_mse_reg[[fs[6]]][, -c(1,8)], 2, summary)
)
desc_reg <- round(desc_reg, 3)
colnames(desc_reg) <- rep(paste0("T", 1:G), 2)
rownames(desc_reg) <- c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum")
desc_reg

print_tab(table_reg,
          name_rownames = "Up to",
          file = paste0(dir_manuscript, "/tables/table_regression.tex"),
          addr = list("pos" = list(-1),
                      "command" = c("\\hline \n & \\multicolumn{6}{c}{Kolmogorov-Smirnov p-value} & \\multicolumn{6}{c}{Wasserstein Distance} \\\\ \\cmidrule{2-7} \\cmidrule{8-13}")
          ),
          caption = "Tests for Differences in the Distributions of Mean Squared Error for the Multivariate Random Forests versus Lasso Regressions.",
          label = "tb:test_mse_regression")

print_tab(desc_reg,
          name_rownames = "Statistic",
          file = paste0(dir_manuscript, "/tables/descriptives_regression.tex"),
          addr = list("pos" = list(-1),
                      "command" = c("\\hline \n & \\multicolumn{6}{c}{Training Predictors up to T1} & \\multicolumn{6}{c}{Training Predictors up to T6} \\\\ \\cmidrule{2-7} \\cmidrule{8-13}")
          ),
          caption = "Descriptives on the Distributions of Mean Squared Error for the Training Predictors up to T1 versus T6 for the Lasso Regressions.",
          label = "tb:desc_mse_regression")
