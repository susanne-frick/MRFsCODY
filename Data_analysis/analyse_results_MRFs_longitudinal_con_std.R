####-------------------- analyse results MRFs longitudinal outcome ------------------####

library(ggplot2)
library(colorspace)
library(gridExtra)
library(colorspace)

mean_mse <- readRDS("Data_analysis/mean_mse_MRFs_con_std.rds")

# observed variances per outcome from the complete data in each of the 16 fit_list_short cases
moments_split <- readRDS("Data_analysis/moments_splits_con_std.rds")

dir_manuscript <- "~/Dokumente/FAIR/Reha/paper/mrfs_cody_man"
dir_som <- "~/Dokumente/FAIR/Reha/paper/SOM"

# re-build design list
fit_list_short <- readRDS("Data_analysis/fit_list_short.rds")
G <- 6

color_columns <- function(df, identifier = "T") {
  df_unlist <- unlist(df[, grep(identifier, colnames(df))])
  df_unlist[is.na(df_unlist)] <- ""
  df[, grep(identifier, colnames(df))] <- matrix(paste(
    paste0("\\cellcolor{",
           ifelse(is.na(df), "white", ifelse(df < .05, "green!30", "white")),
           "}"),
    format(df_unlist, nsmall = 3)
  ), nrow = nrow(df))

  return(df)
}

# cutoffs for R^2 = 1 - MSE/Var(Outcome); with Var(Outcome) = 1
custom_cut <- function(vec) cut(vec,
                                breaks = c(-100, .09, .16, .25, 10),
                                include.lowest = TRUE,
                                labels = c("tiny", "small", "medium", "large"))


print_tab <- function(tab, name_rownames = NULL, file, addr = NULL, sideways = TRUE, ...) {
  # change rownames to first column if rownames title is given
  if(isFALSE(is.null(name_rownames))) {
    cn <- colnames(tab)
    tab <- data.frame(rownames(tab), tab)
    colnames(tab) <- c(name_rownames, cn)
    hl <- c(0, nrow(tab))
  } else {
    hl <- c(-1, 0, nrow(tab))
  }

  print(xtable::xtable(tab, digits=2, ...),
        include.colnames = T, include.rownames = is.null(name_rownames),
        hline.after = hl,
        add.to.row = addr,
        sanitize.rownames.function=function(x){x},
        sanitize.colnames.function = function(x){x},
        sanitize.text.function = function(x){x},
        NA.string = "", table.placement = "htp",
        caption.placement = "top",
        floating.environment = ifelse(sideways, "sidewaystable", "table"),
        file=file)
}


####------------------- basic predictors -----------------------####
f_multi_basic <- which(fit_list_short$outcome == "multi_long" & fit_list_short$predictors == "basic")
head(mean_mse[[f_multi_basic]])
mse_long_basic <- data.frame("id" = rep(mean_mse[[f_multi_basic]]$ids, G),
                             "timepoint" = rep(1:G, each = nrow(mean_mse[[f_multi_basic]])),
                             "mse" = do.call(c, mean_mse[[f_multi_basic]][ , 2:(G + 1)]))
mse_long_basic$timepoint <- factor(mse_long_basic$timepoint)
mse_long_basic$r2 <- 1 - mse_long_basic$mse/moments_split[[f_multi_basic]]$var[mse_long_basic$timepoint]
head(mse_long_basic)


# tried out log scale -> not used in favor of interpretability
plot_basic_long <- ggplot(data = mse_long_basic, aes(y = mse, x=timepoint)) +
# plot_basic_long <- ggplot(data = mse_long_basic, aes(y = -log(mse), x=timepoint)) +
  geom_violin(width = 1.2, show.legend=FALSE, fill = qualitative_hcl(3)[1]) +
  geom_boxplot(width = 0.2, alpha = 0.2, color = "black") +
  ylim(0, 5) +
  # labs(y="log(MSE)", x="Timepoint") +
  labs(y="MSE", x="Timepoint") +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11),
        title=element_text(size = 11)) +
  ggtitle(label = "Multivariate, Pre-test & Grade")


####---------------------- training predictors ----------------------------####

fs_multi_training <- which(fit_list_short$outcome == "multi_long" & fit_list_short$predictors == "training")
fit_list_short[fs_multi_training, ]

## re-structured: until T1,...T6 on the x axis, panel for each outcome
mse_long <- vector("list", length(fs_multi_training))
for (f in fs_multi_training) {
  mse_long[[f - fs_multi_training[1] + 1]] <- data.frame("id" = rep(mean_mse[[f]]$ids, G),
                         "timepoint" = rep(1:G, each = nrow(mean_mse[[f]])),
                         "mse" = do.call(c, mean_mse[[f]][ , 2:(G + 1)]))
  mse_long$timepoint <- factor(mse_long$timepoint)
  print(head(mse_long))
}
lapply(mse_long, nrow)
mse_long_until <- do.call(rbind, mse_long)
mse_long_until$until <- rep(1:G, do.call(c, lapply(mse_long, nrow)))
mse_long_until$until <- paste0("P(T", mse_long_until$until, ")")
mse_long_until$until <- factor(mse_long_until$until)

# mean variances for multi_long (MRFs) and training P(T1) to P(T6)
var_multi_training <- do.call(rbind, lapply(moments_split[fs_multi_training], function(m) m$var))
# MSE = Var(1 - r2)
cutoffs_multi_training <- apply(var_multi_training, 2, function(v) {
  cm <- data.frame(v %*% t(1 - c(.09, .16, .25)))
  colnames(cm) <- c("large", "medium", "small")
  cm$x <- seq(.5, 6, 1)
  cm$xend <- seq(1.5, 6.5, 1)
  cm
  }
)
cutoffs_multi_training

plot_list <- vector("list", length = G)

for (g in 1:G) {
  plot_list[[g]] <- ggplot(data = mse_long_until[mse_long_until$timepoint == g,], aes(y = mse, x=until)) +
    geom_violin(width = 1, show.legend=FALSE, fill = RColorBrewer::brewer.pal(n = 3, name = "Set3")[2]) +
    geom_boxplot(width = 0.2, alpha = 0.2) +
    geom_segment(aes(x = x, xend = xend, y = large, yend = large), dat = cutoffs_multi_training[[g]], color = "azure4", linetype = 3) +
    geom_segment(aes(x = x, xend = xend, y = medium, yend = medium), dat = cutoffs_multi_training[[g]], color = "azure4", linetype = 2) +
    geom_segment(aes(x = x, xend = xend, y = small, yend = small), dat = cutoffs_multi_training[[g]], color = "azure4", linetype = 1) +
    labs(y="MSE", x="Training Predictors") +
    theme(axis.text=element_text(size=11),
          axis.title=element_text(size=11),
          title=element_text(size = 11)) +
    scale_y_continuous(limits = c(0, 4)) + #(0, 6.5)
    ggtitle(label = paste0("O(T", g, ")"))

}
ggsave(marrangeGrob(plot_list, layout_matrix = matrix(1:G, 2, 3, byrow = TRUE), top = NULL),
       file="plots/plots_mse_longitudinal_until.pdf",
       width=30, height=15, units="cm")
file.copy(from = "plots/plots_mse_longitudinal_until.pdf", to = paste0(dir_manuscript, "/figures"), overwrite = TRUE)


# Kolmogorov-Smirnov test for T6, data until T6 vs. T5 to T1
lapply(mean_mse[fs_multi_training], nrow)

# Kolmogorov-Smirnov test and Wasserstein distance for predictors until T1 to T5 vs. T6
# Wasserstein with p = 2 (squared distance)
# for each outcome (T1-T6) separately
table_w_timepoints <- table_p_timepoints <- matrix(NA, 5, 6)
rownames(table_p_timepoints) <- rownames(table_w_timepoints) <- paste0("P(T", 1:(G-1), ")")
colnames(table_p_timepoints) <- colnames(table_w_timepoints) <- paste0("O(T", 1:G, ")")

for(o in 1:6) {
  for(d in 1:5) {
    ks <- ks.test(mean_mse[[fs_multi_training[6]]][, o + 1], mean_mse[[fs_multi_training[d]]][, o + 1])
    table_w_timepoints[d, o] <- transport::wasserstein1d(mean_mse[[fs_multi_training[6]]][, o + 1], mean_mse[[fs_multi_training[d]]][, o + 1], p = 2)
    table_p_timepoints[d, o] <- ks$p.value
  }
}

table_p_timepoints <- round(table_p_timepoints, 2)
table_w_timepoints <- round(table_w_timepoints, 2)
summary(c(table_w_timepoints))

# exemplarily for until T1 versus until T6
desc_timepoints <- rbind(
  cbind(apply(mean_mse[[fs_multi_training[1]]][, -c(1,8)], 2, summary),
        apply(mean_mse[[fs_multi_training[6]]][, -c(1,8)], 2, summary)),
  cbind(apply(1 - t(t(mean_mse[[fs_multi_training[1]]][, -c(1,8)])/moments_split[[fs_multi_training[1]]]$var), 2, function(cl) table(custom_cut(cl))/length(cl)),
        apply(1 - t(t(mean_mse[[fs_multi_training[6]]][, -c(1,8)])/moments_split[[fs_multi_training[6]]]$var), 2, function(cl) table(custom_cut(cl))/length(cl)))[4:1,])
desc_timepoints <- round(desc_timepoints, 2)
colnames(desc_timepoints) <- rep(paste0("O(T", 1:G, ")"), 2)
rownames(desc_timepoints) <- c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum",
                               "\\% large", "\\% medium", "\\% small", "\\% tiny")
desc_timepoints

# across all timepoints for number in text
tb_cutoffs_timepoints <- matrix(NA, G, 4)
colnames(tb_cutoffs_timepoints) <- c("tiny", "small", "medium", "large")
rownames(tb_cutoffs_timepoints) <- paste0("P", 1:G)
for(g in 1:G) tb_cutoffs_timepoints[g,] <-  table(custom_cut(unlist(1 - t(t(mean_mse[[fs_multi_training[g]]][, -c(1,8)])/moments_split[[fs_multi_training[g]]]$var))))
tb_cutoffs_timepoints

sum(tb_cutoffs_timepoints[, "large"])/sum(tb_cutoffs_timepoints)
sum(tb_cutoffs_timepoints[, "tiny"])/sum(tb_cutoffs_timepoints)

print_tab(desc_timepoints,
          name_rownames = "Statistic",
          file = paste0(dir_som, "/tables/descriptives_timepoints.tex"),
          addr = list("pos" = list(-1, 6),
                      "command" = c("\\hline \n & \\multicolumn{6}{c}{P(T1)} & \\multicolumn{6}{c}{P(T6)} \\\\ \\cmidrule{2-7} \\cmidrule{8-13}",
                                    "\\hline \n Explained Var. & \\multicolumn{12}{c}{} \\\\ \\hline \n")
          ),
          caption = "Descriptives on the Distributions of Mean Squared Error for the Training Predictors until T1 P(T1) versus T6 P(T6) for the Multivariate Random Forests.",
          label = "tb:desc_mse_timepoints")

####-------------- basic versus training predictors ----------------------------------####

# Kolmogorov-Smirnov test and Wasserstein distance
table_w_basic <- table_p_basic <- matrix(NA, 6, 6)
rownames(table_p_basic) <- rownames(table_w_basic) <- paste0("P(T", 1:G, ")")
colnames(table_p_basic) <- colnames(table_w_basic) <- paste0("O(T", 1:G, ")")

for(o in 1:6) {
  for(d in 1:6) {
    ks <- ks.test(mean_mse[[f_multi_basic]][, o + 1], mean_mse[[fs_multi_training[d]]][, o + 1])
    table_w_basic[d, o] <- transport::wasserstein1d(mean_mse[[f_multi_basic]][, o + 1], mean_mse[[fs_multi_training[d]]][, o + 1], p = 2)
    table_p_basic[d, o] <- ks$p.value
  }
}
# ties in the training data

table_p_basic <- round(table_p_basic, 2)
table_w_basic <-  round(table_w_basic, 2)
summary(c(table_w_basic))

# descriptives multi basic
desc_basic <- rbind(apply(mean_mse[[f_multi_basic]][, -c(1,8)], 2, summary),
                    apply(1 - t(t(mean_mse[[f_multi_basic]][, -c(1,8)])/moments_split[[f_multi_basic]]$var), 2, function(cl) table(custom_cut(cl))/length(cl))[4:1,])
desc_basic <- round(desc_basic, 2)
colnames(desc_basic) <- paste0("O(T", 1:G, ")")
rownames(desc_basic) <- c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum",
                          "\\% large", "\\% medium", "\\% small", "\\% tiny")
desc_basic

# across all timepoints for number in text
tb_cutoffs_basic <- table(custom_cut(unlist(1 - t(t(mean_mse[[f_multi_basic]][, -c(1,8)])/moments_split[[f_multi_basic]]$var))))
round(tb_cutoffs_basic/sum(tb_cutoffs_basic), 2)

print_tab(desc_basic,
          name_rownames = "Statistic",
          file = paste0(dir_som, "/tables/descriptives_basic.tex"),
          addr = list("pos" = list(-1, 6),
                      "command" = c("\\hline \n",
                                    "\\hline \n Explained Var. & \\multicolumn{6}{c}{} \\\\ \\hline \n")
          ),
          caption = "Descriptives on the Distributions of Mean Squared Error for the Basic Predictors for the Multivariate Random Forests.",
          label = "tb:desc_mse_basic")

#####---------------- joint tables for timepoints and basic vs. training --------------####

table_w_timepoints
table_w_basic

table_w_basic_timepoints <- cbind(table_w_basic, rbind(table_w_timepoints, rep(NA, ncol(table_w_timepoints))))
table_w_basic_timepoints

print_tab(table_w_basic_timepoints,
          name_rownames = "Predictors",
          file = paste0(dir_manuscript, "/tables/table_predictors.tex"),
          addr = list("pos" = list(-1),
                      "command" = c("\\hline \n & \\multicolumn{6}{c}{Basic versus Training} & \\multicolumn{6}{c}{Versus P(T6)} \\\\ \\cmidrule{2-7} \\cmidrule{8-13}")
          ),
          caption = "Wasserstein Distances between the Distributions of Mean Squared Error Based on Different Predictor Sets for the Multivariate Random Forests.",
          label = "tb:wasserstein_mse_predictors")


table_p_timepoints
table_p_basic

table_p_basic_timepoints <- cbind(table_p_basic, rbind(table_p_timepoints, rep(NA, ncol(table_p_timepoints))))
table_p_basic_timepoints <- color_columns(table_p_basic_timepoints)

print_tab(table_p_basic_timepoints,
          name_rownames = "Predictors",
          file = paste0(dir_som, "/tables/table_tests_predictors.tex"),
          addr = list("pos" = list(-1),
                      "command" = c("\\hline \n & \\multicolumn{6}{c}{Basic versus Training} & \\multicolumn{6}{c}{Versus P(T6)} \\\\ \\cmidrule{2-7} \\cmidrule{8-13}")
          ),
          caption = "Kolmogorov-Smirnov p-values for Comparisons of the Distributions of Mean Squared Error Based on Different Predictor Sets for the Multivariate Random Forests.",
          label = "tb:ks_mse_predictors")

####-------------- comparison with univariate sum score, separately for the 6 timepoints --------####

# basic predictors
fs_uni_basic <- which(fit_list_short$outcome == "uni" & fit_list_short$predictors == "basic")
lapply(mean_mse[fs_uni_basic], nrow)

mse_uni_basic <- data.frame("id" = do.call(c, lapply(mean_mse[fs_uni_basic], function(m) m$ids)),
                             "timepoint" = do.call(c, lapply(1:G, function(g, m) rep(g, nrow(m[[g]])), m = mean_mse[fs_uni_basic])),
                             "mse" = do.call(c, lapply(mean_mse[fs_uni_basic], function(m) m[,2])))
mse_uni_basic$timepoint <- factor(mse_uni_basic$timepoint)
mse_uni_basic$r2 <- 1 - mse_uni_basic$mse/do.call(c, lapply(moments_split[fs_uni_basic], function(m) m$var))[mse_uni_basic$timepoint]
head(mse_uni_basic)


# training predictors
# ! until T3 for T4-T6
fs_uni_training3 <- which(fit_list_short$outcome == "uni" & fit_list_short$predictors == "training3")
lapply(mean_mse[fs_uni_training3], nrow)

mse_uni_training <- data.frame("id" = do.call(c, lapply(mean_mse[fs_uni_training3], function(m) m$ids)),
                            "timepoint" = do.call(c, lapply(1:3, function(g, m) rep(g, nrow(m[[g]])), m = mean_mse[fs_uni_training3])),
                            "mse" = do.call(c, lapply(mean_mse[fs_uni_training3], function(m) m[,2])))
mse_uni_training$timepoint <- mse_uni_training$timepoint + 3
mse_uni_training$timepoint <- factor(mse_uni_training$timepoint)
mse_uni_training$r2 <- 1 - mse_uni_training$mse/do.call(c, lapply(moments_split[fs_uni_training3], function(m) m$var))[mse_uni_training$timepoint]
head(mse_uni_training)


f_multi_training3 <- which(fit_list_short$outcome == "multi_long" & fit_list_short$predictors == "training" & fit_list_short$timepoint == 3)
fit_list_short[f_multi_training3, ]

mse_long_T3 <- data.frame("id" = rep(mean_mse[[f_multi_training3]]$ids, G),
                          "timepoint" = rep(1:G, each = nrow(mean_mse[[f_multi_training3]])),
                          "mse" = do.call(c, mean_mse[[f_multi_training3]][ , 2:(G + 1)]))
mse_long_T3 <- mse_long_T3[mse_long_T3$timepoint > 3,]
mse_long_T3$timepoint <- factor(mse_long_T3$timepoint)
mse_long_T3$r2 <- 1 - mse_long_T3$mse/moments_split[[f_multi_training3]]$var[mse_long_T3$timepoint]
head(mse_long_T3)


# Kolmogorov-Smirnov test and Wasserstein distance for predictors until T1 to T5 vs. T6
# Wasserstein with p = 2 (squared distance)
# for each outcome (T1-T6) separately
table_w_uni <- table_p_uni <- matrix(NA, 2, 6)
rownames(table_p_uni) <- rownames(table_w_uni) <- c("Basic", "P(T3)")
colnames(table_p_uni) <- colnames(table_w_uni) <- paste0("O(T", 1:G, ")")

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

table_p_uni <- round(table_p_uni, 2)
table_w_uni <-  round(table_w_uni, 2)
summary(c(table_w_uni))


# basic
desc_uni <- rbind(
  cbind(round(do.call(cbind, tapply(mse_uni_basic$mse, mse_uni_basic$timepoint, summary)), 2),
        round(do.call(cbind, tapply(mse_long_basic$mse, mse_long_basic$timepoint, summary)), 2)
  ),
  cbind(round(do.call(cbind, tapply(mse_uni_basic$r2, mse_uni_basic$timepoint, function(cl) table(custom_cut(cl))/length(cl))), 2),
        round(do.call(cbind, tapply(mse_long_basic$r2, mse_long_basic$timepoint, function(cl) table(custom_cut(cl))/length(cl))), 2)
  )[4:1,],
  #training
  cbind(matrix(NA, 6, 3),
        round(do.call(cbind, tapply(mse_uni_training$mse, mse_uni_training$timepoint, summary)), 2),
        matrix(NA, 6, 3),
        round(do.call(cbind, tapply(mse_long_T3$mse, mse_long_T3$timepoint, summary)), 2)
  ),
  cbind(matrix(NA, 4, 3),
        round(do.call(cbind, tapply(mse_uni_training$r2, mse_uni_training$timepoint, function(cl) table(custom_cut(cl))/length(cl))), 2)[4:1,],
        matrix(NA, 4, 3),
        round(do.call(cbind, tapply(mse_long_T3$r2, mse_long_T3$timepoint, function(cl) table(custom_cut(cl))/length(cl))), 2)[4:1,]
  )
)
desc_uni <- round(desc_uni, 2)
colnames(desc_uni) <- rep(paste0("O(T", 1:G, ")"), 2)
rownames(desc_uni) <- rep(c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum",
                            "\\% large", "\\% medium", "\\% small", "\\% tiny"), 2)
desc_uni

tb_cutoffs_uni <- table(custom_cut(unlist(mse_uni_training$r2)))
round(tb_cutoffs_uni/sum(tb_cutoffs_uni), 2)

print_tab(desc_uni,
          name_rownames = "Statistic",
          file = paste0(dir_som, "/tables/descriptives_univariate.tex"),
          addr = list("pos" = list(-1, 0, 6, 10, 16),
                      "command" = c("\\hline \n & \\multicolumn{6}{c}{Univariate} & \\multicolumn{6}{c}{Multivariate} \\\\ \\cmidrule{2-7} \\cmidrule{8-13}",
                                   "\\hline \n Basic &&&&&&&&&&&& \\\\",
                                   "\\hline \n Explained Var. & \\multicolumn{12}{c}{} \\\\ \\hline \n",
                                   "\\hline \n P(T3) &&&&&&&&&&&& \\\\ \\hline \n ",
                                   "\\hline \n Explained Var. & \\multicolumn{12}{c}{} \\\\ \\hline \n")
          ),
          caption = "Descriptives on the Distributions of Mean Squared Error for the Univariate versus Multivariate Random Forests.",
          label = "tb:desc_mse_univariate")

####-------------------- MRFs versus regression -------------------####

mean_mse_reg <- readRDS("Data_analysis/mean_mse_reg_con_std.rds")

fit_list_short[fs_multi_training, ]

# Kolmogorov-Smirnov test and Wasserstein distance
# Wasserstein with p = 2 (squared distance)
# for each outcome (T1-T6) separately
table_w_reg <- table_p_reg <- matrix(NA, 6, 6)
rownames(table_p_reg) <- rownames(table_w_reg) <- paste0("P(T", 1:G, ")")
colnames(table_p_reg) <- colnames(table_w_reg) <- paste0("O(T", 1:G, ")")

for(o in 1:6) {
  for(d in 1:6) {
    ks <- ks.test(mean_mse[[fs_multi_training[d]]][, o + 1], mean_mse_reg[[fs_multi_training[d]]][, o + 1])
    table_w_reg[d, o] <- transport::wasserstein1d(mean_mse[[fs_multi_training[d]]][, o + 1], mean_mse_reg[[fs_multi_training[d]]][, o + 1], p = 2)
    table_p_reg[d, o] <- ks$p.value
  }
}

table_p_reg <- round(table_p_reg, 2)
table_w_reg <- round(table_w_reg, 2)
summary(c(table_w_reg))

# exemplarily for until T1 versus until T6 (fit_list_short was the same as for MRFs)
desc_reg <- rbind(cbind(apply(mean_mse_reg[[fs_multi_training[1]]][, -c(1,8)], 2, summary),
                        apply(mean_mse_reg[[fs_multi_training[6]]][, -c(1,8)], 2, summary)),
                  cbind(apply(1 - t(t(mean_mse_reg[[fs_multi_training[1]]][, -c(1,8)])/moments_split[[fs_multi_training[1]]]$var), 2, function(cl) table(custom_cut(cl))/length(cl)),
                        apply(1 - t(t(mean_mse_reg[[fs_multi_training[6]]][, -c(1,8)])/moments_split[[fs_multi_training[1]]]$var), 2, function(cl) table(custom_cut(cl))/length(cl)))[4:1,]
)
desc_reg <- round(desc_reg, 2)
colnames(desc_reg) <- rep(paste0("O(T", 1:G, ")"), 2)
rownames(desc_reg) <- c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum",
                        "\\% large", "\\% medium", "\\% small", "\\% tiny")
desc_reg

# across all timepoints for number in text
tb_cutoffs_multi_training_reg <- matrix(NA, G, 4)
colnames(tb_cutoffs_multi_training_reg) <- c("tiny", "small", "medium", "large")
rownames(tb_cutoffs_multi_training_reg) <- paste0("P", 1:G)
for(g in 1:G) tb_cutoffs_multi_training_reg[g,] <-  table(custom_cut(unlist(1 - t(t(mean_mse[[fs_multi_training[g]]][, -c(1,8)])/moments_split[[fs_multi_training[g]]]$var))))
tb_cutoffs_multi_training_reg

sum(tb_cutoffs_multi_training_reg[, "large"])/sum(tb_cutoffs_multi_training_reg)
sum(tb_cutoffs_multi_training_reg[, "tiny"])/sum(tb_cutoffs_multi_training_reg)


print_tab(desc_reg,
          name_rownames = "Statistic",
          file = paste0(dir_som, "/tables/descriptives_regression.tex"),
          addr = list("pos" = list(-1, 6),
                      "command" = c("\\hline \n & \\multicolumn{6}{c}{P(T1)} & \\multicolumn{6}{c}{Ṕ(T6)} \\\\ \\cmidrule{2-7} \\cmidrule{8-13}",
                                    "\\hline \n Explained Var. & \\multicolumn{12}{c}{} \\\\ \\hline \n")
          ),
          caption = "Descriptives on the Distributions of Mean Squared Error for the Training Predictors until T1 P(T1) versus T6 P(T6) for the Lasso Regressions.",
          label = "tb:desc_mse_regression")

####----------------------- joint tables for multi vs. uni vs. regression until T3 ------------------------------------------####

table_w_reg
table_w_uni

table_w_reg_uni <- rbind(table_w_reg, table_w_uni)
rownames(table_w_reg_uni) <- sub("Training T3", "P(T3)", rownames(table_w_reg_uni))

print_tab(table_w_reg_uni,
          name_rownames = "Predictors",
          sideways = FALSE,
          file = paste0(dir_manuscript, "/tables/table_regression_univariate.tex"),
          addr = list("pos" = list(-1, nrow(table_w_reg)),
                      "command" = c("\\hline \n & \\multicolumn{6}{c}{Lasso Regression} \\\\ \\cmidrule{2-7}",
                                    "\\hline \n & \\multicolumn{6}{c}{Univariate Random Forests} \\\\ \\cmidrule{2-7} ")
          ),
          caption = "Wasserstein Distances between the Distributions of Mean Squared Error for the Multivariate Random Forests versus Lasso Regressions and Univariate Random Forests.",
          label = "tb:wasserstein_mse_models")

table_p_reg
table_p_uni

table_p_reg_uni <- rbind(table_p_reg, table_p_uni)
table_p_reg_uni <- color_columns(table_p_reg_uni)

print_tab(table_p_reg_uni,
          name_rownames = "Predictors",
          sideways = FALSE,
          file = paste0(dir_som, "/tables/table_tests_regression_univariate.tex"),
          addr = list("pos" = list(-1, nrow(table_p_reg)),
                      "command" = c("\\hline \n & \\multicolumn{6}{c}{Lasso Regression} \\\\ \\cmidrule{2-7}",
                                    "\\hline \n & \\multicolumn{6}{c}{Univariate Random Forests} \\\\ \\cmidrule{2-7} ")
          ),
          caption = "Kolmogorov-Smirnov p-values for Comparisons of the Distributions of Mean Squared Error for the Multivariate Random Forests versus Lasso Regressions and Univariate Random Forests.",
          label = "tb:ks_mse_models")


####----------------------- plot multi vs. uni vs. regression until T3 ------------------------------------------####

# another version: uni vs. multi next to each other

mse_reg_T3 <- data.frame("id" = rep(mean_mse_reg[[f_multi_training3]]$ids, G),
                          "timepoint" = rep(1:G, each = nrow(mean_mse_reg[[f_multi_training3]])),
                          "mse" = do.call(c, mean_mse_reg[[f_multi_training3]][ , 2:(G + 1)]))
mse_reg_T3 <- mse_reg_T3[mse_reg_T3$timepoint > 3,]
mse_reg_T3$timepoint <- factor(mse_reg_T3$timepoint)
head(mse_reg_T3)

mse_reg_basic <- data.frame("id" = rep(mean_mse_reg[[f_multi_training3]]$ids, G),
                         "timepoint" = rep(1:G, each = nrow(mean_mse_reg[[f_multi_training3]])),
                         "mse" = do.call(c, mean_mse_reg[[f_multi_training3]][ , 2:(G + 1)]))
mse_reg_basic$timepoint <- factor(mse_reg_basic$timepoint)
head(mse_reg_basic)

# P(T1) to P(T3) for O(T4 to O(T6))
fs_uni_training <- which(fit_list_short$outcome == "uni" & fit_list_short$predictors == "training3")
var_uni_training <- do.call(rbind, lapply(moments_split[fs_uni_training], function(m) m$var))
cutoffs_uni_training <- data.frame(var_uni_training %*% t(1 - c(.09, .16, .25)))
colnames(cutoffs_uni_training) <- c("large", "medium", "small")

cutoffs_P3 <- rbind(
  do.call(rbind, lapply(cutoffs_multi_training[4:6], function(cm) cm[3, c("large", "medium", "small")])), #multi
  cutoffs_uni_training,
  do.call(rbind, lapply(cutoffs_multi_training[4:6], function(cm) cm[3, c("large", "medium", "small")])) # same for reg (same datasets)
)
cutoffs_P3$Model <- rep(c("multivariate", "univariate", "lasso"), each = 3)
cutoffs_P3$x <- c(seq(.5, 2.5, by = 1),
                  seq(.5 + 1/3, 2.5 + 1/3, by = 1),
                  seq(.5 + 2/3, 2.5 + 2/3, by = 1))
cutoffs_P3$xend <- c(seq(.5 + 1/3, 2.5 + 1/3, by = 1),
                     seq(.5 + 2/3, 2.5 + 2/3, by = 1),
                     seq(.5 + 3/3, 2.5 + 3/3, by = 1))
cutoffs_P3

# basic
var_uni_basic <- do.call(rbind, lapply(moments_split[fs_uni_basic], function(m) m$var))
cutoffs_uni_basic <- data.frame(var_uni_basic %*% t(1 - c(.09, .16, .25)))
colnames(cutoffs_uni_basic) <- c("large", "medium", "small")

var_multi_basic <- moments_split[[f_multi_basic]]$var
cutoffs_multi_basic <- data.frame(var_multi_basic %*% t(1 - c(.09, .16, .25)))
colnames(cutoffs_multi_basic) <- c("large", "medium", "small")

cutoffs_basic <- rbind(cutoffs_multi_basic, cutoffs_uni_basic, cutoffs_multi_basic) # same for reg (same datasets)
cutoffs_basic$Model <- rep(c("multivariate", "univariate", "lasso"), each = 6)
cutoffs_basic$x <- c(seq(.5, 5.5, by = 1),
                  seq(.5 + 1/3, 5.5 + 1/3, by = 1),
                  seq(.5 + 2/3, 5.5 + 2/3, by = 1))
cutoffs_basic$xend <- c(seq(.5 + 1/3, 5.5 + 1/3, by = 1),
                     seq(.5 + 2/3, 5.5 + 2/3, by = 1),
                     seq(.5 + 3/3, 5.5 + 3/3, by = 1))
cutoffs_basic

mse_long_uni_T3 <- rbind(mse_long_T3[, -c(4)], mse_uni_training[, -c(4)], mse_reg_T3) # without r2
mse_long_uni_T3$timepoint <- paste0("O(T", mse_long_uni_T3$timepoint, ")")
mse_long_uni_T3$Model <- rep(c("multivariate", "univariate", "lasso"),
                               c(nrow(mse_long_T3), nrow(mse_uni_training), nrow(mse_reg_T3)))

mse_long_uni_basic <- rbind(mse_long_basic[, -c(4)], mse_uni_basic[, -c(4)], mse_reg_basic) # without r2
mse_long_uni_basic$timepoint <- paste0("O(T", mse_long_uni_basic$timepoint, ")")
mse_long_uni_basic$Model <- rep(c("multivariate", "univariate", "lasso"),
                               c(nrow(mse_long_basic), nrow(mse_uni_basic), nrow(mse_reg_basic)))

plot_long_uni_T3 <- ggplot(data = mse_long_uni_T3, aes(y = mse, x=timepoint, fill=Model)) +
  geom_violin(show.legend = TRUE, position = position_dodge(0.9)) +
  geom_boxplot(width = 0.2, alpha = 0.2, position = position_dodge(0.9), show.legend = FALSE) +
  geom_segment(aes(x = x, xend = xend, y = large, yend = large), dat = cutoffs_P3, color = "azure4", linetype = 3) +
  geom_segment(aes(x = x, xend = xend, y = medium, yend = medium), dat = cutoffs_P3, color = "azure4", linetype = 2) +
  geom_segment(aes(x = x, xend = xend, y = small, yend = small), dat = cutoffs_P3, color = "azure4", linetype = 1) +
  ylim(0, 4) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(n = 3, name = "Set3")) +
  labs(y="MSE", x="Status Test") +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11),
        title=element_text(size = 11)) +
  ggtitle(label = "P(T3)")

plot_long_uni_basic <- ggplot(data = mse_long_uni_basic, aes(y = mse, x=timepoint, fill=Model)) +
  geom_violin(show.legend = FALSE, position = position_dodge(0.9)) +
  geom_boxplot(width = 0.2, alpha = 0.2, position = position_dodge(0.9), show.legend = FALSE) +
  geom_segment(aes(x = x, xend = xend, y = large, yend = large), dat = cutoffs_basic, color = "azure4", linetype = 3) +
  geom_segment(aes(x = x, xend = xend, y = medium, yend = medium), dat = cutoffs_basic, color = "azure4", linetype = 2) +
  geom_segment(aes(x = x, xend = xend, y = small, yend = small), dat = cutoffs_basic, color = "azure4", linetype = 1) +
  ylim(0, 4) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(n = 3, name = "Set3")) +
  labs(y="MSE", x="Status Test") +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11),
        title=element_text(size = 11)) +
  ggtitle(label = "Basic Predictors")

ggsave(grid.arrange(plot_long_uni_basic, plot_long_uni_T3,
                    nrow = 1, ncol = 2, widths = c(3,2)),
       file = "plots/plot_multi-uni-reg.pdf",
       width = 30, height = 12, units = "cm")
file.copy(from = "plots/plot_multi-uni-reg.pdf", to = paste0(dir_manuscript, "/figures"), overwrite = TRUE)
