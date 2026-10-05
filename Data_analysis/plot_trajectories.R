####----------------------------------- plot predicted trajectories ---------------------####

# load functions
library(colorspace)
library(gam)
devtools::load_all()

recode.df <- function(variable, levels.old, levels.new) {
  levels.new[match(as.character(variable), as.character(levels.old))]
}

# load fit lists
fit_list_short <- readRDS("Data_analysis/fit_list_short.rds")

# load data splits
data_splitted <- readRDS("Data_analysis/data_splitted.rds")

get_f <- function(outcome, predictors, timepoint, fl = fit_list_short) {
  which(fl$outcome == outcome &
          fl$predictors == predictors &
          fl$timepoint == timepoint)
}

G <- 6

tp <- 5 # training up to tp
f <- get_f("multi_long", "training", tp)
# fit_f <- readRDS(paste0("results_MRFs/fit_MRF_f", f, ".rds"))
fit_f <- readRDS(paste0("results_MRFs_con_std/fit_MRF_con_std_f", f, ".rds"))

# obtain predictions for test data each
predictions_f <- vector("list", length(fit_f))
for (s in 1:length(fit_f)) {
  predictions_f[[s]] <- obtain_predictions(data_splitted[[f]][[s]]$test, fit_f[[s]], avg = FALSE)
}
predictions_f <- do.call(rbind, predictions_f)
colnames(predictions_f)[2:7] <- paste0("obs_", colnames(predictions_f)[2:7])
# columns of output: y_obs, x_obs, y_prd
head(predictions_f)

# differences between start and end
predictions_f$diff_T6_T1 <- predictions_f$Score_gold_sum_6 - predictions_f$Score_gold_sum_1
predictions_f$mse <- rowMeans((predictions_f[, grep("^Score_gold_sum", colnames(predictions_f))] -
                                 predictions_f[, grep("obs_Score_gold_sum", colnames(predictions_f))])^2)

plot_spaghetti <- function(dat, plot_cols, color_col, n_color = 4, main = "Status Test",
                           xlab = "Training Day", labels_x = TRUE, yl = NULL,
                           deficiency = NULL,
                           file = NULL) {
  # deficiency: deutan, tritan, protan = simulation functions in colorspace
  dat_plot <- dat[, plot_cols]

  if(length(unique(dat[, color_col])) <= n_color) {
    br <- sort(unique(dat[, color_col]))
    cols <- recode.df(dat[, color_col],
                      br,
                      # sequential_hcl(length(br) + 1, palette = "YlGnBu", rev = TRUE)[-c(1)]
                      rev(RColorBrewer::brewer.pal(min(length(br), 11), name = "RdYlBu"))
                      )
  } else {
    br <- quantile(dat[, color_col], seq(0, 1, 1/n_color))
    cols <- as.character(cut(dat[, color_col],
                             breaks = br,
                             labels = rev(RColorBrewer::brewer.pal(min(length(br) - 1, 11), name = "RdYlBu"))
                             # labels = sequential_hcl(length(br), palette = "YlGnBu", rev = TRUE)[-c(1)])
                             ))
  }

  if(isFALSE(is.null(deficiency))) cols <- deficiency(cols)
  if (is.null(yl)) yl <- range(dat_plot, na.rm = TRUE)

  pdf(file = file, width = 7, height = 6)
  par(mar = c(5,4,4,7)) #default c(5,4,4,2) + 0.1
  plot(1:ncol(dat_plot), rep(0, ncol(dat_plot)),
       ylim = yl,
       ylab = "Score", xlab = xlab, main = main, type = "n", xaxt = "none")
  axis(side = 1, at = 1:ncol(dat_plot), labels = labels_x)
  for(p in 1:nrow(dat_plot)) lines(1:ncol(dat_plot), dat_plot[p,], col = cols[p])
  dev.off()

}


plot_spaghetti(predictions_f,
               plot_cols = grep("obs_Score_gold_sum", colnames(predictions_f)),
               color_col = "obs_Score_gold_sum_1",
               main = "Observed",
               xlab = "Status Test", labels_x = paste0("O(T", 1:G, ")"),
               yl = c(-2.5, 3),
               file = "plots/plot_spaghetti_gold_observed.pdf")

plot_spaghetti(predictions_f,
               plot_cols = grep("obs_Score_gold_sum", colnames(predictions_f)),
               color_col = "Grade",
               main = "Observed",
               file = "plots/plot_spaghetti_gold_observed_coloredByGrade.pdf")

plot_spaghetti(predictions_f,
               plot_cols = grep("^Score_gold_sum", colnames(predictions_f)),
               color_col = "Score_gold_sum_1",
               main = "Predicted",
               xlab = "Status Test",
               yl = c(-2.5, 3), labels_x = paste0("O(T", 1:G, ")"),
               file = "plots/plot_spaghetti_gold_predicted.pdf")

plot_spaghetti(predictions_f,
               plot_cols = grep("^Score_gold_sum", colnames(predictions_f)),
               color_col = "Score_cody_1",
               main = "Predicted - Colored Pre-Test",
               file = "plots/plot_spaghetti_gold_pretest.pdf")

plot_spaghetti(predictions_f,
               plot_cols = grep("^Score_gold_sum", colnames(predictions_f)),
               color_col = "Grade",
               n_color = 3,
               main = "Predicted - Colored Grade",
               file = "plots/plot_spaghetti_gold_grade.pdf")

plot_spaghetti(predictions_f,
               plot_cols = grep("^Score_gold_sum", colnames(predictions_f)),
               color_col = "leveldif_PC2",
               main = "Level Difference PC2",
               file = "plots/plot_spaghetti_gold_leveldifPC2.pdf")

plot_spaghetti(predictions_f,
               plot_cols = grep("^Score_gold_sum", colnames(predictions_f)),
               color_col = "score_PC1",
               main = "Level PC1",
               file = "plots/plot_spaghetti_gold_levelPC1.pdf")

plot_spaghetti(predictions_f,
               plot_cols = grep("^Score_gold_sum", colnames(predictions_f)),
               color_col = "diff_T6_T1",
               main = "Predicted - Colored T6 - T1",
               file = "plots/plot_spaghetti_gold_T6minusT1.pdf")

plot_spaghetti(predictions_f,
               plot_cols = grep("^Score_gold_sum", colnames(predictions_f)),
               color_col = "mse",
               main = "Predicted - Mean Squared Error",
               file = "plots/plot_spaghetti_gold_MSE.pdf")

####------------------------- mean trajectories based on GAM --------------------------####
# GAM by grade

obtain_gam <- function(dat, plot_cols){
  dat_plot <- dat[, plot_cols]
  head(dat_plot)
  dat_long <- data.frame(y = unlist(dat_plot),
                         time = rep(1:ncol(dat_plot), each = nrow(dat_plot)),
                         row.names = NULL)
  head(dat_long)
  fit <- gam(y  ~  s(time), data = dat_long)
  plot(fit, ylim = c(-1,1))
  dispersion  <- summary(fit)$dispersion # standard deviation
  pred <- predict(fit, newdata = data.frame(time = 1:ncol(dat_plot))) # kurve
  return(list(pred = pred, disp = dispersion))
}

pred_grade2 <- obtain_gam(predictions_f[predictions_f$Grade == 2,],
                          plot_cols = grep("obs_Score_gold_sum", colnames(predictions_f)))
pred_grade3 <- obtain_gam(predictions_f[predictions_f$Grade == 3,],
                          plot_cols = grep("obs_Score_gold_sum", colnames(predictions_f)))
pred_grade4 <- obtain_gam(predictions_f[predictions_f$Grade == 4,],
                          plot_cols = grep("obs_Score_gold_sum", colnames(predictions_f)))

cols <- RColorBrewer::brewer.pal(3, "Dark2")
dat_plot <- predictions_f[, grep("obs_Score_gold_sum", colnames(predictions_f))]

pdf(file = "plots/plot_spaghetti_gam_observed.pdf", width = 7, height = 6)
par(mar = c(5,4,4,5), xpd = TRUE)
plot(1:ncol(dat_plot), rep(0, ncol(dat_plot)),
     ylim = range(dat_plot, na.rm = TRUE),
     ylab = "Score", xlab = "Status Test", main = "Observed Trajectories", type = "n",
     xaxt = "n")
axis(side = 1, at = 1:ncol(dat_plot), labels = paste0("O(T", 1:ncol(dat_plot), ")"))
for(p in 1:nrow(dat_plot)) lines(1:ncol(dat_plot), dat_plot[p,], col = "lightgrey")
lines(1:G, pred_grade2$pred, type = "l", col = cols[1])
polygon(x = c(1:G, G:1),
        y = c(pred_grade2$pred - pred_grade2$disp, pred_grade2$pred + pred_grade2$disp),
        col = adjustcolor(cols[1], alpha.f = .2), border = NA)
lines(1:G, pred_grade3$pred, col = cols[2])
polygon(x = c(1:G, G:1),
        y = c(pred_grade3$pred - pred_grade3$disp, pred_grade3$pred + pred_grade3$disp),
        col = adjustcolor(cols[2], alpha.f = .2), border = NA)
lines(1:G, pred_grade4$pred, col = cols[3])
polygon(x = c(1:G, G:1),
        y = c(pred_grade4$pred - pred_grade4$disp, pred_grade4$pred + pred_grade4$disp),
        col = adjustcolor(cols[3], alpha.f = .2), border = NA)
abline(h = 0, col = adjustcolor("black", alpha.f = .8), xpd = FALSE, lty = "dashed")
legend(x = 6.3, y = 3,
       legend = 4:2, title = "Grade",
       fill = rev(cols),
       bty = "n")
dev.off()

####--------------------- regressions -----------------------####

library(glmnet)

res <- readRDS("Data_analysis/results_reg_con_std.rds")

tp <- 5 # training up to tp
f <- get_f("multi_long", "training", tp)
fs <- get_f("multi_long", "training", tp, fl = fit_list)

# obtain predictions for test data each
predictions_reg_f <- vector("list", length(fs))
for (s in 1:length(fs)) {

  dat <- data_splitted[[f]][[s]]$test
  dat <- dat[, grep("Id|^Score_gold_sum", colnames(dat), invert = TRUE)]

  predictions_reg_f[[s]] <- cbind(data_splitted[[f]][[s]]$test,
                                  predict(object = res[[fs[s]]]$fit, s = "lambda.min", newx = as.matrix(dat))[,,1])

}
predictions_reg_f <- do.call(rbind, predictions_reg_f)
colnames(predictions_reg_f)[2:7] <- paste0("obs_", colnames(predictions_reg_f)[2:7])
# columns of output: y_obs, x_obs, y_prd
head(predictions_reg_f)

plot_spaghetti(predictions_reg_f,
               plot_cols = grep("^Score_gold_sum", colnames(predictions_reg_f)),
               color_col = "Score_gold_sum_1",
               main = "Predicted Regression",
               file = "plots/plot_spaghetti_gold_predicted_regression_con_std.pdf")

####--------------------- Status tests without z-standardization -----------------------####

load("Data_preparation/Daten_cody_gold_sub_wide.RData")
head(Daten_wide)
summary(Daten_wide[, grep("^Score_gold_sum", colnames(Daten_wide))])
psych::describe(Daten_wide[, grep("^Score_gold_sum", colnames(Daten_wide))])

plot_spaghetti(Daten_wide,
               plot_cols = grep("^Score_gold_sum", colnames(Daten_wide)),
               color_col = "Score_gold_sum_1",
               main = "Observed Unstandardized",
               file = "plots/plot_spaghetti_gold_observed_unstandardized.pdf")

gold_long <- data.frame("id" = rep(Daten_wide$Id, G),
                        "status_score" = unlist(Daten_wide[, grep("^Score_gold_sum", colnames(Daten_wide))]),
                        "timepoint" = rep(1:G, each = nrow(Daten_wide)))
gold_long$timepoint <- factor(gold_long$timepoint)
head(gold_long)
head(Daten_wide[, c("Id", grep("^Score_gold_sum", colnames(Daten_wide), value  = TRUE))])
library(ggplot2)

plot_observed_unstandardized <- ggplot(data = gold_long, aes(y = status_score, x=timepoint, fill = qualitative_hcl(1))) +
  geom_violin(show.legend = FALSE, position = position_dodge(0.9)) +
  geom_boxplot(width = 0.2, alpha = 0.2, position = position_dodge(0.9), show.legend = FALSE) +
  scale_fill_manual(values=qualitative_hcl(3)) +
  labs(y="Score", x="Repetition") +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11),
        title=element_text(size = 11)) +
  ggtitle(label = "Observed Unstandardized")
ggsave("plots/density_gold_observed_unstandardized.pdf", plot_observed_unstandardized, width = 7, height = 6)

# GAM by grade
# concurrent standardization
Daten_wide[, paste0("con_std_", grep("Score_gold_sum", colnames(Daten_wide), value = TRUE))] <-
  (Daten_wide[, grep("Score_gold_sum", colnames(Daten_wide))] -
        mean(do.call(c, Daten_wide[, grep("Score_gold_sum", colnames(Daten_wide))]), na.rm = TRUE)) /
        sd(do.call(c, Daten_wide[, grep("Score_gold_sum", colnames(Daten_wide))]), na.rm = TRUE)
psych::describe(Daten_wide[, grep("con_std_Score_gold_sum", colnames(Daten_wide))])
psych::describe(do.call(c, Daten_wide[, grep("con_std_Score_gold_sum", colnames(Daten_wide))]))


pred_con_std_grade2 <- obtain_gam(Daten_wide[Daten_wide$Grade == 2,],
                          plot_cols = grep("con_std_Score_gold_sum", colnames(Daten_wide)))
pred_con_std_grade3 <- obtain_gam(Daten_wide[Daten_wide$Grade == 3,],
                          plot_cols = grep("con_std_Score_gold_sum", colnames(Daten_wide)))
pred_con_std_grade4 <- obtain_gam(Daten_wide[Daten_wide$Grade == 4,],
                          plot_cols = grep("con_std_Score_gold_sum", colnames(Daten_wide)))

cols <- RColorBrewer::brewer.pal(3, "Dark2")
dat_plot <- Daten_wide[, grep("con_std_Score_gold_sum", colnames(Daten_wide))]

pdf(file = "plots/plot_spaghetti_gam_observed_concurrent_standardization.pdf", width = 7, height = 6)
par(mar = c(5,4,4,5), xpd = TRUE)
plot(1:ncol(dat_plot), rep(0, ncol(dat_plot)),
     ylim = range(dat_plot, na.rm = TRUE),
     ylab = "Score", xlab = "Timepoint", main = "Observed Trajectories", type = "n")
for(p in 1:nrow(dat_plot)) lines(1:ncol(dat_plot), dat_plot[p,], col = "lightgrey")
lines(1:G, pred_con_std_grade2$pred, type = "l", col = cols[1])
polygon(x = c(1:G, G:1),
        y = c(pred_con_std_grade2$pred - pred_con_std_grade2$disp, pred_con_std_grade2$pred + pred_con_std_grade2$disp),
        col = adjustcolor(cols[1], alpha.f = .2), border = NA)
lines(1:G, pred_con_std_grade3$pred, col = cols[2])
polygon(x = c(1:G, G:1),
        y = c(pred_con_std_grade3$pred - pred_con_std_grade3$disp, pred_con_std_grade3$pred + pred_con_std_grade3$disp),
        col = adjustcolor(cols[2], alpha.f = .2), border = NA)
lines(1:G, pred_con_std_grade4$pred, col = cols[3])
polygon(x = c(1:G, G:1),
        y = c(pred_con_std_grade4$pred - pred_con_std_grade4$disp, pred_con_std_grade4$pred + pred_con_std_grade4$disp),
        col = adjustcolor(cols[3], alpha.f = .2), border = NA)
abline(h = 0 #mean(Daten_wide[, "con_std_Score_gold_sum_1"], na.rm = TRUE) # mean at T1
       , col = adjustcolor("black", alpha.f = .8), xpd = FALSE, lty = "dashed")
legend(x = 6.3, y = 3,
       legend = 4:2, title = "Grade",
       fill = rev(cols),
       bty = "n")
dev.off()
