####------------------ descriptives on levels -----------------------------####

library(colorspace)

recode.df <- function(variable, levels.old, levels.new) {
  levels.new[match(as.character(variable), as.character(levels.old))]
}

# within the tasks, the levels differ in whether they include addition, subtraction, and which numbers (<= 20, 100, 1000)

# load data
load("Data_preparation/Daten_cody_gold_sub_wide.RData")
head(Daten_wide)

# leveldif_PC2
# has fairly consistent high loadings for G01 (Fischfreunde) and G03 (Drachenwaage 1)
# positive: rep2 and rep3
# negative: rep1

apply(Daten_wide[, grep("Score_G01", colnames(Daten_wide))], 2, table)
apply(Daten_wide[, grep("Score_G03", colnames(Daten_wide))], 2, table)
# rather wide distributions, G01 is getting smaller (see also figure)

# G01 (Fischfreunde)
# according to MCT_Levels/01_Fischfreunde [24-04-2023].xlsx
# from Level 5: Zahlenraum 100, from Level 15: Zahlenraum 1000

apply(Daten_wide[, grep("Score_G01", colnames(Daten_wide))], 2, function(r) {
  round(c(sum(r < 5, na.rm = TRUE), sum((r >= 5) & (r < 15), na.rm = TRUE))/length(r), 2)
  })

# first calculation (sum below 10) at level 3
apply(Daten_wide[, grep("Score_G01", colnames(Daten_wide))], 2, function(r) {
  round(c(sum(r < 3, na.rm = TRUE))/length(r), 2)
})

# first subtraction at level 7
apply(Daten_wide[, grep("Score_G01", colnames(Daten_wide))], 2, function(r) {
  round(c(sum(r < 7, na.rm = TRUE))/length(r), 2)
})

# first addition over 10 level 10
apply(Daten_wide[, grep("Score_G01", colnames(Daten_wide))], 2, function(r) {
  round(c(sum(r < 10, na.rm = TRUE))/length(r), 2)
})

apply(Daten_wide[, grep("Level_dif_G01", colnames(Daten_wide))], 2, table)
# G01 has (also according to Sarah Chromik) the largest proportion of negative improvements

# summary (across children) for each timepoint
apply(Daten_wide[, grep("Score_G01", colnames(Daten_wide))], 2, summary, na.rm = TRUE)


# G03 (Drachenwaage I) ####
# according to MCT_Levels/03_Drachenwaage I [24-04-2023].xlsx
# from Level 8: Zahlenraum 100, from Level 21: Zahlenraum 1000

apply(Daten_wide[, grep("Score_G03", colnames(Daten_wide))], 2, function(r) {
  round(c(sum(r < 8, na.rm = TRUE), sum((r >= 8) & (r < 21), na.rm = TRUE))/length(r), 2)
})

# addition from level 3
apply(Daten_wide[, grep("Score_G03", colnames(Daten_wide))], 2, function(r) {
  round(c(sum(r > 3, na.rm = TRUE))/length(r), 2)
})
# almost everyone

# subtraction, over 10 from level 6
apply(Daten_wide[, grep("Score_G03", colnames(Daten_wide))], 2, function(r) {
  round(c(sum(r > 6, na.rm = TRUE))/length(r), 2)
})

# summary (across children) for each timepoint
apply(Daten_wide[, grep("Score_G03", colnames(Daten_wide))], 2, summary, na.rm = TRUE)


# summary (across children) for each timepoint
apply(Daten_wide[, grep("Score_G08", colnames(Daten_wide))], 2, summary, na.rm = TRUE)

####-------------------- Spaghetti plots ------------------------####

# Spaghetti plot
plot_spaghetti <- function(task, main, file, levels, labels, dat = Daten_wide,
                           task_cols  = NULL, ylab = "Level") {
  dat_task <- dat[, grep(task, colnames(dat))]
  cols <- sequential_hcl(max(dat_task[, 1], na.rm = TRUE), palette = "Terrain", rev = TRUE)
  cols_task <- recode.df(dat_task[, 1], sort(unique(dat_task[, 1])), cols)

  if(isFALSE(is.null(task_cols))) {
    dat_cols <- dat[, grep(task_cols, colnames(dat))]
    cols <- sequential_hcl(max(dat_cols[, 1], na.rm = TRUE), palette = "Terrain", rev = TRUE)
    cols_task <- recode.df(dat_cols[, 1], sort(unique(dat_cols[, 1])), cols)
  }

  pdf(file = file, width = 7, height = 6)
  par(mar = c(5,4,4,7)) #default c(5,4,4,2) + 0.1
  plot(1:ncol(dat_task), rep(0, ncol(dat_task)), ylim = range(dat_task, na.rm = TRUE),
       ylab = ylab, xlab = "Repetition", main = main, type = "n")
  for(p in 1:nrow(dat_task)) lines(1:ncol(dat_task), dat_task[p,], col = cols_task[p])
  abline(h = levels, col = "darkgrey")
  axis(side = 4, at = levels, labels = labels,
       las = 1, cex.axis = .8)
  dev.off()
}


plot_spaghetti(task = "Score_G03", main = "G03 Drachenwaage",
               file = "plots/plot_spaghetti_G03.pdf",
               levels = c(6, 8, 21),
               labels = c("subtraction,\nover 10", "number space\n100", "number space\n1000"))

plot_spaghetti(task = "Score_G01", main = "G01 Fischfreunde",
               file = "plots/plot_spaghetti_G01.pdf",
               levels =c(5,7,10),
               labels = c("number space\n100", "subtraction", "addition over\n10"))

plot_spaghetti(task = "Score_G08", main = "G08 Zukunftsblick",
               file = "plots/plot_spaghetti_G08.pdf",
               levels = c(9),
               labels = c("area 1-10"))

plot_spaghetti(task = "Score_G12", main = "G12 Pandawald",
               file = "plots/plot_spaghetti_G12.pdf",
               levels = c(3,5,7),
               labels = c("mixed", "number space 18", "mixed"))

# leveldifs

plot_spaghetti(task = "Level_dif_G01", main = "G01 Fischfreunde",
               task_cols = "Score_G01",
               file = "plots/plot_spaghetti_leveldif_G01.pdf",
               ylab = "Level Difference",
               levels =c(0),
               labels = c(""))

plot_spaghetti(task = "Level_dif_G03", main = "G03 Drachenwaage",
               file = "plots/plot_spaghetti_leveldif_G03.pdf",
               ylab = "Level Difference",
               levels =c(0),
               labels = c(""))
plot_spaghetti(task = "Level_dif_G03", main = "G03 Drachenwaage",
               task_cols = "Score_G03",
               file = "plots/plot_spaghetti_leveldif_G03_startlevel.pdf",
               ylab = "Level Difference",
               levels =c(0),
               labels = c(""))


plot_spaghetti(task = "Level_dif_G08", main = "G08 Zukunftsblick",
               task_cols = "Score_G08",
               file = "plots/plot_spaghetti_leveldif_G08_startlevel.pdf",
               ylab = "Level Difference",
               levels =c(0),
               labels = c(""))

plot_spaghetti(task = "Level_dif_G12", main = "G12 Pandawald",
               task_cols = "Score_G12",
               file = "plots/plot_spaghetti_leveldif_G12_startlevel.pdf",
               ylab = "Level Difference",
               levels =c(0),
               labels = c(""))
