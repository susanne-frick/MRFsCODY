####----------------- correlations between outcomes and predictors ----------------------####

# for simplicity for the whole dataset, whithout splitting

library(corrplot)
devtools::load_all("~/Dokumente/packages/DataAnalysisSimulation")

# load data
load("Data_preparation/Daten_cody_gold_sub_wide.RData")
head(Daten_wide)

load("Feature_engineering/pca_scores.RData")

# load fit lists
fit_list_short <- readRDS("Data_analysis/fit_list_short.rds")

####--------------------- define variables for model types ----------------------####

# predictor variables
G <- 6 #number of gold coin tests

# list with entries for each status test
predictors_basic <- predictors_training <- vector("list", length = G)
for (g in 1:G) {
  predictors_basic[[g]] <- c("Score_cody_1", "Grade", "Grade_cody_1")
  predictors_training[[g]] <- c(predictors_basic[[g]],
                                colnames(pca_scores[[g]])[-c(1)]) # without Id
}

outcomes <- paste0("Score_gold_sum_", 1:G)

cor_data <- vector("list", G)
for(g in 1:G) {
  Daten_wide_pca <- merge(Daten_wide, pca_scores[[g]], by = "Id", all = TRUE)
  data_complete <- na.omit(Daten_wide_pca[, c("Id", outcomes, predictors_training[[g]])])

  cor_data[[g]] <- cor(data_complete[, predictors_training[[g]]], data_complete[, outcomes])
  colnames(cor_data[[g]]) <- paste0("T", 1:G)

  rownames(cor_data[[g]]) <- gsub("leveldif_", "Level Diff. ", rownames(cor_data[[g]]))
  rownames(cor_data[[g]]) <- gsub("score_", "Level ", rownames(cor_data[[g]]))
  rownames(cor_data[[g]]) <- gsub("discrep_", "Discrepancy  ", rownames(cor_data[[g]]))
  rownames(cor_data[[g]]) <- gsub("time_", "Time ", rownames(cor_data[[g]]))
  rownames(cor_data[[g]]) <- gsub("Grade_cody_1", "Pretest x Grade", rownames(cor_data[[g]]))
  rownames(cor_data[[g]]) <- gsub("Score_cody_1", "Pretest", rownames(cor_data[[g]]))
}

cor_data

pdf(file = "plots/corplots.pdf", width = 7, height = 21)
for (g in 1:G) {
  corrplot(cor_data[[g]], method = "square", title = paste0("Training Data Until T", g),
           addCoef.col = ifelse(abs(cor_data[[g]]) > .1, yes = "black", no = "transparent"),
           col = colorRampPalette((c("#880C19", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                                     "#FFFFFF","#D1E5F0", "#92C5DE", "#43A7C3", "#217CAC", "#145F89")) )(200),
           tl.col = "black", tl.srt = 0, tl.offset = 1.2, tl.cex = 1, cl.pos = "n", number.cex = 1, mar = c(0,0,2,0))
}
dev.off()

# separately to be put into figure environments in SOM
for (g in 1:G) {
  pdf(file = paste0("plots/corplot_T", g, ".pdf"), width = 7, height = 0.7*nrow(cor_data[[g]]))
  corrplot(cor_data[[g]], method = "square",
           addCoef.col = ifelse(abs(cor_data[[g]]) > .1, yes = "black", no = "transparent"),
           col = colorRampPalette((c("#880C19", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                                     "#FFFFFF","#D1E5F0", "#92C5DE", "#43A7C3", "#217CAC", "#145F89")) )(200),
           tl.col = "black", tl.srt = 0, tl.offset = 1.2, tl.cex = 1, cl.pos = "n", number.cex = 1, mar = c(0,0,2,0))
  dev.off()
}
