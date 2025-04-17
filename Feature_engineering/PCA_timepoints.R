####-------------------- Principal Components Analyses on Training Data ----------------####

library(psych)

load("Data_preparation/Daten_cody_gold_sub_wide.RData")

dir_manuscript <- "~/Dokumente/FAIR/Reha/paper/MRFs CODY/"

# predictor variables
G <- 6 #number of gold coin tests

# list with entries for each status test / timepoint
predictors_training <- vector("list", length = G)
# last training day before each status test
days_gold <- seq(4, 29, by = 5)
for (g in 1:G) {
  predictors_training[[g]] <- grep(paste0("T", 1:days_gold[g], collapse = "$|"), colnames(Daten_wide), value = TRUE) #better: end all with $

}

# separate score, Diskrepanz, and level_dif
predictors_pca <- lapply(predictors_training, function(pr) {
  list("score" = grep("Score_G", pr, value = TRUE),
       "leveldif" = grep("Level_dif_", pr, value = TRUE),
       "discrep" = grep("Diskrepanz_", pr, value = TRUE),
       "time" = grep("D_Bearbeitungszeit_", pr, value = TRUE))
  })


# PCA: number of components is either predictors/2 or those with sd = 0.1*sd(first component)
# conduct also parallel analyses
set.seed(1705)

pcas <- lapply(predictors_pca, function(pred_t, dat)
  lapply(pred_t, function(pred_var, dat) {
    if(length(pred_var) > 0) {
      dt <- na.omit(dat[, c("Id", pred_var)])
      ids <- as.numeric(dt$Id)
      dt <- dt[, grep("Id", colnames(dt), invert = TRUE)]
      list(
        "pca" = prcomp(dt, center = TRUE, scale. = TRUE, tol = 0.1, rank. = length(pred_var)/2),
        "parallel" = fa.parallel(dt, n.iter = 100, fa = "pc"),
        "ids" = ids
        )
    }
  }, dat = dat), dat = Daten_wide)

# print pca summaries, predictors/2, #components extracted, #components parallel
lapply(pcas, function(pca_t) lapply(pca_t, function(pca_t) {
  print(summary(pca_t$pca))
  prop <- round(cumsum(pca_t$pca$sdev^2 / sum(pca_t$pca$sdev^2)), 3)
  n_pca <- ncol(pca_t$pca$rotation)
  n_parallel <- pca_t$parallel$ncomp
  c(nrow(pca_t$pca$rotation)/2, n_pca, n_parallel, prop[n_pca], prop[n_parallel])
} ))
# go with n_parallel
# time needs much more components

# extract predictions for the first n_parallel components
pca_scores <- lapply(pcas, function(pca_t) lapply(pca_t, function(pca_t) {
  n_parallel <- pca_t$parallel$ncomp
  scores <- data.frame("Id" = pca_t$ids, pca_t$pca$x[, 1:n_parallel])
} ))

# combine scores, leveldif and discrep
pca_scores <- lapply(pca_scores, function(pca_st) {
  nm <- names(pca_st)
  for (v in nm) colnames(pca_st[[v]]) <- c("Id", paste0(v, "_", colnames(pca_st[[v]][-c(1)])))
  scores_merged <- merge(pca_st[[1]], pca_st[[2]], by = "Id", all = TRUE)
  for (v in 3:length(pca_st)) {
    scores_merged <- merge(scores_merged, pca_st[[v]], by = "Id", all = TRUE)
  }
  return(scores_merged)
})

str(pca_scores)
# manually assign name to the case where only one component was extracted
colnames(pca_scores[[1]])[4] <- "leveldif_PC1"

save(pca_scores, file = "Feature_engineering/pca_scores.RData")
save(pcas, file = "Feature_engineering/PCA_results.RData")

# describe selected PCA results

# are the PCs positively correlated with the gold coin tests?
pca_scores_gold <- vector("list", G)
names_gold <- grep("Score_gold_sum", colnames(Daten_wide), value = TRUE)
for(g in 1:G) pca_scores_gold[[g]] <- merge(pca_scores[[g]],
                                            Daten_wide[, c("Id", names_gold)],
                                            by = "Id", all.x = TRUE, all.y = FALSE)
cors <- lapply(pca_scores_gold, cor, use = "pairwise")
cors <- lapply(cors, function(cr, ng) cr[ng, grep("score_PC1|leveldif_PC1|leveldif_PC2|leveldif_PC3|time_PC1$", colnames(cr))],
               ng = names_gold)
cors
# score_PC1: always positive loadings
# leveldif_PC1: always positive
# leveldif_PC2: negative for T3
# time_PC1: positive for T4 and T6

# combine for each timepoint, across variable types
pca_loads <- lapply(pcas, function(pca) {
  names_time <- names(pca)
  do.call(cbind, lapply(1:length(pca), function(v, pc, nm) {
    df <- pc[[v]]$pca$rotation[, 1:pc[[v]]$parallel$ncomp, drop = FALSE]
    colnames(df) <- paste(nm[[v]], colnames(df), sep = "_")
    df
  }, pc = pca, nm = names_time))
    })

round(pca_loads[[1]], 3)
round(pca_loads[[4]][order(rownames(pca_loads[[4]])), c("score_PC1","leveldif_PC2")], 3)
round(pca_loads[[4]][order(pca_loads[[4]][, "score_PC1"]), c("score_PC1","leveldif_PC2")], 3)
round(pca_loads[[4]][order(pca_loads[[4]][, "leveldif_PC2"]), c("score_PC1","leveldif_PC2")], 3)

load("Data_preparation/A.RData")
A <- A[A$Trainingstag < 32, ]
A_ord <- A[order(A$uebung), ]
A_ord$rep <- do.call(c, apply(table(A_ord$x)[unique(A_ord$x)], 1, function(count) 1:count))
A_ord

# try to add Bereich
load("Data_preparation/Daten_cody_gold_sub.RData")
area <- Daten_cody_gold[!duplicated(Daten_cody_gold$Uebung), c("Uebung", "Bereich", "Faehigkeit")]
area <- area[grep("CODY|Gold", area$Uebung, invert = TRUE),]
area <- area[is.na(area$Uebung) == FALSE,]

A_ord <- merge(A_ord, area, by.x = "uebung", by.y = "Uebung", all.x = TRUE)

devtools::load_all("~/Dokumente/packages/DataAnalysisSimulation")

# for figure in revision
A_ord_rev <- A_ord
A_ord_rev$area <- recode.df(A_ord_rev$Bereich,
                        c("Zahlen-Größen-Verknüpfung", "Faktenwissen und Rechnen",  "Teil-Ganzes-Verständnis",
                          "Arbeitsgedächtnis", "Dezimalsystem"),
                        c("number-size-connection", "facts and calculating", "part-whole-understanding",
                         "working memory", "decimal system"))
A_ord_rev$area_short <- recode.df(A_ord_rev$area,
                                  c("number-size-connection", "facts and calculating", "part-whole-understanding",
                                    "working memory", "decimal system"),
                                  c("NSC", "FKC", "WPU",
                                    "WM", "DS"))
A_ord_rev <- A_ord_rev[order(A_ord_rev$Trainingstag), c("Trainingstag", "uebung", "area_short")]
A_ord_rev$uebung_num <- as.numeric(gsub("^G|^G0", "", A_ord_rev$uebung))

table_a <- matrix(NA, 5, 30)
rownames(table_a) <- c("WM", "DS", "FKC", "WPU", "NSC")
colnames(table_a) <- 1:30
table_a

for(i in 1:(nrow(A_ord_rev) - 2)) {
  rw <- A_ord_rev[i,]
  if(is.na(table_a[which(rownames(table_a) == rw$area_short), which(colnames(table_a) == rw$Trainingstag)])) {
    table_a[which(rownames(table_a) == rw$area_short), which(colnames(table_a) == rw$Trainingstag)] <- rw$uebung_num
  } else {
    table_a[which(rownames(table_a) == rw$area_short), which(colnames(table_a) == rw$Trainingstag)] <- paste(
      table_a[which(rownames(table_a) == rw$area_short), which(colnames(table_a) == rw$Trainingstag)],
      rw$uebung_num, sep = ",")

  }
}
table_a
write.csv(table_a, file = "Feature_engineering/table_revision.csv")

prep_loads <- function(loads, cols = c("score_PC1|score_PC2|leveldif_PC1|leveldif_PC2|time_PC1$"), A_o = A_ord, ar = area)  {
  #if(! (cols[2] %in% colnames(loads))) cols[2] <- "leveldif_PC1"
  pca_loads_t <- data.frame(round(loads[, grep(cols, colnames(loads))], 2))
  pca_loads_t$score_rep <- recode.df(rownames(pca_loads_t),
                                      paste0("Score_", A_o$uebung, "_T", A_o$Trainingstag),
                                      paste0(A_o$uebung, "_rep", A_o$rep))
  uebungen_t <- do.call(c, lapply(strsplit(pca_loads_t$score_rep, split = "_"), function(l) l[1]))
  pca_loads_t$area <- ar$Bereich[match(uebungen_t, ar$Uebung)]
  pca_loads_t$ability <- ar$Faehigkeit[match(uebungen_t, ar$Uebung)]
  pca_loads_t
}

pca_loads_timepoints <- lapply(pca_loads, prep_loads)
lapply(pca_loads_timepoints, function(pl) pl[order(pl$score_PC1), -c(2:5)])
lapply(pca_loads_timepoints, function(pl) pl[order(pl$score_PC2), -c(1,3:5)])
lapply(pca_loads_timepoints, function(pl) pl[order(pl$leveldif_PC1), -c(1,2,4:5)])
lapply(pca_loads_timepoints[1], function(pl) pl[order(pl$leveldif_PC1), ])
lapply(pca_loads_timepoints[2:6], function(pl) pl[order(pl$leveldif_PC2), -c(1:3,5)]) # reversed for T3!
# higher loadings for 2nd and 3rd, later also 4,5,6th of G01 and G03
# T5 and T6: negative loadings for later repetitions of G12 and G08

loads4 <- prep_loads(pca_loads[[4]], cols = c("score_PC1|leveldif_PC3")) # reversed for T3!
loads4[order(loads4$leveldif_PC3),]

# time is not that clear
lapply(pca_loads_timepoints, function(pl) pl[order(pl$time_PC1), -c(1:4)])

prep_loads_paper <- function(loads) {
  loads$area <- recode.df(loads$area,
                          c("Zahlen-Größen-Verknüpfung", "Faktenwissen und Rechnen",  "Teil-Ganzes-Verständnis",
                            "Arbeitsgedächtnis", "Dezimalsystem"),
                          c("NSC", "FKC", "WPU", "WM", "DS"))
  # WM = working memory, DS = Decimal
  # System, FKC = Factual Knowledge and Calculating, WPU = whole-part-understanding,
  # NSC = number-size-connection.
  loads$trial <- substr(loads$score_rep, 8, 8)
  loads$score_rep <- substr(loads$score_rep, 1, 3)
  loads

  loads[, grep("PC", colnames(loads))] <- matrix(paste(
    paste0("\\cellcolor{",
           cut(as.matrix(loads[, grep("PC", colnames(loads))]),
               breaks = c(-1, -.5, -.3, -.1, .1, .3, .5, 1),
               labels = c("blue!40", "blue!25", "blue!10", "white", "red!10", "red!25", "red!40"),
               include.lowest = TRUE),
           "}"),
    format(unlist(loads[, grep("PC", colnames(loads))]), nsmall = 2)
  ), nrow = nrow(loads))

  # add column day
  loads$day <- rep(1:(nrow(loads)/2), each = 2)

  # re-order columns
  if("leveldif_PC2" %in% colnames(loads)) {
    loads <- loads[, c("day", "area", "score_rep", "trial", "score_PC1", "score_PC2", "leveldif_PC1", "leveldif_PC2", "time_PC1")]
    colnames(loads) <- c("Day", "Area", "Task", "Trial", "PC1 level", "PC2 level", "PC1 level difference", "PC2 level difference", "PC1 time")
  } else {
    loads <- loads[, c("day", "area", "score_rep", "trial", "score_PC1", "score_PC2", "leveldif_PC1", "time_PC1")]
    colnames(loads) <- c("Day", "Area", "Task", "Trial", "PC1 level", "PC2 level", "PC1 level difference", "PC1 time")
  }
  loads

}

loads_paper <- lapply(pca_loads_timepoints, prep_loads_paper)
loads_paper

# textables
header_T1 <- list()
header_T1$pos <- list(-1)
header_T1$command <- c("\\hline \n Day & Area & Task & Trial & \\multicolumn{2}{c}{Level} & \\multicolumn{1}{c}{Level Difference} & Time \\\\
                    \\hline \\multicolumn{4}{c}{} & PC1 & PC2 & PC1 & PC1 \\\\",
                          "\\hline \\multicolumn{8}{l}{\\small \\textit{Note.} WM = working memory, DS = Decimal System,}\\\\
                    \\multicolumn{8}{l}{\\small FKC = Factual Knowledge and Calculating, WPU = whole-part-understanding,}\\\\
                    \\multicolumn{8}{l}{\\small NSC = number-size-connection.} \\\\")
header_short <- list()
header_short$pos <- list(-1)
header_short$command <- c("\\hline \n Day & Area & Task & Trial & \\multicolumn{2}{c}{Level} & \\multicolumn{2}{c}{Level Difference} & Time \\\\
                    \\hline \\multicolumn{4}{c}{} & PC1 & PC2 & PC1 & PC2 & PC1 \\\\",
                    "\\hline \\multicolumn{9}{l}{\\small \\textit{Note.} WM = working memory, DS = Decimal System,}\\\\
                    \\multicolumn{9}{l}{\\small FKC = Factual Knowledge and Calculating, WPU = whole-part-understanding,}\\\\
                    \\multicolumn{9}{l}{\\small NSC = number-size-connection.} \\\\")
header_long <- list()
header_long$pos <- list(-1)
header_long$command <- c("\\hline \n Day & Area & Task & Trial & \\multicolumn{2}{c}{Level} & \\multicolumn{2}{c}{Level Difference} & Time \\\\
                    \\hline \\multicolumn{4}{c}{} & PC1 & PC2 & PC1 & PC2 & PC1 \\\\ \\hline \\endfirsthead
                    {{\\bfseries \\tablename \\hspace*{0.7pt} \\thetable{} -- continued}} \\\\
                    \\hline \n Day & Area & Task & Trial & \\multicolumn{2}{c}{Level} & \\multicolumn{2}{c}{Level Difference} & Time \\\\
                    \\hline \\multicolumn{4}{c}{} & PC1 & PC2 & PC1 & PC2 & PC1 \\\\ \\hline \\endhead
                    {\\small -- continued --} \\endfoot
                    \\multicolumn{9}{l}{\\small \\textit{Note.} WM = working memory, DS = Decimal System,}\\\\
                    \\multicolumn{9}{l}{\\small FKC = Factual Knowledge and Calculating, WPU = whole-part-understanding,}\\\\
                    \\multicolumn{9}{l}{\\small NSC = number-size-connection.}\\endlastfoot", "")

for(i in 1:6) {
  if(i == 1) {
    header <- header_T1
    fl <- TRUE
    tab <- "tabular"
  } else if(i < 5) {
    header <- header_short
    fl <- TRUE
    tab <- "tabular"
  } else {
    header <- header_long
    fl <- FALSE
    tab <- "longtable"
  }
  header$pos <- list(-1, nrow(loads_paper[[i]]))
  print(xtable::xtable(loads_paper[[i]],
                       digits = c(rep(0, 5), rep(2, ifelse(i == 1, 4, 5))),
                       align = rep("r", ncol(loads_paper[[i]]) + 1),
                       caption = paste0("Loadings for the Principal Components extracted from the Predictors Until T", i),
                       label = paste0("tb:PCA_loadings_T", i)),
        include.colnames = FALSE, include.rownames = FALSE,
        hline.after=c(0),
        sanitize.rownames.function=function(x){x},
        sanitize.colnames.function = function(x){x},
        sanitize.text.function = function(x){x},
        NA.string = "", table.placement = "htp",
        add.to.row = header,
        caption.placement = "top",
        tabular.environment = tab, floating = fl,
        file = paste0(dir_manuscript, "tables/textable_pca_loads_T", i, ".tex"))
  file.copy(from = paste0(dir_manuscript, "tables/textable_pca_loads_T", i, ".tex"),
            to = paste0("~/Dokumente/FAIR/Reha/paper/SOM/tables/textable_pca_loads_T", i, ".tex"),
            overwrite = TRUE)
}
# note: T5 and T6 were modified by hand: put into spacing environment
