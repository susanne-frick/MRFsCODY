####------------------------ analyze selections of lasso regressions ----------------------------####

# setwd("~/Dokumente/FAIR/Reha")
library(glmnet)

res <- readRDS("Data_analysis/results_reg.rds")

# re-build design list
G <- 6
fit_list_short <- rbind(expand.grid("timepoint" = 1:G,
                                    "outcome" = c("multi", "uni", "row", "sub", "add"),
                                    "predictors" = c("basic", "training", "status")),
                        expand.grid("timepoint" = 1:G,
                                    "outcome" = "multi_long",
                                    "predictors" = c("basic", "training")) # multi_long without status predictors
)

outer_folds <- 10
fit_list <- rbind(expand.grid("timepoint" = 1:G,
                              "outcome" = c("multi", "uni", "row", "sub", "add"),
                              "predictors" = c("basic", "training", "status"),
                              "fold" = 1:outer_folds),
                  expand.grid("timepoint" = 1:G,
                              "outcome" = "multi_long",
                              "predictors" = c("basic", "training"),
                              "fold" = 1:outer_folds)
)


# outcome variables
gold_uni <- lapply(1:G, function(g) paste0("Score_gold_sum_", g))
gold_add <- lapply(1:G, function(g) paste0("score_add_", g))
gold_sub <- lapply(1:G, function(g) paste0("score_sub_", g))
gold_row <- lapply(1:G, function(g) paste0("score_row_", g))
gold_multi <- lapply(1:G, function(g) paste0(c("score_row_", "score_sub_", "score_add_"), g))
gold_uni_long <- lapply(1:G, function(g) paste0("Score_gold_sum_", 1:G))
#list, all others vector
outcomes <- list("multi" = gold_multi, "uni" = gold_uni, "add" = gold_add, "sub" = gold_sub, "row" = gold_row,
                 "multi_long" = gold_uni_long)

####----------------------------------------------------------------------------------------------------------####
coef_list <- n_sel <- vector("list", nrow(fit_list_short))

for(f in 1:nrow(fit_list_short)) {

  if((fit_list_short[f, "timepoint"] == 1) & (fit_list_short[f, "predictors"] == "status")) {
    next
  } else {

    fs <- which((fit_list$timepoint == fit_list_short$timepoint[f]) &
                  (fit_list$outcome == fit_list_short$outcome[f]) &
                  fit_list$predictors == fit_list_short$predictors[f])
    outcomes_f <- outcomes[[fit_list_short[[f, "outcome"]]]][[fit_list_short[f, "timepoint"]]]

    if(fit_list_short$outcome[f] %in% c("multi", "multi_long")) {
      coefs <- lapply(res[fs], function(r) do.call(cbind, lapply(coef(r$fit), as.matrix)))
    } else {
      coefs <- lapply(res[fs], function(r) as.matrix(coef(r$fit)))
    }
    # without intercept
    coefs <- lapply(coefs, function(cf) cf[-c(1),,drop = FALSE])
    coef_list[[f]] <- array(do.call(c, coefs), dim = c(dim(coefs[[1]]), length(fs)))

    # count how often a variable was selected
    n_sel[[f]] <- apply(coef_list[[f]], c(1,2), function(cf) sum(cf != 0))
    if(fit_list_short$outcome[f] %in% c("multi", "multi_long")) {
      # selections are identical across multivariate outcomes
      n_sel[[f]] <- n_sel[[f]][,1,drop = FALSE]
      dimnames(n_sel[[f]]) <- list(rownames(as.matrix(res[[fs[1]]]$fit$beta[[1]])), fit_list_short$outcome[f])
      dimnames(coef_list[[f]]) <- list(rownames(as.matrix(res[[fs[1]]]$fit$beta[[1]])), outcomes_f, paste0("fold", 1:10))
    } else {
      dimnames(n_sel[[f]]) <- list(rownames(res[[fs[1]]]$fit$beta), fit_list_short$outcome[f])
      dimnames(coef_list[[f]]) <- list(rownames(res[[fs[1]]]$fit$beta), outcomes_f, paste0("fold", 1:10))
    }

  }
}

get_f <- function(outcome, predictors, timepoint) {
  which(fit_list_short$outcome == outcome &
          fit_list_short$predictors == predictors &
          fit_list_short$timepoint == timepoint)
}

rm_0 <- function(imp, cut = 0) round(imp[rowSums(round(imp[, -c(ncol(imp)), drop = FALSE], 2) <= cut) < (ncol(imp) - 1),], 2)

# examine selected results: basic, T1, multi vs. uni
n_sel[[get_f("multi_long", "basic", 1)]]
n_sel[[get_f("uni", "basic", 1)]]

# training T6
cbind(n_sel[[get_f("multi_long", "training", 6)]],
      n_sel[[get_f("uni", "training", 6)]])

# estimated coefficients
apply(coef_list[[get_f("multi_long", "basic", 1)]], c(1,2), function(cf) mean(cf[cf != 0]))
apply(coef_list[[get_f("multi_long", "basic", 1)]], c(1,2), function(cf) summary(cf[cf != 0]))

imp_multi_T6 <- rm_0(cbind(round(apply(coef_list[[get_f("multi_long", "training", 6)]], c(1,2), function(cf) mean(cf[cf != 0])), 2),
                                  n_sel[[get_f("multi_long", "training", 6)]])[n_sel[[get_f("multi_long", "training", 6)]] > 0,], cut = .02)
imp_multi_T3 <- rm_0(cbind(round(apply(coef_list[[get_f("multi_long", "training", 3)]], c(1,2), function(cf) mean(cf[cf != 0])), 2),
                           n_sel[[get_f("multi_long", "training", 3)]])[n_sel[[get_f("multi_long", "training", 3)]] > 0,], cut = .02)

imp_uni_status_T6 <- rm_0(cbind(round(apply(coef_list[[get_f("uni", "status", 6)]], c(1,2), function(cf) mean(cf[cf != 0])), 2),
                                n_sel[[get_f("uni", "status", 6)]])[n_sel[[get_f("uni", "status", 6)]] > 0, ], cut = .02)


imp_multi_T6 <- imp_multi_T6[order(rownames(imp_multi_T6)), ]
imp_multi_T3 <- imp_multi_T3[order(rownames(imp_multi_T3)), ]
imp_uni_T6 <- imp_uni_T6[order(rownames(imp_uni_T6)), ]

colnames(imp_multi_status_T6) <- c("Rows", "Subtraction", "Addition", " ")
colnames(imp_uni_status_T6) <- c("Mean $\\beta$", "N")

rownames(imp_multi_status_T6) <- gsub("_", " ", rownames(imp_multi_status_T6))
rownames(imp_uni_status_T6) <- gsub("_", " ", rownames(imp_uni_status_T6))

header <- list()
header$pos <- list(-1)
header$command <- c("\\hline \n & \\multicolumn{3}{c}{Mean $\\beta$} & N \\\\ \\cmidrule{2-4}")

print(xtable::xtable(imp_multi_status_T6, digits=c(0,rep(2,3),0)), include.colnames = T, include.rownames=T,
      hline.after=c(0, nrow(imp_multi_status_T6)),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp", add.to.row = header,
      caption.placement = "top", latex.environments = NULL,
      file="presentations/Update_Meeting_March2023/tables/textable_importance_reg_multi-status-T6.tex")

print(xtable::xtable(imp_uni_status_T6, digits=c(0,2,0)), include.colnames = T, include.rownames=T,
      hline.after=c(-1, 0, nrow(imp_uni_status_T6)),
      sanitize.rownames.function=function(x){x}, sanitize.colnames.function = function(x){x},
      sanitize.text.function = function(x){x},
      NA.string = "", table.placement = "htp",
      caption.placement = "top", latex.environments = NULL,
      file="presentations/Update_Meeting_March2023/tables/textable_importance_reg_uni-status-T6.tex")

####--------------------------------------------------------------------------------
# how many predictors were selected in single fits?
apply(coef_list[[get_f("multi", "status", 6)]], 3, function(s) sum(s != 0)/length(s))
