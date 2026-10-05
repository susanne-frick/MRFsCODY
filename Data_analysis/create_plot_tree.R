####-------------- create plot of a single tree (for presentations) ----------------####

library(randomForestSRC)

fit1 <- readRDS("results_MRFs/results_MRF_f1.rds")

tree1 <- get.tree(fit1$fit, tree.id = 1)

# did not work for some reason -> saved from RStudio
pdf("plots/plot_tree.pdf")
plot(tree1)
dev.off()

file.copy(from = "plots/plot_tree.pdf", to = "~/Dokumente/FAIR/Reha/presentations/Kolloquium_Bamberg_June2023/figures",
          overwrite = TRUE)
