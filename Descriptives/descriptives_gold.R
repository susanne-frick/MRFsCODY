####--------------------- descriptives on status tests -------------------------####

library(ggplot2)
library(colorspace)
library(gridExtra)

# load data
load("Data_preparation/Daten_cody_gold_sub_wide.RData")
head(Daten_wide)

Daten_gold <- Daten_wide[, grep("Score_gold_sum", colnames(Daten_wide))]
head(Daten_gold)

# number of status tests
G <- 6

####--------------------------------- plots -----------------------------------####

# distributions for each timepoint separately
for (g in 1:G) {
  hist(Daten_gold[, g], breaks = 50, main = paste("Status Test", g), xlab = "Score")
}

# kernel density estimates together
# data to long format
Daten_gold_long <- data.frame("ID" = rep(1:nrow(Daten_gold), times = G),
                              "Timepoint" = rep(1:G, each = nrow(Daten_gold)),
                              "Score" = unlist(Daten_gold),
                              row.names = NULL)
Daten_gold_long$Timepoint <- factor(Daten_gold_long$Timepoint)
head(Daten_gold_long)

plot_violin <- ggplot(data = Daten_gold_long,
                    aes(y = Score, x = Timepoint, fill = Timepoint)) +
  geom_violin(show.legend=FALSE) +
  scale_fill_manual(values=sequential_hcl(6)) +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11),
        title=element_text(size = 11)) +
  ggtitle(label = "Status Tests")
ggsave(plot_violin,
       file="plots/plot_score_violin.pdf",
       width=15, height=15, units="cm")

####--------------------------------- figures -----------------------------------####

range(Daten_gold, na.rm = TRUE)

apply(Daten_gold, 2, summary)
psych::describe(Daten_gold)

# intraindividual differences: changes to T1
Daten_gold_intra <- Daten_gold[, 2:G] - Daten_gold[, 1]
psych::describe(Daten_gold_intra)
