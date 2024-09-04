####------------ plot densitities of gold coin scores --------------------####

load("Data_preparation/Daten_cody_gold_sub_wide.RData")

# to long format
gold_names <- grep("Score_gold_sum_", colnames(Daten_wide), value = TRUE)
Daten_gold <- Daten_wide[, c("Id", gold_names)]
colnames(Daten_gold) <- c("Id", paste0("gold_", 1:6))
gold_long <- reshape(Daten_gold,
                     direction = "long", varying = colnames(Daten_gold)[-c(1)],
                     sep = "_", idvar = "Id")
head(gold_long)

# density plots
library(ggplot2)
library(colorspace)
library(gridExtra)

gold_long$time <- factor(gold_long$time)
plot_gold_violin <- ggplot(data = gold_long,
                              aes(y = gold, x = time)) +
  geom_violin(color = "gold") +
    labs(y = "Score", x = "Timepoint") +
    theme(axis.text=element_text(size=11),
          axis.title=element_text(size=11))

plot_gold_violin

spaghetti_gold <- ggplot(data = gold_long,
                         aes(y = gold, x = time, group = Id, color = Id)) +
  geom_line() +
  labs(y = "Score", x = "Timepoint") +
  theme(legend.position = "none") +
  stat_summary(aes(group = 1), geom = "line", fun.y = mean,
               color = "black", size = 1)
ggsave(spaghetti_gold, file = "plots/plot_spaghetti_gold.pdf", height = 10, width = 12, units = "cm")
file.copy(from = "plots/plot_spaghetti_gold.pdf", to = "~/Dokumente/FAIR/Reha/presentations/Kolloquium_Bamberg_June2023/figures",
          overwrite = TRUE)
