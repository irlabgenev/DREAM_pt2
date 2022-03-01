require(tidyverse)
require(ggpubr)
require(patchwork)

# common layout
common_layout <- theme_classic(base_size = 10, base_family = "ArialMT") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.key.size = unit(0.7, "lines")
  )

setwd("D:/Projects/dream-rnascope/e20220201_MK-aceto_exposure/")

d <- read_tsv("e20220201_MK-aceto_exposure-mustn1.puncta_per_cell.txt")
p.mustn1.puncta_per_cell <- ggplot(d, aes(x=condition, y=value)) +
  geom_violin(color = "dodgerblue1", fill = "dodgerblue1", alpha = 0.4) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..), alpha = 0.5) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..), alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", alpha = 0.8) +
  stat_compare_means(aes(label = ..p.signif..), 
                     method = "wilcox.test",
                     label.x = 1.5, label.y = 29) +
  common_layout +
  scale_y_continuous(limits = c(0, 30), expand = c(0,0)) +
  ylab("puncta per cell") +
  ggtitle("Mustn1")

wilcox.test(value ~ condition, data = d)


d <- read_tsv("e20220201_MK-aceto_exposure-olfr16.total_int.txt")
p.olfr16.totalint <- ggplot(d, aes(x=condition, y=value)) +
  geom_violin(color = "dodgerblue1", fill = "dodgerblue1", alpha = 0.4) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..), alpha = 0.5) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..), alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", alpha = 0.8) +
  stat_compare_means(aes(label = ..p.signif..), 
                     method = "wilcox.test",
                     label.x = 1.5, label.y = 550) +
  common_layout +
  scale_y_continuous(limits = c(0, 600), expand = c(0,0)) +
  ylab("total fluo. int.") +
  ggtitle("Olfr16")

wilcox.test(value ~ condition, data = d)

d <- read_tsv("e20220201_MK-aceto_exposure-olfr16.pixel_per_cell.txt")
p.olfr16.pixels <- ggplot(d, aes(x=condition, y=value)) +
  geom_violin(color = "dodgerblue1", fill = "dodgerblue1", alpha = 0.4) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..), alpha = 0.5) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..), alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", alpha = 0.8) +
  stat_compare_means(aes(label = ..p.signif..), 
                     method = "wilcox.test",
                     label.x = 1.5, label.y = 1900) +
  common_layout +
  scale_y_continuous(limits = c(0, 2000), expand = c(0,0)) +
  ylab("puncta per cell") +
  ggtitle("Olfr16")

wilcox.test(value ~ condition, data = d)

p <- wrap_plots(p.mustn1.puncta_per_cell,
                p.olfr16.totalint,
                p.olfr16.pixels, nrow = 1)



ggsave(filename = "e20220201_MK-aceto_exposure-measures.pdf", 
       plot = p,
       device = "pdf", 
       units = "cm",
       width = 10, 
       height = 5, 
       useDingbats=FALSE)

