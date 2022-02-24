require(tidyverse)
require(ggpubr)
require(patchwork)
require(lsa)
require(ggtree)
require(ggtreeExtra)
require(ggnewscale)
require(viridis)
require(fields)
require(ggrastr)
require(ggstar)
require(gtools)
require(broom)

# PATH (to be edited)
current.path <- "N:/Novell/__articles in progress__/DREAM transcriptomic adaptations/"
setwd(current.path)

### LAYOUT
source(file.path(current.path, "figures/R scripts/common_layout.r"))

### FUNCTIONS
source(file.path(current.path, "figures/R scripts/Figure2-functions.r"))

### DATASETS
file_names <- c(
  "pairwise_distances.comp" = "DREAM_datasets/RData_20220209-pairwise_distances.comp-centroids-min10cells-PC1-15.rds",
  "or_data" = "DREAM_datasets/RData_20220209-or_data.rds",
  "mouse_or_phylo" = "DREAM_datasets/Mus_musculus.GRCm38.olfrs.clustalo_full.trim.phyml.rooted.20210916.txt"
)
file_names <- sapply(file_names, function(x) file.path(current.path, x))

pairwise_distances.comp <- read_rds(file_names["pairwise_distances.comp"])
or_data <- read_rds(file_names["or_data"])
mouse_or_phylo <- ape::read.tree(file_names["mouse_or_phylo"])

### FIGURES
### Fig 2

# A) OR phylogeny with cluster assignation 

# gene cluster color code
gene_clusters <- filter(or_data, !singleton) %>% 
  pull(gene_cluster) %>% 
  unique() %>%
  sort() %>%
  as.character()
gene_clusters.colors <- wesanderson::wes_palette("Zissou1", 
                                                 length(gene_clusters),
                                                 type = "continuous")
gene_clusters.colors <- as.character(gene_clusters.colors)
names(gene_clusters.colors) <- gene_clusters
gene_clusters.colors <- c(gene_clusters.colors, c("singleton" = "#6E3A07"))

gene_clusters.colors.subset <- subset(gene_clusters.colors,
                                      names(gene_clusters.colors) %in% c("8(chr7)", "5(chr14)"))
  
or.tree.p <- ggtree(mouse_or_phylo, 
                      branch.length = "none",
                      layout = "fan",
                      open.angle = 5,
                      size = 0.25) %<+%
    
  # append metadata
  #filter(or_data, n.osn >= 10) +
  or_data +
  
  # transcriptome cluster identity as colored tile
  geom_fruit(
    geom = geom_tile,
    mapping = aes(fill=transcriptome_cluster),
    width = 3,
    offset = 0.1) +
  scale_fill_manual(values = clusters.fill.osn, na.value = "white") +

  # gene cluster annotation
  new_scale_fill() +
  geom_fruit(
    geom = geom_star,
    mapping = aes(fill = ifelse(gene_cluster %in% c("8(chr7)", "5(chr14)") & n.osn >= 10, 
                                as.character(gene_cluster), NA)),
    color = NA,
    starshape = 23,
    size = 1.5,
    offset = 0.1) +
  scale_fill_manual(values = gene_clusters.colors.subset, na.value = NA) +
  theme(
    legend.position = "none"
  )
  
ggsave(filename = "figures/fig2-or_phylogeny_20220214.pdf",
       plot = or.tree.p,
       device = "pdf", 
       units = "cm",
       width = 6, 
       height = 6, 
       useDingbats=FALSE)

# B) Distribution of pairwise amino acid difference across bins of 
#    transcriptomic distances

# Almost identical sequence definition 
dist.aa.limits <- filter(pairwise_distances.comp, same_class) %>%
  pivot_longer(cols = starts_with("dist.aa"), 
               names_to = "seq.range", 
               values_to = "dist.aa") %>%
  group_by(seq.range) %>%
  summarise(almost_ident.thresh = quantile(dist.aa, probs = 0.01))
dist.aa.limits.names <- dist.aa.limits$seq.range
dist.aa.limits <- pull(dist.aa.limits, almost_ident.thresh)
names(dist.aa.limits) <- dist.aa.limits.names

#dist.aa.limits["dist.aa.miyata.trimmed"] <- 85

# 95th percentile of the distribution of intergenic distance
intergenic.dist <- filter(or_data, gene_cluster != "singleton") %>%
  group_by(gene_cluster) %>%
  arrange(start) %>%
  summarise(x = list(start[-1] - start[-n()])) %>%
  pull(x) %>%
  unlist()
intergenic.dist.mean <- mean(intergenic.dist)
intergenic.dist.95th <- quantile(intergenic.dist, probs = 0.95)
dist.genome.limits <- intergenic.dist.95th
names(dist.genome.limits) <- "dist.genome"

dist.limits <- c(dist.aa.limits, dist.genome.limits)


# definintion of bins of transcriptomic distances values for all pairwise
# distances within clusters
bin.levels = c(as.character(1:6), "7-10")
pairwise_distances.comp.bins <- filter(
  pairwise_distances.comp, same_class, same_cluster, gene_cluster != "singleton") %>%
  mutate(bin.trans.10th = as.numeric(cut(dist.trans, breaks = 10)),
         bin.trans.10th = ifelse(bin.trans.10th >= 7, "7-10", as.character(bin.trans.10th)),
         bin.trans.10th = factor(bin.trans.10th, levels = bin.levels),
         bin.genome.10th = as.numeric(cut(dist.genome, breaks = 10)),
         bin.genome.10th = ifelse(bin.genome.10th >= 7, "7-10", as.character(bin.genome.10th)),
         bin.genome.10th = factor(bin.genome.10th, levels = bin.levels),
         bin.genome.4genes = ceiling(log2(dist.genome)-log2(intergenic.dist.mean*4)),
         bin.genome.4genes = as.character(ifelse(bin.genome.4genes < 0, 0, bin.genome.4genes) + 1),
         bin.genome.4genes = factor(bin.genome.4genes, levels = as.character(1:7)),
         bin.aa.10th = as.numeric(cut(dist.aa.miyata.trimmed, breaks = 10)),
         bin.aa.10th = ifelse(bin.aa.10th >= 7, "7-10", as.character(bin.aa.10th)),
         bin.aa.10th = factor(bin.aa.10th, levels = bin.levels)) %>%
  select(pair, dist.genome, dist.trans, 
         dist.aa.miyata.trimmed, bin.trans.10th, bin.genome.10th, bin.genome.4genes, bin.aa.10th) 

# statistics on pairwise_distances.comp.bins combuted with transcriptomic 
# distance bins
pairwise_distances.comp.transbins.stat <- group_by(pairwise_distances.comp.bins, bin.trans.10th) %>%
  summarise(n = n(),
            almost_ident.prop = length(which(dist.aa.miyata.trimmed < dist.limits["dist.aa.miyata.trimmed"]))/n,
            nearby.prop = length(which(dist.genome < dist.limits["dist.genome"]))/n,
            cor.aa.pearson = cor(dist.trans, dist.aa.miyata.trimmed, method = "pearson"),
            cor.aa.spearman = cor(dist.trans, dist.aa.miyata.trimmed, method = "spearman"),
            cor.genome.pearson = cor(dist.trans, dist.genome, method = "pearson"),
            cor.genome.spearman = cor(dist.trans, dist.genome, method = "spearman"))

pairwise_distances.comp.noconfound <- bind_rows(
  
  # distant & homologous
  filter(pairwise_distances.comp, 
         !same_cluster | gene_cluster == "singleton",
       dist.aa.miyata.trimmed < dist.limits["dist.aa.miyata.trimmed"]) %>%
    mutate(cat = "D-H",
           genome = "distant",
           aa = "homologous"),
  
  # close & non-homologous
  filter(pairwise_distances.comp, 
         same_cluster, gene_cluster != "singleton",
         dist.aa.miyata.trimmed > dist.limits["dist.aa.miyata.trimmed"],
         dist.genome < dist.limits["dist.genome"]) %>%
    mutate(cat = "C-N",
           genome = "close",
           aa = "non-homologous"),
  
  # distant & non-homologous
  filter(pairwise_distances.comp, 
         !same_cluster | gene_cluster == "singleton",
         dist.aa.miyata.trimmed > dist.limits["dist.aa.miyata.trimmed"]) %>%
    mutate(cat = "D-N",
           genome = "distant",
           aa = "non-homologous"),
  
  # close & homologous
  filter(pairwise_distances.comp, 
         same_cluster, gene_cluster != "singleton",
         dist.aa.miyata.trimmed < dist.limits["dist.aa.miyata.trimmed"],
         dist.genome < dist.limits["dist.genome"]) %>%
    mutate(cat = "C-H",
           genome = "close",
           aa = "homologous")
) %>%
  mutate(cat = factor(cat, levels = c("C-H", "C-N", "D-H", "D-N")))

pairwise_distances.comp.noconfound.stat <- select(pairwise_distances.comp.noconfound, cat) %>%
  group_by(cat) %>%
  summarise(n = n())
pairwise_distances.comp.noconfound.tests <- select(pairwise_distances.comp.noconfound, cat, dist.trans) %>%
  ungroup() %>%
  mutate(CHvsCN = cat %in% c("C-H", "C-N"),
         CHvsDH = cat %in% c("C-H", "D-H"),
         DHvsDN = cat %in% c("D-H", "D-N")) %>%
  pivot_longer(cols = c("CHvsCN", "CHvsDH", "DHvsDN"),
               names_to = "comparison") %>%
  filter(value) %>%
  nest(data = c(cat, dist.trans)) %>%
  mutate(wilcox.test = map(data, ~ tidy(wilcox.test(dist.trans ~ cat, data = .x))),
         w = unlist(map(wilcox.test, ~ .x$statistic)),
         p = unlist(map(wilcox.test, ~ .x$p.value)),
         p.adj = p.adjust(p, method = "bonferroni"),
         p.signif = stars.pval(p.adj)) %>%
  select(-data, -wilcox.test)

write_delim(x = pairwise_distances.comp.noconfound.tests,
            file = "figures/Figure2G_stats_20220214.txt",
            delim = "\t")

p.close_distant <- ggplot(pairwise_distances.comp.noconfound, aes(x=cat, y=dist.trans)) +
  geom_violin(color = "dodgerblue1", fill = "dodgerblue1", alpha = 0.4, scale = "width") +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..), alpha = 0.5) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..), alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", alpha = 0.8) +
  stat_compare_means(aes(label = ..p.signif..), 
                     method = "wilcox.test", comparisons = list(
                       c("C-H", "C-N"),
                       c("D-H", "D-N"), 
                       c("C-H", "D-H")),
                     label.y = c(28, 32, 34)) +
  common_layout +
  scale_y_continuous(limits = c(0, 38), expand = c(0,0), breaks = c(10, 20, 30)) +
  ylab("transcriptomic distance")

# 
p.binning.prop <- ggplot(pairwise_distances.comp.transbins.stat, 
                         aes(x=bin.trans.10th)) +
  geom_col(aes(y=almost_ident.prop), fill = "#8D85BE", color = "#8D85BE", width = 0.4, position = position_nudge(-0.25)) +
  stat_cor(aes(x=as.numeric(bin.trans.10th), y=almost_ident.prop), method = "spearman", label.y = 0.3, label.x = 4, color = "#8D85BE") +
  geom_col(aes(y=nearby.prop), fill = "#66C2A4", color = "#66C2A4", width = 0.4, position = position_nudge(+0.25)) +
  stat_cor(aes(x=as.numeric(bin.trans.10th), y=nearby.prop), method = "spearman", label.y = 0.25, label.x = 4, color = "#66C2A4") +
  scale_x_discrete(limits = c(as.character(1:6), "7-10")) +
  scale_y_continuous(limits = c(0, 0.33), expand = c(0,0)) +
  common_layout +
  xlab("transcriptomic distance bin") +
  ylab("proportion of pairs")

dist.trans.cuts <- levels(cut(pairwise_distances.comp.bins$dist.trans, breaks = 10)) %>%
  str_remove(., ".+,") %>%
  str_remove(., "\\]") %>%
  as.numeric()

dist.aa.cuts <- levels(cut(pairwise_distances.comp.bins$dist.aa.miyata.trimmed, breaks = 10)) %>%
  str_remove(., ".+,") %>%
  str_remove(., "\\]") %>%
  as.numeric()

dist.genome.cuts <- levels(cut(pairwise_distances.comp.bins$dist.genome, breaks = 10)) %>%
  str_remove(., ".+,") %>%
  str_remove(., "\\]") %>%
  as.numeric()

dist.genome.cuts.4genes <- intergenic.dist.mean*4*2**c(0, seq_along(unique(as.numeric(pairwise_distances.comp.bins$bin.genome.4genes))))

# binning and violin plots of transcriptomic distance vs genomic distance
p.binning.genome <- ggplot(pairwise_distances.comp.bins, 
                       aes(x=bin.genome.10th, y=dist.trans)) +
  geom_violin(aes(group=bin.genome.10th), color = "dodgerblue1", fill = "dodgerblue1", alpha = 0.4, scale = "width") +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..), alpha = 0.5) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..), alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", alpha = 0.8) +
  scale_y_continuous(limits = c(0, 38), expand = c(0,0)) +
  common_layout +
  xlab("genomic distance (bin #)") +
  ylab("transcriptomic distance")

# scatterplot of the same dataset
p.scatterplot.genome <- ggplot(pairwise_distances.comp.bins, 
                               aes(x=dist.genome, y=dist.trans)) +
  rasterize(geom_point(alpha=0.3, size = 0.2, color="dodgerblue4", shape = 19), dpi = 600) +
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",
                  bins = 10, alpha = 0.4, show.legend = F) +
  #geom_vline(xintercept = dist.genome.cuts, linetype = "dashed") +
  #geom_vline(xintercept = dist.genome.cuts.4genes[1:6], linetype = "dashed", color = "red") +
  #stat_smooth(method = "lm", se = F) +
  #stat_cor(method = "spearman") +
  #stat_smooth(data = filter(pairwise_distances.comp.bins, as.numeric(bin.genome.10th) <= 3),
  #            method = "lm", se = F, color = "red") +
  scale_x_continuous(limits = c(0, 5.1e+6), expand = c(0,0)) +
  common_layout +
  xlab("genomic distance") +
  ylab("transcriptomic distance")


p.binning.genome.2 <- ggplot(pairwise_distances.comp.bins, 
                           aes(x=bin.genome.4genes, y=dist.trans)) +
  geom_violin(aes(group=bin.genome.4genes), color = "dodgerblue1", fill = "dodgerblue1", alpha = 0.4, scale = "width") +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..), alpha = 0.5) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..), alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", alpha = 0.8) +
  scale_y_continuous(limits = c(0, 38), expand = c(0,0)) +
  common_layout +
  xlab("genomic distance (bin #)") +
  ylab("transcriptomic distance")

x <- filter(pairwise_distances.comp.bins, bin.genome.10th %in% c("1", "2", "3")) %>%
  pull(dist.trans)
y <- filter(pairwise_distances.comp.bins, bin.genome.10th %in% c("1", "2", "3")) %>%
  pull(dist.genome)
p.corr.genome <- cor.test(x, y, method = "spearman")
y.aa <- filter(pairwise_distances.comp.bins, bin.genome.10th %in% c("1", "2", "3")) %>%
  pull(dist.aa.miyata.trimmed)
p.corr.genome.but.aa <- cor.test(x, y.aa, method = "spearman")
length(x)

x <- filter(pairwise_distances.comp.bins, bin.aa.10th %in% c("1", "2", "3")) %>%
  pull(dist.trans)
y <- filter(pairwise_distances.comp.bins, bin.aa.10th %in% c("1", "2", "3")) %>%
  pull(dist.aa.miyata.trimmed)
p.corr.aa <- cor.test(x, y, method = "spearman")
y.genome <- filter(pairwise_distances.comp.bins, bin.aa.10th %in% c("1", "2", "3")) %>%
  pull(dist.genome)
p.corr.aa.but.genome <- cor.test(x, y.genome, method = "spearman")
length(x)

x <- filter(pairwise_distances.comp.bins, bin.genome.4genes %in% c("1", "2", "3")) %>%
  pull(dist.trans)
y <- filter(pairwise_distances.comp.bins, bin.genome.4genes %in% c("1", "2", "3")) %>%
  pull(dist.genome)
p.corr.genome.4genes <- cor.test(x, y, method = "spearman")
length(x)

p.binning.aa <- ggplot(pairwise_distances.comp.bins, 
                       aes(x=as.character(bin.aa.10th), y=dist.trans)) +
  geom_violin(aes(group=bin.aa.10th), color = "dodgerblue1", fill = "dodgerblue1", alpha = 0.4, scale = "width") +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..), alpha = 0.5) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..), alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", alpha = 0.8) +
  scale_x_discrete(limits = c(as.character(1:6), "7-10")) +
  scale_y_continuous(limits = c(0, 38), expand = c(0,0)) +
  common_layout +
  xlab("amino acid difference (bin #)") +
  ylab("transcriptomic distance")

# scatterplot of the same dataset (--> supp)
p.scatterplot.aadiff <- ggplot(pairwise_distances.comp.bins, 
                               aes(x=dist.aa.miyata.trimmed, y=dist.trans)) +
  rasterize(geom_point(alpha=0.3, size = 0.2, color="dodgerblue4", shape = 19), dpi = 600) +
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",
                  bins = 10, alpha = 0.4, show.legend = F) +
  #geom_vline(xintercept = dist.aa.cuts, linetype = "dashed") +
  #stat_smooth(method = "lm", se = F) +
  #stat_cor(method = "spearman") +
  #stat_smooth(data = filter(pairwise_distances.comp.bins, as.numeric(bin.aa.10th) <= 3),
  #            method = "lm", se = F, color = "red") +
  scale_x_continuous(limits = c(0, 430), expand = c(0,0)) +
  common_layout +
  xlab("amino acid difference") +
  ylab("transcriptomic distance")


layout <- "
AABB
CCD#
"

p.binning <- wrap_plots(p.binning.genome, p.binning.prop, p.binning.aa, p.close_distant, design = layout)

ggsave(filename = "figures/fig2efg_binning_20220214.pdf", 
       plot = p.binning,
       device = "pdf", 
       units = "cm",
       width = 20, 
       height = 10, 
       useDingbats=FALSE)

p.density.trans <- ggplot(pairwise_distances.comp.bins,
                           aes(x=dist.trans)) +
  stat_density(fill = "dodgerblue1", color = "transparent", size = 1, alpha = 0.4) +
  scale_x_continuous(limits = c(0, 32), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 7.3e-02), expand = c(0,0)) +
  geom_vline(xintercept = dist.trans.cuts, color = "white", size = 1.5) +
  geom_text(data = data.frame(label=paste(c("bin ", rep("", length(dist.trans.cuts)-1)),
                                            seq_along(dist.trans.cuts), sep = ""),
                              x=dist.trans.cuts-(dist.trans.cuts[2]-dist.trans.cuts[1])/2), 
            mapping = aes(x = x, y = 0.01, label = label)) +
  common_layout +
  xlab("transcriptomic distance\n(euclidian distances between PCs)") +
  ylab("density")

p.density.genome.1 <- ggplot(pairwise_distances.comp.bins,
                           aes(x=dist.genome)) +
  stat_density(fill = "#66C2A4", color = "transparent", size = 1) +
  scale_x_continuous(limits = c(0, 5.2e+06), expand = c(0,0), 
                     breaks = (0:5)*1e+06, labels = c(0:5)) +
  scale_y_continuous(limits = c(0, 7.8e-07), expand = c(0,0)) +
  geom_vline(xintercept = dist.genome.cuts, color = "white", size = 1.5) +
  geom_vline(xintercept = dist.limits["dist.genome"], color = "red") +
  common_layout +
  xlab("genomic distance (Kb)") +
  ylab("density")

p.density.genome.2 <- ggplot(pairwise_distances.comp.bins,
                           aes(x=dist.genome)) +
  stat_density(fill = "#66C2A4", color = "transparent", size = 1) +
  scale_x_continuous(limits = c(0, 5.2e+06), expand = c(0,0), 
                     breaks = (0:5)*1e+06, labels = c(0:5)) +
  scale_y_continuous(limits = c(0, 7.8e-07), expand = c(0,0)) +
  geom_vline(xintercept = dist.genome.cuts.4genes, color = "white", size = 1.5) +
  geom_vline(xintercept = dist.limits["dist.genome"], color = "red") +
  common_layout +
  xlab("genomic distance (Kb)") +
  ylab("density")

p.density.aa <- ggplot(pairwise_distances.comp.bins,
                           aes(x=dist.aa.miyata.trimmed)) +
  stat_density(fill = "#8D85BE", color = "transparent", size = 1) +
  scale_x_continuous(limits = c(0, 420), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 8.4e-03), expand = c(0,0)) +
  geom_vline(xintercept = dist.aa.cuts, color = "white", size = 1.5) +
  geom_vline(xintercept = dist.limits["dist.aa.miyata.trimmed"], color = "red") +
  common_layout +
  xlab("amino acid difference (Miyata score)") +
  ylab("density")

layout <- "
AB
FC
GD
#E
"

p.densities <- wrap_plots(p.density.trans,
                          p.density.genome.1,
                          p.density.genome.2, 
                          p.binning.genome.2,
                          p.density.aa,
                          ncol = 1)


p.scatterplots <- wrap_plots(p.scatterplot.aadiff, 
                             p.scatterplot.genome,
                             ncol = 1)

ggsave(filename = "figures/figS4_densities_20220215.pdf", 
       plot = p.densities,
       device = "pdf", 
       units = "cm",
       width = 8, 
       height = 20, 
       useDingbats=FALSE)

ggsave(filename = "figures/figS4_scatterplots_20220214.pdf", 
       plot = p.scatterplots,
       device = "pdf", 
       units = "cm",
       width = 11, 
       height = 15, 
       useDingbats=FALSE)

# P element regulated genes
interesting_clusters <- filter(or_data, !is.na(cluster_annotation)) %>% 
  pull(gene_cluster) %>% 
  unique() %>%
  as.character()
pairwise_distances.clusters <- filter(pairwise_distances.comp, as.character(gene_cluster) %in% interesting_clusters) %>%
  select(pair, gene_cluster, dist.trans, dist.aa.miyata.trimmed, dist.genome) %>%
  mutate(gene_name.1 = str_remove(pair, ":.+"),
         gene_name.2 = str_remove(pair, ".+:"))

# mirror dataset for matrices
pairwise_distances.clusters <- rename(pairwise_distances.clusters,
                                      gene_name.b = gene_name.1,
                                      gene_name.a = gene_name.2) %>%
  rename(gene_name.1 = gene_name.a,
         gene_name.2 = gene_name.b) %>%
  bind_rows(pairwise_distances.clusters)

# add cluster annotations
pairwise_distances.clusters <- left_join(pairwise_distances.clusters, select(or_data, label, start, transcriptome_cluster, cluster_annotation),
                                         by = c("gene_name.1" = "label")) %>%
  left_join(select(or_data, label, start, transcriptome_cluster, cluster_annotation),
            by = c("gene_name.2" = "label"), suffix = c(".1", ".2")) %>%
  group_by(gene_cluster) %>%
  mutate(genome.rank.1 = order(start.1),
         genome.rank.2 = order(start.2))

# gradient genome
d <- data.frame(
  fill = c(min(pairwise_distances.comp$dist.genome, na.rm = T),
           max(pairwise_distances.comp$dist.genome, na.rm = T)),
  x = c(1,1),
  y = c(1,2)
)

ggplot(d, aes(x, y)) +
  geom_tile(aes(fill = fill))

# transcriptomic distance pairwise matrices, ordered with genomic locations
pairwise_distances.clusters.list <- split(pairwise_distances.clusters, pairwise_distances.clusters$gene_cluster)
p.matrices <- mapply(
  d = pairwise_distances.clusters.list,
  gene_cluster = names(pairwise_distances.clusters.list),
  SIMPLIFY = F,
  FUN = function(d, gene_cluster) {
    
    # x axis
    gene_name.x <- select(d, gene_name.1, start.1) %>%
      arrange(start.1) %>%
      distinct() %>%
      pull(gene_name.1)
    
    # y axis
    gene_name.y <- select(d, gene_name.2, start.2) %>%
      arrange(start.2) %>%
      distinct() %>%
      pull(gene_name.2)
    
    p.regulation <- ggplot(d, aes(x=1, y=gene_name.2)) +
      geom_point(aes(color=cluster_annotation.2)) +
      scale_y_discrete(limits = gene_name.y) +
      coord_fixed() +
      common_layout +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      ) +
      ylab("gene name")
    
    p.dist.trans <- ggplot(d, aes(x=gene_name.1, y=gene_name.2)) +
      geom_tile(aes(fill = dist.trans)) +
      scale_x_discrete(limits = gene_name.x) +
      scale_y_discrete(limits = gene_name.y) +
      scale_fill_distiller(limits = c(0,32), direction = -1) +
      coord_fixed() +
      common_layout +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
      )
    
    p.dist.aa <- ggplot(d, aes(x=gene_name.1, y=gene_name.2)) +
      geom_tile(aes(fill = dist.aa.miyata.trimmed)) +
      scale_x_discrete(limits = gene_name.x) +
      scale_y_discrete(limits = gene_name.y) +
      scale_fill_distiller(limits = c(0, 375), palette = "BuPu", direction = -1) +
      coord_fixed() +
      common_layout +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
      )
    
    p <- wrap_plots(p.regulation, p.dist.trans, p.dist.aa, nrow = 1, guides = "collect") +
      plot_annotation(title = gene_cluster)
    
    return(p)
  }
)

for (i in 1:4) {
  filename <- str_replace("figures/fig2-gene_cluster#_regulation_20220214.pdf",
                          "#", as.character(i))
  ggsave(filename = filename,
         plot = p.matrices[[i]],
         device = "pdf", 
         units = "cm",
         width = 20, 
         height = 20, 
         useDingbats=FALSE)
}

### Supplementary figure 2
# OR phylogeny with numbers of OSN, horizontal
or.tree.nosn <- ggtree(mouse_or_phylo, 
                       branch.length = "none",
                       layout = "fan",
                       open.angle = 10,
                       size = 0.25) %<+%
  
  # append metadata
  or_data +
  
  # transcriptome cluster identity as colored tile
  #geom_fruit(
  #  geom = geom_tile,
  #  mapping = aes(fill=transcriptome_cluster),
  #  width = 3,
  #  offset = 0.1) +
  #scale_fill_manual(values = clusters.fill.osn, na.value = "white") +
  
  # transcriptome cluster identity as colored tile
  #new_scale_fill() +
  #geom_fruit(
  #  geom = geom_tile,
  #  mapping = aes(fill=gene_cluster),
  #  width = 3,
  #  offset = 0.1) +
  #scale_fill_manual(values = gene_clusters.colors, na.value = "white") +
  
  # number of OSN expressing a given OR
  geom_fruit(
    geom = geom_bar,
    mapping = aes(y=ID, x=n.osn),
    stat="identity",
    orientation="y",
    pwidth=2,
    axis.params=list(
      axis       = "x",
      text.size  = 1,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 3,
    ),
    grid.params=list()
  ) +
  layout_rectangular()

ggsave(filename = "figures/figs1d-phylo_osn_counts_20220214.pdf",
       plot = or.tree.nosn,
       device = "pdf",
       units = "cm",
       width = 18,
       height = 18,
       useDingbats=FALSE)

# R2 with PC usage
dist_files <- list.files("DREAM_datasets/", pattern = "RData_20220209-pairwise_distances-simple-centroids-min10cells-PC1-")
pc.num <- as.numeric(str_extract(dist_files, "[0-9]+(?=\\.rds)"))
d.dist_files <- data.frame(
  file.name = dist_files,
  pc.num = pc.num
) %>%
  arrange(pc.num)

rhos <- mapply(
  file.name = paste("DREAM_datasets", d.dist_files$file.name, sep = "/"),
  pc.num = d.dist_files$pc.num,
  MoreArgs = list("metadata" = select(pairwise_distances.comp.bins, 
                                      pair, dist.aa.miyata.trimmed, dist.genome,
                                      bin.genome.10th, bin.genome.4genes, bin.aa.10th)),
  SIMPLIFY = F,
  FUN = function(file.name, pc.num, metadata) {
    
    d <- readRDS(file = file.name) %>%
      inner_join(metadata, by = "pair")
    
    rho.genome.10th.all <- cor(
      x = pull(d, dist.trans),
      y = pull(d, dist.genome),
      method = "spearman"
    )
    rho.genome.10th.3firstbins <- cor(
      x = pull(d, dist.trans)[which(d$bin.genome.10th %in% as.character(1:3))],
      y = pull(d, dist.genome)[which(d$bin.genome.10th %in% as.character(1:3))],
      method = "spearman"
    )
    rho.genome.4genes.3firstbins <- cor(
      x = pull(d, dist.trans)[which(d$bin.genome.4genes %in% as.character(1:3))],
      y = pull(d, dist.genome)[which(d$bin.genome.4genes %in% as.character(1:3))],
      method = "spearman"
    )
    rho.aa.all <- cor(
      x = pull(d, dist.trans),
      y = pull(d, dist.aa.miyata.trimmed),
      method = "spearman"
    )
    rho.aa.3firstbins <- cor(
      x = pull(d, dist.trans)[which(d$bin.aa.10th %in% as.character(1:3))],
      y = pull(d, dist.aa.miyata.trimmed)[which(d$bin.aa.10th %in% as.character(1:3))],
      method = "spearman"
    )
    
    d.rho <- data.frame(
      rho = c(rho.genome.10th.all, rho.genome.10th.3firstbins,
              rho.genome.4genes.3firstbins, rho.aa.all, rho.aa.3firstbins),
      distance.y = c("genome", "genome", "genome", "aa", "aa"),
      dataset = c("all", "3 first bins", "3 first bins", "all", "3 first bins"),
      binning = c("genome-10th", "genome-10th", "genome-4th", "aa-10th", "aa-10th"),
      PC.num = rep(pc.num, 5)
    )
    
    return(d.rho)
    
  }
) %>%
  do.call(what = rbind, args = .) %>%
  remove_rownames()

p.rhos <- ggplot(rhos, aes(x = PC.num, y = rho)) +
  geom_point(aes(color = binning, shape = dataset)) +
  scale_y_continuous(limits = c(0, 0.32)) +
  scale_shape_manual(values = c("all" = 19, "3 first bins" = 17)) +
  scale_x_continuous(breaks = 2:15, labels = paste("1", as.character(2:15), sep = "-")) +
  xlab("PCs") +
  common_layout

ggsave(filename = "figures/fig2supp_rho_PCs_20220211.pdf", 
       plot = p.rhos,
       device = "pdf", 
       units = "cm",
       width = 10, 
       height = 5, 
       useDingbats=FALSE)

# add Ns  
d.1 <- group_by(pairwise_distances.comp.bins, first3 = bin.genome.10th %in% c("1", "2", "3")) %>%
  summarise(n = n()) %>%
  mutate(binning = "bin.genome.10th")
d.2 <- group_by(pairwise_distances.comp.bins, first3 = bin.genome.4genes %in% c("1", "2", "3")) %>%
  summarise(n = n()) %>%
  mutate(binning = "bin.genome.4genes")
d.3 <- group_by(pairwise_distances.comp.bins, first3 = bin.aa.10th %in% c("1", "2", "3")) %>%
  summarise(n = n()) %>%
  mutate(binning = "bin.aa.10th")
d <- rbind(d.1, d.2, d.3)

p.ns <- ggplot(d, aes(x=binning, y=n)) +
  geom_col(aes(fill=first3), color = "black", width = 1, show.legend = F) +
  scale_y_continuous(limits = c(0,3750), expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  ylab("num. of comp.") +
  common_layout 

ggsave(filename = "figures/fig2supp_comp_num_20220211.pdf", 
plot = p.ns,
device = "pdf", 
units = "cm",
width = 5, 
height = 5, 
useDingbats=FALSE)
