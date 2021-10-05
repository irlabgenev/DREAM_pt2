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

# PATH (to be edited)
current.path <- "..."
setwd(current.path)

### LAYOUT
source(file.path(current.path, "common_layout.r"))

### FUNCTIONS
source(file.path(current.path, "Figure2-functions.r"))

### DATASETS
file_names <- c(
  "pairwise_distances.comp" = "datasets/RData_20210929-pairwise_distances.comp-centroids-min10cells-PC1-19.rds",
  "or_data" = "datasets/RData_20211003-or_data.rds",
  "osn_GEP.loading" = "datasets/RData_20210929-osn_GEP.loading.rds",
  "mouse_or_phylo" = "datasets/Mus_musculus.GRCm38.olfrs.clustalo_full.trim.phyml.rooted.20210916.txt"
)
file_names <- sapply(file_names, function(x) file.path(current.path, x))

pairwise_distances.comp <- read_rds(file_names["pairwise_distances.comp"])
or_data <- read_rds(file_names["or_data"])
osn_GEP.loading <- read_rds(file_names["osn_GEP.loading"])
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
  
ggsave(filename = "figures/fig2-or_phylogeny.pdf",
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
       !same_cluster,
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
         !same_cluster,
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

p.close_distant <- ggplot(pairwise_distances.comp.noconfound, aes(x=cat, y=dist.trans)) +
  geom_violin(color = "dodgerblue1", fill = "dodgerblue1", alpha = 0.4) +
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

p.binning.violin.aadiff <- ggplot(pairwise_distances.comp.bins, 
                                  aes(x=bin.trans.10th, y=dist.aa.miyata.trimmed)) +
  geom_violin(color = "grey80", fill = "grey80") +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..), alpha = 0.5) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..), alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", alpha = 0.8) +
  stat_compare_means(method = "anova", label.y = 400) +
  scale_y_continuous(limits = c(0, 450), expand = c(0,0)) +
  common_layout +
  xlab("transcriptomic distance bin") +
  ylab("amino acid differences")

p.binning.prop <- ggplot(pairwise_distances.comp.transbins.stat, 
                         aes(x=bin.trans.10th)) +
  geom_col(aes(y=almost_ident.prop), fill = "#8D85BE", color = "#8D85BE", width = 0.4, position = position_nudge(-0.25)) +
  stat_cor(aes(x=as.numeric(bin.trans.10th), y=almost_ident.prop), method = "spearman", label.y = 0.3, label.x = 4, color = "#8D85BE") +
  geom_col(aes(y=nearby.prop), fill = "#66C2A4", color = "#66C2A4", width = 0.4, position = position_nudge(+0.25)) +
  stat_cor(aes(x=as.numeric(bin.trans.10th), y=nearby.prop), method = "spearman", label.y = 0.25, label.x = 4, color = "#66C2A4") +
  scale_x_discrete(limits = c(as.character(1:6), "7-10")) +
  scale_y_continuous(limits = c(0, 0.32), expand = c(0,0)) +
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

p.binning.genome <- ggplot(pairwise_distances.comp.bins, 
                       aes(x=bin.genome.10th, y=dist.trans)) +
  geom_violin(aes(group=bin.genome.10th), color = "dodgerblue1", fill = "dodgerblue1", alpha = 0.4) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..), alpha = 0.5) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..), alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", alpha = 0.8) +
  scale_y_continuous(limits = c(0, 38), expand = c(0,0)) +
  common_layout +
  xlab("genomic distance (bin #)") +
  ylab("transcriptomic distance")

p.binning.genome.2 <- ggplot(pairwise_distances.comp.bins, 
                           aes(x=bin.genome.4genes, y=dist.trans)) +
  geom_violin(aes(group=bin.genome.4genes), color = "dodgerblue1", fill = "dodgerblue1", alpha = 0.4) +
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
length(x)

x <- filter(pairwise_distances.comp.bins, bin.aa.10th %in% c("1", "2", "3")) %>%
  pull(dist.trans)
y <- filter(pairwise_distances.comp.bins, bin.aa.10th %in% c("1", "2", "3")) %>%
  pull(dist.aa.miyata.trimmed)
p.corr.aa <- cor.test(x, y, method = "spearman")
length(x)

x <- filter(pairwise_distances.comp.bins, bin.genome.4genes %in% c("1", "2", "3")) %>%
  pull(dist.trans)
y <- filter(pairwise_distances.comp.bins, bin.genome.4genes %in% c("1", "2", "3")) %>%
  pull(dist.aa.miyata.trimmed)
p.corr.genome.4genes <- cor.test(x, y, method = "spearman")
length(x)

p.binning.aa <- ggplot(pairwise_distances.comp.bins, 
                       aes(x=as.character(bin.aa.10th), y=dist.trans)) +
  geom_violin(aes(group=bin.aa.10th), color = "dodgerblue1", fill = "dodgerblue1", alpha = 0.4) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..), alpha = 0.5) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..), alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", alpha = 0.8) +
  scale_x_discrete(limits = c(as.character(1:6), "7-10")) +
  scale_y_continuous(limits = c(0, 38), expand = c(0,0)) +
  common_layout +
  xlab("amino acid difference (bin #)") +
  ylab("transcriptomic distance")

layout <- "
AABB
CCD#
"

p.binning <- wrap_plots(p.binning.genome, p.binning.prop, p.binning.aa, p.close_distant, design = layout)

ggsave(filename = "figures/fig2_binning.pdf", 
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
                          p.density.genome.1, p.binning.genome,
                          p.density.genome.2, p.binning.genome.2,
                          p.density.aa, p.binning.aa,
                          design = layout)

ggsave(filename = "figures/fig2supp_densities.pdf", 
       plot = p.densities,
       device = "pdf", 
       units = "cm",
       width = 15, 
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
  fill = c(min(pairwise_distances.comp.transbins$dist.genome),
           max(pairwise_distances.comp.transbins$dist.genome)),
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
  filename <- str_replace("figures/fig2-gene_cluster#_regulation.pdf",
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
  geom_fruit(
    geom = geom_tile,
    mapping = aes(fill=transcriptome_cluster),
    width = 3,
    offset = 0.1) +
  scale_fill_manual(values = clusters.fill.osn, na.value = "white") +
  
  # transcriptome cluster identity as colored tile
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(fill=gene_cluster),
    width = 3,
    offset = 0.1) +
  scale_fill_manual(values = gene_clusters.colors, na.value = "white") +
  
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

ggsave(filename = "figures/fig1supp-phylo_osn_counts.pdf",
       plot = ggplotify::as.ggplot(or.tree.nosn, angle=90, scale = 1.5) + coord_fixed(4),
       device = "pdf",
       units = "cm",
       width = 18,
       height = 18,
       useDingbats=FALSE)

