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

# PATH (to be edited)
current.path <- "..."
setwd(current.path)

### FUNCTIONS
source(file.path(current.path, "Figure2-functions.r"))
source(file.path(current.path, "common_layout.r"))

### DATASETS
file_names <- c(
  
  # coordinates of individual functional OR genes and their respective genomic
  # cluster
  "or_genomic_clusters" = "datasets/Mus_musculus.mm10.functional_OR_gene_clusters.20210914.rds",
  
  # functional OR phylogeny
  "or_phylogeny" = "datasets/Mus_musculus.GRCm38.olfrs.clustalo_full.trim.phyml.rooted.20210916.txt",
  
  # All subtrees including OR genes detected in the scRNA-seq
  "subtrees" = "datasets/Mus_musculus.GRCm38.olfrs.clustalo_full.trim.phyml.rooted.subtrees.20210916.txt",
  
  # functional OR class I / class II assignation
  "or_classes" = "datasets/Mus_musculus.mm10.olfr_classes.20210907.rds",
  
  # different metrics of distances between OSN populations and their
  # corresponding OR
  "pairwise_distances" = "datasets/Mus_musculus.mm10.functional_OR_pairwise_distances.20210810.rds",
  
  # OSN transcriptome metadata, after removing cells expressing more than one OR
  "osn_metadata" = "datasets/Mus_musculus.mm10.filtered_mature_OSN_metadata.20210910.rds",
  
  # OSN integrated count matrix
  "osn_count_matrix" = "datasets/Mus_musculus.mm10.filtered_mature_OSN_integrated_count_matrix.20210910.rds",
  
  # OSN PCA
  "osn_PCA" = "datasets/Mus_musculus.mm10.filtered_mature_OSN_PCA.20210920.rds",
  
  # amino-acid distances
  "dist.aa.levenshtein.trimmed" = "datasets/Mus_musculus.GRCm38.olfrs_and_taars.clustalo_full-trimmed-levenshtein.20210920.txt",
  "dist.aa.levenshtein.tm" = "datasets/Mus_musculus.GRCm38.olfrs_and_taars.clustalo_full-transmembrane-levenshtein.20210920.txt",
  "dist.aa.levenshtein.cons" = "datasets/Mus_musculus.GRCm38.olfrs_and_taars.clustalo_full-conserved-levenshtein.20210920.txt",
  "dist.aa.grantham.trimmed" = "datasets/Mus_musculus.GRCm38.olfrs_and_taars.clustalo_full-trimmed-grantham.20210922.txt",
  "dist.aa.miyata.trimmed" = "datasets/Mus_musculus.GRCm38.olfrs_and_taars.clustalo_full-trimmed-withgaps-miyata.20210924.txt",
  "dist.aa.miyata.transmembrane" = "datasets/Mus_musculus.GRCm38.olfrs_and_taars.clustalo_full-transmembrane-miyata-withgaps.20210927.txt",
  "dist.aa.miyata.intra" = "datasets/Mus_musculus.GRCm38.olfrs_and_taars.clustalo_full-intra-miyata-withgaps.202109257.txt",
  "dist.aa.miyata.extra" = "datasets/Mus_musculus.GRCm38.olfrs_and_taars.clustalo_full-extra-miyata-withgaps.20210927.txt",
  
  # GEP loading matrix
  "osn_GEP" = "datasets/cNMF_GEPs_no_olfr.usages.k_29.dt_0_35.consensus.txt"
)

file_names <- sapply(file_names, function(x) file.path(current.path, x))

# OR phylogeny
mouse_or_phylo <- ape::read.tree(file_names["or_phylogeny"])

# OR class information
or_classes <- readRDS(file_names["or_classes"])

# OR gene cluster data
or_cluster <- readRDS(file_names["or_genomic_clusters"]) %>% 
  select(gene_name, cluster_name, chr, start, singleton) %>%
  arrange(singleton, chr, start) %>%
  rename(gene_cluster = cluster_name)

gene_cluster_levels <- select(or_cluster, gene_cluster, chr, singleton) %>%
  mutate(cluster_rank = as.numeric(ifelse(singleton,
                                          0,
                                          str_remove(gene_cluster, "\\(.+")))) %>%
  arrange(singleton, order(gtools::mixedorder(chr)), cluster_rank) %>%
  pull(gene_cluster) %>%
  unique()

or_cluster <- mutate(or_cluster, 
                     gene_cluster = factor(gene_cluster, levels=gene_cluster_levels))

# pooled OR data with OSN statistics
osn_metadata <- read_rds(file_names["osn_metadata"])
osn_ident <- pull(osn_metadata, olfr.ident)
names(osn_ident) <- rownames(osn_metadata)
or_data <- osn_metadata %>%
  rename(label = olfr.ident) %>%
  group_by(label) %>%
  summarise(n.osn = n(),
            transcriptome_cluster = choose.cluster(clusters)) %>%
  inner_join(or_cluster, c("label" = "gene_name")) %>%
  inner_join(or_classes, "label") 

# P-element controlled gene and neighboring cluster
pelement.controlled.genes.strong <- c("Olfr710", "Olfr6", "Olfr711", "Olfr2", "Olfr713", "Olfr714")
pelement.controlled.genes.moderate <- c("Olfr695", "Olfr705", "Olfr706", "Olfr715")
pelement.cluster <- "8(chr7)"
pelement.neighboring_cluster <- "9(chr7)"
helement.controlled.genes <- c("Olfr1507", "Olfr1508", "Olfr1509")
helement.cluster <- "5(chr14)"
another_cluster <- "1(chr19)"
cluster_labels <- c("8(chr7)" = "P cluster",
                    "9(chr7)" = "9(chr7)",
                    "1(chr19)" = "1(chr19)",
                    "5(chr14)" = "H cluster")
gene_labels <- c("Olfr710" = "P-element controlled genes (strong)", 
                 "Olfr6" = "P-element controlled genes (strong)", 
                 "Olfr711" = "P-element controlled genes (strong)", 
                 "Olfr2" = "P-element controlled genes (strong)", 
                 "Olfr713" = "P-element controlled genes (strong)", 
                 "Olfr714" = "P-element controlled genes (strong)",
                 "Olfr695" = "P-element controlled genes (moderate)", 
                 "Olfr705" = "P-element controlled genes (moderate)",
                 "Olfr706" = "P-element controlled genes (moderate)",
                 "Olfr715" = "P-element controlled genes (moderate)",
                 "Olfr1507" = "H-element controlled genes", 
                 "Olfr1508" = "H-element controlled genes",
                 "Olfr1509" = "H-element controlled genes")


or_data <- mutate(or_data, 
                  cluster_annotation = ifelse(gene_cluster %in% names(cluster_labels),
                                              cluster_labels[as.character(gene_cluster)],
                                              NA),
                  cluster_annotation= ifelse(label %in% names(gene_labels),
                                             gene_labels[label],
                                             cluster_annotation))
osn_n <- pull(or_data, n.osn)
names(osn_n) <- pull(or_data, label)

# select another gene cluster with:
# - closely related sequence
# - OSN populations detected with > 10 cells
group_by(or_data, gene_cluster) %>%
  summarise(n = n(),
            pop10 = length(which(n.osn >= 10)),
            prop = pop10/n,
            good = pop10/n == 1) %>%
  filter(n > 10) %>%
  arrange(-prop)

# fetch values of a given set of PCs
#osn_pca <- read.table(file_names["osn_PCA"])[1:7] %>%
osn_pca <- read_rds(file_names["osn_PCA"])[,1:19] %>%  
  as.data.frame(.) %>%
  rownames_to_column("cell.ident") %>%
  mutate(olfr.ident = osn_ident[cell.ident]) %>%
  filter(osn_n[olfr.ident] >= 10)

# split by OSN population
osn_pca.list <- split(osn_pca, osn_pca$olfr.ident)

# calculate centroids for OSN populations
osn_pca.centroids.matrix <- lapply(osn_pca.list, function(mat) {
  mat <- mat[,which(startsWith(colnames(mat), "PC_"))]
  
  # optionally: outlier identification using mean distance of k-nearest
  # neighbors
  find_centroid(mat) }) %>%
  do.call(what = rbind, args = .)

# calculate distances between centroids and make the pairwise_distances dataset
#osn_pca.centroids.matrix <- pivot_wider(osn_pca.centroids,
#                                        names_from = PC, values_from = centroid) %>%
#  column_to_rownames("olfr.ident")
osn_dist <- rdist(as.matrix(osn_pca.centroids.matrix))
rownames(osn_dist) <- rownames(osn_pca.centroids.matrix)
colnames(osn_dist) <- rownames(osn_pca.centroids.matrix)
pairwise_distances <- as.data.frame(osn_dist) %>%
  rownames_to_column("gene_name.1") %>%
  pivot_longer(cols = where(is.numeric), values_to = "dist.trans", names_to = "gene_name.2") %>%
  mutate(pair = unlist(mapply(
    a = gene_name.1,
    b = gene_name.2,
    SIMPLIFY = F,
    FUN = function(a,b) {
      x <- sort(c(a,b))
      return(paste(x[1], x[2], sep = ":"))
    }
  ))) 

# calculated the distance using the shared nearest neighbor method
osn_pca.snn.matrix <- compute_pairwise_shared_neighbors(osn_pca.list)
osn_pca.snn <- as.data.frame(osn_pca.snn.matrix) %>%
  rownames_to_column("gene_name.1") %>%
  pivot_longer(where(is.numeric), names_to = "gene_name.2", values_to = "dist.trans.snn") %>%
  filter(!is.na(dist.trans.snn)) %>%
  mutate(pair = unlist(mapply(
    a = gene_name.1,
    b = gene_name.2,
    SIMPLIFY = F,
    FUN = function(a,b) {
      x <- sort(c(a,b))
      return(paste(x[1], x[2], sep = ":"))
    }
  ))) %>%
  select(pair, dist.trans.snn)
pairwise_distances <- left_join(pairwise_distances, osn_pca.snn, by = "pair") %>%
  mutate(dist.trans.snn = ifelse(is.na(dist.trans.snn), 1, dist.trans.snn))

# amino acid differences
for (dist_name in c("dist.aa.miyata.trimmed", 
                    "dist.aa.miyata.intra", 
                    "dist.aa.miyata.extra",
                    "dist.aa.miyata.transmembrane")) {
  d.x <- read_tsv(file_names[dist_name],
                  col_names = c("gene_name.1", "gene_name.2", dist_name)) %>%
    mutate(pair = unlist(mapply(
      a = gene_name.1,
      b = gene_name.2,
      SIMPLIFY = F,
      FUN = function(a,b) {
        x <- sort(c(a,b))
        return(paste(x[1], x[2], sep = ":"))
      }
    ))) %>%
    select(all_of(c("pair", dist_name)))
  pairwise_distances <- left_join(pairwise_distances, d.x, by = "pair") 
}

# add metadata and genomic distance
pairwise_distances <- left_join(pairwise_distances, 
                                select(or_data, label, gene_cluster, class, 
                                       transcriptome_cluster, chr, start), 
                                by = c("gene_name.1" = "label")) %>%
  rename(class.1 = class,
         gene_cluster.1 = gene_cluster,
         transcriptome_cluster.1 = transcriptome_cluster,
         chr.1 = chr, start.1 = start) %>%
  left_join(select(or_data, label, gene_cluster, class,
                   transcriptome_cluster, chr, start), 
            by = c("gene_name.2" = "label")) %>%
  rename(class.2 = class,
         gene_cluster.2 = gene_cluster,
         transcriptome_cluster.2 = transcriptome_cluster,
         chr.2 = chr, start.2 = start) %>%
  mutate(dist.genome = ifelse(chr.1 == chr.2, abs(start.1-start.2), NA))

#saveRDS(pairwise_distances, "D:/Projects/dream-or_phylo/pairwise_distances.backup20210929.rds")
#pairwise_distances <- read_rds("D:/Projects/dream-or_phylo/pairwise_distances.backup20210929.rds")

# pairwise difference in GEP
# osn_GEP.matrix <- read_tsv(file_names["osn_GEP"]) %>%
#   as.data.frame() %>%
#   column_to_rownames("...1")
# osn_GEP.matrix.rowsums <- rowSums(osn_GEP.matrix)
# osn_GEP <- lapply(1:nrow(osn_GEP.matrix), function(i) {
#   return(osn_GEP.matrix[i,] / osn_GEP.matrix.rowsums[i])
# }) %>% 
#   do.call(rbind, .) %>%
#   as.data.frame() %>%
#   rownames_to_column("cell.ident") %>%
#   mutate(cell.ident = str_replace(cell.ident, "\\.", "-")) %>%
#   mutate(olfr.ident = osn_ident[cell.ident]) %>%
#   pivot_longer(where(is.numeric), names_to = "GEP", values_to = "loading") %>%
#   group_by(olfr.ident, GEP) %>%
#   summarise(mean.loading = mean(loading)) 
# 
# density.GEPloading <- lapply(split(osn_GEP, osn_GEP$GEP), function(d) {
#   p <- ggplot(d, aes(x = mean.loading)) +
#     stat_density() + 
#     common_layout +
#     ggtitle(paste("GEP", unique(d$GEP), sep = ""))
# })
# 
# osn_GEP.loading <- mutate(osn_GEP, GEP = paste("GEP", GEP, sep = "")) %>%
#   pivot_wider(names_from = GEP, values_from = mean.loading) %>%
#   column_to_rownames("olfr.ident")
# 
# saveRDS(osn_GEP.loading, "datasets/RData_20210929-osn_GEP.loading.rds")
# 
# combinations <- expand.grid(row.1 = 1:nrow(osn_GEP.loading), row.1 = 1:nrow(osn_GEP.loading))
# 
# pairwise_delta_GEP <- apply(X = combinations, MARGIN = 1, FUN = function(x) {
#   score.1 <- unname(unlist(osn_GEP.loading[x[1],]))
#   gene_name.1 <- rownames(osn_GEP.loading)[x[1]]
#   score.2 <- unname(unlist(osn_GEP.loading[x[2],]))
#   gene_name.2 <- rownames(osn_GEP.loading)[x[2]]
#   return(data.frame("gene_name.1" = rep(gene_name.1, ncol(osn_GEP.loading)),
#                     "gene_name.2" = rep(gene_name.2, ncol(osn_GEP.loading)),
#                     "GEP.delta" = abs(score.1 - score.2),
#                     "GEP" = colnames(osn_GEP.loading)))
# }) %>% 
#   bind_rows() %>%
#   pivot_wider(names_from = GEP, values_from = GEP.delta)
# 
# pairwise_distances <- left_join(pairwise_distances, pairwise_delta_GEP,
#                                 by = c("gene_name.1", "gene_name.2")) 
# unique comparisons
pairwise_distances.comp <- arrange(pairwise_distances, gene_name.1, gene_name.2) %>%
  
  # remove self-comparisons
  filter(gene_name.1 != gene_name.2) %>%
  
  # identify pairs belonging to the same gene cluster, transcriptome cluster or
  # OR class
  mutate(same_cluster = gene_cluster.1 == gene_cluster.2,
         same_class = class.1 == class.2,
         same_transcriptome_cluster = transcriptome_cluster.1 == transcriptome_cluster.2) %>%
  mutate(class = ifelse(same_class, class.1, NA),
         gene_cluster = ifelse(same_cluster, as.character(gene_cluster.1), NA),
         transcriptome_cluster = ifelse(same_transcriptome_cluster, transcriptome_cluster.1, NA)) %>%
  select(starts_with("GEP") | starts_with("dist.") | starts_with("same_") | all_of(c("pair", "class", "gene_cluster", "transcriptome_cluster"))) %>%
  distinct()

saveRDS(object = pairwise_distances.comp,
        file = "datasets/RData_20210929-pairwise_distances.comp-centroids-min10cells-PC1-19.rds")

saveRDS(object = or_data,
        file = "datasets/RData_20211003-or_data.rds")
