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
current.path <- "N:/Novell/__articles in progress__/DREAM transcriptomic adaptations/"
setwd(current.path)

### FUNCTIONS
source("figures/R scripts/Figure2-functions.r")
source("figures/R scripts/common_layout.r")

### DATASETS
file_names <- c(
  
  # coordinates of individual functional OR genes and their respective genomic
  # cluster
  "or_genomic_clusters" = "DREAM_datasets/Mus_musculus.mm10.functional_OR_gene_clusters.20210914.rds",
  
  # functional OR phylogeny
  "or_phylogeny" = "DREAM_datasets/Mus_musculus.GRCm38.olfrs.clustalo_full.trim.phyml.rooted.20210916.txt",
  
  # All subtrees including OR genes detected in the scRNA-seq
  "subtrees" = "DREAM_datasets/Mus_musculus.GRCm38.olfrs.clustalo_full.trim.phyml.rooted.subtrees.20210916.txt",
  
  # functional OR class I / class II assignation
  "or_classes" = "DREAM_datasets/Mus_musculus.mm10.olfr_classes.20210907.rds",
  
  # different metrics of distances between OSN populations and their
  # corresponding OR
  "pairwise_distances" = "DREAM_datasets/Mus_musculus.mm10.functional_OR_pairwise_distances.20210810.rds",
  
  # OSN transcriptome metadata, after removing cells expressing more than one OR
  "osn_metadata" = "DREAM_datasets/Mus_musculus.mm10.filtered_mature_OSN_metadata.20220209.rds",
  
  # OSN integrated count matrix
  "osn_count_matrix" = "DREAM_datasets/Mus_musculus.mm10.filtered_mature_OSN_integrated_count_matrix.20210910.rds",
  
  # OSN PCA
  "osn_PCA" = "DREAM_datasets/Mus_musculus.mm10.filtered_mature_OSN_PCA.20220209.rds",
  
  # amino-acid distances
  "dist.aa.levenshtein.trimmed" = "DREAM_datasets/Mus_musculus.GRCm38.olfrs_and_taars.clustalo_full-trimmed-levenshtein.20210920.txt",
  "dist.aa.levenshtein.tm" = "DREAM_datasets/Mus_musculus.GRCm38.olfrs_and_taars.clustalo_full-transmembrane-levenshtein.20210920.txt",
  "dist.aa.levenshtein.cons" = "DREAM_datasets/Mus_musculus.GRCm38.olfrs_and_taars.clustalo_full-conserved-levenshtein.20210920.txt",
  "dist.aa.grantham.trimmed" = "DREAM_datasets/Mus_musculus.GRCm38.olfrs_and_taars.clustalo_full-trimmed-grantham.20210922.txt",
  "dist.aa.miyata.trimmed" = "DREAM_datasets/Mus_musculus.GRCm38.olfrs_and_taars.clustalo_full-trimmed-withgaps-miyata.20210924.txt",
  "dist.aa.miyata.transmembrane" = "DREAM_datasets/Mus_musculus.GRCm38.olfrs_and_taars.clustalo_full-transmembrane-miyata-withgaps.20210927.txt",
  "dist.aa.miyata.intra" = "DREAM_datasets/Mus_musculus.GRCm38.olfrs_and_taars.clustalo_full-intra-miyata-withgaps.202109257.txt",
  "dist.aa.miyata.extra" = "DREAM_datasets/Mus_musculus.GRCm38.olfrs_and_taars.clustalo_full-extra-miyata-withgaps.20210927.txt",
  
  # GEP loading matrix
  "osn_GEP" = "DREAM_datasets/cNMF_GEPs_no_olfr.usages.k_29.dt_0_35.consensus.txt"
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
for (i in 1:15) {
  osn_pca <- read_rds(file_names["osn_PCA"])[,1:i] %>%  
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
  
  fname <- paste("DREAM_datasets/RData_20220209-pairwise_distances-simple-centroids-min10cells-PC1-", i, ".rds", sep = "")
  saveRDS(pairwise_distances, fname)
}

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
        file = "DREAM_datasets/RData_20220209-pairwise_distances.comp-centroids-min10cells-PC1-15.rds")

saveRDS(object = or_data,
        file = "DREAM_datasets/RData_20220209-or_data.rds")
