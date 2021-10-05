## Euclidean Distances ##

# load packages
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggpubr)
source("N:/Novell/__articles in progress__/DREAM transcriptomic adaptations/figures/common_layout.r")

# Set working directory
setwd(dir = "N:/Novell/___Projects___/DREAM/experiments/Scrnaseq_moe/20210526_moe_scRNAseq_cellranger_grc38irlab/Script_version_upload/")

# Load the OSN data without OR
osn.no.or.sub.integrated.sct <- readRDS("moe_noexp_scRNAseq10x_cellranger_filter_OSN_NO_OLFR_sct_integrated_clusters.rds")

##################################
#                                #        
#   ED all cells without OR      #
#                                #
##################################

# Compute the euclidean distance 
pca.coord.without.or <- osn.no.or.sub.integrated.sct@reductions$pca@cell.embeddings[,1:19]

dist.mat.without.or <- fields::rdist(x1 = pca.coord.without.or, x2 = NULL)

colnames(dist.mat.without.or) <- rownames(pca.coord.without.or)
rownames(dist.mat.without.or) <- rownames(pca.coord.without.or)

osn.idents <- osn.no.or.sub.integrated.sct$olfr.ident

dist.mat.without.or[1:10,1:10]

dist.mat.without.or[lower.tri(dist.mat.without.or, diag = TRUE)] <- NA

dist.mat.without.or <- melt(data = dist.mat.without.or,
                            na.rm = TRUE,
                            varnames = c("OSN_1", "OSN_2"),
                            value.name = "euclidean_distance")

dist.mat.without.or <- data.table::data.table(dist.mat.without.or)
dist.mat.without.or$OSN_1 <- as.character(x = dist.mat.without.or$OSN_1)
dist.mat.without.or$OSN_2 <- as.character(x = dist.mat.without.or$OSN_2)

dist.mat.without.or$OSN_1_ident <- osn.idents[dist.mat.without.or$OSN_1]
dist.mat.without.or$OSN_2_ident <- osn.idents[dist.mat.without.or$OSN_2]

dist.mat.without.or$osn_pairs <- ifelse(test = dist.mat.without.or$OSN_1_ident == dist.mat.without.or$OSN_2_ident,
                                        yes = "same OSN pop.",
                                        no = "different OSN pop.")

dist.mat.without.or$OSN_1_class <- osn.no.or.sub.integrated.sct$class_first[dist.mat.without.or$OSN_1]
dist.mat.without.or$OSN_2_class <- osn.no.or.sub.integrated.sct$class_first[dist.mat.without.or$OSN_2]

dist.mat.without.or[,
                    class_pairs := paste(sort(x = c(OSN_1_class, OSN_2_class)),
                                         collapse = "_"),
                    by = seq_len(length.out = nrow(x = dist.mat.without.or))]

saveRDS(dist.mat.without.or,"euclidean_dist_matrix_without_or_20210920.rds")


# Permutations

# Create a matrix containing the first 19 PC.
pca.coords <- osn.no.or.sub.integrated.sct@reductions$pca@cell.embeddings[,1:19]
osn.idents <- osn.no.or.sub.integrated.sct$olfr.ident

# Determine a number of permutations
permutations <- paste("permutation_",
                      1:1000,
                      sep = "")
names(x = permutations) <- permutations


permuted.distances <- lapply(X = permutations,
                             FUN = function(permutation,
                                            pca.coords,
                                            osn.idents) {
                               print(x = permutation)
                               # Compute euclidean distances of PCA eigenvectors after permuting OSN identities
                               rownames(x = pca.coords) <- sample(x = rownames(x = pca.coords),
                                                                  size = nrow(x = pca.coords),
                                                                  replace = FALSE)
                               
                               dist.mat <- fields::rdist(x1 = pca.coords,
                                                         x2 = NULL)
                               
                               colnames(dist.mat) <- rownames(pca.coords)
                               rownames(dist.mat) <- rownames(pca.coords)
                               
                               dist.mat[lower.tri(dist.mat, diag = TRUE)] <- NA
                               
                               dist.mat <- melt(data = dist.mat,
                                                na.rm = TRUE,
                                                varnames = c("OSN_1", "OSN_2"),
                                                value.name = "euclidean_distance")
                               
                               dist.mat <- data.table::data.table(dist.mat)
                               dist.mat$OSN_1 <- as.character(x = dist.mat$OSN_1)
                               dist.mat$OSN_2 <- as.character(x = dist.mat$OSN_2)
                               
                               dist.mat$OSN_1_ident <- osn.idents[dist.mat$OSN_1]
                               dist.mat$OSN_2_ident <- osn.idents[dist.mat$OSN_2]
                               
                               dist.mat$osn_pairs <- ifelse(test = dist.mat$OSN_1_ident == dist.mat$OSN_2_ident,
                                                            yes = "same OSN pop.",
                                                            no = "different OSN pop.")
                               
                               dist.mat <- dist.mat[dist.mat$osn_pairs == "same OSN pop.",]
                               
                               return(dist.mat)
                             },
                             pca.coords = pca.coords,
                             osn.idents = osn.idents)

#%%% Checkpoint %%%#
saveRDS(permuted.distances,file = "permuted_distances_without_or.rds")

##################################
#                                #        
#           Plots 1              #
#                                #
##################################

# Read the data
permuted.distances <- read_rds("permuted_distances_without_or.rds")
dist.mat.without.or <- read_rds("euclidean_dist_matrix_without_or.rds")
dist.mat.without.or <- dist.mat.without.or[,1:6]

# Select the first 100 permutations
permuted.distances <- permuted.distances[1:1000]

# Combine the permutated matrix
permuted.distances <- dplyr::bind_rows(permuted.distances)
permuted.distances$osn_pairs <- "permutation"

# Combine the ditance matrix to the permutated distance matrix
dist.mat.without.or.per <- dplyr::bind_rows(dist.mat.without.or, permuted.distances)
dist.mat.without.or.per$osn_pairs <- factor(dist.mat.without.or.per$osn_pairs,
                                            levels = c("same OSN pop.",
                                                       "different OSN pop.",
                                                       "permutation"))

euclidean.dist.plot.1 <- ggplot(data = dist.mat.without.or.per,
                                mapping = aes(y = euclidean_distance,
                                              x = osn_pairs,
                                              fill = osn_pairs,
                                              colour = osn_pairs)) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..)) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..)) +
  stat_summary(geom = "point", fun = "median")+
  geom_violin(alpha = 0.5,scale = "width") +
  scale_fill_manual(values = c("same OSN pop." = "deeppink1",
                               "different OSN pop." = "dodgerblue1",
                               "permutation" = "gray"),
                    limits = c("same OSN pop.", "different OSN pop.", "permutation")) +
  scale_colour_manual(values = c("same OSN pop." = "deeppink1",
                                 "different OSN pop." = "dodgerblue1",
                                 "permutation" = "gray"),
                      limits = c("same OSN pop.", "different OSN pop.", "permutation")) +
  labs(y = "Euclidean distance") + 
  ylim(0,70)+
  common_layout +
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        legend.position = c(0,20),
        legend.background = element_rect(fill='transparent'))

ggsave(euclidean.dist.plot.1,
       filename = "eucl_dist_without_or.pdf",
       width = 2.5,
       height = 2.5)



# Statistical test
same.osn.dist <- dist.mat.without.or.per$euclidean_distance[dist.mat.without.or.per$osn_pairs == "same OSN pop."]
different.osn.dist <- dist.mat.without.or.per$euclidean_distance[dist.mat.without.or.per$osn_pairs == "different OSN pop."]
permutation <- dist.mat.without.or.per$euclidean_distance[dist.mat.without.or.per$osn_pairs == "permutation"]

same_diff <- wilcox.test(same.osn.dist,different.osn.dist)

same_per <- wilcox.test(same.osn.dist,permutation)

p.adjust(p = c(same_diff$p.value,
               same_per$p.value),
         method = "bonferroni")


##################################
#                                #        
#           ED per cluster       #
#                                #
##################################

# Eucledian distance Multiple clusters
pca.coord <- osn.no.or.sub.integrated.sct@reductions$pca@cell.embeddings[,1:19]
pca.coord.split <- as.data.frame(cbind(pca.coord,
                                       as.character(osn.no.or.sub.integrated.sct$clusters)))

pca.coord.split <- split(pca.coord.split,pca.coord.split$V20)  

osn.idents <- osn.no.or.sub.integrated.sct$olfr.ident

eucl.dist.without.or.cluster <- lapply(X = pca.coord.split,FUN = function(pca.c){
  pca.c <- pca.c[,-c(20)]
  dist.mat <- fields::rdist(x1 = pca.c, x2 = NULL)
  colnames(dist.mat) <- rownames(pca.c)
  rownames(dist.mat) <- rownames(pca.c)
  dist.mat[lower.tri(dist.mat, diag = TRUE)] <- NA
  dist.mat <- melt(data = dist.mat,
                   na.rm = TRUE,
                   varnames = c("OSN_1", "OSN_2"),
                   value.name = "euclidean_distance")
  dist.mat <- data.table::data.table(dist.mat)
  dist.mat$OSN_1 <- as.character(x = dist.mat$OSN_1)
  dist.mat$OSN_2 <- as.character(x = dist.mat$OSN_2)
  dist.mat$OSN_1_ident <- osn.idents[dist.mat$OSN_1]
  dist.mat$OSN_2_ident <- osn.idents[dist.mat$OSN_2]
  dist.mat$osn_pairs <- ifelse(test = dist.mat$OSN_1_ident == dist.mat$OSN_2_ident,
                               yes = "same OSN pop.",
                               no = "different OSN pop.")
  dist.mat
}
)

saveRDS(eucl.dist.cluster,"euclidean_dist_matrix_clusters_without_or.rds")


# Permutations per cluster
pca.coord <- osn.no.or.sub.integrated.sct@reductions$pca@cell.embeddings[,1:19]
pca.coord.split <- as.data.frame(cbind(pca.coord,
                                       as.character(osn.no.or.sub.integrated.sct$clusters),
                                       as.character(osn.no.or.sub.integrated.sct$olfr.ident)))

pca.coord.split <- split(pca.coord.split,pca.coord.split$V20) 

permutations_clusters <- lapply(pca.coord.split, FUN = function(pca.coord.cluster){
  # Create a matrix containing the first 16 PC.
  
  osn.idents <- pca.coord.cluster["V21"]
  pca.coord.cluster <- pca.coord.cluster[,-c(20,21)]
  
  
  # Determine a number of permutations
  permutations <- paste("permutation_",
                        1:1000,
                        sep = "")
  names(x = permutations) <- permutations
  
  permuted.distances.cluster <- lapply(X = permutations,
                                       FUN = function(permutation) {
                                         print(x = permutation)
                                         # Compute euclidean distances of PCA eigenvectors after permuting OSN identities
                                         rownames(x = pca.coord.cluster) <- sample(x = rownames(x = pca.coord.cluster),
                                                                                   size = nrow(x = pca.coord.cluster),
                                                                                   replace = FALSE)
                                         
                                         pca.coord.cluster
                                         dist.mat <- fields::rdist(x1 = pca.coord.cluster,
                                                                   x2 = NULL)
                                         
                                         colnames(dist.mat) <- rownames(pca.coord.cluster)
                                         rownames(dist.mat) <- rownames(pca.coord.cluster)
                                         
                                         dist.mat[lower.tri(dist.mat, diag = TRUE)] <- NA
                                         
                                         dist.mat <- melt(data = dist.mat,
                                                          na.rm = TRUE,
                                                          varnames = c("OSN_1", "OSN_2"),
                                                          value.name = "euclidean_distance")
                                         
                                         dist.mat <- data.table::data.table(dist.mat)
                                         dist.mat$OSN_1 <- as.character(x = dist.mat$OSN_1)
                                         dist.mat$OSN_2 <- as.character(x = dist.mat$OSN_2)
                                         
                                         dist.mat$OSN_1_ident <- osn.idents[dist.mat$OSN_1,]
                                         dist.mat$OSN_2_ident <- osn.idents[dist.mat$OSN_2,]
                                         
                                         dist.mat$osn_pairs <- ifelse(test = dist.mat$OSN_1_ident == dist.mat$OSN_2_ident,
                                                                      yes = "same OSN pop.",
                                                                      no = "different OSN pop.")
                                         
                                         dist.mat <- dist.mat[dist.mat$osn_pairs == "same OSN pop.",]
                                         
                                         return(dist.mat)
                                       })
}
)


#%%% Checkpoint %%%#
saveRDS(permutations_clusters,file = "permutation_clusters_without_or.rds")


##################################
#                                #        
#           Plots 2              #
#                                #
##################################

# Read data
permutations_clusters <- readRDS("permutation_clusters_without_or.rds")
dist.mat.clusters <- readRDS("euclidean_dist_matrix_clusters_without_or.rds")


# cluster name
clusters <- names(permutations_clusters)

euclide.dist.plots.clusters <- lapply(clusters, FUN = function(cluster){ 
  
  per.clu <- permutations_clusters[[cluster]]
  per.clu <- per.clu[1:1000]
  dist.mat.clu <- dist.mat.clusters[[cluster]]
  # Combine the permutated matrix
  
  per.clu <- dplyr::bind_rows(per.clu)
  per.clu$osn_pairs <- "permutation"
  
  # Combine the ditance matrix to the permutated distance matrix
  dist.mat.full <- dplyr::bind_rows(dist.mat.clu, per.clu)
  
  dist.mat.full$osn_pairs <- factor(dist.mat.full$osn_pairs,
                                    levels = c("same OSN pop.",
                                               "different OSN pop.",
                                               "permutation"))
  dist.mat.full$cluster <- cluster
  
  p <- ggplot(data = dist.mat.full,
              mapping = aes(y = euclidean_distance,
                            x = osn_pairs,
                            fill = cluster,
                            colour = cluster)) + 
    geom_violin(alpha = 0.5,scale = "width") +
    stat_mean(geom="segment",alpha =0.8, mapping = aes(xend=..x.. + 0.5, yend=..y..)) +
    stat_mean(geom="segment",alpha =0.8, mapping = aes(xend=..x.. - 0.5, yend=..y..)) +
    stat_summary(geom = "point", fun = "median", color= "black")+
    scale_fill_manual(values = clusters.fill.osn)+
    scale_color_manual(values = clusters.fill.osn)+
    
    labs(y = "Euclidean distance") + 
    common_layout +
    theme(legend.title = element_blank(),
          axis.title = element_blank(),
          legend.position = c(0,20),
          legend.background = element_rect(fill='transparent'))+
    ggtitle(cluster)+
    ylim(0,70)
  
  p
}
)


ggsave(plot = euclide.dist.plots.clusters[[10]],
       filename =   "euclidean_dist_clusters_ventral_20210921.pdf",
       width = 2.5,
       height = 2.5)


euclide.dist.plots.clusters <- wrap_plots(euclide.dist.plots.clusters)

ggsave(plot = euclide.dist.plots.clusters,
       filename =   "euclidean_dist_clusters_20210921.pdf",
       width = 10,
       height = 9)





##################################
#                                #        
#   ED all cells with OR         #
#                                #
##################################
osn.sub.integrated.sct <- readRDS("moe_noexp_scRNAseq10x_cellranger_filter_OSN_WITH_OLFR_sct_integrated_clusters.rds")

# Compute the euclidean distance 
pca.coord.with.or <- osn.sub.integrated.sct@reductions$pca@cell.embeddings[,1:18]

dist.mat.with.or <- fields::rdist(x1 = pca.coord.with.or, x2 = NULL)

colnames(dist.mat.with.or) <- rownames(pca.coord.with.or)
rownames(dist.mat.with.or) <- rownames(pca.coord.with.or)

osn.idents <- osn.sub.integrated.sct$olfr.ident

dist.mat.with.or[1:10,1:10]

dist.mat.with.or[lower.tri(dist.mat.with.or, diag = TRUE)] <- NA

dist.mat.with.or <- melt(data = dist.mat.with.or,
                         na.rm = TRUE,
                         varnames = c("OSN_1", "OSN_2"),
                         value.name = "euclidean_distance")

dist.mat.with.or <- data.table::data.table(dist.mat.with.or)
dist.mat.with.or$OSN_1 <- as.character(x = dist.mat.with.or$OSN_1)
dist.mat.with.or$OSN_2 <- as.character(x = dist.mat.with.or$OSN_2)

dist.mat.with.or$OSN_1_ident <- osn.idents[dist.mat.with.or$OSN_1]
dist.mat.with.or$OSN_2_ident <- osn.idents[dist.mat.with.or$OSN_2]

dist.mat.with.or$osn_pairs <- ifelse(test = dist.mat.with.or$OSN_1_ident == dist.mat.with.or$OSN_2_ident,
                                     yes = "same OSN pop.",
                                     no = "different OSN pop.")

dist.mat.with.or$osn_pairs <- factor(dist.mat.with.or$osn_pairs,
                                     levels = c("same OSN pop.","different OSN pop."))


#%%% CHECKPOINT %%%#

# Statistical test
same.osn.dist <- dist.mat.with.or$euclidean_distance[dist.mat.with.or$osn_pairs == "same OSN pop."]
different.osn.dist <- dist.mat.with.or$euclidean_distance[dist.mat.with.or$osn_pairs == "different OSN pop."]

same_diff <- wilcox.test(same.osn.dist,different.osn.dist)


# Plots

euclidean.dist.plot.2 <- ggplot(data = dist.mat.with.or,
                                mapping = aes(y = euclidean_distance,
                                              x = osn_pairs,
                                              fill = osn_pairs,
                                              colour = osn_pairs)) + 
  geom_violin(alpha = 0.5,scale = "width") +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..)) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..)) +
  stat_summary(geom = "point", fun = "median")+
  scale_fill_manual(values = c("same OSN pop." = "deeppink1",
                               "different OSN pop." = "dodgerblue1"),
                    limits = c("same OSN pop.", "different OSN pop.")) +
  scale_colour_manual(values = c("same OSN pop." = "deeppink1",
                                 "different OSN pop." = "dodgerblue1"),
                      limits = c("same OSN pop.", "different OSN pop.")) +
  labs(y = "Euclidean distance") + 
  ylim(0,70)+
  common_layout +
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        legend.position = c(0,20),
        legend.background = element_rect(fill='transparent'))


ggsave(euclidean.dist.plot.2,
       filename = "euclidean.dist.all_WITH_OR.pdf",
       width = 1.7,
       height = 2.5)


# Eucledian distance per clusters
pca.coord <- osn.sub.integrated.sct@reductions$pca@cell.embeddings[,1:18]
pca.coord.split <- as.data.frame(cbind(pca.coord,
                                       as.character(osn.sub.integrated.sct$clusters)))

pca.coord.split <- split(pca.coord.split,pca.coord.split$V19)  

osn.idents <- osn.sub.integrated.sct$olfr.ident

eucl.dist.cluster.with_OR <- lapply(X = pca.coord.split,FUN = function(pca.c){
  pca.c <- pca.c[,-c(19)]
  dist.mat <- fields::rdist(x1 = pca.c, x2 = NULL)
  colnames(dist.mat) <- rownames(pca.c)
  rownames(dist.mat) <- rownames(pca.c)
  dist.mat[lower.tri(dist.mat, diag = TRUE)] <- NA
  dist.mat <- melt(data = dist.mat,
                   na.rm = TRUE,
                   varnames = c("OSN_1", "OSN_2"),
                   value.name = "euclidean_distance")
  dist.mat <- data.table::data.table(dist.mat)
  dist.mat$OSN_1 <- as.character(x = dist.mat$OSN_1)
  dist.mat$OSN_2 <- as.character(x = dist.mat$OSN_2)
  dist.mat$OSN_1_ident <- osn.idents[dist.mat$OSN_1]
  dist.mat$OSN_2_ident <- osn.idents[dist.mat$OSN_2]
  dist.mat$osn_pairs <- ifelse(test = dist.mat$OSN_1_ident == dist.mat$OSN_2_ident,
                               yes = "same OSN pop.",
                               no = "different OSN pop.")
  dist.mat
}
)


# Statistical test
same.osn.dist <- dist.mat.with.or$euclidean_distance[dist.mat.with.or$osn_pairs == "same OSN pop."]
different.osn.dist <- dist.mat.with.or$euclidean_distance[dist.mat.with.or$osn_pairs == "different OSN pop."]

same_diff <- wilcox.test(same.osn.dist,different.osn.dist)



# Ventral cluster plot
eucl.dist.cluster.with_OR$Ventral$osn_pairs <- factor(eucl.dist.cluster.with_OR$Ventral$osn_pairs,
                                                      levels = c("same OSN pop.","different OSN pop."))

eucl.dist.cluster.with_OR.plot <- ggplot(data = eucl.dist.cluster.with_OR$Ventral,
                                         mapping = aes(y = euclidean_distance,
                                                       x = osn_pairs,
                                                       fill = osn_pairs,
                                                       colour = osn_pairs)) + 
  geom_violin(alpha = 0.5,scale = "width",color="black",fill="black", alpha=0.5) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..)) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..)) +
  stat_summary(geom = "point", fun = "median")+
  labs(y = "Euclidean distance") + 
  ylim(0,70)+
  common_layout +
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        legend.position = c(0,20),
        legend.background = element_rect(fill='transparent'))

ggsave(eucl.dist.cluster.with_OR.plot,
       filename = "euclidean.dist.ventral_WITH_OR.pdf",
       width = 1.7,
       height = 2.5)
