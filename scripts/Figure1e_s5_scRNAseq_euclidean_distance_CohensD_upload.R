## Euclidean distances ##

# This script compute euclidean distances and generates plots for Figure 1 and 
# supplementary

# Load packages
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(effsize)
source("common_layout.r")

# Load the OSN data without OR
osn.no.or.sub.integrated.sct <- readRDS("moe_noexp_scRNAseq10x_cellranger_filter_OSN_NO_OLFR_sct_integrated_clusters.rds")


##################################
#                                #        
#   ED all cells without OR      #
#                                #
##################################

# Compute the euclidean distance 
pca.coord.without.or <- osn.no.or.sub.integrated.sct@reductions$pca@cell.embeddings[,1:15]

dist.mat.without.or <- fields::rdist(x1 = pca.coord.without.or, x2 = NULL)

colnames(dist.mat.without.or) <- rownames(pca.coord.without.or)
rownames(dist.mat.without.or) <- rownames(pca.coord.without.or)

osn.idents <- osn.no.or.sub.integrated.sct$olfr.ident

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

saveRDS(dist.mat.without.or,"eucl_dist_matrix_without_or.rds")


## Permutations

# Create a matrix containing the first 15 PC.
pca.coords <- osn.no.or.sub.integrated.sct@reductions$pca@cell.embeddings[,1:15]
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

saveRDS(permuted.distances,file = "eucl_dist_permuations_without_or.rds")

# Read the data
permuted.distances <- read_rds("eucl_dist_permuations_without_or.rds")
dist.mat.without.or <- read_rds("eucl_dist_matrix_without_or.rds")

# Cohen's D
same.osn.dist <- dist.mat.without.or$euclidean_distance[dist.mat.without.or$osn_pairs == "same OSN pop."]
different.osn.dist <- unlist(dist.mat.without.or$euclidean_distance[dist.mat.without.or$osn_pairs == "different OSN pop."])
permutation <- lapply(permuted.distances, '[[', 3)

intra.inter <- cohen.d(same.osn.dist,sample(different.osn.dist))


density.cohend <- lapply(c(1:1000), FUN = function(x){
  
  p <- unlist(permutation[x])
  
  sp <- cohen.d(same.osn.dist,p)
  sp <- sp$estimate
  
  dp <- cohen.d(different.osn.dist,p)
  dp <- dp$estimate
  
  c(sp,dp)
  
})

## Checkpoint ##

saveRDS(density.cohend,file = "cohendsD_density_without_or.rds")

density.cohend <- readRDS("cohendsD_density_without_or.rds")

density.cohend.dframe <- data.frame(t(data.frame(density.cohend)))
a <- data.frame(dist = density.cohend.dframe$X1, type = "intra_perm")
b <- data.frame(dist = density.cohend.dframe$X2, type = "inter_perm")

density.cohend.dframe <- rbind(a,b)
density.cohend.dframe$dist <- density.cohend.dframe$dist*c(-1)

stat.test <- ks.test(density.cohend.dframe$dist[density.cohend.dframe$type == "intra_perm"],
                     density.cohend.dframe$dist[density.cohend.dframe$type == "inter_perm"])

cohend_all_clusters <- ggplot(density.cohend.dframe,aes(y = dist,
                                                        x= type))+ 
  geom_violin(scale = "width")+ 
  geom_hline(yintercept = -1*intra.inter$estimate,
             linetype="dotted",color="red")+
  ylim(-0.3,1.7)+
  ylab("Cohen's D")+
  annotate(geom = "text",label ="***",y=1.7,x=1.5)+
  common_layout+
  theme(axis.title.x = element_blank())


# Combine the permuted matrix
permuted.distances <- dplyr::bind_rows(permuted.distances)
permuted.distances$osn_pairs <- "permutation"


# Combine the distance matrix to the permuted distance matrix
dist.mat.without.or.per <- dplyr::bind_rows(dist.mat.without.or, permuted.distances)
dist.mat.without.or.per$osn_pairs <- factor(dist.mat.without.or.per$osn_pairs,
                                            levels = c("same OSN pop.",
                                                       "different OSN pop.",
                                                       "permutation"))

levels(dist.mat.without.or.per$osn_pairs) <- c("intra.","inter.","perm.")

euclidean.dist.noOR.plot <- ggplot(data = dist.mat.without.or.per,
                                   mapping = aes(y = euclidean_distance,
                                                 x = osn_pairs,
                                                 fill = osn_pairs,
                                                 group = osn_pairs)) +
  geom_violin(alpha = 0.5,scale = "width",color="gray50",fill = "gray50") +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..)) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..)) +
  stat_summary(geom = "point", fun = "median")+
  labs(y = "eucl. dist.") + 
  ylim(0,70)+
  common_layout +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0,20),
        legend.background = element_rect(fill='transparent'))


##################################
#                                #        
#   ED all cells with OR         #
#                                #
##################################

osn.sub.integrated.sct <- readRDS("moe_noexp_scRNAseq10x_cellranger_filter_OSN_WITH_OLFR_sct_integrated_clusters.rds")

# Compute the euclidean distance 
pca.coord.with.or <- osn.sub.integrated.sct@reductions$pca@cell.embeddings[,1:14]

dist.mat.with.or <- fields::rdist(x1 = pca.coord.with.or, x2 = NULL)

colnames(dist.mat.with.or) <- rownames(pca.coord.with.or)
rownames(dist.mat.with.or) <- rownames(pca.coord.with.or)

osn.idents <- osn.sub.integrated.sct$olfr.ident

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



# Cohen's D
same.osn.dist <- dist.mat.with.or$euclidean_distance[dist.mat.with.or$osn_pairs == "same OSN pop."]
different.osn.dist <- dist.mat.with.or$euclidean_distance[dist.mat.with.or$osn_pairs == "different OSN pop."]

same_diff <- cohen.d(same.osn.dist,different.osn.dist)


levels(dist.mat.with.or$osn_pairs) <- c("intra.","inter.")

# Plots
euclidean.dist.or.plot <- ggplot(data = dist.mat.with.or,
                                 mapping = aes(y = euclidean_distance,
                                               x = osn_pairs,
                                               fill = osn_pairs,
                                               group = osn_pairs)) + 
  geom_violin(alpha = 0.5,scale = "width",color="gray50",fill = "gray50") +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..)) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..)) +
  stat_summary(geom = "point", fun = "median")+
  labs(y = "eucl. dist.") + 
  ylim(0,70)+
  common_layout +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0,20),
        legend.background = element_rect(fill='transparent'))



euclidean.dist.plot.all <- wrap_plots(A = euclidean.dist.noOR.plot,
                                      B = cohend_all_clusters,
                                      C = euclidean.dist.or.plot,
                                      design = "AABC")

ggsave(euclidean.dist.plot.all,
       filename = "fig1_E_eucldist_all.pdf",
       width = 6,
       height = 2)

##################################
#                                #        
#           ED per cluster       #
#                                #
##################################

# Eucledian distance Multiple clusters
pca.coord <- osn.no.or.sub.integrated.sct@reductions$pca@cell.embeddings[,1:15]

pca.coord.split <- as.data.frame(cbind(pca.coord,
                                       as.character(osn.no.or.sub.integrated.sct$clusters)))

pca.coord.split <- split(pca.coord.split,pca.coord.split$V16)  

osn.idents <- osn.no.or.sub.integrated.sct$olfr.ident

eucl.dist.without.or.cluster <- lapply(X = pca.coord.split,FUN = function(pca.c){
  pca.c <- pca.c[,-c(16)]
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

saveRDS(eucl.dist.without.or.cluster,"euclidean_dist_matrix_clusters_without_or.rds")


# Permutations per cluster
pca.coord <- osn.no.or.sub.integrated.sct@reductions$pca@cell.embeddings[,1:15]
pca.coord.split <- as.data.frame(cbind(pca.coord,
                                       as.character(osn.no.or.sub.integrated.sct$clusters),
                                       as.character(osn.no.or.sub.integrated.sct$olfr.ident)))

pca.coord.split <- split(pca.coord.split,pca.coord.split$V16) 

permutations_clusters <- lapply(pca.coord.split, FUN = function(pca.coord.cluster){
  # Create a matrix containing the first 16 PC.
  
  osn.idents <- pca.coord.cluster["V17"]
  pca.coord.cluster <- pca.coord.cluster[,-c(16,17)]
  
  
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

# Read data
permutations_clusters <- readRDS("permutation_clusters_without_or.rds")
dist.mat.clusters <- readRDS("euclidean_dist_matrix_clusters_without_or.rds")


# cluster name
clusters <- names(permutations_clusters)
clusters <- clusters[c(9,1,5,3,7,10,2,6,4,8)]

osn_noOR_cohensd_cluster <- lapply(clusters,FUN =  function(cluster){
  data.cluster <- dist.mat.clusters[[cluster]]
  intra <- data.cluster[data.cluster$osn_pairs =="different OSN pop."]$euclidean_distance
  inter <- data.cluster[data.cluster$osn_pairs =="same OSN pop."]$euclidean_distance
  
  cohend.intra.inter <- cohen.d(intra,inter)
  
  data.perm <- permutations_clusters[[cluster]]
  
  cohen <- lapply(c(1:1000), FUN = function(perm){
    
    permut <- data.perm[[perm]]$euclidean_distance
    
    inter.perm <- cohen.d(inter,permut)
    intra.perm <- cohen.d(intra,permut)
    
    coh.tra <- c(intra.perm$estimate)
    coh.ter <- c(inter.perm$estimate)
    
    return(c(coh.tra,coh.ter))
  }) 
  
  density.coeh <- data.frame(t(data.frame(cohen)))
  a <- data.frame(dist = density.coeh$X1, type = "intra_perm")
  b <- data.frame(dist = density.coeh$X2, type = "inter_perm")
  
  density.coeh <- rbind(a,b)
  density.coeh$dist <- density.coeh$dist*c(-1)
  density.coeh$cluster <- cluster
  density.coeh$intra.inter.c <- cohend.intra.inter$estimate
  density.coeh
})


saveRDS(osn_noOR_cohensd_cluster, file = "osn_noOR_cohensd_cluster_1000perm.rds")

osn_noOR_cohensd_cluster <- readRDS("osn_noOR_cohensd_cluster_1000perm.rds")

# Plot cohens D per cluster

osn_noOR_cohensd_cluster.plot <- lapply(osn_noOR_cohensd_cluster, FUN = function(x){
  
  csd <- x$intra.inter.c[1]
  
  x$type <- factor(x$type,levels = "intra_perm","inter_perm")
  
  a <- ggplot(data = x,
              mapping = aes(y = dist,
                            x = type,
                            fill = cluster,
                            color= cluster))+ 
    scale_fill_manual(values = clusters.fill.osn)+
    scale_color_manual(values = clusters.fill.osn)+
    geom_violin(alpha=0.5) + 
    geom_hline(yintercept = csd)+
    ylab("Cohen's D")+
    ylim(-0.3,1.3)+
    common_layout+
    theme(legend.position = "none")
  
})

coehdsd <- lapply(osn_noOR_cohensd_cluster, FUN = function(x){
  
  csd <- x$intra.inter.c[1]
  names(csd) <- x$cluster[1]
  csd
  
})

coehdsd.cluster <- data.frame(c.d. = unlist(coehdsd))
coehdsd.cluster$cluster <- rownames(coehdsd.cluster)

coehdsd.cluster.boxplot <- ggplot(coehdsd.cluster,aes(x = "cluster",
                                                      y = c.d.,
                                                        color=cluster))+ 
  geom_boxplot(color="black")+
  geom_point()+
  geom_hline(yintercept = 1.438092)+
  ylim(-0.3,1.6)+
  common_layout

ggsave(coehdsd.cluster.boxplot,filename = "coehnsd_bocplot.pdf",
       height = 3,width = 3)

Kstest.cluster <- lapply(osn_noOR_cohensd_cluster,  function(x){
  
  intra <- x[x$type =="intra_perm",]$dist
  inter <- x[x$type =="inter_perm",]$dist
  
  ks <- ks.test(intra,inter)
  ks
  
  
}
) 

p.adjust(method = "bonferroni",rep(unlist(Kstest.cluster),10))

euclide.dist.plots.clusters <- lapply(clusters, FUN = function(cluster){ 
  
  per.clu <- permutations_clusters[[cluster]]
  per.clu <- per.clu[1:100]
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
  
  levels(dist.mat.full$osn_pairs) <- c("intra.","inter.","perm.")
  
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

a <- lapply(c(1:10), FUN = function(x){
  plot <- wrap_plots(A=euclide.dist.plots.clusters[[x]],B=osn_noOR_cohensd_cluster.plot[[x]],
                     design = "AAB") 
}
)

euclide.dist.plots.clusters.wrap <- wrap_plots(c(a[1],
                                                 a[6],
                                                 a[2],
                                                 a[7],
                                                 a[3],
                                                 a[8],
                                                 a[4],
                                                 a[9],
                                                 a[5],
                                                 a[10]),ncol = 2)


ggsave(plot = euclide.dist.plots.clusters.wrap,
       filename =   "euclidean_dist_clusters_4.pdf",
       width = 5,
       height = 8)











##################################
#                                #        
#         Ventral cluster        #
#                                #
##################################


# Read data
permutations.clusters <- readRDS("permutation_clusters_without_or.rds")
dist.mat.clusters <- readRDS("euclidean_dist_matrix_clusters_without_or.rds")

permutations.clusters.v <- permutations.clusters$Ventral
dist.mat.clusters.v <- dist.mat.clusters$Ventral

same.osn.dist.v <- dist.mat.clusters.v$euclidean_distance[dist.mat.clusters.v$osn_pairs == "same OSN pop."]
different.osn.dist.v <- unlist(dist.mat.clusters.v$euclidean_distance[dist.mat.clusters.v$osn_pairs == "different OSN pop."])
permutation.v <- lapply(permutations.clusters.v, '[[', 3)


intra.inter.v <- cohen.d(same.osn.dist.v,sample(different.osn.dist.v))

density.cohend.v <- lapply(c(1:1000), FUN = function(x){
  
  p <- unlist(permutation.v[x])
  
  sp <- cohen.d(same.osn.dist.v,p)
  sp <- sp$estimate
  
  dp <- cohen.d(different.osn.dist.v,p)
  dp <- dp$estimate
  
  c(sp,dp)
  
})

## Checkpoint ##

saveRDS(density.cohend.v,file = "density_cohend_ventral_20220214.rds")

density.cohend.v <- readRDS("density_cohend_ventral_20220214.rds")

density.cohend.dframe.v <- data.frame(t(data.frame(density.cohend.v)))
a <- data.frame(dist = density.cohend.dframe.v$X1, type = "intra_perm")
b <- data.frame(dist = density.cohend.dframe.v$X2, type = "inter_perm")

density.cohend.dframe.v <- rbind(a,b)
density.cohend.dframe.v$dist <- density.cohend.dframe.v$dist*-1


stat.test <- ks.test(density.cohend.dframe.v$dist[density.cohend.dframe.v$type == "intra_perm"],
                     density.cohend.dframe.v$dist[density.cohend.dframe.v$type == "inter_perm"])

cohend_all_clusters.v <- ggplot(density.cohend.dframe.v,aes(y = dist,
                                                            x= type))+ 
  geom_violin(scale = "width")+ 
  geom_hline(yintercept = -1*intra.inter.v$estimate,
             linetype="dotted",color="red")+
  ylim(-0.3,1.7)+
  ylab("Cohen's D")+
  annotate(geom = "text",label ="***",y=1.7,x=1.5)+
  common_layout+
  theme(axis.title.x = element_blank())



# Select the first 100 permutations
permutations.clusters.v <- permutations.clusters.v[1:100]


# Combine the permutated matrix
permutations.clusters.v <- dplyr::bind_rows(permutations.clusters.v)
permutations.clusters.v$osn_pairs <- "permutation"


# Combine the ditance matrix to the permutated distance matrix
dist.mat.without.or.per.v <- dplyr::bind_rows(dist.mat.clusters.v, permutations.clusters.v)
dist.mat.without.or.per.v$osn_pairs <- factor(dist.mat.without.or.per.v$osn_pairs,
                                              levels = c("same OSN pop.",
                                                         "different OSN pop.",
                                                         "permutation"))

levels(dist.mat.without.or.per.v$osn_pairs) <- c("intra.","inter.","perm.")

euclidean.dist.noOR.v.plot <- ggplot(data = dist.mat.without.or.per.v,
                                     mapping = aes(y = euclidean_distance,
                                                   x = osn_pairs,
                                                   fill = osn_pairs,
                                                   group = osn_pairs)) +
  geom_violin(alpha = 0.5,scale = "width",color="green",fill = "green") +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..)) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..)) +
  stat_summary(geom = "point", fun = "median")+
  labs(y = "eucl. dist.") + 
  ylim(0,70)+
  common_layout +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0,20),
        legend.background = element_rect(fill='transparent'))


osn.sub.integrated.sct <- readRDS("moe_noexp_scRNAseq10x_cellranger_filter_OSN_WITH_OLFR_sct_integrated_clusters.rds")

# Eucledian distance per clusters
pca.coord <- osn.sub.integrated.sct@reductions$pca@cell.embeddings[,1:14]
pca.coord.split <- as.data.frame(cbind(pca.coord,
                                       as.character(osn.sub.integrated.sct$clusters)))

pca.coord.split <- split(pca.coord.split,pca.coord.split$V15)  

osn.idents <- osn.sub.integrated.sct$olfr.ident

eucl.dist.cluster.with_OR <- lapply(X = pca.coord.split,FUN = function(pca.c){
  pca.c <- pca.c[,-c(15)]
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

eucl.dist.cluster.with_OR.v <- eucl.dist.cluster.with_OR$Ventral

# Statistical test
same.osn.dist <- eucl.dist.cluster.with_OR.v$euclidean_distance[eucl.dist.cluster.with_OR.v$osn_pairs == "same OSN pop."]
different.osn.dist <- eucl.dist.cluster.with_OR.v$euclidean_distance[eucl.dist.cluster.with_OR.v$osn_pairs == "different OSN pop."]

same_diff <- cohen.d(same.osn.dist,different.osn.dist)


# Ventral cluster plot
eucl.dist.cluster.with_OR.v$osn_pairs <- factor(eucl.dist.cluster.with_OR.v$osn_pairs,
                                                levels = c("same OSN pop.","different OSN pop."))

eucl.dist.cluster.with_OR.v.plot <- ggplot(data = eucl.dist.cluster.with_OR.v,
                                           mapping = aes(y = euclidean_distance,
                                                         x = osn_pairs,
                                                         fill = osn_pairs,
                                                         group = osn_pairs)) + 
  geom_violin(alpha = 0.5,scale = "width",color="gray50",fill = "gray50") +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..)) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..)) +
  stat_summary(geom = "point", fun = "median")+
  labs(y = "eucl. dist.") + 
  ylim(0,70)+
  common_layout +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0,20),
        legend.background = element_rect(fill='transparent'))



euclidean.dist.plot.ventral <- wrap_plots(A=euclidean.dist.noOR.v.plot, 
                                          B=cohend_all_clusters.v,
                                          C =eucl.dist.cluster.with_OR.v.plot,
                                          design = "AABC")



ggsave(euclidean.dist.plot.ventral,
       filename = "fig_1_E_eucldist_ventral.pdf",
       width = 6,
       height = 2)







