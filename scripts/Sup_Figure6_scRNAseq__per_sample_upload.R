## Euclidean Distances per sample ##

# load packages
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

# Python package
library(reticulate)
py_install("kneed")
py_install("matplotlib")

## Analysis par sample
osn.no.or.sub.integrated.sct <- readRDS("moe_noexp_scRNAseq10x_cellranger_filter_OSN_NO_OLFR_sct_integrated_clusters.rds")

a <- SplitObject(osn.no.or.sub.integrated.sct,split.by = "orig.ident")

# Filter cells per size (only pop > 3 cells)
a.filt <-lapply(a, FUN = function(samples){
  
  or.pop.size <- as.data.frame(x = table(samples$olfr.ident))
  olfrs.to.keep <- as.character(x = or.pop.size$Var1[or.pop.size$Freq >= 3])
  cells.to.keep <- rownames(x = samples@meta.data)[samples@meta.data$olfr.ident %in% olfrs.to.keep]
  
  samples <- subset(x = samples,
                    cells = cells.to.keep)
  
})


# Normalize, Scale and compute PC per sample


b <- lapply(a.filt, FUN = function(x){
  
  DefaultAssay(x) <- "RNA"
  
  x <- SCTransform(x,variable.features.n = 5000)
  
  
  x <- RunPCA(x)
  
})

kn <- import("kneed")


# Euclidean distance per sample
d <- lapply(b, FUN = function(x){
  
  
  knaa <- kn$KneeLocator(x = seq_len(length.out = length(x = x@reductions$pca@stdev)),
                         y = x@reductions$pca@stdev,
                         S = 1,
                         curve = "convex",
                         direction = "decreasing")
  pc <- knaa$knee
  
  pca.coord.without.or <- x@reductions$pca@cell.embeddings[,1:pc]
  
  dist.mat.without.or <- fields::rdist(x1 = pca.coord.without.or, x2 = NULL)
  
  colnames(dist.mat.without.or) <- rownames(pca.coord.without.or)
  rownames(dist.mat.without.or) <- rownames(pca.coord.without.or)
  
  osn.idents <- x$olfr.ident
  
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
  
  
  dist.mat.without.or
  
})

# Permutated distance per sample
permuation.samples <- lapply(b, FUN = function(samples){
  
  knaa <- kn$KneeLocator(x = seq_len(length.out = length(x = samples@reductions$pca@stdev)),
                         y = samples@reductions$pca@stdev,
                         S = 1,
                         curve = "convex",
                         direction = "decreasing")
  pc <- knaa$knee
  
  pca.coords <- samples@reductions$pca@cell.embeddings[,1:pc]
  osn.idents <- samples$olfr.ident
  
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
  
})

saveRDS(permuation.samples,file = "permutation_per_sample_20220227")
permuation.samples <- readRDS("permutation_per_sample_20220227")


plots.samples <- lapply(c(1:4), FUN = function(s){
  
  # Combine the permuted matrix
  permuted.distances <- dplyr::bind_rows(permuation.samples[s][[1]])
  permuted.distances$osn_pairs <- "permutation"

  #Combine the distance matrix to the permuted distance matrix
  dist.mat.without.or.per <- as.data.frame(dplyr::bind_rows(d[s][[1]], permuted.distances))
  dist.mat.without.or.per$osn_pairs <- factor(dist.mat.without.or.per$osn_pairs,
                                              levels = c("same OSN pop.",
                                                         "different OSN pop.",
                                                         "permutation"))
  
  levels(dist.mat.without.or.per$osn_pairs) <- c("intra.","inter.","perm.")
  
  permuted.distances
  
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
  
  
})

cohenD.sample <- lapply(c(1:4), FUN = function(s){
  # Cohen's D
  same.osn.dist <- d[s][[1]]$euclidean_distance[d[s][[1]]$osn_pairs == "same OSN pop."]
  different.osn.dist <- unlist(d[s][[1]]$euclidean_distance[d[s][[1]]$osn_pairs == "different OSN pop."])
  permutation <- lapply(permuation.samples[s][[1]], '[[', 3)
  
  intra.inter <- cohen.d(same.osn.dist,sample(different.osn.dist))
  
  
  density.cohend <- lapply(c(1:1000), FUN = function(x){
    
    p <- unlist(permutation[x])
    
    sp <- cohen.d(same.osn.dist,p)
    sp <- sp$estimate
    
    dp <- cohen.d(different.osn.dist,p)
    dp <- dp$estimate
    
    c(sp,dp)
    
  })
})

saveRDS(cohenD.sample,"f.cohensD.sct.rds")
cohenD.sample <- readRDS("f.cohensD.sct.rds")

cohenD.plot <- lapply(c(1:4), FUN = function(s){
  
  same.osn.dist <- d[s][[1]]$euclidean_distance[d[s][[1]]$osn_pairs == "same OSN pop."]
  different.osn.dist <- unlist(d[s][[1]]$euclidean_distance[d[s][[1]]$osn_pairs == "different OSN pop."])
  
  intra.inter <- cohen.d(different.osn.dist,same.osn.dist)
  
  density.cohend.dframe <- data.frame(t(data.frame(cohenD.sample[s][[1]])))
  a <- data.frame(dist = density.cohend.dframe$X1, type = "intra_perm")
  b <- data.frame(dist = density.cohend.dframe$X2, type = "inter_perm")
  
  density.cohend.dframe <- rbind(a,b)
  density.cohend.dframe$dist <- density.cohend.dframe$dist*c(-1)
  
  stat.test <- ks.test(density.cohend.dframe$dist[density.cohend.dframe$type == "intra_perm"],
                       density.cohend.dframe$dist[density.cohend.dframe$type == "inter_perm"])
  
  cohend_all_clusters <- ggplot(density.cohend.dframe,aes(y = dist,
                                                          x= type))+ 
    geom_violin(scale = "width")+ 
    geom_hline(yintercept = intra.inter$estimate,
               color="gray")+
    ylim(-0.3,1.7)+
    ylab("Cohen's D")+
    annotate(geom = "text",label ="***",y=1.7,x=1.5)+
    common_layout+
    theme(axis.title.x = element_blank())
  
  
})

final.plot <- lapply(c(1:4), FUN = function(x){
  wrap_plots(plots.samples[x][[1]],cohenD.plot[x][[1]],ncol = 2)
})


f <- wrap_plots(final.plot)



ggsave(plot = f, filename = "Sup_6.pdf",height = 6,width = 10)

