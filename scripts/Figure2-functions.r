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

### FUNCTIONS
find_furthest_point <- function(mat) {
  subcluster.densities <- mapply(
    i = 1:nrow(mat),
    MoreArgs = list("mat" = mat),
    SIMPLIFY = F,
    FUN = function(i, mat) {
      submat <- mat[-i,]
      centroid <- colMeans(submat)
      dist.to.centroid <- apply(submat, 1, function(x) { rdist(rbind(x, centroid))[1,2] }) 
      return(mean(dist.to.centroid))
    }) %>% unlist()
  worse <- which(subcluster.densities == min(subcluster.densities))[1]
  return(worse)
}

remove_furthest_point <- function(mat) {
  return(mat[-find_furthest_point(mat),])
}

find_outliers <- function(mat, sd.max = 1.96) {
  return(which(as.numeric(scale(rowMeans(rdist(mat))) > sd.max)))
}

remove_outliers <- function(mat, sd.max = 1.96) {
  return(mat[-find_outliers(mat),])
}

find_centroid <- function(mat, trim=0, sd.max = NA) {
  mat.x <- mat
  if (trim >= 1) {
    for (i in 1:trim) {
      mat.x <- remove_furthest_point(mat.x)
    }  
  }
  if (!is.na(sd.max)) {
    mat.x <- remove_outliers(mat.x, sd.max=sd.max)
  }
  
  return(colMeans(mat.x))
}

# pairwise distance using k-nearest neighbors identity
# !! OUTLIER FUNCTION TO BE FIXED
compute_pairwise_shared_neighbors <- function(
  
  # Each element of the list contains a set of coordinates for the
  # points of a given population. All possible pairs of population are
  # being compared.
  osn_pca.list,
  
  # Performed the operation on a subset of each population. The
  # provided value (0;1] corresponds to the fraction of each population
  # to be sampled for the calculation.
  random_subsampling = 1,
  
  # Remove outliers prior to the compuation. The provided value [0;1)
  # corresponds to the fraction of each population to be discarded.
  outlier_removal = 0) {
  
  # iterate over gene_name.1
  mat <- mapply(
    gene_name.1 = names(osn_pca.list),
    i = seq_along(names(osn_pca.list)),
    MoreArgs = list("osn_pca.list" = osn_pca.list),
    SIMPLIFY = F,
    FUN = function(gene_name.1, osn_pca.list, i) {
      
      # iterate over gene_name.2
      rows <- mapply(
        gene_name.2 = names(osn_pca.list),
        j = seq_along(names(osn_pca.list)),
        MoreArgs = list("gene_name.1" = gene_name.1, 
                        "osn_pca.list" = osn_pca.list, 
                        "i" = i),
        SIMPLIFY = T,
        FUN = function(gene_name.1, gene_name.2, osn_pca.list, i, j) {
          
          # do not compare with self, and not reversely
          if (j <= i) {
            return(NA)
          } else {
            
            # get PC values of gene 1
            mat.1 <- osn_pca.list[[gene_name.1]]
            
            # subsampling
            if (random_subsampling < 1) {
              n.1 <- round(random_subsampling*nrow(mat.1))
              mat.1 <- mat.1[sample(1:nrow(mat.1), n.1, replace = F),]
            }
            
            # outlier removal
            if (outlier_removal > 0) {
              outliers.n.1 <- round(outlier_removal*nrow(mat.1))
              if (nrow(mat.1) - outliers.n.1 > 0 & outliers.n.1 < nrow(mat.1)) {
                n.1 <- 1 - round(outlier_removal*nrow(mat.1))
                for (iteration in 1:(1-n.1)) {
                  outlier.i <- find_outlier(mat.1[,grep("^PC_", colnames(mat.1))])
                  mat.1 <- mat.1[-outlier.i,]
                }  
              }
            }
            
            cell.ids.1 <- mat.1[,"cell.ident"]
            n.1 <- length(cell.ids.1)
            
            # get PC values of gene 2
            mat.2 <- osn_pca.list[[gene_name.2]]
            
            # subsampling
            if (random_subsampling < 1) {
              n.2 <- round(random_subsampling*nrow(mat.2))
              mat.2 <- mat.2[sample(1:nrow(mat.2), n.2, replace = F),]
            }
            
            # outlier removal
            if (outlier_removal > 0) {
              outliers.n.2 <- round(outlier_removal*nrow(mat.2))
              if (nrow(mat.2) - outliers.n.2 > 0 & outliers.n.2 < nrow(mat.2)) {
                n.2 <- 1 - round(outlier_removal*nrow(mat.2))
                for (iteration in 1:(1-n.2)) {
                  outlier.i <- find_outlier(mat.2[,grep("^PC_", colnames(mat.2))])
                  mat.2 <- mat.2[-outlier.i,]
                }  
              }
            }
            
            cell.ids.2 <- mat.2[,"cell.ident"]
            n.2 <- length(cell.ids.2)
            
            # merge matrices
            mat <- rbind(mat.1, mat.2)
            olfr.ident <- mat[,"olfr.ident"]
            cell.ident <- c(cell.ids.1, cell.ids.2)
            names(olfr.ident) <- cell.ident
            
            # keep only cells with PC values
            mat <- mat[,which(startsWith(x = colnames(mat), prefix = "PC_"))]
            
            # calculate pairwise distance
            dist.mat <- rdist(mat)
            rownames(dist.mat) <- cell.ident
            colnames(dist.mat) <- cell.ident
            
            # search k-neighbors
            k <- max(c(n.1, n.2)) - 1
            nearest.neighbors <- FNN::get.knn(dist.mat, k)$nn.index
            #colnames(nearest.neighbors) <- paste("nn", 1:ncol(nearest.neighbors), sep = "_")
            rownames(nearest.neighbors) <- cell.ident
            
            # for each cell, get the fraction of nearest neighbor that are
            # of the same identity
            nearest.neighbors.1 <- nearest.neighbors[which(olfr.ident[cell.ident] == gene_name.1),]
            cell.scores.1 <- apply(nearest.neighbors.1, 1, function(nearest.neighbor) {
              same <- length(which(olfr.ident[cell.ident[nearest.neighbor[1:(n.1-1)]]] == gene_name.1))
              return(same/(n.1-1))
            })
            nearest.neighbors.2 <- nearest.neighbors[which(olfr.ident[cell.ident] == gene_name.2),]
            cell.scores.2 <- apply(nearest.neighbors.2, 1, function(nearest.neighbor) {
              same <- length(which(olfr.ident[cell.ident[nearest.neighbor[1:(n.2-1)]]] == gene_name.2))
              return(same/(n.2-1))
            })
            return(mean(c(cell.scores.1, cell.scores.2)))  
          }
        }
      )
      return(rows)
    })
  return(do.call(rbind, mat))
}
