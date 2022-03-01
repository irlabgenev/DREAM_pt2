## Control for PCs choice with OR ##

library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(Hmisc)
library(stringr)
library(reshape2)
library(effsize)

source("common_layout.r")

# Read data

osn.or.integrated.sct <- readRDS(file = "moe_noexp_scRNAseq10x_cellranger_filter_OSN_WITH_OLFR_sct_integrated_clusters.rds")


coehnsd.pcs.with.or <- lapply(seq(2,30,2), FUN = function(x){
  
  print(paste("PC_",x))
  # Compute the euclidean distance 
  pca.coord.without.or <- osn.or.integrated.sct@reductions$pca@cell.embeddings[,1:x]
  
  dist.mat.without.or <- fields::rdist(x1 = pca.coord.without.or, x2 = NULL)
  
  colnames(dist.mat.without.or) <- rownames(pca.coord.without.or)
  rownames(dist.mat.without.or) <- rownames(pca.coord.without.or)
  
  osn.olfr <- osn.or.integrated.sct$olfr.ident
  
  
  dist.mat.without.or[lower.tri(dist.mat.without.or, diag = TRUE)] <- NA
  
  dist.mat.without.or <- melt(data = dist.mat.without.or,
                              na.rm = TRUE,
                              varnames = c("OSN_1", "OSN_2"),
                              value.name = "euclidean_distance")
  
  head(dist.mat.without.or)
  
  dist.mat.without.or <- data.table::data.table(dist.mat.without.or)
  dist.mat.without.or$OSN_1 <- as.character(x = dist.mat.without.or$OSN_1)
  dist.mat.without.or$OSN_2 <- as.character(x = dist.mat.without.or$OSN_2)
  
  
  dist.mat.without.or$OSN_1_ident <- osn.olfr[dist.mat.without.or$OSN_1]
  dist.mat.without.or$OSN_2_ident <- osn.olfr[dist.mat.without.or$OSN_2]
  
  
  
  dist.mat.without.or$osn_pairs <- ifelse(test = dist.mat.without.or$OSN_1_ident == dist.mat.without.or$OSN_2_ident,
                                          yes = "intra OSN pop.",
                                          no = "different OSN pop.")
  
  
  # Within 
  
  inter <- sample(dist.mat.without.or[dist.mat.without.or$osn_pairs =="different OSN pop.",]$euclidean_distance)
  intra <- sample(dist.mat.without.or[dist.mat.without.or$osn_pairs =="intra OSN pop.",]$euclidean_distance)
  
  dist.plot <- rbind(data.frame(dist=intra, type = "intra",pc = x),data.frame(dist=inter, type = "inter",pc = x))
  
  
  dist.plot
  
}

)


## Checkpoint

#saveRDS(coehnsd.pcs.with.or,file = "test.multiple.pc_withOR_2.rds")

#coehnsd.pcs.with.or <- readRDS("test.multiple.pc.rds")

coh <- lapply(coehnsd.pcs.with.or, FUN = function(x){
  d2 <- x[x[2] =="intra",]$dist
  d1 <- x[x[2] =="inter",]$dist
  
  a <- cohen.d(d1,d2)
  a$estimate
  
})

coh.plot <- data.frame("Cohen" = unlist(coh), pc = seq(2,30,2))

cd.plot.pc <- ggplot(coh.plot,aes(x=pc,y=Cohen))+ geom_point() +
  geom_vline(xintercept = 14,color="red",linetype="dotted")+
  scale_x_continuous(breaks = seq(2,30,2),limits = c(2,30))+
  ylab("Cohen's D")+
  xlab("PCs")+
  common_layout

ggsave(cd.plot.pc,
       filename = "Sup_3e,f,g.pdf",
       height = 6,
       width = 6)
