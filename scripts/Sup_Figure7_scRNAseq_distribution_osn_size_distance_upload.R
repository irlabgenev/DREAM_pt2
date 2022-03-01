## Integrastion control ##

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
osn.no.or.sub.integrated.sct <- readRDS(file = "moe_noexp_scRNAseq10x_cellranger_filter_OSN_NO_OLFR_sct_integrated_clusters.rds")

# Number of OSN per OR pop.
osn.size <- data.frame(table(osn.no.or.sub.integrated.sct$olfr.ident))
colnames(osn.size) <- c("OR","count")

osn.size$bin <-cut2(osn.size$count,g = 6)

table(osn.size$bin)

osn.in.name <- setNames(as.character(osn.size$bin), osn.size$OR)

# Extract digit
regexp <- "[[:digit:]]+"
bins.loc <- c(str_extract(levels(osn.size$bin), regexp))
bins.loc <- as.numeric(bins.loc[-1])

osn.size$size <- ifelse(test = osn.size$bin %in% c("[ 3, 5)","[21,83]"),
                        yes= ifelse(test = osn.size$bin %in% c("[ 3, 5)"),
                                    yes = "small",
                                    no = "large"),
                        no = "no")



dist.bin <- ggplot(osn.size,aes(x=osn.size$count))+ geom_density(fill="lightblue",color="white")+
  geom_vline(xintercept = bins.loc,color="white",size=1) +
  common_layout + 
  xlab("OSN pop. size") + 
  annotate(geom = "text",label="small",x = 2.5,y=0.005)+
  annotate(geom = "text",label="large",x = 25,y=0.005)+
  scale_x_continuous(expand = c(0, 3)) +
  scale_y_continuous(expand = c(0, 0)) 




osn.no.or.sub.integrated.sct$bin <- osn.in.name[osn.no.or.sub.integrated.sct$olfr.ident]

table(osn.no.or.sub.integrated.sct$bin)

subset(osn.no.or.sub.integrated.sct, bin %in% c("[ 3, 5)","[21,83]"))

osn.no.or.sub.integrated.sct.sub <- subset(osn.no.or.sub.integrated.sct, bin %in% c("[ 3, 5)","[21,83]"))


# Compute the euclidean distance 
pca.coord.without.or <- osn.no.or.sub.integrated.sct.sub@reductions$pca@cell.embeddings[,1:15]

dist.mat.without.or <- fields::rdist(x1 = pca.coord.without.or, x2 = NULL)

colnames(dist.mat.without.or) <- rownames(pca.coord.without.or)
rownames(dist.mat.without.or) <- rownames(pca.coord.without.or)

osn.bin <- osn.no.or.sub.integrated.sct.sub$bin
osn.olfr <- osn.no.or.sub.integrated.sct.sub$olfr.ident

dist.mat.without.or[1:10,1:10]

dist.mat.without.or[lower.tri(dist.mat.without.or, diag = TRUE)] <- NA

dist.mat.without.or <- melt(data = dist.mat.without.or,
                            na.rm = TRUE,
                            varnames = c("OSN_1", "OSN_2"),
                            value.name = "euclidean_distance")

head(dist.mat.without.or)

dist.mat.without.or <- data.table::data.table(dist.mat.without.or)
dist.mat.without.or$OSN_1 <- as.character(x = dist.mat.without.or$OSN_1)
dist.mat.without.or$OSN_2 <- as.character(x = dist.mat.without.or$OSN_2)

dist.mat.without.or$OSN_1_bin <- osn.bin[dist.mat.without.or$OSN_1]
dist.mat.without.or$OSN_2_bin <- osn.bin[dist.mat.without.or$OSN_2]

dist.mat.without.or$OSN_1_ident <- osn.olfr[dist.mat.without.or$OSN_1]
dist.mat.without.or$OSN_2_ident <- osn.olfr[dist.mat.without.or$OSN_2]

dist.mat.without.or$bin_pairs <- ifelse(test = dist.mat.without.or$OSN_1_bin == dist.mat.without.or$OSN_2_bin,
                                        yes = ifelse(test = dist.mat.without.or$OSN_1_bin == "[21,83]",
                                                     yes = "large-large",
                                                     no  = "small-small"),
                                        no = "small-large")

dist.mat.without.or$osn_pairs <- ifelse(test = dist.mat.without.or$OSN_1_ident == dist.mat.without.or$OSN_2_ident,
                                        yes = "intra OSN pop.",
                                        no = "different OSN pop.")



## Plots and subsampling


# Subsample 1/3 cells from each pop and compute the Cohen's D between small and large
subsample <- lapply(c(1:10000), FUN = function(x){
  
  print(paste("Subsample_",x))
  
  data.olfr.s.l <- data.frame(olfr.ident=osn.no.or.sub.integrated.sct.sub$olfr.ident)
  data.olfr.s.l.split <- split(data.olfr.s.l,data.olfr.s.l$olfr.ident)
  
  cells.names <- lapply(names(data.olfr.s.l.split), FUN = function(x){
    a <- rownames(data.olfr.s.l.split[x][[1]])
    a <- sample(a,(1/3)*length(a))
    a
    
  })
  
  cells.names <- unlist(cells.names)
  
  data.subsample <- dist.mat.without.or[dist.mat.without.or$OSN_1 %in% cells.names & 
                                          dist.mat.without.or$OSN_2 %in% cells.names,]
  
  ss <- data.subsample[data.subsample$bin_pairs =="small-small" &
                       data.subsample$osn_pairs == "different OSN pop.",]$euclidean_distance
  
  ll<- data.subsample[data.subsample$bin_pairs =="large-large" &
                      data.subsample$osn_pairs == "different OSN pop.",]$euclidean_distance
  
  intra <- dist.mat.without.or[dist.mat.without.or$osn_pairs =="intra OSN pop.",]$euclidean_distance
  
  ss.ll<-cohen.d(ss,ll)
  ss.inta <- cohen.d(ss,intra)
  ll.intra <- cohen.d(ll,intra)
  
  c(ss.ll$estimate,ss.inta$estimate,ll.intra$estimate)

}
)

saveRDS(subsample,"lagre_small_subsample_corrected.rds")

subsample <- readRDS("lagre_small_subsample_corrected.rds")

subsample.cohenD <- data.frame(t(data.frame(subsample)))
a <- data.frame(dist = subsample.cohenD$X1, type = "small-large")
b <- data.frame(dist = subsample.cohenD$X2, type = "small-intra")
c <- data.frame(dist = subsample.cohenD$X3, type = "large-intra")

subsample.cohenD.dframe <- rbind(a,b,c)

stat.test <- ks.test(subsample.cohenD.dframe[subsample.cohenD.dframe$type == "small-intra",]$dist,
                     subsample.cohenD.dframe[subsample.cohenD.dframe$type == "large-intra",]$dist)



subsample.cohenD.plot <- ggplot(subsample.cohenD.dframe,aes(y = dist,
                                                            x= type))+ 
  geom_violin(scale = "width")+ 
  ylim(-0.3,2)+
  ylab("Cohen's D")+
  annotate(geom = "text",label ="***",y=1.7,x=1.5)+
  common_layout+
  theme(axis.title.x = element_blank())


dist.mat.without.or.inter <- dist.mat.without.or[dist.mat.without.or$osn_pairs =="different OSN pop."]


violin.plot.bin <- ggplot(dist.mat.without.or.inter,aes(x=bin_pairs,
                                                  y = euclidean_distance))+
  geom_violin(fill = "lightblue",color="lightblue") + 
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..)) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..)) +
  stat_summary(geom = "point", fun = "median")+
  common_layout

violin.plot.dist <- ggplot(dist.mat.without.or,aes(x=osn_pairs,y = euclidean_distance))+
  geom_violin(color = "gray", fill = "gray",alpha = 0.5)+
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..)) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..)) +
  stat_summary(geom = "point", fun = "median")+
  common_layout


large.small.dist.final <- wrap_plots(A = dist.bin,  
                                     B = violin.plot.dist,
                                     C = violin.plot.bin,
                                     D = subsample.cohenD.plot,
                                     design = "AAA
                     BCD")


ggsave(large.small.dist.final,filename = "Sup_7.pdf",
       height = 8,width = 12)



