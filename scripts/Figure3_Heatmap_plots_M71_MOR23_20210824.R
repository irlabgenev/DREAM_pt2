#HEATMAP Plots

library(ggplot2)
library(gplots)
library(reshape2);
library(RColorBrewer);
library(pheatmap)
library(viridis)

resM71mutated.FC <- resM71mutatedLFC
resMOR23mutated.FC <- resMOR23mutatedLFC
#Prepare the data to plot the heatmaps ordered by FC

#M71

data.M71 <- resM71mutated.FC[resM71mutated.FC$color_flag != "grey",]

data.M71<- data.M71[,c(1,3,7:12)]
rownames(data.M71) <- data.M71$GeneName
data.M71$GeneName <- NULL
a <- data.M71[,1, drop=T]
data.M71$log2FoldChange <- NULL
data.log.M71 <- apply(X = data.M71, MARGIN = 2, log1p)

data.log.M71 <- cbind(data.log.M71,a)
data.log.M71 <- as.data.frame(data.log.M71)

data.log.M71 <- data.log.M71[order(data.log.M71$a, decreasing = T),]
data.log.M71$a <- NULL

longData.M71<- data.matrix(data.log.M71)

#MOR23

data.Mor23 <- resMOR23mutated.FC[resMOR23mutated.FC$color_flag != "grey",]

data.Mor23<- data.Mor23[,c(1,3,7:12)]
rownames(data.Mor23) <- data.Mor23$GeneName
data.Mor23$GeneName <- NULL
b <- data.Mor23[,1, drop=T]
data.Mor23$log2FoldChange <- NULL

data.log.Mor23 <- apply(X = data.Mor23, MARGIN = 2, log1p)

data.log.Mor23 <- cbind(data.log.Mor23,b)
data.log.Mor23 <- as.data.frame(data.log.Mor23)
data.log.Mor23 <- data.log.Mor23[order(data.log.Mor23$b, decreasing = T),]
data.log.Mor23$b <- NULL

longData.Mor23<- data.matrix(data.log.Mor23)



#Scaling the data
longdata_norm_Mor23 <- scale(x = t(x = longData.Mor23),
                              scale = TRUE,
                              center = TRUE)


longdata_norm_M71 <- scale(x = t(x = longData.M71),
                              scale = TRUE,
                              center = TRUE)


#Preparing the annotations - M71

#annotation_row.M71 = data.frame(
 # State = factor(rep(c("Upregulated", "Downregulated"), c(367,434)))
#)

annotation_col.M71 = data.frame(
  Condition = factor(rep(c("Exposed", "Control"), each = 3)))

rownames(annotation_col.M71) = c("Ace1", "Ace2","Ace3", "F1","F2", "F3")

ann_colors.M71 = list(
  Condition = c(Exposed = "red", Control = "blue")
)


#Preparing the annotations - MOR23


#annotation_row.Mor23 = data.frame(
 # State = factor(rep(c("Upregulated", "Downregulated"), c(665,484)))
#)

annotation_col.Mor23 = data.frame(
  Condition = factor(rep(c("Exposed", "Control"), each = 3)))

rownames(annotation_col.Mor23) = c("L1", "L3","L5","D1", "D3","D5")

ann_colors.Mor23 = list(
  Condition = c(Exposed = "red", Control = "blue")
)




#Basic HEATMAPS

c <-pheatmap(t(longdata_norm_Mor23), 
         annotation_col = annotation_col.Mor23,
         annotation_colors = ann_colors.Mor23,
         color = viridis(10) ,
         show_rownames = F,
         cluster_rows = F,
         cluster_cols = F)


a <- pheatmap(t(longdata_norm_M71),
         annotation_col = annotation_col.M71,
         annotation_colors = ann_colors.M71,
         color = viridis(10),
         show_rownames = F,
         cluster_rows = F,
         cluster_cols = F)








