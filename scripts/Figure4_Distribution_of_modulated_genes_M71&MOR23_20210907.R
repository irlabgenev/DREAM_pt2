#Required packages
require(DESeq2)
require(edgeR)
require(ggplot2)
require(ggrepel)
require(dplyr)
require(apeglm)
require(ggbeeswarm)



MOR23.d <- resMOR23mutatedLFC[resMOR23mutatedLFC$state_flag == "Downreg",]
M71.d <- resM71mutatedLFC[resM71mutatedLFC$state_flag == "Downreg",]




MOR23.df <- resMOR23mutatedLFC

M71.df <- resM71mutatedLFC




common.MOR23 <- MOR23.df[MOR23.df$GeneName %in% M71.df$GeneName,]
Common.M71 <- M71.df[M71.df$GeneName %in% common.MOR23$GeneName,]
M71.FC <- Common.M71[,c(3,16)]
MOR23.M71 <- cbind(common.MOR23,M71.FC)
MOR23.M71 <- MOR23.M71[,c(1,3,16,18:19)]
colnames(MOR23.M71) <- c("Genename", "Mor23.log2FoldChange","Mor23.state_flag","M71.log2FoldChange","M71.state_flag")



MOR23.M71$color <-  ifelse(MOR23.M71$Mor23.log2FoldChange < 0 & MOR23.M71$M71.log2FoldChange <0 & MOR23.M71$Mor23.state_flag == "Downreg" & MOR23.M71$M71.state_flag == "Downreg", "blue",
                      ifelse(MOR23.M71$Mor23.log2FoldChange < 0 & MOR23.M71$M71.log2FoldChange <0 & MOR23.M71$Mor23.state_flag == "Downreg" & MOR23.M71$M71.state_flag != "Downreg", "cyan",
                             ifelse(MOR23.M71$Mor23.log2FoldChange < 0 & MOR23.M71$M71.log2FoldChange <0 & MOR23.M71$Mor23.state_flag != "Downreg" & MOR23.M71$M71.state_flag == "Downreg", "cyan",
                                    ifelse(MOR23.M71$Mor23.log2FoldChange < 0 & MOR23.M71$M71.log2FoldChange >0 & MOR23.M71$Mor23.state_flag != "n.s." & MOR23.M71$M71.state_flag == "n.s.", "cyan",
                                           ifelse(MOR23.M71$Mor23.log2FoldChange < 0 & MOR23.M71$M71.log2FoldChange > 0 & MOR23.M71$Mor23.state_flag =="n.s." & MOR23.M71$M71.state_flag != "n.s.", "tomato",
                                                  ifelse(MOR23.M71$Mor23.log2FoldChange > 0 & MOR23.M71$M71.log2FoldChange < 0 & MOR23.M71$Mor23.state_flag != "n.s." & MOR23.M71$M71.state_flag == "n.s.", "tomato",
                                                         ifelse(MOR23.M71$Mor23.log2FoldChange > 0 & MOR23.M71$M71.log2FoldChange < 0 & MOR23.M71$Mor23.state_flag =="n.s." & MOR23.M71$M71.state_flag != "n.s.", "cyan",
                                                                ifelse(MOR23.M71$Mor23.log2FoldChange > 0 & MOR23.M71$M71.log2FoldChange > 0 & MOR23.M71$Mor23.state_flag == "Upreg" & MOR23.M71$M71.state_flag == "Upreg", "red",
                                                                       ifelse(MOR23.M71$Mor23.log2FoldChange > 0 & MOR23.M71$M71.log2FoldChange > 0 & MOR23.M71$Mor23.state_flag == "Upreg" & MOR23.M71$M71.state_flag != "Upreg", "tomato",
                                                                              ifelse(MOR23.M71$Mor23.log2FoldChange > 0 & MOR23.M71$M71.log2FoldChange > 0 & MOR23.M71$Mor23.state_flag != "Upreg" & MOR23.M71$M71.state_flag == "Upreg", "tomato",
                                                                                     "grey"))))))))))



remove.non.significant <- MOR23.M71[MOR23.M71$color != "grey",]

##Function to plot

 p2 <- ggplot(data = remove.non.significant,
             aes(Mor23.log2FoldChange,
                 M71.log2FoldChange)) +
  geom_point(size = 0.5) +
  geom_point(aes(colour = color),
             size = 1) +
  scale_color_manual(values=c("blue" = "blue",
                              "cyan" = "cyan",
                              "tomato" = "tomato",
                              "red" = "red",
                              "grey" = "grey"),
                     guide=FALSE)  +
  labs(x = "M71 FoldChange (log2)", 
       y = "MOR23 FoldChange (log2)")  +
  geom_vline(xintercept = 0, colour = "black") + 
  geom_abline(intercept = 0, colour = "blue") +
  geom_hline(yintercept = 0, colour = "black")+
  #xlim(-5,10) +
  #ylim(-5,10)+
  theme(aspect.ratio = 1) +
  theme_classic()
  #ggtitle("Distribution of the downregulated genes from M71 & MOR23")



p2$data$color <- factor(p2$data$color, levels = c("grey","tomato","cyan","lightblue","red","blue"))
p2$data <- p2$data[order(p2$data$color),]
 
p2




#old colours

MOR23.M71$color <-  ifelse(MOR23.M71$Mor23.log2FoldChange < 0 & MOR23.M71$M71.log2FoldChange <0 & MOR23.M71$Mor23.state_flag == "Downreg" & MOR23.M71$M71.state_flag == "Downreg", "blue",
                           ifelse(MOR23.M71$Mor23.log2FoldChange < 0 & MOR23.M71$M71.log2FoldChange <0 & MOR23.M71$Mor23.state_flag == "Downreg" & MOR23.M71$M71.state_flag != "Downreg", "lightblue",
                                  ifelse(MOR23.M71$Mor23.log2FoldChange < 0 & MOR23.M71$M71.log2FoldChange <0 & MOR23.M71$Mor23.state_flag != "Downreg" & MOR23.M71$M71.state_flag == "Downreg", "cyan",
                                         ifelse(MOR23.M71$Mor23.log2FoldChange < 0 & MOR23.M71$M71.log2FoldChange >0 & MOR23.M71$Mor23.state_flag != "n.s." & MOR23.M71$M71.state_flag == "n.s.", "gold",
                                                ifelse(MOR23.M71$Mor23.log2FoldChange < 0 & MOR23.M71$M71.log2FoldChange > 0 & MOR23.M71$Mor23.state_flag =="n.s." & MOR23.M71$M71.state_flag != "n.s.", "darkgoldenrod",
                                                       ifelse(MOR23.M71$Mor23.log2FoldChange > 0 & MOR23.M71$M71.log2FoldChange < 0 & MOR23.M71$Mor23.state_flag != "n.s." & MOR23.M71$M71.state_flag == "n.s.", "green",
                                                              ifelse(MOR23.M71$Mor23.log2FoldChange > 0 & MOR23.M71$M71.log2FoldChange < 0 & MOR23.M71$Mor23.state_flag =="n.s." & MOR23.M71$M71.state_flag != "n.s.", "seagreen",
                                                                     ifelse(MOR23.M71$Mor23.log2FoldChange > 0 & MOR23.M71$M71.log2FoldChange > 0 & MOR23.M71$Mor23.state_flag == "Upreg" & MOR23.M71$M71.state_flag == "Upreg", "red",
                                                                            ifelse(MOR23.M71$Mor23.log2FoldChange > 0 & MOR23.M71$M71.log2FoldChange > 0 & MOR23.M71$Mor23.state_flag == "Upreg" & MOR23.M71$M71.state_flag != "Upreg", "tomato",
                                                                                   ifelse(MOR23.M71$Mor23.log2FoldChange > 0 & MOR23.M71$M71.log2FoldChange > 0 & MOR23.M71$Mor23.state_flag != "Upreg" & MOR23.M71$M71.state_flag == "Upreg", "orangered",
                                                                                          "grey"))))))))))

