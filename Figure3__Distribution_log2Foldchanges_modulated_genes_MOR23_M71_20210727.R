library(ggplot2)

Mor23.df <- resMOR23mutatedLFC
M71.df <- resM71mutatedLFC


M71.corrected.down <- M71.df[M71.df$state_flag == "Downreg",]
M71.corrected.up <- M71.df[M71.df$state_flag == "Upreg",]

MOR23.corrected.down <- Mor23.df[Mor23.df$state_flag == "Downreg",]
MOR23.corrected.up <- Mor23.df[Mor23.df$state_flag == "Upreg",]


max.x.val.MOR23 <- ceiling(max(abs(x = c(MOR23.corrected.down$log2FoldChange,
                                         MOR23.corrected.up$log2FoldChange))))


max.x.val.M71 <- ceiling(max(abs(x = c(M71.corrected.down$log2FoldChange,
                                   M71.corrected.up$log2FoldChange))))



#MOR23.Down

ggplot(data = MOR23.corrected.down,
       aes(x = abs(x = log2FoldChange))) + 
  geom_histogram( colour="black", fill="white")+
  labs(x = "Abs. FoldChange (log2)", y= "Number of genes")+
  scale_x_continuous(expand =c(0, 0), limits = c(0, max.x.val.MOR23)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()



#MOR23 Up

ggplot(data=MOR23.corrected.up,
       aes(x = abs(x = log2FoldChange))) + 
  geom_histogram( colour="black", fill="white")+
  labs(x = "Abs. FoldChange(log2)", y= "Number of genes")+
  scale_x_continuous(expand =c(0, 0), limits = c(0, max.x.val.MOR23)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()


#M71.Down

ggplot(data=M71.corrected.down,
       aes(x = abs(x = log2FoldChange))) + 
  geom_histogram(colour="black", fill="white")+
  labs(x = "Abs. FoldChange(log2)", y= "Number of genes")+
  scale_x_continuous(limits = c(0, max.x.val.M71)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()

#M71.up

ggplot(data=M71.corrected.up,
       aes(x = abs(x = log2FoldChange))) + 
  geom_histogram(colour="black", fill="white") +
  labs(x = "Abs. FoldChange(log2)", y= "Number of genes")+
  scale_x_continuous(limits = c(0, max.x.val.M71)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()






