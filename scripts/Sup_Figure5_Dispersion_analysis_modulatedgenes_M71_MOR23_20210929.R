MOR23.df <- resMOR23mutatedLFC
M71.df <- resM71mutatedLFC

MOR23.df$variance.control <- apply(MOR23.df[c("D1", "D3", "D5")], 1, var)
MOR23.df$variance.exposed <- apply(MOR23.df[c("L1", "L3", "L5")], 1, var)

MOR23.df$mean.control <- apply(MOR23.df[c("D1","D3", "D5")], 1, mean)
MOR23.df$mean.exposed <- apply(MOR23.df[c("L1", "L3", "L5")], 1, mean)

MOR23.df$control<- log1p(MOR23.df[,18]/MOR23.df[,20])
MOR23.df$exposed<- log1p(MOR23.df[,19]/MOR23.df[,21])


dispersion.Mor23 <- MOR23.df[,c(1,16,22:23)]

dispersion.Mor23.down <- dispersion.Mor23[dispersion.Mor23$state_flag == "Downreg",]
rownames(dispersion.Mor23.down) <- dispersion.Mor23.down$GeneName
dispersion.Mor23.down$GeneName <- NULL
dispersion.Mor23.down$state_flag <- NULL
dispersion.Mor23.down <- reshape2::melt(dispersion.Mor23.down, id.vars = NULL)
dispersion.Mor23.down$identity <- rep("Mor23.d", each = 356)


dispersion.Mor23.up <- dispersion.Mor23[dispersion.Mor23$state_flag == "Upreg",]
rownames(dispersion.Mor23.up) <- dispersion.Mor23.up$GeneName
dispersion.Mor23.up$GeneName <- NULL
dispersion.Mor23.up$state_flag <- NULL
dispersion.Mor23.up <- reshape2::melt(dispersion.Mor23.up, id.vars = NULL)

dispersion.Mor23.up$identity <- rep("Mor23.up", each = 419)

dispersion.Mor23.all <- rbind(dispersion.Mor23.down,dispersion.Mor23.up)





#M71  

M71.df$variance.control <- apply(M71.df[c("F1", "F2", "F3")], 1, var)
M71.df$variance.exposed <- apply(M71.df[c("Ace1", "Ace2", "Ace3")], 1, var)

M71.df$mean.control <- apply(M71.df[c("F1", "F2", "F3")], 1, mean)
M71.df$mean.exposed <- apply(M71.df[c("Ace1", "Ace2", "Ace3")], 1, mean)

M71.df$control<- log1p(M71.df[,18]/M71.df[,20])
M71.df$exposed<- log1p(M71.df[,19]/M71.df[,21])


dispersion.M71 <- M71.df[,c(1,16,22:23)]

dispersion.M71.down <- dispersion.M71[dispersion.M71$state_flag == "Downreg",]
rownames(dispersion.M71.down) <- dispersion.M71.down$GeneName
dispersion.M71.down$GeneName <- NULL
dispersion.M71.down$state_flag <- NULL
dispersion.M71.down <- reshape2::melt(dispersion.M71.down, id.vars = NULL)
dispersion.M71.down$identity <- rep("M71.d", each= 752)

dispersion.M71.up <- dispersion.M71[dispersion.M71$state_flag == "Upreg",]
rownames(dispersion.M71.up) <- dispersion.M71.up$GeneName
dispersion.M71.up$GeneName <- NULL
dispersion.M71.up$state_flag <- NULL
dispersion.M71.up <- reshape2::melt(dispersion.M71.up, id.vars = NULL)
dispersion.M71.up$identity <- rep("M71.up", each = 645)



dispersion.d<- rbind(dispersion.M71.down,dispersion.Mor23.down)
dispersion.up <-  rbind(dispersion.M71.up,dispersion.Mor23.up)

dispersion.all <- rbind(dispersion.d,dispersion.up)




p8 <- ggplot(dispersion.all, aes(x = variable, y = value)) + 
  geom_violin(colour="black", fill="white")  +
  facet_grid(cols = vars(identity)) +
  stat_summary(fun.y = mean, geom = "crossbar", shape = 10, size = 1, color = "blue", aes(group = 1)) +
  stat_summary(fun.y = median, geom = "point", size = 2, color = "black") +
  labs(x = "Condition", y = "Dispersion") + 
  ylim(0,12)+
  theme_classic()

p8$data$identity <- factor(p8$data$identity, levels = c("Mor23.d","M71.d","Mor23.up","M71.up"))
p8$data <- p8$data[order(p8$data$identity),]

p8



