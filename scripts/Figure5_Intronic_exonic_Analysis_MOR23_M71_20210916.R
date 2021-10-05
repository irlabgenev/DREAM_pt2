#Required packages
require(DESeq2)
require(edgeR)
require(ggplot2)
require(ggrepel)
require(dplyr)
require(apeglm)
require(ggbeeswarm)
library(DESeq2)


#Input exonic files
setwd("N:/Novell/Luis Flores/mapping/")

counts.list.wf<- readRDS(file = "FACSseq_M71_20170405_MOR23_20170112_new_exon_counts_per_sample-remapping20210310_nonRNAs.rds")
M71 <- counts.list.wf$FACSseq_M71_20170405_exon_new_nonRNAs_20210917
MOR23 <- counts.list.wf$FACSseq_MOR23_20170112_exon_new_nonRNAs_20210917
exp_table_M71 <- "sample_listM71_1.txt"
exp_tableMOR23 <- "sample_listMOR23_1.txt"

'%ni%' <- Negate('%in%')



s2cMOR23 <- read.table(exp_tableMOR23,
                       header = TRUE,
                       stringsAsFactors = FALSE)
rownames(s2cMOR23) <- s2cMOR23$sample


s2cM71 <- read.table(exp_table_M71,
                     header = TRUE,
                     stringsAsFactors = FALSE)
rownames(s2cM71) <- s2cM71$sample


Olfr151 <- M71[grep("Olfr151", rownames(M71)),]
Olfr151 <- Olfr151[1,]

Olfr16 <- MOR23[grep("Olfr16", rownames(MOR23)),]
Olfr16 <- Olfr16[1,]

Olfr.genes <- MOR23[grep("Olfr", rownames(MOR23)),]

M71 <- M71[rownames(M71) %ni% rownames(Olfr.genes),]
MOR23 <- MOR23[rownames(MOR23) %ni% rownames(Olfr.genes),]

M71 <- rbind(M71,Olfr151)
MOR23 <- rbind(MOR23, Olfr16)

M71 <- M71[order(rownames(M71)),]
MOR23 <- MOR23[order(rownames(MOR23)),]

filterthresholdM71<- GetCountThreshold_generic(counts=M71,
                                               return.plot = TRUE)
filterthresholdMOR23 <- GetCountThreshold_generic(counts=MOR23,
                                                  return.plot = TRUE)
genes.to.keep.M71 <- rowSums(x = M71) != 0
counts.filtered.M71 <- M71[genes.to.keep.M71,]

# Normalize the data to determine threshold of gene expression
norm.counts.M71 <- (t(x = t(x = counts.filtered.M71) / colSums(x = counts.filtered.M71)) * 1e6)

genes.to.keep.MOR23 <- rowSums(x = MOR23) != 0
counts.filtered.MOR23 <- MOR23[genes.to.keep.MOR23,]

# Normalize the data to determine threshold of gene expression
norm.counts.MOR23 <- (t(x = t(x =  counts.filtered.MOR23) / colSums(x = counts.filtered.MOR23)) * 1e6)





filteredM71 <- FilterCountsMatrix(counts=norm.counts.M71,
                                  condition.1.names = s2cM71$sample[1:3],
                                  condition.2.names = s2cM71$sample[4:6],
                                  count.threshold = filterthresholdM71$count.threshold)

filteredMOR23 <- FilterCountsMatrix(counts=norm.counts.MOR23,
                                    condition.1.names = s2cMOR23$sample[1:3],
                                    condition.2.names = s2cMOR23$sample[4:6],
                                    use.pattern = FALSE,
                                    condition.1.pattern =,
                                    condition.2.pattern =,
                                    count.threshold = filterthresholdMOR23$count.threshold)

genes.filtered.MOR23 <- rownames(filteredMOR23$filtered.counts)
genes.filtered.M71 <- rownames(filteredM71$filtered.counts)

MOR23 <- MOR23[rownames(MOR23) %in% genes.filtered.MOR23,]
M71 <- M71[rownames(M71) %in% genes.filtered.M71,]

ddsMOR23 <- ImportCountsMatrixDESeq2(counts = MOR23, 
                                     colData = s2cMOR23,
                                     condition.levels = c("DMSO",
                                                          "Lyral"))

ddsM71 <- ImportCountsMatrixDESeq2(counts = M71, 
                                   colData = s2cM71,
                                   condition.levels = c("DMSO",
                                                        "Aceto"))
resM71 <- RunDESeq2(dds = ddsM71, 
                    #                    normGenes = filteredM71$genes, 
                    #   genestokeep = filteredM71$genes,
                    colcontrolsamples = s2cM71$sample[1:3],
                    colexposedsamples = s2cM71$sample[4:6],
                    coef = "condition_Aceto_vs_DMSO")

resMOR23 <- RunDESeq2(dds = ddsMOR23, 
                      #                    normGenes = filteredM71$genes, 
                      #   genestokeep = filteredM71$genes,
                      colcontrolsamples = s2cMOR23$sample[1:3],
                      colexposedsamples = s2cMOR23$sample[4:6],
                      coef = "condition_Lyral_vs_DMSO")


#Mutate and highlight significant/ns genes

resM71mutatedLFC.exon <- MutateCountsdfLFC(res = resM71$resLFCMerged,
                                           threshold.value = 0.1)

resMOR23mutatedLFC.exon <- MutateCountsdfLFC(res = resMOR23$resLFCMerged,
                                             threshold.value = 0.1)



#Scatterplots

#M71
x <-GetScatterplot(Mutated.df = resM71mutatedLFC.exon,
                   OlfrOfInterest = "Olfr151",
                   GeneOfInterest = "")

x$p2$data$state_flag <- factor(x$p2$data$state_flag, levels = c("n.s.","Downreg","Upreg"))
x$p2$data <- x$p2$data[order(x$p2$data$state_flag),]

x


#MOR23

y <-GetScatterplot(Mutated.df = resMOR23mutatedLFC.exon,
                   OlfrOfInterest = "Olfr16",
                   GeneOfInterest = "")


y$p2$data$state_flag <- factor(y$p2$data$state_flag, levels = c("n.s.","Downreg","Upreg"))
y$p2$data <- y$p2$data[order(y$p2$data$state_flag),]

y


#Input Intronic files


setwd("N:/Novell/Luis Flores/mapping/")

Intron.count.list<- readRDS(file = "FACSseq_M71_20170405_MOR23_20170112_new_intron_counts_per_sample-remapping20210310_nonRNAs.rds")
Intron.count.list.M71 <- Intron.count.list$FACSseq_M71_20170405_intron_new_nonRNAs_20210917
Intron.count.list.MOR23 <- Intron.count.list$FACSseq_MOR23_20170112_intron_new_nonRNAs_20210917
exp_table_M71 <- "sample_listM71_1.txt"
exp_tableMOR23 <- "sample_listMOR23_1.txt"


'%ni%' <- Negate('%in%')

s2cMOR23 <- read.table(exp_tableMOR23,
                       header = TRUE,
                       stringsAsFactors = FALSE)
rownames(s2cMOR23) <- s2cMOR23$sample


s2cM71 <- read.table(exp_table_M71,
                     header = TRUE,
                     stringsAsFactors = FALSE)
rownames(s2cM71) <- s2cM71$sample







intron.M71 <- Intron.count.list.M71[rownames(Intron.count.list.M71) %in% rownames(M71), ]
intron.MOR23 <- Intron.count.list.MOR23[rownames(Intron.count.list.MOR23) %in%  rownames(MOR23),]

M71.forintrons <- M71[rownames(M71) %in% rownames(intron.M71),]
MOR23.forintrons <- MOR23[rownames(MOR23) %in% rownames(intron.MOR23),]

M71.forintrons <- M71.forintrons[order(rownames(M71.forintrons)),]
MOR23.forintrons <- MOR23.forintrons[order(rownames(MOR23.forintrons)),]

intron.M71 <- (intron.M71-M71.forintrons)
intron.MOR23 <- (intron.MOR23-MOR23.forintrons)



intron.M71 <- apply(intron.M71, 2 , function(x) ifelse(x < 0, 0, x))
intron.MOR23 <- apply(intron.MOR23, 2 , function(x) ifelse(x < 0, 0, x))




intron.ddsMOR23 <- ImportCountsMatrixDESeq2(counts = intron.MOR23, 
                                     colData = s2cMOR23,
                                     condition.levels = c("DMSO",
                                                          "Lyral"))

intron.ddsM71 <- ImportCountsMatrixDESeq2(counts = intron.M71, 
                                   colData = s2cM71,
                                   condition.levels = c("DMSO",
                                                        "Aceto"))
intron.resM71 <- RunDESeq2(dds = intron.ddsM71, 
                    #                    normGenes = filteredM71$genes, 
                    #   genestokeep = filteredM71$genes,
                    colcontrolsamples = s2cM71$sample[1:3],
                    colexposedsamples = s2cM71$sample[4:6],
                    coef = "condition_Aceto_vs_DMSO")

intron.resMOR23 <- RunDESeq2(dds = intron.ddsMOR23, 
                      #                    normGenes = filteredM71$genes, 
                      #   genestokeep = filteredM71$genes,
                      colcontrolsamples = s2cMOR23$sample[1:3],
                      colexposedsamples = s2cMOR23$sample[4:6],
                      coef = "condition_Lyral_vs_DMSO")




#Mutate and highlight significant/ns genes

intron.resM71mutatedLFC <- MutateCountsdfLFC(res = intron.resM71$resLFCMerged,
                                      threshold.value = 0.1)

intron.resMOR23mutatedLFC <- MutateCountsdfLFC(res = intron.resMOR23$resLFCMerged,
                                        threshold.value = 0.1)






M71.analysis <- resM71mutatedLFC.exon
M71.intronic <-  intron.resM71mutatedLFC

MOR23.intronic <- intron.resMOR23mutatedLFC
MOR23.analysis <- resMOR23mutatedLFC.exon



# M71
M71.exon <- M71.analysis[M71.analysis$GeneName %in% M71.intronic$GeneName,]
M71.intron <- M71.intronic[M71.intronic$GeneName %in% M71.analysis$GeneName,]
M71.FC.intron <- M71.intron[,c(3,16)]
M71.complete <- cbind(M71.exon,M71.FC.intron)
M71.complete <- M71.complete[,c(1,3,16,18:19)]
colnames(M71.complete) <- c("Genename", "exon.FC","exon.state_flag","intron.FC","intron.state_flag")


M71.complete$color <-  ifelse(M71.complete$exon.state_flag == "Downreg" & M71.complete$intron.state_flag == "Downreg","blue",
                              ifelse(M71.complete$exon.state_flag == "n.s." & M71.complete$intron.state_flag == "Downreg", "cyan",
                                     ifelse(M71.complete$exon.state_flag == "Downreg" & M71.complete$intron.state_flag == "n.s.", "cyan",
                                            ifelse(M71.complete$exon.state_flag == "Upreg" & M71.complete$intron.state_flag == "Upreg", "red",
                                                   ifelse(M71.complete$exon.state_flag == "Upreg" & M71.complete$intron.state_flag == "n.s.", "tomato",
                                                          ifelse(M71.complete$exon.state_flag == "n.s." & M71.complete$intron.state_flag == "Upreg", "tomato",
                                                                 "grey"))))))

p1 <- ggplot(data = M71.complete,
             aes(exon.FC, intron.FC)) +
  geom_point(size = 0.5) +
  geom_point(aes(colour = color),
             size = 1) +
  scale_color_manual(values=c("blue" = "blue",
                              "tomato" = "tomato",
                              "cyan" = "cyan",
                              "grey" = "grey",
                              "red" = "red"),
                     guide=FALSE) +
  labs(x = "Exonic rel. FoldChange (log2)", 
       y = "Intronic rel. FoldChange (log2)")  +
  #geom_text_repel(data = MOR23.complete[MOR23.complete$Genename %in% genes.toplot,],
   #               aes(label = Genename)) + 
  #geom_abline(intercept = 0,
  #           colour = "blue") +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #xlim(-12,12)+
  #ylim(-12,12)+
  theme(aspect.ratio = 1) +
  theme_classic()



p1$data$color <- factor(p1$data$color, levels =  c("grey","tomato","cyan","red","blue"))
p1$data <- p1$data[order(p1$data$color),]

p1





#MOR23
MOR23.exon <- MOR23.analysis[MOR23.analysis$GeneName %in% MOR23.intronic$GeneName,]
MOR23.intron <- MOR23.intronic[MOR23.intronic$GeneName %in% MOR23.analysis$GeneName,]
MOR23.FC.intron <- MOR23.intron[,c(3,16)]
MOR23.complete <- cbind(MOR23.exon,MOR23.FC.intron)
MOR23.complete <- MOR23.complete[,c(1,3,16,18:19)]
colnames(MOR23.complete) <- c("Genename", "exon.FC","exon.state_flag","intron.FC","intron.state_flag")


MOR23.complete$color <-  ifelse(MOR23.complete$exon.state_flag == "Downreg" & MOR23.complete$intron.state_flag == "Downreg","blue",
                                ifelse(MOR23.complete$exon.state_flag == "n.s." & MOR23.complete$intron.state_flag == "Downreg", "cyan",
                                       ifelse(MOR23.complete$exon.state_flag == "Downreg" & MOR23.complete$intron.state_flag == "n.s.", "cyan",
                                              ifelse(MOR23.complete$exon.state_flag == "Upreg" & MOR23.complete$intron.state_flag == "Upreg", "red",
                                                     ifelse(MOR23.complete$exon.state_flag == "Upreg" & MOR23.complete$intron.state_flag == "n.s.", "tomato",
                                                            ifelse(MOR23.complete$exon.state_flag == "n.s." & MOR23.complete$intron.state_flag == "Upreg", "tomato",
                                                                   "grey"))))))

##Function to plot

p2 <- ggplot(data = MOR23.complete,
             aes(exon.FC, intron.FC)) +
  geom_point(size = 0.5) +
  geom_point(aes(colour = color),
             size = 1) +
  scale_color_manual(values=c("blue" = "blue",
                              "tomato" = "tomato",
                              "cyan" = "cyan",
                              "grey" = "grey",
                              "red" = "red"),
                     guide=FALSE)  +
  labs(x = "Exonic rel. FoldChange (log2)", 
       y = "Intronic rel. FoldChange (log2)")  +
  geom_abline(intercept = 0,
              colour = "blue") +  
  #geom_text_repel(data = MOR23.complete[MOR23.complete$Genename %in% genes.toplot,],
  #                                   aes(label = Genename)) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  #xlim(-12,12)+
  #ylim(-12,12)+
  theme(aspect.ratio = 1) +
  theme_classic()
#ggtitle("Distribution of the downregulated genes from M71 & MOR23")





p2$data$color <- factor(p2$data$color, levels = c("grey","tomato","cyan","red","blue"))
p2$data <- p2$data[order(p2$data$color),]

p2

genes.toplot <- c("Ebf4","Lingo2","Tom1","Cidea","S100a5","Srxn1")

M71.complete[M71.complete$Genename %in% genes.toplot,]




























