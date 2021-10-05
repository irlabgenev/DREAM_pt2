#Required packages
require(DESeq2)
require(edgeR)
require(ggplot2)
require(ggrepel)
require(dplyr)
require(apeglm)
require(ggbeeswarm)
library(DESeq2)


#Input functions
#source(file = "N:/Novell/___Projects___/DREAM/experiments/Scrnaseq_moe/MOE_scrna_seurat/DESeq2_helper_functions.R")

#Input files
setwd("N:/Novell/Luis Flores/mapping/")

counts.list.wf<- readRDS(file = "FACSseq_M71_20170405_MOR23_20170112_gene_counts_per_sample-remapping20210310.rds")
counts.list.wf.M71 <- counts.list.wf$FACSseq_M71_20170405
counts.list.wf.MOR23 <- counts.list.wf$FACSseq_MOR23_20170112
exp_table_M71 <- "sample_listM71_1.txt"
exp_tableMOR23 <- "sample_listMOR23_1.txt"

'%ni%' <- Negate('%in%')
spikesM71 <- counts.list.wf.M71[grep("ERCC", rownames(counts.list.wf.M71)),]
M71 <- counts.list.wf.M71[rownames(counts.list.wf.M71) %ni% rownames(spikesM71),]


spikesMor23 <- counts.list.wf.MOR23[grep("ERCC", rownames(counts.list.wf.MOR23)),]
MOR23 <- counts.list.wf.MOR23[rownames(counts.list.wf.MOR23) %ni% rownames(spikesMor23),]



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

#Data visualization
#PCA
GetPCAplot(dds = resM71$dds,
           intgroup = "condition")

GetPCAplot(dds = resMOR23$dds,
           intgroup = "condition")


#Mutate and highlight significant/ns genes

resM71mutatedLFC <- MutateCountsdfLFC(res = resM71$resLFCMerged,
                                      threshold.value = 0.1)

resMOR23mutatedLFC <- MutateCountsdfLFC(res = resMOR23$resLFCMerged,
                                        threshold.value = 0.1)



#Scatterplots

#M71
x <-GetScatterplot(Mutated.df = resM71mutatedLFC,
                   OlfrOfInterest = "Olfr151",
                   GeneOfInterest = "")

x$p2$data$state_flag <- factor(x$p2$data$state_flag, levels = c("n.s.","Downreg","Upreg"))
x$p2$data <- x$p2$data[order(x$p2$data$state_flag),]

x


#MOR23

y <-GetScatterplot(Mutated.df = resMOR23mutatedLFC,
                   OlfrOfInterest = "Olfr16",
                   GeneOfInterest = "")


y$p2$data$state_flag <- factor(y$p2$data$state_flag, levels = c("n.s.","Downreg","Upreg"))
y$p2$data <- y$p2$data[order(y$p2$data$state_flag),]

y





#VolcanoPlot
GetVolcanoplot(Mutated.df = M71,
               threshold = 0.1,
               inf.limit = -11,
               sup.limit = 11,
               OlfrOfInterest = "Olfr151",
               GenesOfInterest = "Mustn1")

GetVolcanoplot(Mutated.df = Mor23.df,
               threshold = 0.1,
               inf.limit = -11,
               sup.limit =11,
               OlfrOfInterest = "Olfr16",
               GenesOfInterest = "Mustn1")



###FUNCTIONS you have to run to perform the Analysis

peakfinder <- function(d){ 
  dh <- hist(d, plot = FALSE) 
  ins <- dh[["counts"]]
  nbins <- length(ins) 
  ss <- which(rank(ins) %in% seq(from = nbins-2, to = nbins)) ## pick the top 3 intensities 
  dh[["mids"]][ss]
}


GetCountThreshold_generic <- function(counts,
                                      scale.factor = 1e6,
                                      return.plot = TRUE){
  # Filter all genes not expressed in the dataset
  genes.to.keep <- rowSums(x = counts) != 0
  counts.filtered <- counts[genes.to.keep,]
  
  # Normalize the data to determine threshold of gene expression
  norm.counts <- (t(x = t(x = counts.filtered) / colSums(x = counts.filtered)) * scale.factor)
  log.mean.norm.counts <- log(x = rowMeans(x = norm.counts) + 1)
  
  # Determine threshold of gene expression from the density distribution
  density.data <- density(x = log.mean.norm.counts)
  local.minima <- density.data$x[which(x = diff(x = sign(x = diff(x = density.data$y))) == 2) + 1][1]
  count.threshold <- expm1(x = round(x = local.minima))
  
  if (return.plot) {
    
    plot.density <- ggplot2::ggplot(data = as.data.frame(x = log.mean.norm.counts),
                                    mapping = ggplot2::aes(x = log.mean.norm.counts,
                                                           y = ..density..)) +
      ggplot2::geom_density() + 
      ggplot2::geom_vline(xintercept = log(x = count.threshold + 1)) +
      ggplot2::ggtitle(label = "Density distribution of mean expression values")
    
    return(list(counts.filtered = counts.filtered,
                genes.to.keep = genes.to.keep,
                norm.counts = norm.counts,
                log.mean.norm.counts = log.mean.norm.counts,
                density.data = density.data,
                count.threshold = count.threshold,
                plot.density = plot.density))
  } else {
    return(list(counts.filtered = counts.filtered,
                genes.to.keep = genes.to.keep,
                norm.counts = norm.counts,
                log.mean.norm.counts = log.mean.norm.counts,
                density.data = density.data,
                count.threshold = count.threshold))
  }
}


#Filter a count matrix based on the sum of the count values of a condition's triplicate (n=3) being superior
#to the count.threshold previously determined






FilterCountsMatrix <- function(counts,
                               condition.1.names,
                               condition.2.names,
                               use.pattern = FALSE,
                               condition.1.pattern,
                               condition.2.pattern,
                               count.threshold){
  conditions <- colnames(x = counts)
  
  if (use.pattern) {
    if (missing(x = condition.1.pattern) | missing(x = condition.2.pattern)) {
      stop("One of 'condition.1.pattern' or 'condition.2.pattern' is missing:\n
           Please define both patterns to proceed")
    }
    condition.1.names <- grep(pattern = condition.1.pattern,
                              x = conditions,
                              value = TRUE)
    condition.2.names <- grep(pattern = condition.2.pattern,
                              x = conditions,
                              value = TRUE)
  }
  
  rows.keep <- apply(X = counts,
                     MARGIN = 1,
                     FUN = function(counts.matrix,
                                    condition.1,
                                    condition.2) { 
                       sum(counts.matrix[condition.1] > count.threshold) >=3 | 
                         sum(counts.matrix[condition.2] > count.threshold) >=3},
                     condition.1 = condition.1.names,
                     condition.2 = condition.2.names)
  
  #We filter out non-expressed genes, by requiring more than "filter_treshold" reads in at least 3 samples for each gene
  filtered.counts <- counts[rows.keep,]
  filtered.counts <- data.frame(apply(X = filtered.counts,
                                      MARGIN = 2,
                                      FUN = function(count.matrix) {
                                        storage.mode(count.matrix) <- 'integer'; count.matrix}))
  
  counts.rownames <- rownames(x = filtered.counts)
  genes <- counts.rownames
  
  return(list(filtered.counts = filtered.counts, 
              genes = genes))}





#Make a dds object from filtered counts 
ImportCountsMatrixDESeq2 <- function(counts, 
                                     colData,
                                     design.formula = "~ condition",
                                     condition.levels){
  if (!all(colnames(x = counts) == rownames(x = colData))) {
    counts <- counts[, rownames(x = colData)]
  }
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design = as.formula(object = design.formula))
  dds$condition <- factor(x = dds$condition, 
                          levels = condition.levels)
  return(dds)
}


RunDESeq2 <- function(dds, 
                      #                      normGenes, 
                      #genestokeep, 
                      colcontrolsamples, 
                      colexposedsamples,
                      coef){
  # option1: normalize with normalization genes
  
  # dds <- DESeq2::estimateSizeFactors(object = dds, 
  #                                 controlGenes=match(normGenes,
  #                                                   rownames(dds)))
  #  dds <- DESeq2::estimateDispersions(object = dds)
  # dds <- DESeq2::nbinomWaldTest(object = dds, 
  #                              betaPrior = FALSE)
  # option2: or dont...
  dds <- DESeq(dds)
  
  # running Differential Expression Analysis
  res <- results(object = dds, 
                 cooksCutoff=TRUE,
                 independentFiltering=TRUE )
  #res <- res[genestokeep,]
  res <- res[order(res$log2FoldChange),]
  CountsPerSample <- counts(dds, 
                            normalized=TRUE)
  resMerged <- merge(res, 
                     CountsPerSample,
                     by=0,
                     col.names=TRUE)
  colnames(resMerged)[1] <- "GeneName"
  rownames(resMerged) <-  res$GeneName
  #Add log2 value for DMSO and Lyral BaseMeans
  resMerged$control_log2mean <- log2(rowMeans(resMerged[colcontrolsamples]+1))
  resMerged$exposed_log2mean <- log2(rowMeans(resMerged[colexposedsamples]+1)) 
  #Perform Log fold change shrinkage for visualization and ranking
  resLFC <- lfcShrink(dds = dds,
                      coef =coef,
                      type = "apeglm")
  resLFCMerged <- merge(resLFC, 
                        CountsPerSample,
                        by=0,
                        col.names=TRUE)
  colnames(resLFCMerged)[1] <- "GeneName"
  rownames(resLFCMerged) <-  resLFCMerged$GeneName
  resLFCMerged$control_log2mean <- log2(rowMeans(resLFCMerged[colcontrolsamples]+1))
  resLFCMerged$exposed_log2mean <- log2(rowMeans(resLFCMerged[colexposedsamples]+1)) 
  
  return(list(dds = dds,
              res = res,
              resMerged = resMerged,
              resLFC =resLFC,
              resLFCMerged = resLFCMerged))
}

MutateCountsdfLFC <- function(res, 
                              threshold.value){
  row.has.na <- apply(res,
                      1,
                      function(x){
                        any(is.na(x))
                      })
  sum(row.has.na) 
  df <- res[!row.has.na,]
  #Highlight rows with a padj<0.05
  df1 <- mutate(df,
                significance=ifelse(df$padj<threshold.value & df$log2FoldChange <(- 0.5),
                                    paste0("padj<",
                                           threshold.value),
                                    ifelse(df$padj<threshold.value & df$log2FoldChange > 0.5,
                                           paste0("padj<",
                                                  threshold.value),
                                           "Not significant")))
  df1$significance <- factor(df1$significance,
                             levels=c(paste0("padj<",
                                             threshold.value),
                                      "Not significant"))
  df1$state_flag <- ifelse(df1$log2FoldChange > 0.5 & df1$significance %in% paste0("padj<", threshold.value),
                           "Upreg",
                           ifelse(df1$log2FoldChange < 0.5 & df1$significance %in% paste0("padj<", threshold.value),
                                  "Downreg",
                                  "n.s."))
  df1$color_flag <- ifelse(df1$log2FoldChange > 0.5 & df1$significance %in% paste0("padj<", threshold.value),
                           "tomato",
                           ifelse(df1$log2FoldChange < 0.5 & df1$significance %in% paste0("padj<", threshold.value),
                                  "blue",
                                  "grey"))
  return(ResMutated = df1)
}



#PCA plot
GetPCAplot <- function(dds,
                       intgroup){
  rld <- rlog(object = dds)
  pca.data <- DESeq2::plotPCA(object = rld, 
                              ntop = nrow(dds), 
                              returnData = TRUE,
                              intgroup = intgroup) 
  
  percentVar <- round(100*attr(pca.data, "percentVar"))
  
  p1 <- ggplot(pca.data, aes(PC1, 
                             PC2,
                             color = condition)) +
    geom_point(size=3) +
    geom_text_repel(aes(label = name)) +
    xlab(paste0("PC1: ",
                percentVar[1],
                "% variance")) +
    ylab(paste0("PC2: ",
                percentVar[2],
                "% variance")) +
    coord_fixed() +
    theme_bw()
  return(p1)
}




#Scatter plot of DEseq2-normalized counts
GetScatterplot <- function(Mutated.df, 
                           GeneOfInterest,
                           OlfrOfInterest) {
  #discard row with na (pvalue) and switch to long format 
  max.rna <- ceiling(max(Mutated.df$control_log2mean,
                         Mutated.df$exposed_log2mean))
  totalupreg <- length(which(x = Mutated.df$state_flag == "Upreg") == T)
  totaldowreg <- length(which(x = Mutated.df$state_flag == "Downreg") == T)
  totalNS <- length(which(x = Mutated.df$state_flag == "n.s.") == T)
  
  p2 <- ggplot(data = Mutated.df,
               aes(exposed_log2mean,
                   control_log2mean)) +
    geom_point(size = 0.5) +
    geom_point(aes(colour = state_flag),
               size = 1) +
    scale_color_manual(values=c("Upreg" = "tomato",
                                "Downreg" = "blue",
                                "n.s." = "grey"),
                       guide=FALSE) +
    geom_text_repel(data = Mutated.df[Mutated.df$GeneName == OlfrOfInterest,],
                    aes(label = GeneName))  +
    labs(x = "Odorant-exposed mean norm. counts (log2)", 
         y = "DMSO-exposed mean norm. counts (log2)") + 
    expand_limits(x = c(0, max.rna),
                  y = c(0, max.rna)) +
    theme(aspect.ratio = 1) +
    theme_classic()
  return(list(totalupreg = totalupreg,
              totaldowreg = totaldowreg,
              totalNS = totalNS,
              p2 = p2))
}

#Volcano plot of DEseq2-normalized counts

GetVolcanoplot <- function(Mutated.df,
                           threshold,
                           inf.limit,
                           sup.limit,
                           GenesOfInterest,
                           OlfrOfInterest) {
  
  plot1 <- ggplot(data = Mutated.df,
                  mapping = aes(log2FoldChange,
                                -log10(padj))) +
    geom_hline(yintercept = -log10(0.1),
               col = "green3",
               lty = 2) +
    geom_point(mapping = aes(x = log2FoldChange,
                             y = -log10(padj),
                             colour= state_flag)) +
    scale_color_manual(values = c("Upreg" = "tomato",
                                  "Downreg" = "blue",
                                  "n.s." = "grey")) +
    
    geom_text_repel(data = subset(Mutated.df,
                                  padj < threshold & (log2FoldChange > sup.limit | log2FoldChange < (inf.limit))),
                    aes(label = GeneName),
                    size = 3,
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines")) +
    geom_text_repel(data = Mutated.df[Mutated.df$GeneName == OlfrOfInterest,],
                    aes(label = GeneName),
                    size = 3,
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines")) +
    geom_text_repel(data = Mutated.df[Mutated.df$GeneName == GenesOfInterest,],
                    aes(label = GeneName),
                    size = 3,
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines")) +
    labs(x = "Foldchange(log2)",
         y = "-log10 adjusted pvalue") +
    theme_classic()
  return(plot1 = plot1)
}





