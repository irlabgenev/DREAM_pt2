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
setwd("N:/Novell/Luis Flores/DREAM_ethyiso_paper1")

ethyl.data.exon <- readRDS(file = "RNAseq_EthylIso_20170405_remappedLFH_20210925_exons_gtf_onlywith_proteincodingIRLAB_annotations.rds")
exp_table_ethyl <- "sample_listEthyl_RNAseq.txt"

'%ni%' <- Negate('%in%')




s2cethyl <- read.table(exp_table_ethyl,
                         header = TRUE,
                         stringsAsFactors = FALSE)
rownames(s2cethyl) <- s2cethyl$sample




filterthresholdethyl.exon<- GetCountThreshold_generic(counts=ethyl.data.exon,
                                                   return.plot = TRUE)


genes.to.keep.ethyl <- rowSums(x = ethyl.data.exon) != 0
counts.filtered.ethyl <- ethyl.data.exon[genes.to.keep.ethyl,]

# Normalize the data to determine threshold of gene expression
norm.counts.ethyl <- (t(x = t(x = counts.filtered.ethyl) / colSums(x = counts.filtered.ethyl)) * 1e6)




filtered.ethyl <- FilterCountsMatrix(counts=norm.counts.ethyl,
                                       condition.1.names = s2cethyl$sample[1:3],
                                       condition.2.names = s2cethyl$sample[4:6],
                                       count.threshold = filterthresholdethyl.exon$count.threshold)



genes.filtered.ethyl <- rownames(filtered.ethyl$filtered.counts)


ethyl.data.exon <- ethyl.data.exon[rownames(ethyl.data.exon) %in% genes.filtered.ethyl,]

ethyl.data.exon <- ethyl.data.exon[order(rownames(ethyl.data.exon)),]

ddsethyl.exon <- ImportCountsMatrixDESeq2(counts = ethyl.data.exon, 
                                       colData = s2cethyl,
                                       condition.levels = c("DMSO",
                                                            "Ethyl"))

resethyl.exon <- RunDESeq2(dds = ddsethyl.exon, 
                        #                    normGenes = filteredM71$genes, 
                        #   genestokeep = filteredM71$genes,
                        colcontrolsamples = s2cethyl$sample[1:3],
                        colexposedsamples = s2cethyl$sample[4:6],
                        coef = "condition_Ethyl_vs_DMSO")


#Data visualization
#PCA
GetPCAplot(dds = resethyl.exon$dds,
           intgroup = "condition")




#Mutate and highlight significant/ns genes

resethylmutatedLFC.exon <- MutateCountsdfLFC.ethyl(res = resethyl.exon$resLFCMerged,
                                          threshold.value = 0.1)





##Input the intron data

setwd("N:/Novell/Luis Flores/DREAM_ethyiso_paper1")

ethyl.data.intron <- readRDS(file = "RNAseq_EthylIso_20170405_remappedLFH_20210925_introns_gtf_onlywith_proteincodingIRLAB_annotations.rds")



ethyl.data.intron <- ethyl.data.intron[rownames(ethyl.data.intron) %in% rownames(ethyl.data.exon),]



ethyl.data.intron <- ethyl.data.intron[order(rownames(ethyl.data.intron)),]

ethyl.data.exon <- ethyl.data.exon[order(rownames(ethyl.data.exon)),]

ethyl.data.intron <- (ethyl.data.intron - ethyl.data.exon)

ethyl.data.intron <- apply(ethyl.data.intron, 2 , function(x) ifelse(x < 0, 0, x))


ddsethyl.intron <- ImportCountsMatrixDESeq2(counts = ethyl.data.intron, 
                                            colData = s2cethyl,
                                            condition.levels = c("DMSO",
                                                                 "Ethyl"))

resethyl.intron <- RunDESeq2(dds = ddsethyl.intron, 
                             #                    normGenes = filteredM71$genes, 
                             #   genestokeep = filteredM71$genes,
                             colcontrolsamples = s2cethyl$sample[1:3],
                             colexposedsamples = s2cethyl$sample[4:6],
                             coef = "condition_Ethyl_vs_DMSO")

intron.resethylmutatedLFC <- MutateCountsdfLFC.ethyl(res = resethyl.intron$resLFCMerged,
                                                     threshold.value = 0.1)


olfrs.ethyl <- c("Olfr741","Olfr1364","Olfr15","Olfr60","Olfr169","Olfr166","Olfr167")

Olfrs.exon <- resethylmutatedLFC.exon[resethylmutatedLFC.exon$GeneName %in% olfrs.ethyl, ]

Olfrs.intron <- intron.resethylmutatedLFC[intron.resethylmutatedLFC$GeneName %in% olfrs.ethyl,]

intron.downreg <- Olfrs.intron[,c(7:12)]
exon.downreg <- Olfrs.exon[,c(7:12)]

#Input this output in the statistical downregulation script

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

MutateCountsdfLFC.ethyl<- function(res, 
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
                significance=ifelse(df$padj<threshold.value,
                                    paste0("padj<",
                                           threshold.value), 
                                    "Not significant"))
  df1$significance <- factor(df1$significance,
                             levels=c(paste0("padj<",
                                             threshold.value),
                                      "Not significant"))
  df1$state_flag <- ifelse(df1$log2FoldChange > 0 & df1$significance %in% paste0("padj<", threshold.value),
                           "Upreg",
                           ifelse(df1$log2FoldChange < 0 & df1$significance %in% paste0("padj<", threshold.value),
                                  "Downreg",
                                  "n.s."))
  df1$color_flag <- ifelse(df1$log2FoldChange > 0 & df1$significance %in% paste0("padj<", threshold.value),
                           "tomato",
                           ifelse(df1$log2FoldChange < 0 & df1$significance %in% paste0("padj<", threshold.value),
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
    geom_text_repel(data = Mutated.df[Mutated.df$GeneName == GeneOfInterest,],
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





