#Required packages
require(DESeq2)
require(edgeR)
require(ggplot2)
require(ggrepel)
require(dplyr)
require(apeglm)
require(ggbeeswarm)

library(DESeq2)


#Input files
setwd("N:/Novell/Luis Flores/mapping/")

counts.list.wf<- readRDS(file = "FACSseq_M71_20170405_MOR23_20170112_gene_counts_per_sample-remapping20210310.rds")
counts.list.wf.M71 <- counts.list.wf$FACSseq_M71_20170405
counts.list.wf.MOR23 <- counts.list.wf$FACSseq_MOR23_20170112

'%ni%' <- Negate('%in%')


spikesM71 <- counts.list.wf.M71[grep("ERCC", rownames(counts.list.wf.M71)),]
M71 <- counts.list.wf.M71[rownames(counts.list.wf.M71) %ni% rownames(spikesM71),]


spikesMor23 <- counts.list.wf.MOR23[grep("ERCC", rownames(counts.list.wf.MOR23)),]
MOR23 <- counts.list.wf.MOR23[rownames(counts.list.wf.MOR23) %ni% rownames(spikesMor23),]

All <- cbind(MOR23,M71)

Olfr151 <- All[grep("Olfr151", rownames(All)),]
Olfr151 <- Olfr151[1,]

Olfr16 <- All[grep(("Olfr16"), rownames(All)),]
Olfr16 <- Olfr16[1,]


Olfr.genes <- All[grep(("Olfr"), rownames(All)),]


All <- All[rownames(All) %ni% rownames(Olfr.genes),]




All <- rbind(All, Olfr151)
All <- rbind(All, Olfr16)


exp_table_all <- "sample_list.txt"

s2call <- read.table(exp_table_all,
                     header = TRUE,
                     stringsAsFactors = FALSE)
rownames(s2call) <- s2call$short_name



filterthresholdall <- GetCountThreshold_generic(counts=All,
                                                return.plot = TRUE)
genes.to.keep.all <- rowSums(x = All) != 0
counts.filtered.all <- All[genes.to.keep.all,]

# Normalize the data to determine threshold of gene expression
norm.counts.all <- (t(x = t(x = counts.filtered.all) / colSums(x = counts.filtered.all)) * 1e6)


filteredAll<- FilterCountsMatrixAll_D1(counts=norm.counts.all,
                                    condition.1.names = s2call$short_name[1:3],
                                    condition.2.names = s2call$short_name[4:6],
                                    condition.3.names = s2call$short_name[7:9],
                                    condition.4.names = s2call$short_name[10:12],
                                    count.threshold = filterthresholdall$count.threshold)

genesall.tokeep <- rownames(filteredAll$filtered.counts)

All <- All[rownames(All) %in% genesall.tokeep, ]

ddsall <- ImportCountsMatrixDESeq2(counts = All, 
                                   colData = s2call,
                                   condition.levels = c("DMSOA",
                                                        "Aceto",
                                                        "DMSOL",
                                                        "Lyral"))

GetPCAplot(dds = ddsall,
           intgroup = "condition")




#Functions to run 

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




FilterCountsMatrixAll_D1 <- function(counts,
                                  condition.1.names,
                                  condition.2.names,
                                  condition.3.names,
                                  condition.4.names,
                                  count.threshold){
  conditions <- colnames(x = counts)
  
  
  rows.keep <- apply(X = counts,
                     MARGIN = 1,
                     FUN = function(counts.matrix,
                                    condition.1,
                                    condition.2,
                                    condition.3,
                                    condition.4) { 
                       sum(counts.matrix[condition.1] > count.threshold) >=3|  
                         sum(counts.matrix[condition.2] > count.threshold) >=3|
                         sum(counts.matrix[condition.3] > count.threshold) >=3| 
                         sum(counts.matrix[condition.4] > count.threshold) >=3},
                     condition.1 = condition.1.names,
                     condition.2 = condition.2.names,
                     condition.3 = condition.3.names,
                     condition.4 = condition.4.names)
  
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
    theme_classic()
  return(p1)
}
