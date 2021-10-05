
setwd("N:/Novell/Luis Flores/mapping/")

##Input genes and its expression to test upregulation

exon.M71.upr <- read.csv("exon_M71_up.csv")
intron.M71.upr <- read.csv("intron_M71_up.csv")

exon.MOR23.upr <- read.csv("exon_MOR23_up.csv")
intron.MOR23.upr <- read.csv("intron_MOR23_up.csv")

# Import data set

d <- exon.MOR23.up



# Define Control and Exposed mice
control.mice <- c("D1","D3","D5")
exposed.mice <- c("L1","L3","L5")

# Define gene order
gene <- c("Cidea","Hagh","Nr4a1","S100a5","Srxn1")

#Values for Boxplot Rq values output from qPCR analysis
# slots: Define the number of columns in the boxplot
slots <- 5

#install.packages("beeswarm")
library(beeswarm)

#Create a dataset dn containing only Rq values normed by DMSO
averageDMSO <- colMeans(d[1:3,])

dn <- data.frame(t(apply(d[4:6,], 1, function(x) x/averageDMSO)))

# change levels of longd$gene
dn <- dn[gene]

# -----------------------------------------------------------------------------
#########################################################################################################################

#########################################################################################################################
#### Functions ####

# Melt data frame by gene and name varible and value names
melt_by_gene <- function(x, variable.name = "Mouse", value.name = "rqCT") {
  require(reshape2)
  x$Gene <- rownames(x)
  y <- melt(x, variable.name = variable.name, value.name = value.name)
  return(y)
}

# Assign asterisk to significance levels
asterisk <- function(x) {
  ifelse(x <=  0.001, "***",
         ifelse( x > 0.001 & x <= 0.01, "**",
                 ifelse( x > 0.01 & x <= 0.1, "*", "")))
}


#########################################################################################################################

#########################################################################################################################
#### Prepare data set for subsequent statistical analysis ####

# Transpose data so that genes are in rownames
td <- data.frame(t(d))

# Define conditions
Condition <- c("Control", "Exposed")

# Melt transposed data by Gene
df <- melt_by_gene(td, variable.name = "Mouse", value.name = "rqCT")

# Define df$Condition for each mouse
df$Condition <- ifelse(df$Mouse %in% control.mice, Condition[1],
                       ifelse(df$Mouse %in% exposed.mice, Condition[2], NA))

# Factorize df$Condition for statistical analysis
df$Condition <- factor(df$Condition, levels = Condition)

# Split data set by Gene in a large list
# Each list element will be a seperate data frame containing all the data for a given gene
df.list <- split(df, df$Gene)
#########################################################################################################################

#########################################################################################################################
#### Statistical analysis ####

# Test for normality using Shapiro-Wilk test
# if p-value < 0.05, the normality of the data set is rejected

# For controls
shapiro.p.value.Control <- unlist(lapply(df.list, function(x) { shapiro.test(x$rqCT[x$Condition == "Control"])$p.value }), recursive = F)
shapiro.p.value.Control <- shapiro.p.value.Control[gene]
shapiro.p.value.Control <- p.adjust(shapiro.p.value.Control, method="fdr", n= length(shapiro.p.value.Control))
shapiro.p.value.Control[shapiro.p.value.Control <= 0.1]

# For exposed
shapiro.p.value.Exposed <- unlist(lapply(df.list, function(x) { shapiro.test(x$rqCT[x$Condition == "Exposed"])$p.value }), recursive = F)
shapiro.p.value.Exposed <- shapiro.p.value.Exposed[gene]
shapiro.p.value.Exposed <- p.adjust(shapiro.p.value.Exposed, method="fdr", n= length(shapiro.p.value.Exposed))
shapiro.p.value.Exposed[shapiro.p.value.Exposed <= 0.1]

# Non parametric test :
# independent 2-group Mann-Whitney U Test 
# To be implemented if null hypothesis of shapiro.wilk test is rejected;
# that is, p-value <= 0.05
wilcoxon.p.value <- unlist(lapply(df.list, function(x) { wilcox.test(data = x, rqCT ~ Condition)$p.value }), recursive = F)
wilcoxon.p.value <- wilcoxon.p.value[gene]
wilcoxon.p.adj <- p.adjust(wilcoxon.p.value, method="fdr", n= length(wilcoxon.p.value))
wilcoxon.p.adj[wilcoxon.p.adj <= 0.1]

# Parametric test :
# Welsh two-sample t.test
# To be implemented if null hypothesis of shapiro.wilk test is not rejected;
# that is, p-value > 0.05
t.p.value <- unlist(lapply(df.list, function(x) { t.test(data = x, rqCT ~ Condition)$p.value }), recursive = F)
t.p.value <- t.p.value[gene]
t.p.adj <- p.adjust(t.p.value, method="fdr", n= length(t.p.value))
t.p.adj[t.p.adj <= 0.1]

# Pairwise t.test using lm (USE THIS ONE)
# Extract p.value of the model
lm.p.value <- unlist(lapply(df.list, function(x) { summary(lm(data = x, rqCT ~ Condition))$coeff[2,4] }), recursive = F)
# ...
# order p.values by the order of significantly downregulated genes from the RNAseq data
lm.p.value <- lm.p.value[gene]
#...
# Adjust p.value for multiple testing using 'FDR' correction
lm.p.adj <- p.adjust(lm.p.value, method="fdr", n= length(lm.p.value))
# ...
# Select genes having a significant p.adj.value
lm.p.adj[lm.p.adj <= 0.1]

# Run individual models as in 
# model <- lm(data = list$gene, rqCT ~Condition)
# plot(model) plots diagnostics of the model
#old_par <- par(mfrow=c(2,2))
#plot(model)
#par(old_par)

# Plot model residuals, to test if they are normally distributed
#hist(resid(model))

# Plot individual values
library(ggplot2)
ggplot(df.list$Nr4a1, aes(Condition, rqCT)) +
  geom_boxplot(aes(color=Condition)) +
  geom_point()

#########################################################################################################################

#########################################################################################################################
#########################################################################################################################

#########################################################################################################################

# -----------------------------------------------------------------------------


# colors
group_colors <- c("royalblue", "deeppink", "black")

# assign group number 
make_genegroup <- function(genenames) {
  genegroup <- c()
  group_i <- 1
  for (x in genenames) {
    if (startsWith(x, "X")) {
      group_i <- group_i + 1
    } else {
      genegroup <- append(genegroup, group_i)
    }
  }
  genegroup <- data.frame(genegroup, row.names = genenames[!startsWith(genenames, "X")])
  return(genegroup)
}

# convert to long format
longd <- reshape(dn, 
                 varying = colnames(dn), times = colnames(dn), 
                 v.names = "foldchange", timevar = "gene",
                 direction = "long", new.row.names = 1:(nrow(dn)*length(dn)))

# discard id column
longd <- longd[ , !(names(longd) == "id")]

# discard Xs
longd <- longd[!startsWith(longd$gene, "X"),]

# discard rows with NA
longd <- longd[!is.na(longd$foldchange),]

############################
# add genegroup
genegroup <- make_genegroup(colnames(dn))
genelabel_colors <- ifelse(rownames(genegroup) %in% names(lm.p.adj[lm.p.adj <= 0.1]), group_colors[1], group_colors[2])

longd$group <- genegroup[longd$gene,]

# add gaps
gap_level <- seq(1.5, length.out = length(unique(genegroup$genegroup))-1)
remaining_slots <- slots - (length(genegroup$genegroup)+length(gap_level))
if (remaining_slots > 0) {
  gap_level <- append(gap_level, seq(gap_level[-1]+1,10.5,1))
} else if (remaining_slots < 0) {
  stop("Error: Not enough slots, please increase slots variable")
}
gap_names <- paste(rep("_", length(gap_level)), seq(1, length.out=length(gap_level)), sep="")
gaps <- data.frame(gap_names, rep(4, length(gap_level)), gap_level)
colnames(gaps) <- c("gene", "foldchange", "group")

# merge
longd <- rbind(longd, gaps)

# group as factor
longd$group=ordered(longd$group)

# genes as factor
gaps = data.frame(gaps$group, row.names=gaps$gene)
colnames(gaps) <- c("genegroup")
genegroup <- rbind(genegroup, gaps)
longd$gene <- factor(longd$gene, levels = row.names(genegroup)[order(genegroup$genegroup)])



ymin <- 0
ymax <- 4

yaxis_tick <- log(c(1, seq(2, 20, 2)), 2)
yaxis_labels <- ifelse(is.element(c(1, seq(2, 20, 2)), c(1, 2, 4, 8, 16)), as.character(c(1, seq(2, 20, 2))), "")
plot(ymin+1, ymax+1, ylim=c(ymin, ymax), xlim=c(0.5, slots+0.5), type='n', bty='n', xaxt="n", yaxt="n", xlab="", ylab="Rel. fold change")

# gene labels
xaxis_tick <- which(genegroup[order(genegroup$genegroup),] == round(genegroup[order(genegroup$genegroup),]))
xaxis_labels <- grep("^_", rownames(genegroup)[order(genegroup$genegroup)], value=T, invert=T)

# grey zone
#polygon(x=c(-1, -1, slots+2, slots+2), y=c(-0.322, 0.322, 0.322, -0.322), col="grey20", border=NA)

# middle line
abline(h=0, lty=2)

# x axis
axis(1, at=xaxis_tick, labels=F)
mtext(xaxis_labels, side=1, at=xaxis_tick, las=3, col=genelabel_colors, line=1)

# y axis
axis(4.32, at=yaxis_tick, labels=yaxis_labels, las=1, side = 2)

P <- boxplot(log(foldchange,2) ~ gene, data=longd, add=T, yaxt="n", xaxt="n",
             # median parameters    
             medcol="white",
             
             # staple and whisker parameters
             whisklty=0, staplelty=0, outlty=0, 
             
             # box parameters
             boxwex=0.45, boxfill=c("white"),
             
             # outliers are not drawn
             outline=F,
             
             las=1)


# add mean bars
mean_foldchange <- aggregate(longd$foldchange, list(longd$gene), mean)
colnames(mean_foldchange) <- c("gene", "foldchange")
mean_colors <- ifelse(rownames(genegroup) %in% names(lm.p.adj[lm.p.adj <= 0.1]), group_colors[1], group_colors[2])
segments(1:length(P$n)+0.3, log(mean_foldchange$foldchange,2), 1:length(P$n)-0.3, log(mean_foldchange$foldchange,2), lwd=5, lend=2, col=mean_colors)

# add vertical bar between groups
higher_group <- max(genegroup[genegroup$genegroup == round(genegroup$genegroup),])
n_group_separation <- length(gaps[gaps$genegroup < higher_group,])

abline(v = which(startsWith(P$names, "_"))[1:n_group_separation], lty="19", lwd=2)

# add points
beeswarm(log(foldchange,2) ~ gene, spacing=0.3, data=longd, pch=21, cex=1, add=T, bg="white")

# add asterisk representing statistical significance (following the lm computed padj values)
significance <- asterisk(lm.p.adj)

# add "" as many times as there are samples -1 to the 'significance' vector for each gene
# this is important because each element in the 'significance' vector corresponds 
# to one sample on the graph and not to a single gene
signif <- as.character(sapply(significance, function(x) { c(x, rep("", times = length(longd$gene[longd$gene %in% longd$gene[1]])-1))}))

# define y limit for asterisk display
# -3 corresponds to 0.125 on the graph
y <- 0.1

# add text to display level of significance for each gene on the graph
text(x=longd$gene, y = y, signif, cex=1.5)

