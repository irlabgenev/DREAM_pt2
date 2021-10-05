
library(BiocStyle)
library(rmarkdown)
library(geneplotter)
library(ggplot2)
library(plyr)
library(LSD)
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(stringr)
library(topGO)
library(genefilter)
library(biomaRt)
library(dplyr)
library(EDASeq)
library(fdrtool)
library(org.Mm.eg.db)
library('biomaRt')


setwd("N:/Novell/___Projects___/DREAM/experiments/Scrnaseq_moe/MOE_scrna_seurat/")


MOR23.df <- resMOR23mutatedLFC
M71.df <- resM71mutatedLFC

Mor23.topgo <- MOR23.df[MOR23.df$GeneName %in% M71.df$GeneName,]
M71.topgo <- M71.df[M71.df$GeneName %in% MOR23.df$GeneName,]

Mor23.topgo <- Mor23.topgo[,2, drop= F]
M71.topgo <- M71.topgo[,2,drop=F]

all.topgo <- cbind(Mor23.topgo,M71.topgo)
colnames(all.topgo) <- c("Mor23", "M71")

all.topgo$baseMean <- apply(all.topgo, 1, mean)
all.topgo$baseMean <- ceiling(all.topgo$baseMean)


Mor23.topgo.d <- MOR23.df[MOR23.df$state_flag == "Downreg",]
Mor23.topgo.up <- MOR23.df[MOR23.df$state_flag == "Upreg",]

M71.topgo.d  <- M71.df[M71.df$state_flag == "Downreg",]
M71.topgo.up  <- M71.df[M71.df$state_flag == "Upreg",]

common.d <- intersect(M71.topgo.d$GeneName,Mor23.topgo.d$GeneName)
common.up <- intersect(M71.topgo.up$GeneName,Mor23.topgo.up$GeneName)

comm.mod <- union(common.d,common.up)


#Input functions
Interested_Genes <-  comm.mod
Allgens <- all.topgo
Allgens$GeneName <- rownames(Allgens)


#### GOseq gene ontologies analysis ####
#Get Gene official symbol into GeneID

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- rownames(Allgens)

G_list <- getBM(filters= "external_gene_name", attributes= c("ensembl_gene_id","external_gene_name"),values=genes,mart= mart)

resLvl_ID <- merge(Allgens,G_list,by.x="GeneName",by.y="external_gene_name")

rownames(resLvl_ID) <-  resLvl_ID$ensembl_gene_id

resLvl_ID$ensembl_gene_id <- NULL

#We extract the significant upregulated genes (padj<0.1) and their annotation
a <- Interested_Genes
sigGenes <- subset(resLvl_ID, GeneName %in% a) 
sigGenes <- rownames(sigGenes)

anno <- AnnotationDbi::select(org.Mm.eg.db, 
                              keys=rownames(resLvl_ID), 
                              columns=c("SYMBOL","SYMBOL", "GENENAME"),
                              keytype="ENSEMBL")

anSig <- as.data.frame(subset(anno, ENSEMBL %in% sigGenes))

#We first get average gene expressions for each of the genes and then find non DE genes 
#that show a similar expression as the DE genes. These genes are then our background

#Get basemean from our whole dataset
overallBaseMean <- as.matrix(resLvl_ID[, "baseMean", drop = F])

#Match our significantly upregulated genes in these Basemeans
sig_idx <- match(anSig$ENSEMBL, rownames(overallBaseMean))

#Create the background dataset
backG <- c()

for(i in sig_idx){
  ind <- genefinder(overallBaseMean, i, 50, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
  
}

backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]

#We now remove DE genes from background and the get the total number of genes in the background
backG <- setdiff(backG,  anSig$ENSEMBL)
length(backG)
multidensity( list( 
  all= log2(resLvl_ID[,"baseMean"]) ,
  foreground =log2(resLvl_ID[anSig$ENSEMBL, "baseMean"]), 
  background =log2(resLvl_ID[backG, "baseMean"])), 
  xlab="log2 mean normalized counts", main = "Matching for enrichment analysis")


#Running topGO
#We first create a factor alg which indicates for every gene in our universe (union of background and DE-genes),
#whether it is differentially expressed or not. It has the ENSEMBL IDs of the genes in our universe as names 
#and contains 1 if the gene is DE and 0 otherwise.

#Create the 3 GO categories
onts = c( "MF", "BP", "CC" )

geneIDs = rownames(overallBaseMean)

inUniverse = geneIDs %in% c(anSig$ENSEMBL,  backG) 


inSelection =  geneIDs %in% anSig$ENSEMBL 

alg <- factor( as.integer( inSelection[inUniverse] ) )

names(alg) <- geneIDs[inUniverse]

tab = as.list(onts)

names(tab) = onts

for(i in 1:3){
  
  ## prepare data
  tgd <- new( "topGOdata", ontology= "MF", allGenes = alg, nodeSize=5,
              annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl" )
  

## run tests
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  
  resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
  
  #resultweight <- runTest(tgd, algorithm = "", statistic = "Fischer")
  
  ## look at results
  tab[[i]] <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                        Fisher.classic = resultTopGO.classic,
                        orderBy = "Fisher.classic" , topNodes = 1048)
  
}


for(i in 1:3){
  
  ## prepare data
  tgd <- new( "topGOdata", ontology= "CC", allGenes = alg, nodeSize=5,
              annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl" )
  
  
  ## run tests
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  
  resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
  
  #resultweight <- runTest(tgd, algorithm = "", statistic = "Fischer")
  
  ## look at results
  tab[[i]] <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                        Fisher.classic = resultTopGO.classic,
                        orderBy = "Fisher.classic" , topNodes = 841)
  
}


for(i in 1:3){
  
  ## prepare data
  tgd <- new( "topGOdata", ontology= "BP", allGenes = alg, nodeSize=5,
              annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl" )
  
  
  ## run tests
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  
  resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
  
  #resultweight <- runTest(tgd, algorithm = "", statistic = "Fischer")
  
  ## look at results
  tab[[i]] <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                        Fisher.classic = resultTopGO.classic,
                        orderBy = "Fisher.classic" , topNodes = 5730)
  
}







MF_topgo <- tab$MF
BP_topgo <- tab$BP
CC_topgo <- tab$CC
MF_topgo$ontology <- rep("MF", each = 1048)
BP_topgo$ontology <- rep("BP", each = 5730)
CC_topgo$ontology <- rep("CC", each = 841)


topgo.data <- rbind(MF_topgo,BP_topgo,CC_topgo)

topgo.data$p.adj <- p.adjust(topgo.data$Fisher.classic, method = "BH")

topgo.data

topgo.terms.data <-topgo.data[,c(1:2,9:10)]

colnames(topgo.terms.data) <- c("GOTerm", "Description", "GOTermType", "p.adj.val")



topgo.terms.data$GOTerm <- as.character(x = topgo.terms.data$GOTerm)
topgo.terms.data$GOTermType <- factor(x = topgo.terms.data$GOTermType,
                                      levels = c("BP", "MF", "CC"))
topgo.terms.data <- topgo.terms.data[order(topgo.terms.data$GOTermType, topgo.terms.data$GOTerm),]
topgo.terms.data$GOTerm <- factor(x = topgo.terms.data$GOTerm,
                                  levels = topgo.terms.data$GOTerm)



ggplot(data = topgo.terms.data,
       mapping = aes(x = GOTerm,
                     y = -log10(as.numeric(p.adj.val)),
                     colour = GOTermType)) +
  geom_point() +
  geom_text_repel(data = topgo.terms.data[topgo.terms.data$p.adj.val <= 0.05 ,],
                  mapping = aes(x = GOTerm,
                                y = -log10(as.numeric(p.adj.val)),
                                label = Description),
                  colour = "black") +
  geom_hline(yintercept = -log10(0.1),
             lty = 2)+
  ylab("-log10(p.value)") + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())



