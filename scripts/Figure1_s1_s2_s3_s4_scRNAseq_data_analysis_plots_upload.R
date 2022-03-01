## MOE and mature OSN plots ##

# This script generates plots for Figure 1 and supplementary

# Load packages
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(patchwork)
library(reshape2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(colorRamps)
library(fsbrain)
library(igraph)
source("common_layout.r")
source("center.palette.R") # from https://github.com/AmyOlex/RNASeqBits/blob/master/R/center.palette.R

# Python package
library(reticulate)
py_install("kneed")
py_install("matplotlib")

kn <- import("kneed")


##################################
#                                #        
#     MOE pre-filter plots       #
#                                #
##################################

moe.sct.integrated <- readRDS(file = "moe_noexp_scRNAseq10x_cellranger_nofilter.rds")

# Compute clusters
moe.sct.integrated <- FindNeighbors(object = moe.sct.integrated, 
                                    dims = 1:9)

moe.sct.integrated <- FindClusters(object = moe.sct.integrated,
                                   resolution = 1.1)

moe.sct.integrated <- RunUMAP(object = moe.sct.integrated, dims = 1:9 )

clusters.moe.no.filter <- c("0" ="mature OSN (7 clusters)",
                            "1" ="Mv2 (1 cluster)",
                            "2" ="mature OSN (7 clusters)" ,
                            "3" ="mature OSN (7 clusters)",
                            "4" ="mature OSN (7 clusters)",
                            "5" ="mature OSN (7 clusters)",
                            "6" ="mature OSN (7 clusters)",
                            "7" ="mature OSN (7 clusters)",
                            "8" ="developing OSN (4 clusters)",
                            "9" ="developing OSN (4 clusters)",
                            "10"="HBC1 + HBC2 (1 cluster)",
                            "11"="Immune cells (7 clusters)",
                            "12"="high % mito. genes (4 clusters)",
                            "13"="Immune cells (7 clusters)",
                            "14"="high % mito. genes (4 clusters)",
                            "15"="Sus (1 cluster)",
                            "16"="high % mito. genes (4 clusters)",
                            "17"="high % mito. genes (4 clusters)",
                            "18"="Immune cells (7 clusters)",
                            "19"="Immune cells (7 clusters)",
                            "20"="developing OSN (4 clusters)",
                            "21"="GBC",
                            "22"="Mv1 (1 cluster)",
                            "23"="Blood cells (2 clusters)",
                            "24"="Immune cells (7 clusters)",
                            "25"="Immune cells (7 clusters)",
                            "26"="Blood cells (2 clusters)",
                            "27"="Immune cells (7 clusters)")

moe.sct.integrated <- RenameIdents(object = moe.sct.integrated,clusters.moe.no.filter)

Clusters.to.keep.1 <- unique(clusters.moe.no.filter)
Clusters.to.keep.1 <- Clusters.to.keep.1[-6]

moe.sct.integrated.plot <- subset(x = moe.sct.integrated, 
                                  idents = Clusters.to.keep.1)

# Plot
no.filt.moe.dimplot <- DimPlot(object = moe.sct.integrated.plot,
                               reduction = "umap",
                               label = TRUE, cols = clusters.fill.moe.no.filter)  +
  common_layout + 
  theme(legend.position = "none",
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title = element_blank(),) + 
  xlim(-15, 15) +
  ylim(-15, 15)+
  coord_fixed(ratio=1)

S1.feature <- FeaturePlot(moe.sct.integrated.plot,features = c("Ptprc","Gypa"),
                          cols = c("lightskyblue1","darkorange1"),
                          min.cutoff = 0, 
                          ncol = 1,
                          order = T)

for(i in 1:2) {
  S1.feature[[i]] <- S1.feature[[i]]+xlim(-15, 15) +
    ylim(-15, 15)+
    coord_fixed(ratio=1)
}

S1 <- wrap_plots(A=no.filt.moe.dimplot,B=S1.feature,design = "AAB
                                                              AAB")
ggsave(filename = "Sup_1_A_moe.pdf",
       plot = S1,
       device = "pdf",
       width = 12,height = 6)


##################################
#                                #        
#           MOE plots            #
#                                #
##################################

# Load the MOE data
moe.filter.sct.integrated <- readRDS(file = "moe_noexp_scRNAseq10x_cellranger_filter_sct_integrated_clusters.rds")

x <- moe.filter.sct.integrated@reductions$umap@cell.embeddings


moe.filter.sct.integrated@reductions$umap@cell.embeddings[,2] <- ifelse(x[,2] < 0,
                                                                        yes = abs(x[,2]),
                                                                        no = -x[,2])
# Dimplot MOE
moe.dimplot <- DimPlot(object = moe.filter.sct.integrated,
                       reduction = "umap",
                       label = TRUE,
                       cols = clusters.fill.moe)  +
  common_layout + 
  theme(legend.position = "none",
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title = element_blank(),) + 
  xlim(-15, 15) +
  ylim(-15, 15)+
  coord_fixed(ratio=1)

ggsave(filename = "Fig_1_B_moe.pdf",
       plot = moe.dimplot,
       device = "pdf",
       width = 5,height = 5)

# Violin plot MOE
DefaultAssay(object = moe.filter.sct.integrated) <- "RNA"
moe.filter.sct.integrated <- NormalizeData(object = moe.filter.sct.integrated)

# Genes to plot
genes <- c("Ascl1",
           "Kit",
           "Neurod1",
           "Gng8",
           "Gap43",
           "Omp",
           "Adcy3",
           "Trp63",
           "Krt5",
           "Trpm5",
           "Ascl3",
           "Cftr",
           "Cyp2g1",
           "Cyp1a2")

# Extract the normalized count martix add the cluster identity of the cells
matrix.moe <- as.data.frame(x = t(x = moe.filter.sct.integrated@assays$RNA@data[genes,]))
matrix.moe$clusters <- factor(x = moe.filter.sct.integrated@active.ident,levels = c("Sus",
                                                                                    "Mv2",
                                                                                    "Mv1",
                                                                                    "HBC2",
                                                                                    "HBC1",
                                                                                    "mOSN",
                                                                                    "iOSN2",
                                                                                    "iOSN1",
                                                                                    "INP",
                                                                                    "GBC"))
# Violin plot MOE
moe.vlnplot <- lapply(X = colnames(matrix.moe)[-15],
                      FUN = function(gene){
                        a <- ggplot(data = matrix.moe,
                                    mapping = aes(x = matrix.moe[,gene],
                                                  y = clusters))+ 
                          geom_violin(scale = "width",
                                      aes(color = clusters,
                                          fill = clusters))+
                          common_layout + 
                          theme(legend.position = "none",
                                axis.title = element_blank())+
                          ggtitle(gene) +
                          scale_fill_manual(values = clusters.fill.moe)+
                          scale_color_manual(values = clusters.fill.moe)
                        
                        
                        if (gene != "Ascl1"){
                          a <- a + theme(axis.text.y = element_blank())
                          return(a)
                        }
                        else {return(a)}
                      }
)

moe.vlnplot <- wrap_plots(moe.vlnplot,nrow = 1)

sup.fig.1.moe <- wrap_plots(moe.dimplot,moe.vlnplot)

ggsave(filename = "Sup_2_A_moe.pdf",
       plot = sup.fig.1.moe,device = "pdf",
       width = 13,height = 4)

# Featureplot MOE
moe.featureplot <- FeaturePlot(moe.filter.sct.integrated,
                               features = c("Omp",
                                            "Adcy3"),
                               order = T,
                               ncol = 1,
                               cols = c("lightskyblue1","darkorange1"),
                               pt.size = 0.001) 

for(i in 1:2) {
  moe.featureplot[[i]] <- moe.featureplot[[i]] + common_layout+
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title = element_blank())+
    xlim(-15, 15) +
    ylim(-15, 15)+
    coord_fixed(ratio=1)}

ggsave(filename = "Fig_1_B_featureplot.pdf",
       plot = moe.featureplot,
       device = "pdf",
       width = 5,height = 10)


##################################
#                                #        
#           OSN plots            #
#                                #
##################################

# Read the OSN data 
osn.no.or.sub.integrated.sct <- readRDS("moe_noexp_scRNAseq10x_cellranger_filter_OSN_NO_OLFR_sct_integrated_clusters.rds")

# Dimplot OSN
osn.dimplot <- DimPlot(object = osn.no.or.sub.integrated.sct,
                       reduction = "umap",
                       label = TRUE, cols = clusters.fill.osn)  +
  common_layout + 
  theme(legend.position = "none",
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title = element_blank()) +
  coord_fixed(ratio = 1)

ggsave(filename = "Fig_1_C_osn.pdf",
       plot = osn.dimplot,
       device = "pdf",
       width = 5,height = 5)

# VlnPlot OSN
DefaultAssay(osn.no.or.sub.integrated.sct) <-"RNA"
osn.no.or.sub.integrated.sct <- NormalizeData(osn.no.or.sub.integrated.sct)

# Select the genes to plot
genes <- c("Nfix",
           "Ncam2",
           "Nqo1",
           "Acsm4",
           "Calb2",
           "Dlg2",
           "Cd36",
           "Tshz1",
           "Cd55")

# Extract the normalized count martix and add the cluster identity 
matrix.olfr <- as.data.frame(x = t(x = osn.no.or.sub.integrated.sct@assays$RNA@data[genes,]))
matrix.olfr$clusters <- osn.no.or.sub.integrated.sct$clusters

# Determine the order of the clusters
matrix.olfr$clusters <- factor(osn.no.or.sub.integrated.sct$clusters,levels = rev(c("Ventral",
                                                                                    "Calb2+ V.",
                                                                                    "Dlg2+ V.",
                                                                                    "Cd36+ V.",
                                                                                    "Cd55+ V.",
                                                                                    "Dorsal",
                                                                                    "Calb2+ D.",
                                                                                    "Dlg2+ D.",
                                                                                    "Cd36+ D.",
                                                                                    "Cd55+ D.")))

# Violin plot
osn.vlnplot <- lapply(X = colnames(matrix.olfr)[-10],
                      FUN = function(gene){
                        gene
                        a <- ggplot(matrix.olfr,aes(x = matrix.olfr[,gene],
                                                    y = clusters))+ 
                          geom_violin(scale = "width",aes(color=clusters,fill=clusters))+
                          common_layout + 
                          theme(legend.position = "none",
                                axis.title = element_blank())+
                          ggtitle(gene) +
                          scale_fill_manual(values = clusters.fill.osn)+
                          scale_color_manual(values = clusters.fill.osn)
                        
                        
                        if (gene != "Nfix"){
                          a <- a + theme(axis.text.y = element_blank())
                          return(a)
                        }
                        else {return(a)}
                      }
)

osn.vlnplot <- wrap_plots(osn.vlnplot,nrow = 1)

supl.fig.1.osn <- wrap_plots(osn.dimplot,osn.vlnplot)

ggsave(filename = "Sup_1_B_osn.pdf",
       plot = supl.fig.1.osn,device = "pdf",
       width = 10,height = 4)

# Feature plot
osn.featureplot <- FeaturePlot(osn.no.or.sub.integrated.sct,
                               features = c("Nfix",
                                            "Nqo1"),
                               order = T,
                               ncol = 1,
                               cols = c("lightskyblue1","darkorange1"),
                               pt.size = 0.001)+
  coord_fixed(ratio=1)

for(i in 1:2) {
  osn.featureplot[[i]] <- osn.featureplot[[i]] + common_layout+
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title = element_blank())+
    xlim(-15, 15) +
    ylim(-15, 15)+
    coord_fixed(ratio=1)}


ggsave(filename = "Fig_1_C_featureplot.pdf",
       plot = osn.featureplot,
       device = "pdf",
       width = 5,height = 10)


# OSN population UMAP
olfr.to.plot.list.f.plot <-  c("Olfr677",
                               "Olfr549",
                               "Olfr552",
                               "Olfr599",
                               "Olfr592",
                               "Olfr354",
                               "Olfr553",
                               "Olfr750",
                               "Olfr56",
                               "Olfr728",
                               "Olfr1183",
                               "Olfr1195")

# Create a TRUE/FALSE variable in meta.data for each Olfr
olfr.to.plot <- lapply(olfr.to.plot.list.f.plot, function(olfr){
  ident_olfr <- as.data.frame(osn.no.or.sub.integrated.sct@meta.data$olfr.ident == olfr)
  colnames(ident_olfr) <- olfr
  ident_olfr
})
olfr.to.plot <- do.call(olfr.to.plot,what = cbind)

umap_olfr <- osn.no.or.sub.integrated.sct@reductions$umap@cell.embeddings

umap_olfr <- cbind(umap_olfr,
                   olfr.to.plot,
                   osn.no.or.sub.integrated.sct$clusters)

cluster_for_olfr <- lapply(olfr.to.plot.list.f.plot, FUN = function(olfr){
  x <- paste(olfr,"_cluster")
  a <- as.data.frame(ifelse(umap_olfr[,olfr] ==TRUE, 
                            yes = as.character(umap_olfr$`osn.no.or.sub.integrated.sct$clusters`),
                            no = "no"))
  colnames(a) <-  paste(olfr,"cluster",sep = "_")
  a
})

cluster_for_olfr <- do.call(cluster_for_olfr,what = cbind)

umap_olfr <- cbind(umap_olfr,
                   cluster_for_olfr)

# OSN population UMAP plot
olfr.clusters <- lapply(colnames(cluster_for_olfr), FUN = function(olfr){
  n.cells <- sum(umap_olfr[,olfr]!= "no")
  
  ggplot(data = umap_olfr, mapping = aes(x = UMAP_1, y = UMAP_2))+
    
    geom_point(data = umap_olfr[umap_olfr[,olfr] == "no",],
               color = "gray80",size=0.2)+
    
    geom_point(data = umap_olfr[umap_olfr[,olfr] != "no",],
               aes(color = umap_olfr[umap_olfr[,olfr] != "no",][,olfr]),
               size=0.2) + 
    scale_color_manual(values = clusters.fill.osn) +
    annotate("text", x = -5, y = -5, label = n.cells)+
    common_layout +
    theme(legend.position = "none",
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title = element_blank())+
    coord_fixed(ratio = 1) + 
    ggtitle(label = strsplit(x = olfr,split = "_")[[1]][1])
  
})

olfr.clusters <- wrap_plots(olfr.clusters,ncol = 3)

ggsave(filename = "Fig_1_D_osn.pdf",
       plot = olfr.clusters,
       device = "pdf",
       width = 4.5,
       height = 5)


# OSN population violin plot 
olfr.vln <- c(
  # Cd55+ D.
  "Olfr592", 
  "Olfr608",
  
  
  # Cd55+ V.
  "Olfr354",
  "Olfr1440",
  "Olfr1431",
  
  
  # Dlg2+ D.
  "Olfr599", 
  "Olfr419",
  "Olfr65",
  
  
  # Dlg2+ V.
  "Olfr56",
  "Olfr1018",
  "Olfr129", 
  
  
  # Cd36+ D.
  "Olfr553",
  "Olfr810", 
  "Olfr665",
  
  # Cd36+ V.
  "Olfr750",
  "Olfr776", 
  "Olfr806",
  
  # Calb2+ D.
  "Olfr677", 
  "Olfr1039", 
  "Olfr411",
  
  # Calb2+ V.
  "Olfr728",
  "Olfr1490",
  "Olfr1507",
  
  # Dorsal
  "Olfr549",
  "Olfr552",
  "Olfr1361",
  
  # Ventral
  "Olfr1183",
  "Olfr1220",
  "Olfr1228",
  "Olfr1195"
)


# Create a TRUE/FALSE variable in meta.data for each Olfr
olfr.to.plot <- lapply(olfr.vln, function(olfr){
  ident_olfr <- as.data.frame(osn.no.or.sub.integrated.sct@meta.data$olfr.ident == olfr)
  colnames(ident_olfr) <- olfr
  ident_olfr
})

olfr.to.plot <- do.call(olfr.to.plot,what = cbind)

osn.no.or.sub.integrated.sct@meta.data <- cbind(osn.no.or.sub.integrated.sct@meta.data,
                                                olfr.to.plot)

osn.no.or.sub.integrated.sct$olfrs_to_plot <- ifelse(osn.no.or.sub.integrated.sct$olfr.ident %in% olfr.vln,
                                                     yes =osn.no.or.sub.integrated.sct$olfr.ident,
                                                     no ="no")

# Genes to plot
genes.1 <- c("S100a5","Stoml3","Tppp3","Pcdh8","Ptn","Cidea","Adam28","Pcdh10")

# Genes for supplementary Fig 1 C
genes.2 <- c("S100a5","Dlg2","Calb2","Cd55","Cd36","Nfix","Nqo1")

olfr_cells <- colnames(osn.no.or.sub.integrated.sct)[which(osn.no.or.sub.integrated.sct$olfr.ident %in% olfr.vln)]

olfr.data <- as.data.frame(osn.no.or.sub.integrated.sct@assays$RNA@data)
olfr.data <- as.data.frame(t(olfr.data[genes.2,]))
olfr.data <- olfr.data[olfr_cells,]
olfr.data$olfr.ident <- osn.no.or.sub.integrated.sct$olfr.ident[which(rownames(osn.no.or.sub.integrated.sct@meta.data) %in%rownames(olfr.data))]
olfr.data$cluster <- osn.no.or.sub.integrated.sct$clusters[which(rownames(osn.no.or.sub.integrated.sct@meta.data) %in%rownames(olfr.data))]
genes.position <- as.vector(1:(length(colnames(olfr.data))-2))

# Determine the olfr cluster origin
olfr.data.split <- split(x = olfr.data,
                         f = olfr.data$olfr.ident)
cluster.id.olfr <- lapply(X = olfr.data.split, 
                          FUN = function(x){
                            b <- choose.cluster(clusters = x$cluster)
                            b
                          }
)

olfr.data$cluster<- unname(apply(X = olfr.data,
                                 MARGIN = 1,
                                 FUN =  function(x){unlist(cluster.id.olfr[x["olfr.ident"]])}))

# Plot 
olfr.vlnplot <- lapply(X = genes.position,
                       FUN = function(gene.pos){
                         gene.pos
                         a <- ggplot(olfr.data,aes(y = reorder(olfr.ident, 
                                                               S100a5,
                                                               FUN = mean),
                                                   x=olfr.data[,gene.pos]))+ 
                           geom_violin(scale = "width",aes(color=cluster,fill=cluster))+
                           common_layout + 
                           theme(legend.position = "none",
                                 axis.title = element_blank())+
                           ggtitle(colnames(olfr.data)[gene.pos]) +
                           scale_fill_manual(values = clusters.fill.osn)+
                           scale_color_manual(values = clusters.fill.osn)
                         
                         
                         if (gene.pos>1){
                           a <- a + theme(axis.text.y = element_blank())
                           return(a)
                         }
                         else {return(a)}
                       }
)

olfr.vlnplot <- wrap_plots(olfr.vlnplot,nrow = 1)

ggsave(filename = "Sup_2_C.pdf",
       plot = olfr.vlnplot,device = "pdf",
       width = 6.5,height = 5)


##################################
#                                #        
#           Heatmap              #
#                                #
##################################

# Select the OSN populations for the heatmap
osn.split.by.cluster.metadata <- split(osn.no.or.sub.integrated.sct@meta.data, osn.no.or.sub.integrated.sct@meta.data$clusters)

top.1.olfr.per.cluster <- lapply(osn.split.by.cluster.metadata, FUN = function(x){
  a <- names(head(sort(table(x$olfr.ident),decreasing = T),1))
  
})
top.1.olfr.per.cluster <- unname(unlist(top.1.olfr.per.cluster))

big.olfr.ventral <- c("Olfr1348","Olfr239","Olfr1297","Olfr1211","Olfr731","Olfr536")

# Heatmap plots
heat_map_plot <- function(olfrs,cluster){
  
  
  # Create a variable in meta.data for each Olfr
  olfr.to.plot <- lapply(olfrs, function(olfr){
    ident_olfr <- as.data.frame(osn.no.or.sub.integrated.sct@meta.data$olfr.ident == olfr)
    colnames(ident_olfr) <- olfr
    ident_olfr
  })
  olfr.to.plot <- do.call(olfr.to.plot,what = cbind)
  osn.no.or.sub.integrated.sct@meta.data <- cbind(osn.no.or.sub.integrated.sct@meta.data,
                                                  olfr.to.plot)
  
  # Create a varaible in meta.data with the olfr identity of th cells but only for the chosen olfr
  osn.no.or.sub.integrated.sct$olfrs_to_plot <- ifelse(osn.no.or.sub.integrated.sct$olfr.ident %in% olfrs,
                                                       yes =osn.no.or.sub.integrated.sct$olfr.ident,
                                                       no ="no")
  
  # Create a variable
  olfr.cells.sub <-  osn.no.or.sub.integrated.sct
  
  if (cluster == "Ventral") {
    olfr.cells.sub <- subset(osn.no.or.sub.integrated.sct,idents =  "Ventral") 
  }
  
  olfr_cells <- colnames(olfr.cells.sub)[which(olfr.cells.sub$olfr.ident %in% olfrs)]
  
  olfr.cells.sub <- SetIdent(olfr.cells.sub,value = olfr.cells.sub$olfrs_to_plot)
  
  DefaultAssay(olfr.cells.sub) <-"RNA"
  
  olfr.cells.sub[['integrated']] <- NULL
  olfr.cells.sub[['SCT']] <- NULL
  
  olfr.cells.sub <- NormalizeData(olfr.cells.sub)
  olfr.cells.sub <- ScaleData(olfr.cells.sub)
  
  olfr_markers <- FindAllMarkers(object = olfr.cells.sub,
                                 min.pct = 0.7,
                                 only.pos = TRUE)
  
  olfr_markers_small <- split(x = olfr_markers,
                              f = olfr_markers$cluster)
  
  
  olfr_markers_small <- mapply(markers = olfr_markers_small,
                               olfr.ident = names(x = olfr_markers_small),
                               FUN = function(markers,
                                              olfr.ident) {
                                 markers.sub <- markers[markers$p_val_adj < 0.05,]$gene
                                 names(x = markers.sub) <- rep(x = olfr.ident,
                                                               times = length(x = markers.sub))
                                 return(markers.sub)}, 
                               SIMPLIFY = FALSE,
                               USE.NAMES = TRUE)
  
  olfr_markers_small <- unlist(x = unname(obj = olfr_markers_small))
  olfr_markers_small <- olfr_markers_small[order(match(x = names(x = olfr_markers_small), 
                                                       table = olfrs))]
  olfr_markers_small <- olfr_markers_small[-which(names(olfr_markers_small) =="no")]
  
  
  markers.data <- data.frame(markers = olfr_markers_small,
                             olfr.ident = names(x = olfr_markers_small))
  duplicated.genes <- unique(x = markers.data$markers[duplicated(x = markers.data$markers)])
  markers.data$duplicated <- markers.data$markers %in% duplicated.genes
  table(markers.data$olfr.ident, markers.data$duplicated)
  
  # olfr_markers_small <- olfr_markers_small[!duplicated(x = olfr_markers_small)] # Trick to keep names of vectors (not done with the unique() function)
  genes <- markers.data$markers
  
  # Heatmap
  matrix.olfr <- as.data.frame(x = t(x = olfr.cells.sub@assays$RNA@scale.data))
  
  matrix.olfr$olfr.ident <- olfr.cells.sub$olfr.ident
  matrix.olfr$cluster <- olfr.cells.sub$clusters
  matrix.olfr <- matrix.olfr[olfr_cells,]
  
  col_names <- matrix.olfr[c("olfr.ident","cluster")]
  col_names <- col_names[order(match(x = col_names$olfr.ident, table = olfrs)),]
  gaps <- table(col_names$olfr.ident)
  gaps <- gaps[olfrs]
  
  for (i in 1:length(gaps)){
    if (i > 1){ gaps[i] <- gaps[i-1] + gaps[i]}}
  
  matrix.olfr <- matrix.olfr[order(match(x = matrix.olfr$olfr.ident, table = olfrs)), genes]
  
  matrix.olfr.clipped <- clip.data(data = matrix.olfr,
                                   lower = 0.025,
                                   upper = 0.975)
  
  palette.colors <- center.palette(data = matrix.olfr.clipped,
                                   palette_length = 100,
                                   color1 = "lightskyblue1",
                                   color2 = "darkorange1")
  
  # Plot the heatmap
  heat_map <- pheatmap(mat = matrix.olfr.clipped,
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       annotation_row = col_names[-2],
                       show_row = FALSE,
                       show_rownames = FALSE,
                       show_colnames = T,
                       treeheight_col = 0,
                       color = palette.colors$colors,
                       breaks = palette.colors$breaks,
                       gaps_row = gaps)
}

heatmap.inter <- heat_map_plot(top.1.olfr.per.cluster,"all")

heatmap.intra <- heat_map_plot(big.olfr.ventral,"Ventral")

ggsave(filename = "Fig_1_G_heatmap_inter.pdf",
       plot = heatmap.inter,device = "pdf",
       width = 8,
       height = 4)

ggsave(filename = "Fig_1_H_heatmap_intra.pdf",
       plot = heatmap.intra,device = "pdf",
       width = 8,
       height = 4)

# Comparison between with and without OR genes datasets

# NMI score
osn.sub.integrated.sct <- readRDS("moe_noexp_scRNAseq10x_cellranger_filter_OSN_WITH_OLFR_sct_integrated_clusters.rds")
osn.no.or.sub.integrated.sct <- readRDS("moe_noexp_scRNAseq10x_cellranger_filter_OSN_NO_OLFR_sct_integrated_clusters.rds")

clustering.with.OR <- data.frame(colnames(osn.sub.integrated.sct),osn.sub.integrated.sct$clusters)

clustering.no.OR <- data.frame(colnames(osn.no.or.sub.integrated.sct),osn.no.or.sub.integrated.sct$clusters)

compare(clustering.with.OR$osn.sub.integrated.sct.clusters,clustering.no.OR$osn.no.or.sub.integrated.sct.clusters,method = "nmi")

# Dimplot 
osn.dimplot.or <- DimPlot(object = osn.sub.integrated.sct,
                          label = TRUE, cols = clusters.fill.osn)  +
  common_layout + 
  theme(legend.position = "none",
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title = element_blank()) +
  coord_fixed(ratio = 1)

# Elbowplot without OR
elbow.no.or <- osn.no.or.sub.integrated.sct@reductions$pca@stdev
elbow.no.or <- data.frame(PC=c(1:50),standard_deviation=elbow.no.or)

elbow.no.or.plot <- ggplot(elbow.no.or,aes(x=PC,y=standard_deviation))+
  geom_point() + 
  geom_vline(xintercept = 15,color="red")+
  scale_x_continuous(n.breaks = 15,minor_breaks = c(2))+
  ylim(1,7)+
  ylab("standard deviation")+
  common_layout

# Find PC elbow without OR
knee.noOR <- kn$KneeLocator(x = seq_len(length.out = length(x = osn.no.or.sub.integrated.sct@reductions$pca@stdev)),
                            y = osn.no.or.sub.integrated.sct@reductions$pca@stdev,
                            S = 1,
                            curve = "convex",
                            direction = "decreasing")


knee.noOR.frame <- data.frame(y.norm = knee.noOR$y_normalized, 
                              PC = c(1:50),
                              y.diff = knee.noOR$y_difference)

# Elbow
knee.noOR.plot <- ggplot(knee.noOR.frame)+ 
  geom_line(aes(x=PC,
                y=y.norm)) + 
  geom_line(aes(x=PC,
                y=y.diff),color="red") + 
  scale_x_continuous(n.breaks = 10,expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  geom_vline(xintercept = knee.noOR$elbow,color = "red")+
  common_layout


elbow.noOR.plot <- wrap_plots(elbow.no.or.plot,knee.noOR.plot)

# Find PC elbow with OR
elbow.or <- osn.sub.integrated.sct@reductions$pca@stdev
elbow.or <- data.frame(PC=c(1:50),standard_deviation=elbow.or)

elbow.or <- ggplot(elbow.or,aes(x=PC,y=standard_deviation))+
  geom_point() + 
  geom_vline(xintercept = 14,color="red")+
  scale_x_continuous(n.breaks = 15)+
  ylim(1,7)+
  ylab("standard deviation")+
  common_layout


# Elbow
knee.OR <- kn$KneeLocator(x = seq_len(length.out = length(x = osn.sub.integrated.sct@reductions$pca@stdev)),
                          y = osn.sub.integrated.sct@reductions$pca@stdev,
                          S = 1,
                          curve = "convex",
                          direction = "decreasing")

knee.OR.frame <- data.frame(y.norm = knee.OR$y_normalized, 
                            PC = c(1:50),
                            y.diff = knee.OR$y_difference,
                            x.diff = knee.OR$x_difference)


knee.OR.plot <- ggplot(knee.OR.frame)+ 
  geom_line(aes(x=PC,
                y=y.norm)) + 
  geom_line(aes(x=PC,
                y=y.diff),color="red") + 
  scale_x_continuous(n.breaks = 10,expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  geom_vline(xintercept = knee.OR$elbow,color = "red")+
  common_layout


elbow.OR.plot <- wrap_plots(elbow.or,knee.OR.plot)

final_elbow <- wrap_plots(elbow.noOR.plot,elbow.OR.plot, ncol = 1)


ggsave(filename = "Elbows_withoutOR.pdf",
       plot = final_elbow,
       device = "pdf",
       width = 8,height =8)


# Create a variable in meta.data for each Olfr
olfr.to.plot <- lapply(olfr.to.plot.list.f.plot, function(olfr){
  ident_olfr <- as.data.frame(osn.sub.integrated.sct@meta.data$olfr.ident == olfr)
  colnames(ident_olfr) <- olfr
  ident_olfr
})
olfr.to.plot <- do.call(olfr.to.plot,what = cbind)

umap_olfr <- osn.sub.integrated.sct@reductions$umap@cell.embeddings

umap_olfr <- cbind(umap_olfr,
                   olfr.to.plot,
                   osn.sub.integrated.sct$clusters)

cluster_for_olfr <- lapply(olfr.to.plot.list.f.plot, FUN = function(olfr){
  x <- paste(olfr,"_cluster")
  a <- as.data.frame(ifelse(umap_olfr[,olfr] ==TRUE, 
                            yes = as.character(umap_olfr$`osn.sub.integrated.sct$clusters`),
                            no = "no"))
  colnames(a) <-  paste(olfr,"cluster",sep = "_")
  a
})

cluster_for_olfr <- do.call(cluster_for_olfr,what = cbind)

umap_olfr <- cbind(umap_olfr,
                   cluster_for_olfr)

# OSN clustering compariosn between data with and without OR genes
olfr.clusters.or <- lapply(colnames(cluster_for_olfr), FUN = function(olfr){
  ggplot(data = umap_olfr, mapping = aes(x = UMAP_1, y = UMAP_2))+
    
    geom_point(data = umap_olfr[umap_olfr[,olfr] == "no",],
               color = "gray80",size=0.2)+
    
    geom_point(data = umap_olfr[umap_olfr[,olfr] != "no",],
               aes(color = umap_olfr[umap_olfr[,olfr] != "no",][,olfr]),
               size=0.2) + 
    scale_color_manual(values = clusters.fill.osn) +
    common_layout +
    theme(legend.position = "none",
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title = element_blank())+
    coord_fixed(ratio = 1) + 
    ggtitle(label = strsplit(x = olfr,split = "_")[[1]][1])
})


# Create a TRUE/FALSE variable in meta.data for each Olfr
olfr.to.plot <- lapply(olfr.to.plot.list.f.plot, function(olfr){
  ident_olfr <- as.data.frame(osn.no.or.sub.integrated.sct@meta.data$olfr.ident == olfr)
  colnames(ident_olfr) <- olfr
  ident_olfr
})
olfr.to.plot <- do.call(olfr.to.plot,what = cbind)

umap_olfr <- osn.no.or.sub.integrated.sct@reductions$umap@cell.embeddings

umap_olfr <- cbind(umap_olfr,
                   olfr.to.plot,
                   osn.no.or.sub.integrated.sct$clusters)

cluster_for_olfr <- lapply(olfr.to.plot.list.f.plot, FUN = function(olfr){
  x <- paste(olfr,"_cluster")
  a <- as.data.frame(ifelse(umap_olfr[,olfr] ==TRUE, 
                            yes = as.character(umap_olfr$`osn.no.or.sub.integrated.sct$clusters`),
                            no = "no"))
  colnames(a) <-  paste(olfr,"cluster",sep = "_")
  a
})

cluster_for_olfr <- do.call(cluster_for_olfr,what = cbind)

umap_olfr <- cbind(umap_olfr,
                   cluster_for_olfr)

# OSN population UMAP plot
olfr.clusters <- lapply(colnames(cluster_for_olfr), FUN = function(olfr){
  n.cells <- sum(umap_olfr[,olfr]!= "no")
  
  ggplot(data = umap_olfr, mapping = aes(x = UMAP_1, y = UMAP_2))+
    
    geom_point(data = umap_olfr[umap_olfr[,olfr] == "no",],
               color = "gray80",size=0.2)+
    
    geom_point(data = umap_olfr[umap_olfr[,olfr] != "no",],
               aes(color = umap_olfr[umap_olfr[,olfr] != "no",][,olfr]),
               size=0.2) + 
    scale_color_manual(values = clusters.fill.osn) +
    common_layout +
    theme(legend.position = "none",
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title = element_blank())+
    coord_fixed(ratio = 1) + 
    ggtitle(label = strsplit(x = olfr,split = "_")[[1]][1])
  
})

# Combine  with and without OR umap plots
nbr.OR <- c(1:length(olfr.clusters))

umap.plots.or <- wrap_plots(lapply(nbr.OR, function(x){
  plot <- wrap_plots(olfr.clusters[[x]],olfr.clusters.or[[x]])
  plot
}
),ncol=2
)

sup.fig.2 <- wrap_plots(A = osn.dimplot,
                        B = osn.dimplot.or,
                        C = umap.plots.or,
                        design = "AB
                                  CC
                                  CC
                                  CC
                                  CC")

ggsave(filename = "Sup_4_umap.pdf",
       plot = sup.fig.2,
       device = "pdf",
       width = 10,height = 10)
