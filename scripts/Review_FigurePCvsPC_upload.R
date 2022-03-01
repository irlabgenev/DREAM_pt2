
# load packages
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggrastr)
source("N:/Novell/__articles in progress__/DREAM transcriptomic adaptations/figures/R scripts/common_layout.r")

# Set working directory
setwd(dir = "N:/Novell/___Projects___/DREAM/experiments/Scrnaseq_moe/20210526_moe_scRNAseq_cellranger_grc38irlab/Script_version_upload/")

# Load the OSN data without OR
osn.no.or.sub.integrated.sct <- readRDS("__review_RDS/moe_noexp_scRNAseq10x_cellranger_filter_OSN_WITH_OLFR_sct_integrated_clusters.rds")


data.pca <- as.data.frame(osn.no.or.sub.integrated.sct@reductions$pca@cell.embeddings)
data.pca$cluster <- osn.no.or.sub.integrated.sct@meta.data$clusters


PCvsPC.plot <- lapply(c(1:14), function(x){
  plots <- lapply(c(1:14), function(y){
    ggplot(data.pca,aes(x= data.pca[,x], y =data.pca[,y],color=cluster))+ 
      geom_point_rast(raster.dpi = 300)+
      scale_color_manual(values = clusters.fill.osn)+
      common_layout+
      theme(legend.position = "none", axis.title = element_blank())
  }
  )
  plots <- wrap_plots(plots,ncol = 1)
}
)


PCvsPC.plot.ras <- wrap_plots(PCvsPC.plot,ncol = 14)


ggsave(PCvsPC.plot.ras,device = "pdf", filename = "PCvsPC.plot.withOR.pdf",height = 30,width = 30)
