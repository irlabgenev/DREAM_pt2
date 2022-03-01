## Main olfactory epithelium data processing ##

# This script is used to filter and cluster identification/annotation of the MOE and mature OSN. 
# It generates Seurat objects used for plots and further analysis.

##################################
#                                #        
#   MOE analysis: no filtering   #
#                                #
##################################

# Load packages
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(english)

# Python package
library(reticulate)
py_install("kneed")
py_install("matplotlib")

#Custom file
source("common_layout.r")

# Read 10X count matrix
moe1 <- Read10X(data.dir = "s1_filtered_feature_bc_matrix/")
moe1 <- CreateSeuratObject(counts = moe1)

moe2 <- Read10X(data.dir = "s2_filtered_feature_bc_matrix/")
moe2 <- CreateSeuratObject(counts = moe2)

moe3 <- Read10X(data.dir = "s3_filtered_feature_bc_matrix/")
moe3 <- CreateSeuratObject(counts = moe3)

moe4 <- Read10X(data.dir = "s4_filtered_feature_bc_matrix/")
moe4 <- CreateSeuratObject(counts = moe4)

# Merge the samples
moe <- merge(x = moe1,
             y = c(moe2, moe3, moe4))

# Update names of orig.ident
moe@meta.data$orig.ident <- gsub(pattern = ".*_",
                                 replacement = "moe_",
                                 x = rownames(moe@meta.data))

# Measure percentage of counts of mitochondrial genes 
moe[["percent.mt"]] <- PercentageFeatureSet(object = moe,
                                            pattern = "^mt-" )

# Create a variable with nUMI/gene for futur normalization
moe@meta.data$nUMI.per.gene <- moe@meta.data$nCount_RNA/moe@meta.data$nFeature_RNA

# Integrate the data by following the Seurat integration SCTransform normalized data integration pipeline
# https://satijalab.org/seurat/articles/integration_introduction.html)

moe.filt <- moe

moe.filt.list <- SplitObject(moe.filt, 
                             split.by = "orig.ident")

moe.filt.list <- lapply(X = moe.filt.list, 
                        FUN = SCTransform)

features <- SelectIntegrationFeatures(object.list = moe.filt.list, 
                                      nfeatures = 5000)

moe.filt.list <- PrepSCTIntegration(object.list = moe.filt.list,
                                    anchor.features = features)

momoe.anchors <- FindIntegrationAnchors(object.list = moe.filt.list, 
                                      normalization.method = "SCT", 
                                      anchor.features = features)

moe.integrated.sct <- IntegrateData(anchorset = moe.anchors, 
                                    normalization.method = "SCT")

# Compute PCA
moe.integrated.sct <- RunPCA(moe.integrated.sct, 
                             features = VariableFeatures(object = moe.integrated.sct))

ElbowPlot(moe.integrated.sct,ndims = 50)

#%%% Checkpoint %%%#

saveRDS(moe.integrated.sct,"moe_noexp_scRNAseq10x_cellranger_nofilter.rds")


##################################
#                                #        
#   MOE analysis: filtering      #
#                                #
##################################

# read MOE_scRNAseq_no_filter_integ_sct.rds
moe.integrated.sct <- readRDS(file = "moe_noexp_scRNAseq10x_cellranger_nofilter.rds")
DefaultAssay(moe.integrated.sct) <- "integrated"

kn <- import("kneed")
plt <- import(module = "matplotlib.pyplot")
knaa <- kn$KneeLocator(x = seq_len(length.out = length(x = moe.integrated.sct@reductions$pca@stdev)),
                       y = moe.integrated.sct@reductions$pca@stdev,
                       S = 1,
                       curve = "convex",
                       direction = "decreasing")
knaa$all_elbows

# Compute clusters
moe.integrated.sct <- FindNeighbors(object = moe.integrated.sct, 
                                    dims = 1:9)

moe.integrated.sct <- FindClusters(object = moe.integrated.sct,
                                   resolution = 1.1)

moe.integrated.sct <- RunUMAP(object = moe.integrated.sct, dims = 1:9 ,label=T)

# Select clusters to remove:
Immune.cells <- c(11,13,18,19,24,25,27) # based on a the expression level of Igkc, Cd52, Cybb, Ctss,Tyrobp
Apoptotic.cells <- c(16,17,12,14) # based on a high percent.mt and low nFeatures
Blood.cells <- c(23,26) # based on a the expression level of Hbb-a1 and Hbb-a2

Remove.cells <- c(Immune.cells,
                  Apoptotic.cells,
                  Blood.cells)

# Create a vector containing the clusters to keep
Clusters.to.keep <- as.vector(as.numeric(unique(levels(moe.integrated.sct$integrated_snn_res.1.1))))
Clusters.to.keep <- c(Clusters.to.keep[!Clusters.to.keep %in% Remove.cells ])

# Create a new seurat object only containing the Clusters.to.keep
moe.filter <- subset(x = moe.integrated.sct, 
                     idents = Clusters.to.keep) 

# Filter moe.integrated.sct.filter by removing cells with:
moe.filter <- subset(x = moe.filter,
                     subset = nFeature_RNA > 1000 &  # less than 1000 genes expressed
                       percent.mt < 10) # more than 10% mitochondrial genes

# Integrate the data by following the Seurat integration SCTransform normalized data integration pipeline
# https://satijalab.org/seurat/articles/integration_introduction.html)

moe.filt.list <- SplitObject(moe.filter, 
                             split.by = "orig.ident")

moe.filt.list <- lapply(X = moe.filt.list, 
                        FUN = SCTransform)

features <- SelectIntegrationFeatures(object.list = moe.filt.list, 
                                      nfeatures = 5000)

moe.filt.list <- PrepSCTIntegration(object.list = moe.filt.list,
                                    anchor.features = features)

moe.anchors <- FindIntegrationAnchors(object.list = moe.filt.list, 
                                      normalization.method = "SCT", 
                                      anchor.features = features)

moe.filter.sct.integrated <- IntegrateData(anchorset = moe.anchors, 
                                           normalization.method = "SCT")

# Compute the PCA
moe.filter.sct.integrated <- RunPCA(object = moe.filter.sct.integrated)


# Clustering 
DefaultAssay(object = moe.filter.sct.integrated) <- "integrated"

moe.filter.sct.integrated <- FindNeighbors(object = moe.filter.sct.integrated, 
                                           dims = 1:13)

moe.filter.sct.integrated <- FindClusters(object = moe.filter.sct.integrated,
                                          resolution = 0.3)

moe.filter.sct.integrated <- RunUMAP(object = moe.filter.sct.integrated, 
                                     dims = 1:13)


#  Determine clusters indentity (1:15 PC, 0.3 resolution)
moe.filter.sct.integrated <- RenameIdents(object = moe.filter.sct.integrated,
                                          "0" = "mOSN", # Markers Omp, Cnga2, Gnga13
                                          "1" = "mOSN", # Markers Omp, Cnga2, Gnga13
                                          "2" = "mOSN", # Markers Omp, Cnga2, Gnga13
                                          "3" = "Mv2", # Markers Ascl3, Cftr
                                          "4" = "mOSN", ## Markers Omp, Cnga2, Gnga13
                                          "5" = "iOSN2", # Markers Gap43, Omp
                                          "6" = "iOSN1", #  Markers Gap43, Gng8
                                          "7" = "mOSN", # Markers Trp63, Krt5
                                          "8" = "HBC2", # Markers Omp, Cnga2, Gnga13
                                          "9" = "Sus", # Markers Cyp2g1, Cyp1a2
                                          "10" = "INP", # # Markers Neurod1, Lhx2
                                          "11" = "HBC1", # Markers Trp63, Krt5
                                          "12" = "GBC", # Markers Ascl1, Kit
                                          "13" = "HBC2", # Markers Omp, Cnga2, Gnga13
                                          "14" = "Mv1") # Markers Ascl1, Kit
my_levels <- c("Sus",
               "Mv1",
               "Mv2",
               "HBC2",
               "HBC1",
               "mOSN",
               "iOSN2",
               "iOSN1",
               "INP",
               "GBC")

moe.filter.sct.integrated@active.ident <- factor(x = moe.filter.sct.integrated@active.ident, levels = my_levels)


#%%% Checkpoint %%%#

saveRDS(object = moe.filter.sct.integrated,file = "moe_noexp_scRNAseq10x_cellranger_filter_sct_integrated_clusters.rds")



#####################
#                   #        
#  mOSN: filtering  #
#                   #
#####################

# Read the data
moe.data <- readRDS(file = "moe_noexp_scRNAseq10x_cellranger_filter_sct_integrated_clusters.rds")

# Normalize de data
DefaultAssay(moe.data) <- "RNA"
moe.data <- NormalizeData(object = moe.data)

# Remove integrated and SCT assays from the Seurat object
moe.data[['integrated']] <- NULL
moe.data[['SCT']] <- NULL

# Subset mOSN 
osn.data <- subset(moe.data,
                   idents = c("mOSN"))

# Keep genes that are expressed in the dataset
genes.to.keep <- rowSums(x = osn.data) > 0
osn.data <- osn.data[genes.to.keep,]

# Read table containing all Olfr genes identified by Joel
olfr.genes <- readRDS(file = "Mus_musculus_functional_ORs_and_class.rds")

# Remove from olfr.genes, genes not found in the seurat object
olfr.genes <- olfr.genes[olfr.genes$label %in% rownames(x = osn.data),]

# Create a variable containing all Olfr genes
other.olfr.genes <- grep(pattern = "Olfr",
                         x = rownames(x = osn.data),
                         value = TRUE)

# Create a variable containing other olfactory genes
other.olfactory.genes <- grep(pattern = "Gucy2d|Gucy1b2|Taar|Ms4",
                              x = rownames(x = osn.data),
                              value = TRUE)

# Make the union of the different olfactory genes
olfr.gene.names <- union(x = olfr.genes$label,
                         y = c(other.olfr.genes,
                               other.olfactory.genes))

write.table(olfr.gene.names,file = "olfr_gene_name.txt")

# Create a vector containing the normalized counts from the osn.data seurat object
matrix.olfr <- as.data.frame(x = t(x = osn.data@assays$RNA@data))

# Remove non-olfactory genes from the matrix
matrix.olfr <- matrix.olfr[,olfr.gene.names]

# Compute the number of Olfr per cell
osn.data$n.olfr.per.osn <- rowSums(x = matrix.olfr > 0)

# Determine the Olfr gene expressed in each OSN
olfr.max <- t(as.data.frame(apply(X = matrix.olfr, MARGIN = 1, FUN = function (x) {
  olfr.max.lvl <-  max(x)
  olfr.max <- if (olfr.max.lvl == 0) {
    "no olfr"
  } else {
    olfr.gene.names[which(x == max(x))][1]
  }
  return((c(olfr.max,
            as.numeric(x = olfr.max.lvl))))
})))

# Append olfr.max and olfr.second to the Seurat object meta data 
if (all(names(x = olfr.max) == rownames(x = osn.data@meta.data))) {
  osn.data@meta.data$olfr.ident <- olfr.max[,1]
  osn.data@meta.data$olfr.lvl <- as.numeric(olfr.max[,2])
}

# Remove cells that do not express an Olfr gene
osn.no.or <- which(x = osn.data$olfr.ident != "no olfr")
osn.data <- subset(x = osn.data,
                   cells = osn.no.or)


# Detect how many Olfrs are expressed per OSN
olfr.genes.data <- as.matrix(x = osn.data@assays$RNA@data)
olfr.genes.data <- olfr.genes.data[olfr.gene.names,]

n.olfr.per.osn <- colSums(x = olfr.genes.data > 0)

table(n.olfr.per.osn) / length(x = n.olfr.per.osn) * 100

olfr.genes.data.melt <- reshape2::melt(data = olfr.genes.data,
                                       varnames = c("gene",
                                                    "cell"),
                                       value.name = "norm.UMI")
olfr.genes.data.melt$gene <- as.character(x = olfr.genes.data.melt$gene)
olfr.genes.data.melt$cell <- as.character(x = olfr.genes.data.melt$cell)

olfr.genes.data.l <- split(x = olfr.genes.data.melt,
                           f = olfr.genes.data.melt$cell)

olfr.genes.data.l <- lapply(X = olfr.genes.data.l,
                            FUN = function(osn.indiv.data,
                                           n.olfrs.to.keep) {
                              osn.indiv.data <- osn.indiv.data[order(osn.indiv.data$norm.UMI, decreasing = TRUE),]
                              osn.indiv.data <- osn.indiv.data[1:n.olfrs.to.keep, ]
                              olfr.ranks <- as.character(x = as.english(x = 1:n.olfrs.to.keep))
                              osn.indiv.data$olfr.rank <- factor(x = olfr.ranks,
                                                                 levels = olfr.ranks)
                              return(osn.indiv.data)
                            },
                            n.olfrs.to.keep = max(n.olfr.per.osn))

olfr.genes.data <- dplyr::bind_rows(olfr.genes.data.l)
olfr.genes.data <- olfr.genes.data[olfr.genes.data$norm.UMI > 0,]
olfr.genes.data$olfr.num <- n.olfr.per.osn[olfr.genes.data$cell]

# Find the median absolute deviation of the first Olfr expression
first.norm.count <- olfr.genes.data[olfr.genes.data$olfr.rank == "one",]$norm.UMI
med.first <- median(x = first.norm.count)
mad.first <- mad(x = first.norm.count, center = med.first)


y <- olfr.genes.data[!duplicated(olfr.genes.data$cell),]
x <- table(y$olfr.num)
x[1]*100 /(sum(x))

# Filter cells with a lower expression than 3 time the median absolute deviation
filter.val <- med.first - 3*mad.first

cells.to.remove1 <- olfr.genes.data$cell[olfr.genes.data$olfr.rank != "one" & 
                                           olfr.genes.data$norm.UMI > filter.val]

cells.to.remove2 <- olfr.genes.data$cell[olfr.genes.data$olfr.rank == "one" & 
                                           olfr.genes.data$norm.UMI <= filter.val]

cells.to.remove <- c(cells.to.remove1, cells.to.remove2)

osn.data <- osn.data[,!colnames(osn.data) %in% cells.to.remove]


# PLOT

rank.olfr.plot <- ggplot(data = olfr.genes.data, aes(x = olfr.rank, y = norm.UMI)) + 
  geom_violin(fill="royalblue")+
  geom_abline(intercept = filter.val,slope = 0,color="red")+
  common_layout 

ggsave(plot = rank.olfr.plot,
       filename = "olfr.rank.pdf",
       device = "pdf",
       width = 7,height = 4)

y <- olfr.genes.data[!duplicated(olfr.genes.data$cell),]
y <- as.data.frame(table(y$olfr.num))

ggplot(data = y ,aes(x = y$Var1, y = y$Freq))+ geom_boxplot()

# Determine the class of the Olfr gene expressed in each OSN
osn.idents <- osn.data$olfr.ident
genes.data <- merge(x = olfr.genes,
                    y = data.frame(label = unique(x = osn.idents)),
                    by = "label",
                    all = TRUE)
rownames(x = genes.data) <- genes.data$label

genes.data$class <- ifelse(test = is.na(x = genes.data$class),
                           yes = "none",
                           no = genes.data$class)

olfr.classes <- genes.data$class
names(x = olfr.classes) <- genes.data$label
osn.data$class_first <- olfr.classes[osn.data$olfr.ident]

# Get names of Olfr genes expressed in at least 3 cells
or.pop.size <- as.data.frame(x = table(osn.data$olfr.ident))
olfrs.to.keep <- as.character(x = or.pop.size$Var1[or.pop.size$Freq >= 3])
cells.to.keep <- rownames(x = osn.data@meta.data)[osn.data@meta.data$olfr.ident %in% olfrs.to.keep]

# Subset object to include only cells.to.keep
osn.data.sub <- subset(x = osn.data,
                       cells = cells.to.keep)

#%%% Checkpoint %%%#

saveRDS(object = osn.data.sub,file = "moe_noexp_scRNAseq10x_cellranger_filter_OSN.rds")


###############################
#                             #        
#   OSN Without OR analysis   #
#                             #
###############################

# Read the data
osn.data.sub <- readRDS(file = "moe_noexp_scRNAseq10x_cellranger_filter_OSN.rds")


# Remove Olfr genes
olfr.gene.names <- read.table(file = "olfr_gene_name.txt",header = TRUE)
osn.no.or.data.sub <- osn.data.sub[!rownames(x = osn.data.sub) %in% olfr.gene.names$x]

# Clustering
DefaultAssay(object = osn.no.or.data.sub) <- "RNA"

osn.no.or.sub.filt.list <- SplitObject(object = osn.no.or.data.sub, 
                                       split.by = "orig.ident")

osn.no.or.sub.filt.list <- lapply(X = osn.no.or.sub.filt.list, 
                                  FUN = function(x){SCTransform(x,
                                                                vars.to.regress = "percent.mt")})

features.no.or <- SelectIntegrationFeatures(object.list = osn.no.or.sub.filt.list, 
                                            nfeatures = 5000)

osn.no.or.sub.filt.list <- PrepSCTIntegration(object.list = osn.no.or.sub.filt.list,
                                              anchor.features = features.no.or)

osn.no.or.sub.anchors <- FindIntegrationAnchors(object.list = osn.no.or.sub.filt.list, 
                                                normalization.method = "SCT", 
                                                anchor.features = features.no.or)

osn.no.or.sub.integrated.sct <- IntegrateData(anchorset = osn.no.or.sub.anchors, 
                                              normalization.method = "SCT")


## Find PCA
DefaultAssay(osn.no.or.sub.integrated.sct) <- "integrated"

osn.no.or.sub.integrated.sct <-  RunPCA(object = osn.no.or.sub.integrated.sct)


# Select PCA
kn <- import("kneed")
plt <- import(module = "matplotlib.pyplot")
knaa <- kn$KneeLocator(x = seq_len(length.out = length(x = osn.no.or.sub.integrated.sct@reductions$pca@stdev)),
                       y = osn.no.or.sub.integrated.sct@reductions$pca@stdev,
                       S = 1,
                       curve = "convex",
                       direction = "decreasing")
knaa$elbow

# Normalizationa and integration
DefaultAssay(object = osn.no.or.sub.integrated.sct) <- "integrated"

osn.no.or.sub.integrated.sct <- FindNeighbors(object = osn.no.or.sub.integrated.sct, 
                                              dims = 1:15)

osn.no.or.sub.integrated.sct <- FindClusters(object = osn.no.or.sub.integrated.sct,
                                             resolution = 1.6)

osn.no.or.sub.integrated.sct <- RunUMAP(object = osn.no.or.sub.integrated.sct, 
                                        dims = 1:15)

DimPlot(osn.no.or.sub.integrated.sct,label = T,pt.size = 1.5) + ggtitle("NO OR _ PC15")

# Clusters identity
osn.no.or.sub.integrated.sct <- RenameIdents(object = osn.no.or.sub.integrated.sct,
                                             "0" = "Dorsal", 
                                             "1" = "Ventral",  
                                             "2" = "Ventral", 
                                             "3" = "Ventral",
                                             "4" = "Dlg2+ V.", 
                                             "5" = "Ventral", 
                                             "6" = "Ventral",
                                             "7" = "Cd36+ V.", 
                                             "8" = "Dlg2+ D.", 
                                             "9" = "Calb2+ D.", 
                                             "10" = "Calb2+ V.", 
                                             "11" = "Ventral", 
                                             "12" = "Ventral", 
                                             "13" = "Ventral",
                                             "14" = "Cd55+ V.", 
                                             "15" = "Dorsal", 
                                             "16" = "Ventral", 
                                             "17" = "Cd36+ D.", 
                                             "18" = "Dorsal", 
                                             "19" = "Dlg2+ V.", 
                                             "20" = "Cd55+ D.",
                                             "21" = "Calb2+ V.", 
                                             "22" = "Dorsal")
# Determine Idents
osn.no.or.sub.integrated.sct$clusters <- Idents(object = osn.no.or.sub.integrated.sct)

# Determine Idents level
my_levels <- rev(x = c("Ventral",
                       "Calb2+ V.",
                       "Dlg2+ V.",
                       "Cd36+ V.",
                       "Cd55+ V.",
                       "Dorsal",
                       "Calb2+ D.",
                       "Dlg2+ D.",
                       "Cd36+ D.",
                       "Cd55+ D."))

# Create a new variable in meta.data containing the cluster identity
osn.no.or.sub.integrated.sct@active.ident <- factor(x = osn.no.or.sub.integrated.sct@active.ident, levels = my_levels)

#%%%%% Check Point %%%%%#

saveRDS(object = osn.no.or.sub.integrated.sct,file = "moe_noexp_scRNAseq10x_cellranger_filter_OSN_NO_OLFR_sct_integrated_clusters.rds")

osn.no.or.sub.integrated.sct <- readRDS(file = "moe_noexp_scRNAseq10x_cellranger_filter_OSN_NO_OLFR_sct_integrated_clusters.rds")


###############################
#                             #        
#   OSN With OR analysis      #
#                             #
###############################

# Subset object to include only cells.to.keep
osn.data.sub <- readRDS(file = "moe_noexp_scRNAseq10x_cellranger_filter_OSN.rds")

# Normalizationa and integration
DefaultAssay(osn.data.sub) <- "RNA"

osn.sub.filt.list <- SplitObject(osn.data.sub, 
                                 split.by = "orig.ident")

osn.sub.filt.list <- lapply(X = osn.sub.filt.list, 
                            FUN = function(x){SCTransform(x,vars.to.regress = "percent.mt")})

features.no.or <- SelectIntegrationFeatures(object.list = osn.sub.filt.list, 
                                            nfeatures = 5000)

osn.sub.filt.list <- PrepSCTIntegration(object.list = osn.sub.filt.list,
                                        anchor.features = features.no.or)

osn.sub.anchors <- FindIntegrationAnchors(object.list = osn.sub.filt.list, 
                                          normalization.method = "SCT", 
                                          anchor.features = features.no.or)

osn.sub.integrated.sct <- IntegrateData(anchorset = osn.sub.anchors, 
                                        normalization.method = "SCT")

# Find PCA
DefaultAssay(osn.sub.integrated.sct) <- "integrated"

osn.sub.integrated.sct <-  RunPCA(osn.sub.integrated.sct)

# Select PCA

kn <- import("kneed")
plt <- import(module = "matplotlib.pyplot")


knaa <- kn$KneeLocator(x = seq_len(length.out = length(x = osn.sub.integrated.sct@reductions$pca@stdev)),
                       y = osn.sub.integrated.sct@reductions$pca@stdev,
                       S = 1,
                       curve = "convex",
                       direction = "decreasing")
knaa$knee

kneed.res$plot_knee_normalized(figsize = c(10,10)/2.54)
plt$xlim(0, 1)
plt$ylim(0, 1)
plt$savefig("figure_osn.pdf")

# Clustering

osn.sub.integrated.sct <- FindNeighbors(object = osn.sub.integrated.sct, 
                                        dims = 1:14)

osn.sub.integrated.sct <- FindClusters(object = osn.sub.integrated.sct,
                                       resolution = 2.1)

osn.sub.integrated.sct <- RunUMAP(object = osn.sub.integrated.sct, 
                                  dims = 1:14)


DimPlot(osn.sub.integrated.sct,pt.size = 1,label = T) + ggtitle("OR")

ElbowPlot(osn.sub.integrated.sct,ndims = 50) + ggtitle("OR")

FeaturePlot(osn.sub.integrated.sct,features = c("Nfix","Dlg2","Cd55","Cd36"))

# Determine Clusters
osn.sub.integrated.sct <- RenameIdents(object = osn.sub.integrated.sct,
                                       "0" = "Dorsal", 
                                       "1" = "Ventral", 
                                       "2" = "Ventral", 
                                       "3" = "Ventral", 
                                       "4" = "Ventral", 
                                       "5" = "Cd36+ V.", 
                                       "6" = "Dlg2+ V.",
                                       "7" = "Calb2+ D.", 
                                       "8" = "Dlg2+ D.", 
                                       "9" = "Ventral", 
                                       "10" = "Calb2+ V.", 
                                       "11" = "Cd55+ V.", 
                                       "12" = "Ventral", 
                                       "13" = "Dorsal",
                                       "14" = "Dorsal", 
                                       "15" = "Ventral", 
                                       "16" = "Ventral", 
                                       "17" = "Ventral", 
                                       "18" = "Ventral", 
                                       "19" = "Cd36+ D.",
                                       "20" = "Ventral",
                                       "21" = "Dlg2+ V.",
                                       "22" = "Ventral",
                                       "23" = "Cd55+ D.", 
                                       "24" = "Dlg2+ V.", 
                                       "25" = "Ventral", 
                                       "26" = "Calb2+ V.",
                                       "27" = "Dorsal",
                                       "28" = "Calb2+ V.")

# Determine Idents
osn.sub.integrated.sct$clusters <- Idents(osn.sub.integrated.sct)

# Determine Idents level
my_levels <- rev(c("Ventral",
                   "Calb2+ V.",
                   "Dlg2+ V.",
                   "Cd36+ V.",
                   "Cd55+ V.",
                   "Dorsal",
                   "Calb2+ D.",
                   "Dlg2+ D.",
                   "Cd36+ D.",
                   "Cd55+ D."))

# Create a new variable in meta.data containing the cluster identity
osn.sub.integrated.sct@active.ident <- factor(x = osn.sub.integrated.sct@active.ident, levels = my_levels)

DimPlot(osn.sub.integrated.sct,label = T)

#%%%%% Check Point %%%%%#

saveRDS(object = osn.sub.integrated.sct,file = "moe_noexp_scRNAseq10x_cellranger_filter_OSN_WITH_OLFR_sct_integrated_clusters.rds")



