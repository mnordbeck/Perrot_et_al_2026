################################################################################
# Title:        Perrot et al. (2025) - snRNAseq analysis
# Subtitle:     Doublet detection via integration and clustering
# Author:       Mathias Nordbeck (Paris Brain Institute)
# Date:         2025-03-13
# Version:      1.0
# Description:  Integrates samples using Harmony and performs dimensionality 
#               reduction and clustering to identify potential doublets. 
#               Broad CNS populations are first identified using marker gene 
#               expression.
#               Selected populations are then subclustered independently to 
#               identify putative doublet clusters based on co-expression of 
#               markers from distinct cell populations. 
#               Putative doublets identified at both the broad-population and 
#               subcluster levels are removed from the dataset.
################################################################################


################################################################################
# SETUP
################################################################################

# required libraries
library(Seurat)
library(dplyr)
library(harmony)

# set directories
dir_root  <- file.path("/path_to/snRNAseq_mouse_SC_2025_Perrot_et_al")
dir_data  <- file.path(dir_root,"dat")
dir_out   <- file.path(dir_root,"output","QC","Doublet_Detection_Clustering")

# create output directory if it does not already exist 
if (!dir.exists(dir_out)) {
       dir.create(dir_out, recursive = TRUE)}

# read wrapper functions
source(file.path(dir_root,"scripts","Perrot_et_al_snRNAseq_wrapper_functions.R"))

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq.rds")
data <- readRDS(input_file)

# set sample idents
Idents(data) <- as.character(data$Sample_ID)


################################################################################
# Data processing and clustering
################################################################################

# set parameters
n_pcs              <- 30     # n PCs used for integration and clustering
cluster_resolution <- 2      # resolution depth for clustering

# run clustering wrapper function
data <- clustering(
  data, 
  n_pcs = n_pcs,
  resolution = cluster_resolution,
  vars.to.regress = NULL)


################################################################################
# UMAP
################################################################################

message("Plot UMAP.") 
# generate UMAP plot 
Idents(data) <- data$clusters_harmony # created in wrapper function

plt <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.3) 
ggsave(filename = file.path(dir_out, "UMAP_doublet_detection.svg"), 
       plot = plt, width = 8, height = 8)


################################################################################
# Generate marker expression plot
################################################################################

# marker genes used to assign broad CNS populations
markers <- c(
       # oligodendrocytes
       "Mog", "Aspa", "Mbp", "Plp1",

       # neurons
       "Rbfox1", "Syt1", "Grin2b",
       
       # immune cells / microglia
       "Ptprc", "Cx3cr1", "C1qa",
       
       # astrocytes
       "Slc1a2", "Aqp4", "Gfap",
       
       # oligodendrocyte precursor cells (OPC)
       "Vcan", "Pdgfra",  "Cspg4",
       
       # fibroblasts / meninges
       "Col3a1", "Itgbl1",
       
       # pericytes
       "Pdgfrb", "Rgs5",
       
       # ciliated glia (e.g. ependymal cells)
       "Dnah12", "Cfap107",
       
       # endothelial cells
       "Pecam1", "Flt1",

       # Schwann cells
       "Mpz")


message("Generating broad marker expression plot.")

plt <- FeaturePlot(
       data, 
       features = markers, 
       reduction = umap_reduction, 
       ncol = 5)

ggsave(filename = file.path(dir_out, "QC_feature_expression.svg"), 
       plot = plt, width = 15, height = 12)


################################################################################
# Broad CNS annotation
################################################################################

# Cluster identities identified from marker expression. 
# These assignments are specific to this dataset and should be reviewed 
# together with the corresponding marker-expression plots generated above.

oligodendrocytes <- c("0", "6", "1", "2", "16", "3", "4", "10", "12", "13") 
immune           <- c("5", "27", "39", "75") 
opc              <- c("9", "78", "48", "80", "60") 
astrocytes       <- c("86", "52", "8", "37", "82", "76")
vascular         <- c("44", "72")

# Broad-population clusters identified as potential doublets
broad_doublets   <- c("52", "82", "68", "72", "39", "48", "60")


# Initialize all clusters as neurons 
data$annotation_CNS <- "Neurons"

# Assign broad cell-type annotations 
data$annotation_CNS[data$clusters_harmony %in% oligodendrocytes] <- "Oligodendrocytes" 
data$annotation_CNS[data$clusters_harmony %in% opc] <- "OPC"
data$annotation_CNS[data$clusters_harmony %in% immune] <- "Immune"
data$annotation_CNS[data$clusters_harmony %in% astrocytes] <- "Astrocytes"
data$annotation_CNS[data$clusters_harmony %in% vascular] <- "Vascular"

# Mark broad-population doublets last so that doublet annotation takes 
# precedence over the broad cell-type annotation.
data$annotation_CNS[data$clusters_harmony %in% broad_doublets] <- "doublet"

# Set identities to broad CNS annotations
Idents(data) <- data$annotation_CNS



################################################################################
# Subclustering and doublet annotation
################################################################################

# initialize dataframe to store cell-level doublet annotations
doublet_meta <- data.frame(
       Cell_ID = character(),
       annotation_sub = character(),
       stringsAsFactors = FALSE)


################################################################################
# Microglia and immune cells
################################################################################

# set parameters
population         <- "Immune"     # name target population
n_pcs              <- 30           # n PCs used for integration and clustering
cluster_resolution <- 0.5          # resolution depth for clustering

# subset dataset
data_sub <- subset(data, subset = annotation_CNS %in% population)

# run clustering wrapper function
data_sub <- clustering(
  data_sub,
  n_pcs = n_pcs,
  resolution = cluster_resolution,
  vars.to.regress = NULL)

# generate UMAP plot 
Idents(data_sub) <- data_sub$clusters_harmony # created in wrapper function

plt <- DimPlot(data_sub, reduction = "umap", label = TRUE, pt.size = 0.3) 
ggsave(filename = file.path(dir_out, paste0("UMAP_",population,"_doublet_detection.svg")), 
       plot = plt, width = 8, height = 8)

# plot broad marker expression in the subclustered population
plt <- FeaturePlot(
       data_sub, 
       features = markers,
       reduction = "umap",
       ncol = 5)

ggsave(filename = file.path(dir_out, 
       paste0("QC_feature_expression_",population,".svg")), 
       plot = plt, width = 15, height = 12)

# subcluster annotation 
immune_subclusters <- c("8", "7", "10", "9", "11")
immune_doublets <- c("4", "5")

data_sub$annotation_sub <- "Microglia"
data_sub$annotation_sub[data_sub$clusters_harmony_sub %in% 
       immune_subclusters] <- "Immune" # e.g. monocytes, BLY, TLY
data_sub$annotation_sub[data_sub$clusters_harmony_sub %in% 
       immune_doublets] <- "doublet"

# store cell-level annotation 
doublet_meta <- rbind(doublet_meta, 
       data_sub@meta.data[, c("Cell_ID", "annotation_sub")])


################################################################################
# Oligodendrocyte precursor cells (OPC)
################################################################################

# set parameters
population         <- "OPC"     # name target population
n_pcs              <- 30           # n PCs used for integration and clustering
cluster_resolution <- 0.5          # resolution depth for clustering

# subset dataset
data_sub <- subset(data, subset = annotation_CNS %in% population)

# run clustering wrapper function
data_sub <- clustering(
  data_sub,
  n_pcs = n_pcs,
  resolution = cluster_resolution,
  vars.to.regress = NULL)

# generate UMAP plot 
Idents(data_sub) <- data_sub$clusters_harmony # created in wrapper function

plt <- DimPlot(data_sub, reduction = "umap", label = TRUE, pt.size = 0.3) 
ggsave(filename = file.path(dir_out, paste0("UMAP_",population,"_doublet_detection.svg")), 
       plot = plt, width = 8, height = 8)

# plot broad marker expression in the subclustered population
plt <- FeaturePlot(
       data_sub, 
       features = markers,
       reduction = "umap",
       ncol = 5)

ggsave(filename = file.path(dir_out, 
       paste0("QC_feature_expression_",population,".svg")), 
       plot = plt, width = 15, height = 12)

# subcluster annotation
opc_doublets <- c("3", "9", "11")

data_sub$annotation_sub <- "OPC"
data_sub$annotation_sub[data_sub$clusters_harmony_sub %in% 
       opc_doublets] <- "doublet"

# store cell-level annotation 
doublet_meta <- rbind(doublet_meta, 
       data_sub@meta.data[, c("Cell_ID", "annotation_sub")])


################################################################################
# Astrocytes
################################################################################

# set parameters
population         <- "Astrocytes" # name target population
n_pcs              <- 30           # n PCs used for integration and clustering
cluster_resolution <- 0.5          # resolution depth for clustering

# subset dataset
data_sub <- subset(data, subset = annotation_CNS %in% population)

# run clustering wrapper function
data_sub <- clustering(
  data_sub,
  n_pcs = n_pcs,
  resolution = cluster_resolution,
  vars.to.regress = NULL)

# generate UMAP plot 
Idents(data_sub) <- data_sub$clusters_harmony # created in wrapper function

plt <- DimPlot(data_sub, reduction = "umap", label = TRUE, pt.size = 0.3) 
ggsave(filename = file.path(dir_out, paste0("UMAP_",population,"_doublet_detection.svg")), 
       plot = plt, width = 8, height = 8)

# plot broad marker expression in the subclustered population
plt <- FeaturePlot(
       data_sub, 
       features = markers,
       reduction = "umap",
       ncol = 5)

ggsave(filename = file.path(dir_out, 
       paste0("QC_feature_expression_",population,".svg")), 
       plot = plt, width = 15, height = 12)

# subcluster annotation
astrocyte_doublets <- c("4", "9", "7")

data_sub$annotation_sub <- "Astrocytes"
data_sub$annotation_sub[data_sub$clusters_harmony_sub %in% 
       astrocytes_doublets] <- "doublet"

# store cell-level annotation 
doublet_meta <- rbind(doublet_meta, 
       data_sub@meta.data[, c("Cell_ID", "annotation_sub")])


################################################################################
# Oligodendrocytes
################################################################################

# set parameters
population         <- "Oligodendrocytes"     # name target population
n_pcs              <- 30           # n PCs used for integration and clustering
cluster_resolution <- 0.5          # resolution depth for clustering

# subset dataset
data_sub <- subset(data, subset = annotation_CNS %in% population)

# run clustering wrapper function
data_sub <- clustering(
  data_sub,
  n_pcs = n_pcs,
  resolution = cluster_resolution,
  vars.to.regress = NULL)

# generate UMAP plot 
Idents(data_sub) <- data_sub$clusters_harmony # created in wrapper function

plt <- DimPlot(data_sub, reduction = "umap", label = TRUE, pt.size = 0.3) 
ggsave(filename = file.path(dir_out, paste0("UMAP_",population,"_doublet_detection.svg")), 
       plot = plt, width = 8, height = 8)

# plot broad marker expression in the subclustered population
plt <- FeaturePlot(
       data_sub, 
       features = markers,
       reduction = "umap",
       ncol = 5)

ggsave(filename = file.path(dir_out, 
       paste0("QC_feature_expression_",population,".svg")), 
       plot = plt, width = 15, height = 12)

# subcluster annotation
oligodendrocyte_doublets <- c("7")

data_sub$annotation_sub <- "Oligodendrocytes"
data_sub$annotation_sub[data_sub$clusters_harmony_sub %in% 
       oligodendrocyte_doublets] <- "doublet"

# store cell-level annotation 
doublet_meta <- rbind(doublet_meta, 
       data_sub@meta.data[, c("Cell_ID", "annotation_sub")])


################################################################################
# Neurons
################################################################################

# set parameters
population         <- "Neurons"    # name target population
n_pcs              <- 30           # n PCs used for integration and clustering
cluster_resolution <- 0.5          # resolution depth for clustering

# subset dataset
data_sub <- subset(data, subset = annotation_CNS %in% population)

# run clustering wrapper function
data_sub <- clustering(
  data_sub,
  n_pcs = n_pcs,
  resolution = cluster_resolution,
  vars.to.regress = NULL)

# generate UMAP plot 
Idents(data_sub) <- data_sub$clusters_harmony # created in wrapper function

plt <- DimPlot(data_sub, reduction = "umap", label = TRUE, pt.size = 0.3) 
ggsave(filename = file.path(dir_out, paste0("UMAP_",population,"_doublet_detection.svg")), 
       plot = plt, width = 8, height = 8)

# plot broad marker expression in the subclustered population
plt <- FeaturePlot(
       data_sub, 
       features = markers,
       reduction = "umap",
       ncol = 5)

ggsave(filename = file.path(dir_out, 
       paste0("QC_feature_expression_",population,".svg")), 
       plot = plt, width = 15, height = 12)

# subcluster annotation
neuron_doublets <- c("1", "26", "37", "50", "52", "58", "29", "44", "49", 
       "45", "15", "34")

data_sub$annotation_sub <- "Neurons"
data_sub$annotation_sub[data_sub$clusters_harmony_sub %in% 
       neuron_doublets] <- "doublet"

# store cell-level annotation 
doublet_meta <- rbind(doublet_meta, 
       data_sub@meta.data[, c("Cell_ID", "annotation_sub")])


################################################################################
# Vascular cells
################################################################################

# set parameters
population         <- "Vascular"   # name target population
n_pcs              <- 30           # n PCs used for integration and clustering
cluster_resolution <- 0.5          # resolution depth for clustering

# subset dataset
data_sub <- subset(data, subset = annotation_CNS %in% population)

# run clustering wrapper function
data_sub <- clustering(
  data_sub,
  n_pcs = n_pcs,
  resolution = cluster_resolution,
  vars.to.regress = NULL)

# generate UMAP plot 
Idents(data_sub) <- data_sub$clusters_harmony # created in wrapper function

plt <- DimPlot(data_sub, reduction = "umap", label = TRUE, pt.size = 0.3) 
ggsave(filename = file.path(dir_out, paste0("UMAP_",population,"_doublet_detection.svg")), 
       plot = plt, width = 8, height = 8)

# plot broad marker expression in the subclustered population
plt <- FeaturePlot(
       data_sub, 
       features = markers,
       reduction = "umap",
       ncol = 5)

ggsave(filename = file.path(dir_out, 
       paste0("QC_feature_expression_",population,".svg")), 
       plot = plt, width = 15, height = 12)

# subcluster annotation
vascular_doublets <- c("0", "6", "10", "8")

data_sub$annotation_sub <- "Vascular"
data_sub$annotation_sub[data_sub$clusters_harmony_sub %in% 
       vascular_doublets] <- "doublet"

# store cell-level annotation 
doublet_meta <- rbind(doublet_meta, 
       data_sub@meta.data[, c("Cell_ID", "annotation_sub")])


################################################################################
# Save table of doublet annotations
################################################################################

message("Saving subcluster doublet annotations.") 

write.csv(doublet_meta, 
       file = file.path(dir_out, "Subclustering_doublet_annotation.csv"),
       row.names = FALSE)


################################################################################
# Remove doublets from full dataset
################################################################################

message("Adding doublet annotations to the full dataset.")

# match and add subcluster annotations to the full dataset using Cell_ID. 
match_idx <- match(colnames(data), doublet_meta$Cell_ID) 
data$annotation_sub <- doublet_meta$annotation_sub[match_idx]


message("Removing putative doublets.")

# remove cells assigned as broad-population doublets 
data <- subset(data, subset = annotation_CNS != "doublet") 

# remove cells assigned as subcluster-level doublets 
data <- subset(data, subset = annotation_sub != "doublet")

message("Remaining nuclei: ", ncol(data))


################################################################################
# Save singlet dataset
################################################################################

message("Saving singlet dataset.")

# save dataset
output_file <- file.path(dir_data, "Perrot_et_al_snRNAseq.rds")
saveRDS(data, file = output_file)


################################################################################
# End of script
################################################################################

message("Summary:") 
message("Nuclei remaining: ", ncol(data))
message("Features remaining: ", nrow(data))

message("End of script.")
