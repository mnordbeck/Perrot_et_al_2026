################################################################################
# Title:        Perrot et al. (2025) - snRNAseq analysis
# Subtitle:     Data processing, integration, clustering, and CNS annotation
# Author:       Mathias Nordbeck (Paris Brain Institute)
# Date:         2025-04-11
# Version:      2.0
# Description:  Processes QC-filtered singlet nuclei using Harmony integration. 
#               The workflow includes: 
#               - Sample-wise normalization and scaling
#               - PCA, Harmony integration, clustering, and UMAP reduction
#               - Marker gene identification 
#               - Broad CNS population annotation 
#               - Integration of supplementary sample metadata
################################################################################


################################################################################
# SETUP
################################################################################

# Load required libraries
library(Seurat) 
library(harmony) 
library(ggplot2) 
library(org.Mm.eg.db) 
library(AnnotationDbi) 
library(openxlsx)

# set directories
dir_root  <- file.path("/path_to/snRNAseq_mouse_SC_2025_Perrot_et_al")
dir_data  <- file.path(dir_root,"data")
dir_out   <- file.path(dir_root,"output/Clustering_CNS")

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
n_pcs              <- 20     # n PCs used for integration and clustering
cluster_resolution <- 2      # resolution depth for clustering

# run clustering wrapper function
data <- clustering(
  data, 
  n_pcs = n_pcs,
  resolution = cluster_resolution,
  vars.to.regress = "nCount_RNA")


################################################################################
# Find marker genes
################################################################################

# find markers genes for clusters
Idents(data) <- data$clusters_harmony # created in wrapper function
markers <- runMarkers(seurat_obj = data_sub)

# save results table
output_file <- paste0("Markers_CNS_pc_",n_pcs,"_res_",cluster_resolution,".xlsx")
openxlsx::write.xlsx(markers, file.path(dir_out, output_file))


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
# together with corresponding marker-expression plot and table generated above.

oligodendrocyte_clusters <- c("28", "30", "12", "15", "17", "18", "19", "14", 
  "36", "23", "9", "4", "6", "11", "2", "3", "0")

astrocyte_clusters <- c("7", "42") 

vascular_clusters <- c("46", "57", "59") 

opc_clusters <- c("8", "53") 

immune_clusters <- c("1", "34", "54") 

unknown_clusters <- c("56", "58")


# initialize annotation as neurons 
data$annotation_CNS <- "Neurons"

# assign broad CNS annotations
data$annotation_CNS[data$clusters_harmony %in% 
  oligodendrocyte_clusters] <- "Oligodendrocytes"

data$annotation_CNS[data$clusters_harmony %in% 
  astrocyte_clusters] <- "Astrocytes"

data$annotation_CNS[data$clusters_harmony %in% 
  vascular_clusters] <- "Vascular"

data$annotation_CNS[data$clusters_harmony %in% 
  opc_clusters] <- "OPC/COP"

data$annotation_CNS[data$clusters_harmony %in% 
  immune_clusters] <- "Immune" 

data$annotation_CNS[data$clusters_harmony %in%
  unknown_clusters] <- "unknown"


# removing clusters of unknown CNS identity
data <- subset(data, subset = annotation_CNS != "unknown")


message("Nuclei annotation and removal of unknown clusters: ", ncol(data))


################################################################################
# Find marker genes for broad CNS annotations
################################################################################

# find markers genes for clusters
Idents(data) <- data$annotation_CNS
markers <- runMarkers(seurat_obj = data_sub)

# save results table
output_file <- paste0("Markers_annotation_CNS.xlsx")
openxlsx::write.xlsx(markers, file.path(dir_out, output_file))


################################################################################
# Add supplementary sample metadata
################################################################################

message("Loading supplementary sample metadata.")

meta_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_meta.csv")
meta <- read.csv(meta_file, sep = ";", stringsAsFactors = FALSE)


message("Adding supplementary sample metadata.") 

# preserve the existing cell order 
cell_ids <- rownames(data@meta.data) 

# add sample-level metadata 
data@meta.data <- data@meta.data |> left_join(meta, by = "Sample_ID")

# restore cell IDs as row names
rownames(data@meta.data) <- cell_ids


################################################################################
# Save dataset
################################################################################

message("Saving dataset.")

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
