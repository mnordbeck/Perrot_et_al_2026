################################################################################
# Title:        Perrot et al. (2025) - snRNAseq analysis
# Subtitle:     Refine oligodendrocyte annotations
# Author:       Mathias Nordbeck (Paris Brain Institute)
# Date:         2025-06-17
# Version:      2.0
# Description:  This script performs subclustering of selected CNS populations 
#               using sample-wise normalization, variable feature selection, PCA,
#               Harmony integration, graph-based clustering and UMAP visualization.
#               Marker genes are identified using Wilcoxon tests and pairwise 
#               differential expression analyses are performed using MAST.
#
#               Populations analyzed:
#               - Oligodendrocytes and OPC
################################################################################


################################################################################
# SETUP
################################################################################

# Load required libraries
library(Seurat)
library(dplyr)
library(harmony)
library(ggprism)
library(ggplot2)
library(org.Mm.eg.db)
library(openxlsx)

# set directories
dir_root  <- file.path("/path_to/snRNAseq_mouse_SC_2025_Perrot_et_al")
dir_data  <- file.path(dir_root,"data")
dir_out   <- file.path(dir_root,"output/Clustering_sub")

# create output directory if it does not already exist 
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = TRUE)}

# read wrapper functions
source(file.path(dir_root,"scripts","Perrot_et_al_snRNAseq_wrapper_functions.R"))


################################################################################
# Refined oligodendrocytes and progenitors annotation
################################################################################

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq.rds")
data <- readRDS(input_file)

# set sample idents
Idents(data) <- as.character(data$Sample_ID)

# set parameters
population         <- "Oligo_OPC"     # name target population
n_pcs              <- 40              # n PCs used for integration and clustering
cluster_resolution <- 0.3             # resolution depth for clustering

# subset dataset
data_sub <- subset(data, subset = annotation_CNS %in% c("OPC/COP","Oligodendrocytes"))

# run clustering wrapper function
data_sub <- clustering(
  data_sub,
  n_pcs = n_pcs,
  resolution = cluster_resolution,
  vars.to.regress = "nCount_RNA")



# find markers genes for clusters
Idents(data_sub) <- data_sub$clusters_harmony
markers <- runMarkers(seurat_obj = data_sub)

# save results table
output_file <- paste0("Markers_",population,"_pc_",n_pcs,"_res_",cluster_resolution,".xlsx")
openxlsx::write.xlsx(markers, file.path(dir_out, output_file))



# assign sub-populations and transcriptional states
data_sub$annotation_sub <- data_sub$clusters_harmony # created in wrapper function

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("5")] <- "OPC" # precursor cells

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("8")] <- "COP" # committing precursors

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("7")] <- "MFOL" # newly forming oligos

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("0")] <- "MOL_Il33" # mature oligo

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("4","2")] <- "MOL_Klk6" # mature oligo

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("9")] <- "MOL_C4b" # mature oligo; DA oligo?

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("3")] <- "MOL_other" # mature oligo

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("6","1")] <- "MOL_Grm7" # mature oligo



# find markers genes for annotated populations
Idents(data_sub) <- data_sub$annotation_sub
markers <- runMarkers(seurat_obj = data_sub)

# save results table
output_file <- paste0("Markers_",population,"_annotation_sub.xlsx")
openxlsx::write.xlsx(markers, file.path(dir_out, output_file))



# save dataset
output_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_subset_Oligo_OPC.rds")
saveRDS(data, file = output_file)


################################################################################
# DEG analysis: oligodendrocytes and progenitors
################################################################################

# calculate DEG: MS vs HD (mature oligodendrocytes only)
population <- "Oligodendrocytes"
cond1 <- "MS_slow"
cond2 <- "HD"

data_deg <- subset(data_sub,
                   subset = annotation_sub %in% c("OPC","COP"), 
                   invert = TRUE)

Idents(data_deg) <- data_deg$Condition_SSEP
deg <- runDEG(data_deg, ident1 = cond1, ident2 = cond2)

# save results table
output_file <- paste0("DEG_",population,"_",cond1,"_vs_",cond2,".xlsx")
openxlsx::write.xlsx(deg, file.path(dir_out, output_file))



# calculate DEG: MS vs HD (oligodendrocyte progenitors only)
population <- "OPC"
cond1 <- "MS_slow"
cond2 <- "HD"

data_deg <- subset(data_sub,
                   subset = annotation_sub == "OPC")

Idents(data_deg) <- data_deg$Condition_SSEP
deg <- runDEG(data_deg, ident1 = cond1, ident2 = cond2)

# save results table
output_file <- paste0("DEG_",population,"_",cond1,"_vs_",cond2,".xlsx")
openxlsx::write.xlsx(deg, file.path(dir_out, output_file))


################################################################################
# End of script
################################################################################

message("End of script.")