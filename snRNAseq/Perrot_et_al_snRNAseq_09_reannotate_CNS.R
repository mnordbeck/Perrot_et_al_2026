################################################################################
# Title:        Perrot et al. (2025) - snRNAseq analysis
# Subtitle:     Refine microglia and oligodendrocyte annotations
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
#               - Microglia and immune cells
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
# Re-assign refined annotations to global dataset
################################################################################

# read global dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq.rds")
data <- readRDS(input_file)


# refine astrocyte subpopulations from initial CNS clustering
data$annotation_sub <- data$annotation_CNS
data$annotation_sub[data$clusters_harmony %in% c("7")] <- "Astro_Slc7a10"
data$annotation_sub[data$clusters_harmony %in% c("42")] <- "Astro_Gfap"


# recover refined OPC and oligodendrocytes annotations
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_subset_Oligo_OPC.rds")
data_sub <- readRDS(input_file)

data@meta.data <- data@meta.data %>%
  left_join(data_sub@meta.data %>% dplyr::select(Cell_ID, annotation_sub),
            by = "Cell_ID", suffix = c(".old", ".new")) %>%
  mutate(annotation_sub = coalesce(annotation_sub.new, annotation_sub.old))


# recover refined microglia and immune cell annotations
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_subset_Immune.rds")
data_sub <- readRDS(input_file)

data@meta.data <- data@meta.data %>%
  left_join(data_sub@meta.data %>% dplyr::select(Cell_ID, annotation_sub),
            by = "Cell_ID", suffix = c(".old", ".new")) %>%
  mutate(annotation_sub = coalesce(annotation_sub.new, annotation_sub.old))


# select and sort meta data columns
data@meta.data <- data@meta.data %>% dplyr::select(
  c("Cell_ID", "Sample_ID", "Donor_Sex", "Donor_Age", "Donor_Treatment", 
    "Condition", "Condition_SSEP", "SSEP_velocity", 
    "nCount_RNA", "nFeature_RNA",
    "percent.mt", "percent.rb", "percent.hb",
    "annotation_CNS", "annotation_sub"))

# save dataset
Idents(data) <- data$annotation_sub
output_file <- file.path(dir_data, "Perrot_et_al_snRNAseq.rds")
saveRDS(data, file = output_file)


################################################################################
# End of script
################################################################################

message("End of script.")