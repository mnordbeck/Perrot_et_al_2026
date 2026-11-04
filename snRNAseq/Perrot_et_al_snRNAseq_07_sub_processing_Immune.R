################################################################################
# Title:        Perrot et al. (2025) - snRNAseq analysis
# Subtitle:     Refine microglia annotations
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
# Refined microglia and immune cell annotation
################################################################################

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq.rds")
data <- readRDS(input_file)

# set sample idents
Idents(data) <- as.character(data$Sample_ID)

# set parameters
population         <- "Immune"     # name target population
n_pcs              <- 15           # n PCs used for integration and clustering
cluster_resolution <- 1            # resolution depth for clustering

# subset dataset
data_sub <- subset(data, subset = annotation_CNS %in% population)

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
                             c("1","2","3","4")] <- "MG_hom" # microglia homeostatic

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("0","5")] <- "MG_Apoe" # microglia DAM

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("14","10")] <- "MG_Mki67" # microglia proliferating

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("7")] <- "DAM_a" # microglia DAM

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("8")] <- "DAM_b" # microglia DAM

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("6")] <- "IRM" # microglia interferon-response

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("13")] <- "MONO" # monocytes

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("9")] <- "MAC" # macrophages / border-associated

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("11")] <- "BLY" # unclear if human/mouse

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("16")] <- "TLY" # unclear if human/mouse

data_sub$annotation_sub[data_sub[["clusters_harmony"]] %in% 
                             c("12")] <- "NEUT" # neutrophils

# group microglia DAM signatures together
data_sub$annotation_sub[data_sub$annotation_sub %in% 
                             c("DAM_a","DAM_b", "MG_Apoe")] <- "DAM"

# remove unknown / outlier nuclei
data_sub <- subset(data_sub, 
                   subset = data_sub$clusters_harmony %in% c("15"),
                   invert = TRUE)



# find markers genes for annotated populations
Idents(data_sub) <- data_sub$annotation_sub
markers <- runMarkers(seurat_obj = data_sub)

# save results table
output_file <- paste0("Markers_",population,"_annotation_sub.xlsx")
openxlsx::write.xlsx(markers, file.path(dir_out, output_file))



# save dataset
output_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_subset_Immune.rds")
saveRDS(data, file = output_file)


################################################################################
# DEG analysis: microglia
################################################################################

# calculate DEG: DAM vs MG_hom
cond1 <- "DAM"
cond2 <- "MG_hom"
Idents(data_sub) <- data_sub$annotation_sub
deg <- runDEG(data_sub, ident1 = cond1, ident2 = cond2)

# save results table
output_file <- paste0("DEG_",population,"_",cond1,"_vs_",cond2,".xlsx")
openxlsx::write.xlsx(deg, file.path(dir_out, output_file))



# calculate DEG: IRM vs MG_hom
cond1 <- "IRM"
cond2 <- "MG_hom"
Idents(data_sub) <- data_sub$annotation_sub
deg <- runDEG(data_sub, ident1 = cond1, ident2 = cond2)

# save results table
output_file <- paste0("DEG_",population,"_",cond1,"_vs_",cond2,".xlsx")
openxlsx::write.xlsx(deg, file.path(dir_out, output_file))




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
# Clean Seurat object for publication
################################################################################

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_subset_Immune.rds")
data <- readRDS(input_file)

# set sample idents
Idents(data) <- as.character(data$Sample_ID)

# ids not to publish (confidential donor identifyers)
data$SMS_ID <- data$Sample_ID
data$SMS_annotation <- data$annotation_sub

# re-organize meta data columns
data$Cell_ID <- data$cell_id
data$Sample_ID <- data$Paper_ID
data$Donor_Sex <- data$Sex
data$Donor_Age <- data$Age
data$Donor_Treatment <- data$Treatment
data$annotation_sub <- data$annotation_Perrot_et_al

data@meta.data <- data@meta.data %>% dplyr::select(
  c("Cell_ID", "Sample_ID", "Donor_Sex", "Donor_Age", "Donor_Treatment", 
    "Condition", "Condition_SSEP", "SSEP_velocity", 
    "nCount_RNA", "nFeature_RNA",
    "percent.mt", "percent.rb", "percent.hb",
    "annotation_CNS", "annotation_sub", "SMS_ID", "SMS_annotation"))

# set sample idents
Idents(data) <- as.character(data$Sample_ID)

# clean umap
data@reductions$umap <- data@reductions$umap.harmony
data@reductions$umap.harmony <- NULL

# project name
data@project.name <- "Perrot_et_al"

# extract meta.data
meta <- data@meta.data
write.csv(meta, file.path(dir_data,"Perrot_et_al_snRNAseq_meta.csv"))

umap <- data@reductions$umap
write.csv(meta, file.path(dir_data,"Perrot_et_al_snRNAseq_umap_CNS.csv"))


output_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_subset_Immune.rds")
saveRDS(data, input_file)

output_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_subset_Oligo_OPC.rds")
saveRDS(data, input_file)


################################################################################
# End of script
################################################################################

message("End of script.")