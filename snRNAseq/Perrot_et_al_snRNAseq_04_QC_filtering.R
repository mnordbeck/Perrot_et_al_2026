################################################################################
# Title:        Perrot et al. (2025) - snRNAseq analysis
# Subtitle:     Quality control and filtering of nuclei
# Author:       Mathias Nordbeck (Paris Brain Institute)
# Date:         2025-03-13
# Version:      1.0
# Description:  Performs quality control (QC) and filtering of nuclei.
#               QC includes: 
#               - Filtering genes on minimum expression across nuclei 
#               - Filtering nuclei on detected genes and total counts 
#               - Filtering nuclei on mitochondrial and ribosomal RNA content 
#               - Removing nuclei with unusually high gene or count numbers as a 
#                 proxy for potential doublets.
################################################################################


################################################################################
# SETUP
################################################################################

# required libraries
library(Seurat)
library(parallel)

# set directories
dir_root  <- file.path("/path_to/snRNAseq_mouse_SC_2025_Perrot_et_al")
dir_data  <- file.path(dir_root,"data")
dir_out   <- paste0(dir_root,"output/QC/")

# create output directory if it does not already exist 
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = TRUE)}

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq.rds")
data <- readRDS(input_file)


################################################################################
# Calculate QC parameters
################################################################################

message("Calculating QC parameters.")

# calculate mitochondrial gene percentage
data <- PercentageFeatureSet(data, "^mt-", col.name = "percent.mt")

# calculate ribosomal gene percentage
data <- PercentageFeatureSet(data, "^Rt[sl]", col.name = "percent.rb")

# calculate hemoglobin gene percentage
data <- PercentageFeatureSet(data, "^Hb[^(p)]", col.name = "percent.hb")


################################################################################
# Define thresholds for filtering
################################################################################

message("Setting QC filter thresholds.")

# minimum number of nuclei expressing a gene for the gene to be retained 
filter_exp_min <- 5 

# minimum number of detected genes and counts required per nucleus 
filter_genes_min <- 300 
filter_counts_min <- 500 

# maximum number of detected genes and counts allowed per nucleus 
# high values may indicate potential doublets or multiplets 
filter_genes_max <- 7000 
filter_counts_max <- 30000 

# maximum mitochondrial and ribosomal content allowed per nucleus 
filter_mt_max <- 0.10 
filter_rb_max <- 0.16


################################################################################
# Reduce genetic content
################################################################################

message("Filtering genes expressed in fewer than ", filter_exp_min, " nuclei.")

# count features before filtering
nb_features_before <- nrow(data[["RNA"]])

# extract raw count matrix
exp_mtx <- data[["RNA"]]$counts

# retain genes expressed in at least filter_exp_min nuclei
selected_features <- rownames(exp_mtx)[
  Matrix::rowSums(exp_mtx > 0) >= filter_exp_min ]

# filter features
data <- subset(data, features = selected_features)

# count features after filtering
nb_features_after <- nrow(data[["RNA"]])

# save feature-filter table
feature_filter_table <- data.frame(
  Filter = c( "F0_Original", paste0("F1_feature_exp_min_",filter_exp_min,"_nuclei")), 
  n_features = c( nb_features_before, nb_features_after))

write.csv(feature_filter_table, 
  file.path(dir_out, "QC_filter_features.csv"), 
  row.names = FALSE)


################################################################################ 
# INITIALIZE QC TABLE 
###############################################################################

message("Initializing QC summary.") 

# Function to count nuclei per sample while retaining samples with zero nuclei 
count_nuclei <- function(seurat_object, sample_names) {
  counts <- table(factor(seurat_object$Sample_ID,
    levels = sample_names)) 
  
  data.frame(Sample = names(counts),
    Nuclei = as.integer(counts), 
    row.names = NULL)
} 

# Store initial number of nuclei 
QC_table <- count_nuclei(data, samples)
names(QC_table)[2] <- "F0_Original"


################################################################################
# Filter nuclei based on gene and count threholds
################################################################################

message("Filtering nuclei with unusually low gene or count numbers.")

# filter nuclei with too few detected genes 
data <- subset(data, subset = nFeature_RNA >= filter_genes_min)

# filter nuclei with too few counts 
data <- subset(data, subset = nCount_RNA >= filter_counts_min)

# record nuclei retained after low-quality filtering 
QC_table_append <- count_nuclei(data, samples)
names(QC_table_append)[2] <- "F1_Minimum_Genes_Counts"

QC_table <- merge(QC_table, QC_table_append, 
  by = "Sample", all = TRUE, sort = FALSE )



message("Filtering nuclei with unusually high gene or count numbers.")

# remove nuclei with unusually high numbers of detected genes 
data <- subset(data, subset = nFeature_RNA <= filter_genes_max)

# remove nuclei with unusually high total counts 
data <- subset( data, subset = nCount_RNA <= filter_counts_max)

# record nuclei retained after high-count/high-feature filtering 
QC_table_append <- count_nuclei(data, samples)
names(QC_table_append)[2] <- "F2_Potential_Doublets" 

QC_table <- merge(QC_table, QC_table_append, 
  by = "Sample", all = TRUE, sort = FALSE)



message("Filtering nuclei based on mitochondrial and ribosomal content.")

# filter nuclei with excessive mitochondrial content 
data <- subset(data, subset = percent.mt <= filter_mt_max) 

# filter nuclei with excessive ribosomal content 
data <- subset(data, subset = percent.rb <= filter_rb_max) 

# record nuclei retained after mitochondrial/ribosomal filtering 
QC_table_append <- count_nuclei(data, samples)
names(QC_table_append)[2] <- "F3_Mito_Ribo" 

QC_table <- merge(QC_table, QC_table_append, 
  by = "Sample", all = TRUE, sort = FALSE)


message("Saving QC summary tables.")

write.csv(QC_table, file.path(dir_out, "QC_filter_nuclei.csv"), row.names = FALSE)


################################################################################
# Save QC-filtered dataset
################################################################################

message("Saving QC-filtered Seurat object.")

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
