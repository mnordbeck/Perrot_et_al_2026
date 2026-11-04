################################################################################
# Title:        Perrot et al. (2025) - snRNAseq analysis
# Subtitle:     Initialize Seurat object from Cellbender outputs
# Author:       Mathias Nordbeck (Paris Brain Institute)
# Date:         2025-03-11
# Version:      1.0
# Description:  Reads CellBender-processed count matrices, initializes one 
#               Seurat object per sample, merges all samples into a single 
#               Seurat object, and saves the resulting object as an RDS file.
################################################################################


################################################################################
# Setup
################################################################################

# required libraries
library(Seurat)
library(hdf5r)
library(scCustomize)
library(dplyr)

# set directories
dir_cellbender <- file.path("/path_to/snRNAseq_mouse_SC_2025_Perrot_et_al/cellbender_half_learning_rate")
dir_root       <- file.path("/path_to/snRNAseq_mouse_SC_2025_Perrot_et_al")
dir_data       <- file.path(dir_root,"data")

# create output directory if it does not already exist 
if (!dir.exists(dir_data)) {
  dir.create(dir_data, recursive = TRUE)}


################################################################################
# Cellbender to Seurat
################################################################################

# define samples to include in the analysis 
samples <- c("NG",
             "HD11",
             "HD12",
             "MS20",
             "MS23",
             "MS24",
             "MS25")

# initialize an empty list to store Seurat objects
data.list <- list()

# loop through each samples
print("Start importing cellbender processed samples")
for (sample in samples) {

  message("Processing sample: ", sample)

  # define CellBender input file
  input_file <- file.path(dir_cellbender, sample, paste0(sample, "_cellbender_filtered.h5") )

  # read the cellbender filtered count matrix
  sample_counts <- Read_CellBender_h5_Mat(input_file)
  
  # add sample identifier to cell barcodes to ensure unique cell IDs
  colnames(sample_counts) <- paste(sample, colnames(sample_counts), sep = "_")

  # initialize Seurat object
  sample_data <- CreateSeuratObject(counts = sample_counts, 
                                    project = "Perrot_et_al", 
                                    min.cells = 0,
                                    min.features = 0)
  
  # add a unique cell identifier to the metadata
  sample_data$Cell_ID <- colnames(sample_data)

  # store Seurat object in the list
  data.list[[sample]] <- sample_data
}


################################################################################
# Merge samples
################################################################################

message("Starting merge of samples.")

# Merge all Seurat objects into a single object
data <- merge(x = data.list[[1]], 
              y = c(data.list[-1]), 
              project = "Perrot_et_al")


message("Joining Seurat data layers.")

# Join data layers
data <- JoinLayers(data)


message("Updating metadata.") 

# Preserve the original Seurat sample identifier under a descriptive name 
data$Sample_ID <- data$orig.ident 

# Remove the original identifier column 
data$orig.ident <- NULL


################################################################################
# Save output
################################################################################

message("Saving data to RDS.") 

# Save merged Seurat object to RDS file
output_file <- file.path(dir_data, "Perrot_et_al_snRNAseq.rds")
saveRDS(data, file = output_file)


################################################################################
# End of Script
################################################################################

message("End of script")
