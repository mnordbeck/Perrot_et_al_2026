################################################################################
# Title:        Perrot et al. (2025) - snRNAseq analysis
# Subtitle:     Extract per sample cell abundance from annotations
# Author:       Mathias Nordbeck (Paris Brain Institute)
# Date:         2025-06-17
# Version:      1.0
# Description:  
################################################################################


################################################################################
# SETUP
################################################################################

# Load required libraries
library(Seurat)
library(dplyr)
library(reshape2)

# set directories
dir_root  <- file.path("/path_to/snRNAseq_mouse_SC_2025_Perrot_et_al")
dir_data  <- file.path(dir_root,"data")
dir_out   <- file.path(dir_root,"output/Abundance")

# create output directory if it does not already exist 
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = TRUE)}

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq.rds")
data <- readRDS(input_file)

# set ident
Idents(data) <- as.character(data$Sample_ID)


################################################################################
# Extract cell abundance from annotations
################################################################################

annotation_columns <- c("annotation_CNS", "annotation_sub")

for(annotation in annotation_columns){

  message(paste0("Extract abundance from: ",annotation))

  # get absolute numbers of cell abundance per Sample_ID
  abundance_df <- table(
      data@meta.data[[annotation]],
      data$Sample_ID) %>% 
    as.data.frame.matrix()

  # convert absolute numbers in percentages
  for (i in 1:ncol(abundance_df)) {
    percentage <- (abundance_df[i]/sum(abundance_df[i])) * 100
    abundance_df[i] <- percentage[1]
  }
  abundance_df$Cell <- rownames(abundance_df)
  abundance_df <- melt(abundance_df, id.vars = "Cell")
  abundance_df <- abundance_df %>% dplyr::rename("Sample_ID"="variable",
                                                "Percentage"="value")

  # annotate SSEP velocity and conditions
  annotations <- data@meta.data %>% 
    dplyr::select(c("Sample_ID", "Condition_SSEP", "SSEP_velocity")) %>% 
    distinct(Sample_ID, Condition_SSEP, SSEP_velocity)

  abundance_df <- dplyr::left_join(
    x = abundance_df,
    y = annotations,
    by = "Sample_ID")

  # factorize groups to compare (here: Condition_SSEP)
  abundance_df$Condition_SSEP <- as.character(abundance_df$Condition_SSEP)
  abundance_df$Condition_SSEP <- factor(abundance_df$Condition_SSEP,
    levels = c("NG","HD","MS_fast","MS_slow"))

  # save table
  output_file <- paste0("Perrot_et_al_abundance_",annotation,".csv")
  write.csv(abundance_df, file.path(dir_out,output_file))

}


################################################################################
# End of script
################################################################################

message("End of script")
