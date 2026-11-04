################################################################################
# Title:        Perrot et al. (2025) - snRNAseq analysis
# Subtitle:     Extract count matrix and meta data for GEO upload
# Author:       Mathias Nordbeck (Paris Brain Institute)
# Date:         2026-09-30
# Version:      1.0
# Description:  This script extracts the raw cellbender adjusted count matrix
#               and per-nuclei metadatafrom the Seurat object for upload the 
#               the GEO repository.
################################################################################


################################################################################
# Setup
################################################################################

# libraries
library(Seurat)
library(Matrix)
library(dplyr)

# set directories
dir_root       <- file.path("/path_to/snRNAseq_mouse_SC_2025_Perrot_et_al/ScienceAdvances")
dir_data       <- file.path(dir_root,"data")


################################################################################
# Extract cellbender adjusted count matrix
################################################################################

# read dataset
input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq.rds")
data <- readRDS(input_file)

# extract cellbender adjusted count matrix
counts <- LayerData(data, assay = "RNA", layer = "counts")
writeMM(counts, file.path(dir_data,"matrix.mtx"))

# extract barcodes
barcodes <- colnames(counts)
writeLines(barcodes, file.path(dir_data,"barcodes.tsv"))

# extract features
features <- rownames(counts)
feature_df <- data.frame(
  gene_id = features,
  gene_name = features,
  feature_type = "Gene Expression")
write.table(
  feature_df,
  file.path(dir_data, "features.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE)

# compress the individual files
R.utils::gzip(file.path(dir_data,"matrix.mtx"), overwrite = TRUE)
R.utils::gzip(file.path(dir_data,"barcodes.tsv"), overwrite = TRUE)
R.utils::gzip(file.path(dir_data,"features.tsv"), overwrite = TRUE)

# put them together in one archive
tar(
  file.path(dir_data,"Perrot_et_al_feature_bc_matrix.tar.gz"),
  files = c(
    file.path(dir_data,"matrix.mtx.gz"),
    file.path(dir_data,"barcodes.tsv.gz"),
    file.path(dir_data,"features.tsv.gz")
  ),
  compression = "gzip"
)


################################################################################
# Extract per-nuclei metadata (incl. UMAP embeddings)
################################################################################

# annotate UMAP embeddings in per-nuclei metadata
umap <- Embeddings(data, reduction = "umap")
data[["UMAP_1_CNS"]] <- umap[, 1]
data[["UMAP_2_CNS"]] <- umap[, 2]

input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_subset_Immune.rds")
data_sub <- readRDS(input_file)
umap <- Embeddings(data_sub, reduction = "umap")
data@meta.data[rownames(umap), "UMAP_1_Immune"] <- umap[, 1]
data@meta.data[rownames(umap), "UMAP_2_Immune"] <- umap[, 2]

input_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_subset_Oligo_OPC.rds")
data_sub <- readRDS(input_file)
umap <- Embeddings(data_sub, reduction = "umap")
data@meta.data[rownames(umap), "UMAP_1_Oligo"] <- umap[, 1]
data@meta.data[rownames(umap), "UMAP_2_Oligo"] <- umap[, 2]

# extract per-nuclei metadata
meta <- data@meta.data

write.table(
  meta,
  gzfile(file.path(dir_data, "Perrot_et_al_nuclei_metadata.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

################################################################################
# End of script
################################################################################
