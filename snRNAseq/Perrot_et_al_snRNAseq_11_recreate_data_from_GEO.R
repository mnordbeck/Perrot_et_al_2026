################################################################################
# Title:        Perrot et al. (2025) - snRNAseq analysis
# Subtitle:     Recreate Seurat object from GEO repository
# Author:       Mathias Nordbeck (Paris Brain Institute)
# Date:         2026-09-30
# Version:      1.0
# Description:  This script creates the Seurat object from the GEO repository.
################################################################################


################################################################################
# Setup
################################################################################

# libraries
library(Seurat)
library(Matrix)

# set directories
dir_root       <- file.path("/path_to/snRNAseq_mouse_SC_2025_Perrot_et_al/ScienceAdvances")
dir_data       <- file.path(dir_root,"data")


# FASTQ files, per-nuclei meta data and the processed Seurat object can be found
# in the Gene Expression Omnibus (GEO) repository under the accession number:
# GSE351304.

# Download the following files from the GEO repository into dir_data:
# 1) Perrot_et_al_feature_bc_matrix.tar.gz
# 2) Perrot_et_al_nuclei_metadata.csv

# Extract the Perrot_et_al_feature_bc_matrix.tar.gz archive.
# The following files should be present in dir_data:
# 1) barcodes.tsv.gz
# 2) features.tsv.gz
# 3) matrix.mtx.gz


################################################################################
# Recreate Seurat object from cellbender raw count matrix
################################################################################

# read the count matrix
counts <- Read10X(
  data.dir = file.path(dir_data))

# create Seurat object
data <- CreateSeuratObject(
  counts = counts,
  project = "Perrot_et_al")


################################################################################
# Recover per-nuclei metadata
################################################################################

# read metadata
meta <- read.table(
  gzfile(file.path(dir_data,"Perrot_et_al_nuclei_metadata.tsv.gz")),
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE)

# add metadata to Seurat object
meta <- meta[colnames(data), , drop = FALSE]
data <- AddMetaData(
  data,
  metadata = meta)


################################################################################
# Recover UMAP embeddings
################################################################################

# recreate UMAP embeddings for Seurat object
# CNS UMAP
cns_umap <- as.matrix(
  data@meta.data[, c("UMAP_1_CNS", "UMAP_2_CNS")])
colnames(cns_umap) <- c("UMAPCNS_1", "UMAPCNS_2")
data[["umap_CNS"]] <- CreateDimReducObject(
  embeddings = cns_umap,
  key = "UMAPCNS_",
  assay = DefaultAssay(data))

# Immune UMAP
immune_umap <- as.matrix(
  data@meta.data[, c("UMAP_1_Immune", "UMAP_2_Immune")])
colnames(immune_umap) <- c("UMAPImmune_1", "UMAPImmune_2")
data[["umap_Immune"]] <- CreateDimReducObject(
  embeddings = immune_umap,
  key = "UMAPImmune_",
  assay = DefaultAssay(data))

# Oligo/OPC UMAP
oligo_umap <- as.matrix(
  data@meta.data[, c("UMAP_1_Oligo", "UMAP_2_Oligo")])
colnames(oligo_umap) <- c("UMAPOligo_1", "UMAPOligo_2")
data[["umap_Oligo"]] <- CreateDimReducObject(
  embeddings = oligo_umap,
  key = "UMAPOligo_",
  assay = DefaultAssay(data))


################################################################################
# Normalize count matrix
################################################################################

# normalize count matrix
data <- NormalizeData(data)


################################################################################
# Save data
################################################################################

# Save to RDS
output_file <- file.path(dir_data, "Perrot_et_al_snRNAseq_processed.rds")
saveRDS(data, file = output_file)


################################################################################
# End of script
################################################################################

# At this point, a functional Seurat object has been reconstructed.


# Quality-control parameters and nuclei annotations (e.g. annotation_CNS and 
# annotation_sub) are available in the metadata.

# The CNS, Immune, and Oligo/OPC UMAP embeddings are available both as metadata
# columns and as Seurat dimensionality reductions.

# The CellBender-adjusted count matrix has been used to recreate the RNA assay, 
# was normalized, and can be used for downstream analyses, including differential
# gene expression analysis.