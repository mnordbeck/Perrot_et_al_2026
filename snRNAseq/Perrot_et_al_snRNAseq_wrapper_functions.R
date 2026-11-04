################################################################################
# Title:        Perrot et al. (2025) - snRNAseq analysis
# Subtitle:     Helper functions
# Author:       Mathias Nordbeck (Paris Brain Institute)
# Date:         2025-06-17
# Version:      2.0
# Description:  This script contains helper/wrapper functions, facilitating
#               readability in the main analysis scripts.
#
#               Functions included:
#               - clustering(): main processing and clusterin of Seurat objects
#               - runMarkers(): identify marker genes using xilcoxon rank sum
#               - runDEG(): calculate differentially expressed genes using MAST
################################################################################


################################################################################
# Define function to process and cluster Seurat objects
################################################################################

# wrapper function for processing a seurat object
clustering <- function(
  seurat_obj, 
  resolution = 0.5,
  vars.to.regress = NULL,
  n_pcs = n_pcs,
  cluster_name = "clusters_harmony") {
  
  print("Split RNA assay by Sample_ID.")
  # split RNA data by sample ID for integration
  seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], 
                               f = seurat_obj@meta.data$Sample_ID)

  message("Normalizing data.")
  # normalize data using log normalization (scaling factor: 10,000)
  seurat_obj <- NormalizeData(seurat_obj, 
                        normalization.method = "LogNormalize", 
                        scale.factor = 10000)

  message("Finding variable features.")
  # identify the 2,000 most variable features
  seurat_obj <- FindVariableFeatures(seurat_obj, 
                              selection.method = "vst", 
                              nfeatures = 2000)

  message("Scaling data.")
  # scale data
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = vars.to.regress)
  

  message("Calculating principal components.")
  # perform PCA on variable features
  seurat_obj <- RunPCA(seurat_obj, 
    features = VariableFeatures(object = seurat_obj))

  message("Using ", n_pcs, " principal components for downstream analysis.")
  

  message("Harmony integration.")
  # integrate sample-specific layers using Harmony
  seurat_obj <- IntegrateLayers(
        object = seurat_obj,
        method = HarmonyIntegration,
        orig.reduction = "pca",
        new.reduction = "harmony",
        verbose = TRUE)

  message("Finding nearest neighbors.")
  # Find nearest neighbors based on Harmony reduction
  seurat_obj <- FindNeighbors(
        seurat_obj,
        dims = 1:n_pcs,
        reduction = "harmony")

  message("Running clustering.")
  # clustering on predefined resolution
  seurat_obj <- FindClusters(
        seurat_obj,
        resolution = cluster_resolution,
        cluster.name = cluster_name,
        verbose = TRUE)

  message("Calculating UMAP.")
  # calculate umap
  seurat_obj <- RunUMAP(
        seurat_obj, 
        reduction = "harmony",
        reduction.name = "umap",
        dims = 1:n_pcs,
        verbose = TRUE)
  
  # Join layers
  seurat_obj <- JoinLayers(seurat_obj)

  # set identities to Harmony clusters 
  Idents(seurat_obj) <- seurat_obj[[cluster_name]]
    
  return(seurat_obj)
}



################################################################################
# Define function to calculate marker genes
################################################################################

# load libraries
library(org.Mm.eg.db)
library(openxlsx)

# wrapper function for identify marker genes
runMarkers <- function(
  seurat_obj,
  min.pct = 0.25,
  logfc.threshold = 0.5,
  test.use = "wilcox",
  only.pos = FALSE){

  message("Finding markers.")

  # differential expression analysis
  de <- FindAllMarkers(
      seurat_obj, 
      min.pct = min.pct,
      logfc.threshold = logfc.threshold,
      test.use = test.use,
      only.pos = only.pos)

  # add gene annotations from org.Mm.eg.db
  gene_name <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = de$gene,
    keytype = "SYMBOL",
    column = "GENENAME",
    multiVals = "first")
  
  gene_type <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = de$gene,
    keytype = "SYMBOL",
    column = "GENETYPE",
    multiVals = "first") 
  
  de$GeneName <- unname(gene_name) 
  de$GeneType <- unname(gene_type) 
  de$FeatureID <- de$gene 
      
  de$DEG <- c("ns")
  de$DEG[de$p_val_adj <= 0.05 & de$avg_log2FC >= logfc.threshold] <- "UP"
  de$DEG[de$p_val_adj <= 0.05 & de$avg_log2FC <= -logfc.threshold] <- "DOWN"
      
  # select and sort columns 
  de <- de[, c("cluster", "FeatureID", "GeneName", "GeneType", "p_val", 
    "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "DEG")]
  
  return(de)
}


################################################################################
# Define function to calculate differentially expressed genes
################################################################################

# wrapper function for identify differentially expressed genes
runDEG <- function(
  seurat_obj,
  ident1, 
  ident2,
  logfc.threshold = 0.25,
  test.use = "MAST",
  only.pos = FALSE){

  message("Calculating DEG.")
      
  # differential expression analysis
  de <- FindMarkers(
      seurat_obj, 
      ident.1 = ident1, 
      ident.2 = ident2,
      test.use = test.use,
      features = rownames(seurat_obj))

  # add gene annotations from org.Mm.eg.db
  gene_name <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = de$gene,
    keytype = "SYMBOL",
    column = "GENENAME",
    multiVals = "first")
  
  gene_type <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = de$gene,
    keytype = "SYMBOL",
    column = "GENETYPE",
    multiVals = "first") 
  
  de$GeneName <- unname(gene_name) 
  de$GeneType <- unname(gene_type) 
  de$FeatureID <- de$gene
      
  de$DEG <- c("ns")
  de$DEG[de$p_val_adj <= 0.05 & de$avg_log2FC >= logfc.threshold] <- "UP"
  de$DEG[de$p_val_adj <= 0.05 & de$avg_log2FC <= -logfc.threshold] <- "DOWN"
      
  # select and sort columns 
  de <- de[, c("cluster", "FeatureID", "GeneName", "GeneType", "p_val", 
    "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "DEG")]

  return(de)
}

