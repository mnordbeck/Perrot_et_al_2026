# Single nuclei RNA sequencing 

## Overview

This repository contains the R and shell scripts used for the analysis of single
nucleus RNA sequencing data. It includes the following scripts:
1) Align FASTQ files against the GRCm39-2024-A reference genome 
([Perrot_et_al_snRNAseq_01_cellranger.sh](Perrot_et_al_snRNAseq_01_cellranger.sh))
2) Remove ambient RNA from feature-barcode-count matrices 
([Perrot_et_al_snRNAseq_02_cellbender.sh](Perrot_et_al_snRNAseq_02_cellbender.sh))
3) Initiate a Seurat object 
([Perrot_et_al_snRNAseq_03_cellbender_to_seurat.R](Perrot_et_al_snRNAseq_03_cellbender_to_seurat.R)), 
perform QC filtering ([Perrot_et_al_snRNAseq_04_QC_filtering.R](Perrot_et_al_snRNAseq_04_QC_filtering.R)), and doublet 
removal ([Perrot_et_al_snRNAseq_05_doublet_removal.R](Perrot_et_al_snRNAseq_05_doublet_removal.R)).
4) Main processing and clustering of the global dataset 
([Perrot_et_al_snRNAseq_06_CNS_processing.R](Perrot_et_al_snRNAseq_06_CNS_processing.R)).
5) Sub-clustering and processing of the microglia-immune population 
([Perrot_et_al_snRNAseq_07_sub_processing_Immune.R](Perrot_et_al_snRNAseq_07_sub_processing_Immune.R)) and oligodendrocyte lineage ([Perrot_et_al_snRNAseq_08_sub_processing_Oligo.R](Perrot_et_al_snRNAseq_08_sub_processing_Oligo.R)), as well as re-annoation of the
global dataset ([Perrot_et_al_snRNAseq_09_reannotate_CNS.R](Perrot_et_al_snRNAseq_09_reannotate_CNS.R))
6) Extraction ([Perrot_et_al_snRNAseq_10_extract_data_for_GEO.R](Perrot_et_al_snRNAseq_10_extract_data_for_GEO.R)) and recreation 
([Perrot_et_al_snRNAseq_11_recreate_data_from_GEO.R](Perrot_et_al_snRNAseq_11_recreate_data_from_GEO.R)) of the dataset uploaded to 
the Gene Expression Omnibus repository under the accession number GSE351304.
7) Creating the plots that were used for the manuscript 
([Perrot_et_al_snRNAseq_figures.R](Perrot_et_al_snRNAseq_figures.R)).
8) Wrapper functions used for clustering, marker gene identification, and 
differential gene expression analysis ([Perrot_et_al_snRNAseq_wrapper_functions.R](Perrot_et_al_snRNAseq_wrapper_functions.R)).

## Library preparation and sequencing 

Libraries were constructed using the GEM-X Universal 3′ Gene Expression v4 kit 
(10x Genomics) following the manufacturer’s instructions. Library quality and 
fragment size distribution were assessed using a TapeStation 4200 system 
(Agilent Technologies). Quantification of libraries was performed using a Qubit 
fluorometer (Thermo Scientific). Sequencing was carried out on a NovaSeq X Plus 
platform (Illumina). Eight libraries were loaded onto a single flow cell with a 
sequencing capacity of 10 billion reads. Paired-end sequencing was performed 
according to Illumina’s recommendations. 

## Pre-processing and alignment 

FASTQ files were processed using the 10X Genomics CellRanger (v.9.0.0). 
Gene-barcode matrices were generated with the count function, including the 
--include-introns flag and aligning to the GRCm39-2024-A reference mouse genome. 

## Ambient RNA correction

Ambient RNA contamination was corrected for each sample independently using 
CellBender (v.0.3.2). With the CellRanger raw_feature_bc_matrix.h5 files as 
input, the remove-background function was used with half the learning rate 
(learning-rate 0.00005), while other parameters kept at default.

## Quality control and filtering

Nuclei were filtered uniformly across samples on a minimum detection of 300 
features and 500 counts, and by excluding nuclei with > 7,000 features and > 
30,000 counts. Nuclei with mitochondria content > 0.1 % and ribosomal 
content > 0.16 % were removed. Genes expressed in fewer than 5 nuclei were 
excluded. Doublets were identified based on co-expression of canonical marker 
genes in an initial clustering and subsequently removed (see below). After 
filtering, the dataset consisted of seven samples (one non-grafted but 
LPC-lesioned, two HD LY-grafted, and four MS LY-grafted) containing a total 
of 94k nuclei and 31k genes.

## Doublet removal, data processing, and clustering

The nuclei filtered and ambient RNA corrected count matrices were analyzed using 
the Seurat (v.5.3.0) pipeline in R (v.4.3.3). 

All samples were imported to generate a single Seurat object, with the “RNA” 
assay split by sample. Each sample was independently normalized and scaled 
using NormalizeData(), FindVariableFeatures(), and 
ScaleData(vars.to.regress=“nCount_RNA”) with default parameters. 

Principal component analysis (PCA) was performed on variable features 
using RunPCA(), followed by sample integration using Harmony (v.1.2.3) via 
IntegrateLayers(method = HarmonyIntegration, orig.reduction="pca", 
new.reduction="harmony"). A uniform manifold approximation and projection 
(UMAP) was generated using RunUMAP(dims=30, reduction="harmony") for 
visualization. 

To identify and remove potential doublets, an initial clustering analysis was 
performed on all nuclei to define broad cell types, followed by subclustering of 
each annotated cell fraction. Clusters were inspected for co-expression of 
canonical marker genes from distinct cell types, indicative of doublet events. 
Nuclei assigned to such mixed-marker clusters were removed and the dataset 
re-processed. 

Broad cell-type clusters were identified using FindNeighbors(dims=20, 
reduction="harmony") and FindClusters(resolution=2). Clusters were annotated 
based on canonical marker genes: neurons (Rbfox3, Syn1, Grin2b), astrocytes 
(Slc1a2, Aqp4, Gfap), vascular cells including endothelial cells, pericytes, 
and fibroblasts (Col3a1, Flt1, Rgs5), immune cells including microglia and 
peripheral mononuclear cells (Ptprc, C1q, Cx3cr1), OPC/COP (Pdgfra, Cspg4, 
Vcan), and oligodendrocytes (Plp1, Mbp, Mog). 

To further resolve transcriptional states and subpopulations, the immune cell 
fraction and the OPC/COP plus oligodendrocyte fraction were re-integrated and 
clustered as described above. Immune cells were subclustered on the first 15 
dimensions at a resolution of 1, while OPC/COP and oligodendrocytes were 
subclustered on 40 dimensions at a resolution of 0.3.

## Cluster composition and correlation with lesion severity

To assess whether cluster abundance was associated with lesion severity, 
proportions were correlated with the speed of conduction through the spinal 
cord lesioned tissue using simple linear regression. Cell abundance was 
regressed against the speed of conduction, and the strength and significance 
of associations were reported using the coefficient of determination (R²) and 
corresponding p-values. 

Data visualization and statistics were performed in R using the ggplot2 
(v.3.5.2) and ggpmisc (v.0.6.1) package. For each cluster, scatter plots were 
generated with jitter applied to points to reduce overplotting, and regression 
lines with 95% confidence intervals overlaid (e.g. geom_smooth(method = "lm", 
se = TRUE); stat_poly_eq(method=”lm”)).

## Marker gene identification

For marker gene identification, sample layers were joined via JoinLayers(). 
Markers were calculated using FindAllMarkers(min.pct=0.25, logfc.threshold=0.5, 
only.pos=TRUE, test.use=”wilcox”). Only protein-coding genes were considered 
as markers.

## Differential gene expression (DEG) analysis 

For DEG analysis, sample layers were joined via JoinLayers(). DEG was calculated 
using MAST (v.1.32.0) within the FindMarker(min.pct=0.25, logfc.threshold=0.25, 
test.use=”MAST”) function. 

For each cell population, of interest, nuclei were subsetted according to 
cluster identity. Comparisons were made between experimental conditions 
(e.g. MS_slow vs. HD) or transcriptional states (e.g. DAM vs MG_hom). Genes 
with an absolute value of average log2 fold change ≥ 0.25 and a false discovery 
rate (FDR)-adjusted p < 0.05 were considered differentially expressed. 
