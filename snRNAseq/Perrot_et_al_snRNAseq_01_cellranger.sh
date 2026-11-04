#!/bin/bash
#SBATCH --job-name=CellRanger
#SBATCH --cpus-per-task=6
#SBATCH --partition=medium
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=/path_to/snRNAseq_mouse_SC_2025_Perrot_et_al/log/%A_%a.cellranger.log
#SBATCH --array=0-6%4

################################################################################
# Title:        Perrot et al. - snRNAseq analysis
# Subtitle:     CellRanger
# Author:       Mathias Nordbeck (Paris Brain Institute)
# Date:         2024-12-15
# Version:      1.0
# Description:  Aligns FASTQ files generated using the Chromium single cell 3' 
#               v4 (poly A) against the GRCm39-2024-A reference genome and 
#               generates feature-barcode-count matrices using cellranger-9.0.0.
#               This script is configured for single-nucleus RNA sequencing 
#               (snRNA-seq) and therefore includes introns.
################################################################################

set -euo pipefail


################################################################################
# Samples
################################################################################

# list of samples names
SAMPLES=(
NG
HD11
HD12
MS20
MS23
MS24
MS25
)

# select sample corresponding to the SLURM array task
sample=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}


################################################################################
# Directories
################################################################################

# set output directory
DIR_OUT="path_to/snRNAseq_mouse_SC_2025_Perrot_et_al/cellranger_GRCm39"

# set directory of fastq files
DIR_FASTQ="path_to/snRNAseq_mouse_SC_2025_Perrot_et_al/FASTQ"

# choose the reference genome (here: GRCm39-2024-A)
DIR_REF="path_to/Tools/cellranger/refdata-gex-GRCm39-2024-A"

# set directory to cellranger installation
DIR_CELLRANGER="path_to/Tools/cellranger/cellranger-9.0.0/cellranger"


################################################################################
# CellRanger
################################################################################

cd "$DIR_OUT"

echo "CellRanger start." 
echo "Sample: ${sample}" 
echo "Array task: ${SLURM_ARRAY_TASK_ID}"

# run CellRanger count (include introns for nuclei seq)
"$DIR_CELLRANGER" count \
    --id=$sample \
    --sample=$sample \
    --include-introns=true \
    --transcriptome=$DIR_REF \
    --fastqs=$DIR_FASTQ \
    --create-bam=false

echo "CellRanger completed."
