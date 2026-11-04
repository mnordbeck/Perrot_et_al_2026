#!/bin/bash
#SBATCH --job-name=CellBender
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu-ampere
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --output=/path_to/snRNAseq_mouse_SC_2025_Perrot_et_al/log/%A_%a.cellbender.log
#SBATCH --array=0-7%4

################################################################################
# Title:        Perrot et al. - snRNAseq analysis
# Subtitle:     CellBender
# Author:       Mathias Nordbeck (Paris Brain Institute)
# Date:         2025-06-12
# Version:      1.0
# Description:  Uses cellbender-0.3.2 to remove techical artifacs, such as
#               ambient RNA, from feature-barcode-count matrices generated 
#               by CellRanger.
################################################################################

set -euo pipefail


################################################################################
# Samples
################################################################################

# list of samples to process
SAMPLES=(
NG
HD11
HD12
MS20
MS23
MS24
MS25
)

# extract sample of current array
sample=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}


################################################################################
# Directories
################################################################################

DIR_OUTPUT="/path_to/snRNAseq_mouse_SC_2025_Perrot_et_al/cellbender_half_learning_rate/${sample}" 
INPUT_FILE="/path_to/snRNAseq_mouse_SC_2025_Perrot_et_al/cellranger_GRCm39/${sample}/outs/raw_feature_bc_matrix.h5" 
OUTPUT_FILE="${sample}_cellbender.h5"

mkdir -p "$DIR_OUTPUT" 

cd "$DIR_OUTPUT"


################################################################################
# CellBender
################################################################################

# load CUDA
module load CUDA

# activate cellbender environment
source /path_to/Tools/minoforge3_linux/miniforge/bin/activate env_cellbender

echo "CellBender start." 
echo "Sample: ${sample}" 
echo "Array task: ${SLURM_ARRAY_TASK_ID}"

# run cellbender
# NOTE: The learning rate is set to 0.00005, which is half of the default rate
cellbender remove-background \
      --cuda \
      --input "$INPUT_FILE" \
      --output "$OUTPUT_FILE" \
      --learning-rate 0.00005 
      
echo "CellBender completed."
