#!/bin/bash

#SBATCH --job-name=pyscenic
#SBATCH --partition=bigmem2
#SBATCH --time=1:00:00
#SBATCH --mem=100GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

module load java
module load singularity

INPUT_DIR=/scRNA_analyses/2_main_analyses/pyscenic



time srun singularity run -B $INPUT_DIR:$INPUT_DIR \
aertslab-pyscenic-0.11.2.sif pyscenic aucell \
CMP_b_humanBM.loom \
scenic_cluster_CMP_b_regulons.csv \
-o scenic_cluster_CMP_b_auc_mtx.csv \
--num_workers 8
