#!/bin/bash

#SBATCH --job-name=pyscenic
#SBATCH --partition=bigmem2
#SBATCH --time=00:30:00
#SBATCH --mem=100GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

module load java
module load singularity

INPUT_DIR=/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis12_label_transfer_emmreml_edited/pyscenic



time srun singularity run -B $INPUT_DIR:$INPUT_DIR \
aertslab-pyscenic-0.11.2.sif pyscenic aucell \
GMP_b_humanBM.loom \
scenic_cluster_GMP_b_regulons.csv \
-o scenic_cluster_GMP_b_auc_mtx.csv \
--num_workers 8
