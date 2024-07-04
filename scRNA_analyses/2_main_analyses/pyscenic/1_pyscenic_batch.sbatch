#!/bin/bash

#SBATCH --job-name=pyscenic
#SBATCH --partition=bigmem2
#SBATCH --time=18:00:00
#SBATCH --mem=200GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

module load java
module load singularity

INPUT_DIR=/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis12_label_transfer_emmreml_edited/pyscenic

time srun singularity run -B $INPUT_DIR:$INPUT_DIR aertslab-pyscenic-0.11.2.sif pyscenic grn --num_workers 12 -o $2 $1 hs_hgnc_tfs.txt
