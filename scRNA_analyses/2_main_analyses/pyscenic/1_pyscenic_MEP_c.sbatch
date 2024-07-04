#!/bin/bash

#SBATCH --job-name=pyscenic
#SBATCH --partition=bigmem2
#SBATCH --time=4:00:00
#SBATCH --mem=200GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

module load java
module load singularity

INPUT_DIR=/scRNA_analyses/2_main_analyses/pyscenic

singularity run -B $INPUT_DIR:$INPUT_DIR aertslab-pyscenic-0.11.2.sif pyscenic grn --num_workers 12 -o scenic_cluster_MEP_c.tsv MEP_c_humanBM.loom hs_hgnc_tfs.txt
