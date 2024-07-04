#!/bin/bash

#SBATCH --job-name=pyscenic
#SBATCH --partition=bigmem2
#SBATCH --time=2:00:00
#SBATCH --mem=100GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9

module load java
module load singularity

INPUT_DIR=/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis12_label_transfer_emmreml_edited/pyscenic

time srun singularity run -B $INPUT_DIR:$INPUT_DIR aertslab-pyscenic-0.11.2.sif pyscenic ctx scenic_cluster_GMP_a.tsv \
hg19-tss-centered-10kb-7species.mc9nr.feather \
hg19-500bp-upstream-7species.mc9nr.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname GMP_a_humanBM.loom \
--mode "dask_multiprocessing" \
--output scenic_cluster_GMP_a_regulons.csv \
--num_workers 9