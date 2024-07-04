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

INPUT_DIR=/scRNA_analyses/2_main_analyses/pyscenic

time srun singularity run -B $INPUT_DIR:$INPUT_DIR aertslab-pyscenic-0.11.2.sif pyscenic ctx scenic_cluster_MEP_a.tsv \
hg19-tss-centered-10kb-7species.mc9nr.feather \
hg19-500bp-upstream-7species.mc9nr.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname MEP_a_humanBM.loom \
--mode "dask_multiprocessing" \
--output scenic_cluster_MEP_a_regulons.csv \
--num_workers 9
