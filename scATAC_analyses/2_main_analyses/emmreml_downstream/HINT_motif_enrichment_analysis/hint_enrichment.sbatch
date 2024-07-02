#!/bin/bash

#SBATCH --partition=broadwl
#SBATCH --job-name=motif_matching
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB


module load bowtie2
module load STAR
module load java
module load samtools
module load fastqc
module load bedtools
module unload python
module load python/cpython-3.7.0

INPUT_PATH=/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis1_Clusters_harmony/scHINT
OUTPUT_PATH=/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJET/BM_CD34_scATAC/Rprojects/ArchR/analysis1_Clusters_harmony/scHINT/Footprints
PEAK_PATH=/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis1_Clusters_harmony/emmreml_Jedited_downstream/get_DA_peaks

rgt-motifanalysis enrichment --organism=hg38 --output-location=./hint_Enrichment ./hint_input_bed_files/homer_background_peaks.bed ./hint_input_bed_files/$1