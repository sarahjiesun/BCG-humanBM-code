#!/bin/bash

#SBATCH --partition=broadwl
#SBATCH --job-name=overlap
#SBATCH --time=0:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB


module load java
module load samtools
module load bedtools  

INPUT_PATH1=/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/emmreml_downstream/HINT_motif_enrichment_analysis/hint_input_bed_files 
INPUT_PATH2=/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/emmreml_downstream/HINT_motif_enrichment_analysis/bedtools_input
OUTPUT_PATH=/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/emmreml_downstream/HINT_motif_enrichment_analysis/bedtools_output/peaks_with_motifs

bedtools intersect -wo -a $INPUT_PATH1/$1 -b $INPUT_PATH2/$2 > $OUTPUT_PATH/$3