#!/bin/bash


#SBATCH --partition=broadwl
#SBATCH --job-name=annotate_peaks_background
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40GB


module load bowtie2
module load STAR
module load java
module load samtools
module load fastqc
module load bedtools

PATH=$PATH:/project2/lbarreiro/users/Sarah/ChIP8thru13_combined/homer_analysis/homer_software/bin/
INPUT_PATH=/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/emmreml_downstream/HINT_motif_enrichment_analysis/hint_input_bed_files
OUTPUT_PATH=/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/emmreml_downstream/HOMER_annotate/homer_annotate_output

/project2/lbarreiro/users/Sarah/ChIP8thru13_combined/homer_analysis/homer_software/bin/annotatePeaks.pl $INPUT_PATH/homer_background_peaks.bed hg38 -go $OUTPUT_PATH/background_homer_annotate_GO > $OUTPUT_PATH/background_homer_annotate_out.txt
