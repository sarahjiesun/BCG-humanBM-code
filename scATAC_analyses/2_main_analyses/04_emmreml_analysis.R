#R version 4.1.0
#install.packages("EMMREML")
library(ggplot2)
library(statmod)
library(RColorBrewer)
library(edgeR)
library(devtools)
library(cluster)
library(ggrepel)
library(EMMREML)
library(qvalue)
library(dplyr)


# setup ------------------------------------------------------------
setwd("/scATAC_analyses/2_main_analyses")
OUT_DIR <- paste0("0X_analysis_test_filter_thresholds/")
dir.create(OUT_DIR)
IN_DIR <- paste0("0X_test_filter_thresholds/res_full/")


clusters <- c("C1", "C3", "C4", "C5", "C6", "C7","C9", "C10", "C12", "C14",  "C17",  "C21", "C22", "C24","C23", "C15", "C18")
clusters_renamed <- c("unknown1", "CMP1", "CMP2", "HSC1", "HSC2", "MEP1","Pre-BNK", "GMP1","CLP", "MEP2", "CMP3",  "CMP4", "CMP5",  "GMP3", "GMP2", "MEP3", "MEP4")


for(i in 1:length(clusters)){
  option <- c(1, 1.25, 1.5, 1.75, 2, 2.25, 2.5)
  num_original_peaks <- vector(length=length(option))
  num_original_DA_peaks <- vector(length=length(option))
  num_DA_peaks <- vector(length=length(option))
  num_peaks <- vector(length=length(option))
  overlap_DA_peaks <- vector(length=length(option))
  overlap_DA_peaks_percent <- vector(length=length(option))
  overlap_tested_peaks <- vector(length=length(option))
  overlap_tested_peaks_percent <- vector(length=length(option))
  overlap_origDA_peaks_tested <- vector(length=length(option))
  
  original_data <- read.csv(file=paste0("emmreml_betas_Jedited/",clusters_renamed[i],"_res_full.csv"))
  original_sig_peaks <- original_data$X[which(original_data$BH_corrected < 0.1)]
  original_tested_peaks <- original_data$X
  for(k in 1:length(option)){
    test_data <- read.csv(file=paste0(IN_DIR,clusters[i],"_res_full_option",option[k],".csv"))
    test_sig_peaks <- test_data$X[which(test_data$BH_corrected < 0.1)]
    test_tested_peaks <- test_data$X
    
    num_original_peaks[k] <- length(original_tested_peaks)
    num_original_DA_peaks[k] <- length(original_sig_peaks)
    num_DA_peaks[k] <- length(test_sig_peaks)
    num_peaks[k] <- length(test_tested_peaks)
    overlap_DA_peaks[k] <- length(intersect(original_sig_peaks, test_sig_peaks))
    overlap_DA_peaks_percent[k] <- overlap_DA_peaks[k]/length(original_sig_peaks)
    overlap_tested_peaks[k] <- length(intersect(original_tested_peaks, test_tested_peaks))
    overlap_tested_peaks_percent[k] <- overlap_tested_peaks[k]/length(original_tested_peaks)
    overlap_origDA_peaks_tested[k] <- length(intersect(original_sig_peaks, test_tested_peaks))/length(original_sig_peaks)
    
    
  }
  
  df <- cbind(option,num_original_peaks, num_peaks, num_original_DA_peaks, num_DA_peaks, overlap_DA_peaks, overlap_DA_peaks_percent,
                overlap_tested_peaks, overlap_tested_peaks_percent, overlap_origDA_peaks_tested)
  colnames(df) <- c('threshold', "#original_peaks", "#new_peaks", "#original_DA", "#new_DA", "#overlap_DA", "%overlap_DA", "#overlap_tested", "%overlap_tested", "%original_DA_tested")
  write.csv(df,file=paste0(OUT_DIR, clusters_renamed[i],"_data.csv"))
}


## GMP3 p-values for new DA peaks, original DA peaks

GMP3_original_data <- read.csv(file=paste0("../analysis1_Clusters_harmony/emmreml_betas_Jedited/GMP3_res_full.csv"))
GMP3_new_data <- read.csv(file=paste0(IN_DIR,"C24_res_full_option1.csv"))

GMP3_original_pvals <- GMP3_original_data$p_value_vaccinationBCG.timepoint[which(GMP3_original_data$BH_corrected < 0.1)]
GMP3_new_pvals <- GMP3_new_data$p_value_groupBCG.timepoint[which(GMP3_new_data$BH_corrected < 0.1)]

overlap_peaks <- intersect(GMP3_original_data$X[which(GMP3_original_data$BH_corrected < 0.1)], GMP3_new_data$X[which(GMP3_new_data$BH_corrected < 0.1)])
overlap_pvals <- GMP3_new_data$p_value_groupBCG.timepoint[which(GMP3_new_data$X %in% overlap_peaks)]
nonoverlap_peaks <- setdiff(GMP3_new_data$X[which(GMP3_new_data$BH_corrected < 0.1)],GMP3_original_data$X[which(GMP3_original_data$BH_corrected < 0.1)] )
nonoverlap_pvals <- GMP3_new_data$p_value_groupBCG.timepoint[which(GMP3_new_data$X %in% nonoverlap_peaks)]

df <- data.frame(
  group <- c(rep("original", length(GMP3_original_pvals)), rep("new", length(GMP3_new_pvals)), rep("overlap", length(overlap_pvals)), rep("nonoverlap", length(nonoverlap_pvals))),
  vals <- c(GMP3_original_pvals, GMP3_new_pvals, overlap_pvals, nonoverlap_pvals)
)
colnames(df) <- c("group", "vals")

file_name = paste0(OUT_DIR,"GMP3_pvals_comparison_boxplot.tiff")
tiff(file_name)
print(ggplot(df, aes(x=group, y=vals)) + geom_boxplot(aes(fill=group)) + geom_jitter()) 
dev.off()



## GMP3 p-adj for new DA peaks, original DA peaks

GMP3_original_data <- read.csv(file=paste0("../analysis1_Clusters_harmony/emmreml_betas_Jedited/GMP3_res_full.csv"))
GMP3_new_data <- read.csv(file=paste0(IN_DIR,"C24_res_full_option1.csv"))

GMP3_original_pvals <- GMP3_original_data$BH_corrected[which(GMP3_original_data$BH_corrected < 0.1)]
GMP3_new_pvals <- GMP3_new_data$BH_corrected[which(GMP3_new_data$BH_corrected < 0.1)]

overlap_peaks <- intersect(GMP3_original_data$X[which(GMP3_original_data$BH_corrected < 0.1)], GMP3_new_data$X[which(GMP3_new_data$BH_corrected < 0.1)])
write.table(overlap_peaks, file=paste0(OUT_DIR,"GMP3_1_overlap_peaks.txt"))
overlap_pvals <- GMP3_new_data$BH_corrected[which(GMP3_new_data$X %in% overlap_peaks)]
nonoverlap_peaks <- setdiff(GMP3_new_data$X[which(GMP3_new_data$BH_corrected < 0.1)],GMP3_original_data$X[which(GMP3_original_data$BH_corrected < 0.1)] )
nonoverlap_pvals <- GMP3_new_data$BH_corrected[which(GMP3_new_data$X %in% nonoverlap_peaks)]

df <- data.frame(
  group <- c(rep("original", length(GMP3_original_pvals)), rep("new", length(GMP3_new_pvals)), rep("overlap", length(overlap_pvals)), rep("nonoverlap", length(nonoverlap_pvals))),
  vals <- c(GMP3_original_pvals, GMP3_new_pvals, overlap_pvals, nonoverlap_pvals)
)
colnames(df) <- c("group", "vals")

file_name = paste0(OUT_DIR,"GMP3_padj_comparison_boxplot.tiff")
tiff(file_name)
print(ggplot(df, aes(x=group, y=vals)) + geom_boxplot(aes(fill=group))) + geom_jitter()
dev.off()



## GMP3 betas for new DA peaks, original DA peaks

GMP3_original_data <- read.csv(file=paste0("../analysis1_Clusters_harmony/emmreml_betas_Jedited/GMP3_res_full.csv"))
GMP3_new_data <- read.csv(file=paste0(IN_DIR,"C24_res_full_option1.csv"))

GMP3_original_betas <- GMP3_original_data$beta_vaccinationBCG.timepoint[which(GMP3_original_data$BH_corrected < 0.1)]
GMP3_new_betas <- GMP3_new_data$beta_groupBCG.timepoint[which(GMP3_new_data$BH_corrected < 0.1)]

overlap_peaks <- intersect(GMP3_original_data$X[which(GMP3_original_data$BH_corrected < 0.1)], GMP3_new_data$X[which(GMP3_new_data$BH_corrected < 0.1)])
overlap_betas <- GMP3_new_data$beta_groupBCG.timepoint[which(GMP3_new_data$X %in% overlap_peaks)]

df <- data.frame(
  group <- c(rep("original", length(GMP3_original_betas)), rep("new", length(GMP3_new_betas)), rep("overlap", length(overlap_betas))),
  vals <- c(GMP3_original_betas, GMP3_new_betas, overlap_betas)
)
colnames(df) <- c("group", "vals")

file_name = paste0(OUT_DIR,"GMP3_betas_comparison_boxplot.tiff")
tiff(file_name)
print(ggplot(df, aes(x=group, y=vals)) + geom_boxplot(aes(fill=group))) + geom_jitter()
dev.off()
