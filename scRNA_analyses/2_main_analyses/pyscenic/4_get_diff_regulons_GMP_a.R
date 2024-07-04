library(SingleCellExperiment)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(reticulate)
library(tidyr)
library(reshape2)
library(plyr)
library(grid)
library(gridExtra)
library(data.table)
library(SCENIC)
library(SCopeLoomR)
library(doRNG)
library(ggrepel)

setwd("/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis12_label_transfer_emmreml_edited/pyscenic")
OUT_DIR <- "get_diff_regulons_outs/"
dir.create(OUT_DIR)

##read in filter 20 list to see which samples to exclude
filter_20_list <- readRDS(file=paste0("../label_transfer/filter_20_list.rds"))
names(filter_20_list) <- c("HSC_a", "HSC_b", "CMP_a", "CMP_b", "CMP_c", "GMP_b1", "GMP_b2", "MEP_a1", "MEP_a2", "MEP_a3", "mixed_a", "mixed_b", "PreBNK", "unknown1", "unknown2")


## Get data from sce object:
seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"

#Analysis including all donors
auc_mat <- read.csv(file=paste0("scenic_cluster_GMP_a_auc_mtx.csv"))
df <- seurat_obj@meta.data
samples <- unique(df$UNIQUE_ID)
for(i in 1:length(samples))
{
  temp <- row.names(subset(df, df$UNIQUE_ID %in% samples[i]))
  auc_subset <- subset(auc_mat, auc_mat$Cell %in% temp)
  row.names(auc_subset) <- auc_subset$Cell
  auc_subset <- auc_subset[,2:length(colnames(auc_subset))]
  means <- colMeans(auc_subset)
  if(i ==1)
  {
    donor_data <- means
  }
  else
  {
    donor_data <- cbind(donor_data, means)
  }
}
colnames(donor_data) <- samples  #donor data is a data frame with mean regulon activity for cells of each sample

meta_data <- read.csv(file="../all_samples_meta_data.csv")



Tm3_all <- c("SNG-LB-SS-1S-RS-S1-CD34neg_S1_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S2-CD34neg_S2_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S5-CD34neg_S4_R1_001_Tm3_BCG",
             "SNG-LB-SS-1S-RS-S6-CD34neg_S5_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S8-CD34neg_S7_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S9-CD34neg_S8_R1_001_Tm3_BCG",
             "SNG-LB-SS-1S-RS-S10-CD34neg_S9_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S12-CD34neg_S11_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S14-CD34neg_S13_R1_001_Tm3_BCG",
             "SNG-LB-SS-1S-RS-S15-CD34neg_S14_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S16-CD34neg_S15_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S17-CD34neg_S16_R1_001_Tm3_BCG",
             "SNG-LB-SS-1S-RS-S19-CD34neg_S18_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S20-CD34neg_S19_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S21-CD34neg_S20_R1_001_Tm3_BCG", 
             "SNG-LB-SS-1S-RS-S3-CD34neg_S3_R1_001_Tm3_CTL", "SNG-LB-SS-1S-RS-S7-CD34neg_S6_R1_001_Tm3_CTL", "SNG-LB-SS-1S-RS-S11-CD34neg_S10_R1_001_Tm3_CTL",
             "SNG-LB-SS-1S-RS-S13-CD34neg_S12_R1_001_Tm3_CTL", "SNG-LB-SS-1S-RS-S18-CD34neg_S17_R1_001_Tm3_CTL")


Td0_all <- c("SNG-LB-SS-1S-RS-S1-CD34neg_S1_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S2-CD34neg_S2_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S5-CD34neg_S4_R1_001_Td0_BCG",
             "SNG-LB-SS-1S-RS-S6-CD34neg_S5_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S8-CD34neg_S7_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S9-CD34neg_S8_R1_001_Td0_BCG",
             "SNG-LB-SS-1S-RS-S10-CD34neg_S9_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S12-CD34neg_S11_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S14-CD34neg_S13_R1_001_Td0_BCG",
             "SNG-LB-SS-1S-RS-S15-CD34neg_S14_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S16-CD34neg_S15_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S17-CD34neg_S16_R1_001_Td0_BCG",
             "SNG-LB-SS-1S-RS-S19-CD34neg_S18_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S20-CD34neg_S19_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S21-CD34neg_S20_R1_001_Td0_BCG", 
             "SNG-LB-SS-1S-RS-S3-CD34neg_S3_R1_001_Td0_CTL", "SNG-LB-SS-1S-RS-S7-CD34neg_S6_R1_001_Td0_CTL", "SNG-LB-SS-1S-RS-S11-CD34neg_S10_R1_001_Td0_CTL",
             "SNG-LB-SS-1S-RS-S13-CD34neg_S12_R1_001_Td0_CTL", "SNG-LB-SS-1S-RS-S18-CD34neg_S17_R1_001_Td0_CTL")



logfc_mat <- data.frame()
count <- 0
name_short <- vector()
name_count <- 1
for(i in 1:length(Tm3_all))
{
  name <- Tm3_all[i]
  if(name %in% colnames(donor_data))
  {
    d <- name
    Td0_pair_name <- Td0_all[which(Tm3_all == d)]
    if(Td0_pair_name %in% colnames(donor_data) && is.na(donor_data[,Td0_pair_name])==FALSE)
    {
      Td0_col <- donor_data[,Td0_pair_name]
      Tm3_col <- donor_data[,d]
      log2fc <- Tm3_col - Td0_col
      name_short[name_count] <- subset(meta_data, meta_data$Sample %in% as.character(d))$donor
      name_count <- name_count+1
      if(count==0)
      {
        logfc_mat <- log2fc
      }
      else
      {
        logfc_mat <- cbind(logfc_mat, log2fc)
      }
      
      count <- count + 1
    }
    
  }
  
}

colnames(logfc_mat) <- name_short
row.names(logfc_mat) <- substr(row.names(logfc_mat),1,nchar(row.names(logfc_mat))-3)
#logfc_mat is the difference in mean regulon activity for each donor
write.csv(logfc_mat, file=paste0(OUT_DIR, "GMP_a_logfc_mat_no_outliers.csv"))


pval <- vector(length=length(row.names(logfc_mat)))
pval_wilcox <- vector(length=length(row.names(logfc_mat)))
plot_list <- list()
tf_name <- vector()
plot_count <- 1
for(i in 1:length(row.names(logfc_mat)))
{
  values <- logfc_mat[i, ]
  BCG_donors <- unique(meta_data$donor[which(meta_data$vaccination == "BCG")])
  CTL_donors <- unique(meta_data$donor[which(meta_data$vaccination == "CTL")])
  CTL_vals <- values[which(names(values) %in% CTL_donors)]
  BCG_vals <- values[which(names(values) %in% BCG_donors)]
  
  res <- t.test(BCG_vals, CTL_vals)
  pval[i] <- res$p.value
  
  res_wilcox <- wilcox.test(BCG_vals, CTL_vals)
  pval_wilcox[i] <- res_wilcox$p.value
  
}
padj_wilcox <- p.adjust(pval_wilcox, method='BH')
wilcox_data <- data.frame(cbind(pval_wilcox, padj_wilcox))
row.names(wilcox_data) <- row.names(logfc_mat)
write.csv(wilcox_data, file=paste0(OUT_DIR, "GMP_a_wilcox_data_no_outliers.csv")) ##RENAME GMP c1




#plots for specific regulons 

regulons <- c("EGR1","EGR3", "ERG", "KLF5", "KLF6", "KLF10", "ELF2", "ELK3", "ERF", "ETV6", "KLF3")
colors_list <- c("darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "#E15759", "darkgrey", "darkgrey", "darkgrey")
plot_list <- list()

for(i in 1:length(regulons))
{
  if(regulons[i] %in% row.names(logfc_mat))
  {
    values <- logfc_mat[regulons[i], ]
    BCG_donors <- unique(meta_data$donor[which(meta_data$vaccination == "BCG")])
    CTL_donors <- unique(meta_data$donor[which(meta_data$vaccination == "CTL")])
    CTL_vals <- values[which(names(values) %in% CTL_donors)]
    BCG_vals <- values[which(names(values) %in% BCG_donors)]
    
    res <- wilcox.test(BCG_vals, CTL_vals)
    pval <- res$p.value
    
    donors <- c(names(values)[which(names(values) %in% BCG_donors)], names(values)[which(names(values) %in% CTL_donors)])
    lab <- rep("no",length(donors))
    lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
    
    df_plot <- data.frame(
      vals <- c(BCG_vals, CTL_vals),
      group <- c(rep("BCG", length(BCG_vals)), rep("CTL", length(CTL_vals))),
      d <- donors,
      l <- lab
    )
    colnames(df_plot) <- c("vals", "group", "d", "l")
    
    #Make plot where values are shifter relative to controls
    subset <- subset(df_plot, df_plot$group == "CTL")
    subset_median <- median(subset$vals)
    
    new_vals <- df_plot$vals - subset_median
    df_plot_new <- data.frame(
      val_corrected <- new_vals,
      group <- df_plot$group,
      d <- df_plot$d,
      l <- df_plot$l
    )
    colnames(df_plot_new) <- c("vals", "group", "d", "l")
    
    p <- ggplot(data=df_plot_new, aes(x=group, y=vals, fill=group, label=d)) + geom_boxplot(aes(fill=group), alpha=0.8)  + geom_jitter(shape=21, size=2, stroke=1, alpha=1, position=position_jitter(0.2), aes(fill=group))+ theme_classic()+
      labs(x="", y="Normalized TF Activity change", fill="Group", title=paste0(regulons[i]," in GMPa\npval = ", round(pval, digits=2))) + scale_fill_manual(values=c(colors_list[i],"lightgrey")) +
      theme(plot.title=element_text(size=10)) +  theme(legend.position = "none") + theme(axis.title=element_text (size=9))
    
    plot_list[[i]] <- p
  }
  
}


for(k in 1:length(plot_list)) 
{
  file_name = paste0(OUT_DIR, regulons[k], "_regulon_GMP_a.tiff")
  tiff(file_name, units="in", width=2, height=2.5, res=300)
  print(plot_list[[k]])
  dev.off()
  print(k)
  
}




#plots for specific regulons broad TF group names 

regulons <- c( "KLF3", "KLF10", "ELK3")
broad_names <- c( "KLF2/3/6", "KLF/SP_2", "ETS_3")
colors_list <- c("darkgrey", "darkgrey",  "#E15759")
plot_list <- list()

for(i in 1:length(regulons))
{
  values <- logfc_mat[regulons[i], ]
  BCG_donors <- unique(meta_data$donor[which(meta_data$vaccination == "BCG")])
  CTL_donors <- unique(meta_data$donor[which(meta_data$vaccination == "CTL")])
  CTL_vals <- values[which(names(values) %in% CTL_donors)]
  BCG_vals <- values[which(names(values) %in% BCG_donors)]
  
  res <- wilcox.test(BCG_vals, CTL_vals)
  pval <- res$p.value
  
  donors <- c(names(values)[which(names(values) %in% BCG_donors)], names(values)[which(names(values) %in% CTL_donors)])
  lab <- rep("no",length(donors))
  lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
  
  df_plot <- data.frame(
    vals <- c(BCG_vals, CTL_vals),
    group <- c(rep("BCG", length(BCG_vals)), rep("CTL", length(CTL_vals))),
    d <- donors,
    l <- lab
  )
  colnames(df_plot) <- c("vals", "group", "d", "l")
  
  #Make plot where values are shifter relative to controls
  subset <- subset(df_plot, df_plot$group == "CTL")
  subset_median <- median(subset$vals)
  
  new_vals <- df_plot$vals - subset_median
  df_plot_new <- data.frame(
    val_corrected <- new_vals,
    group <- df_plot$group,
    d <- df_plot$d,
    l <- df_plot$l
  )
  colnames(df_plot_new) <- c("vals", "group", "d", "l")
  
  p <- ggplot(data=df_plot_new, aes(x=group, y=vals, fill=group, label=d)) + geom_boxplot(aes(fill=group), alpha=0.8)  + geom_jitter(shape=21, size=2, stroke=1, alpha=1, position=position_jitter(0.2), aes(fill=group))+ theme_classic()+
    labs(x="", y="Normalized TF Activity change", fill="Group", title=paste0(broad_names[i]," in GMP\npval = ", round(pval, digits=2))) + scale_fill_manual(values=c(colors_list[i],"lightgrey")) +
    theme(plot.title=element_text(size=12)) +  theme(legend.position = "none") + theme(axis.title=element_text (size=10))
  
  plot_list[[i]] <- p
}


for(k in 1:length(plot_list)) 
{
  file_name = paste0(OUT_DIR, regulons[k], "_broad_name_GMP_a_excluded.tiff")
  tiff(file_name, units="in", width=2, height=2.5, res=300)
  print(plot_list[[k]])
  dev.off()
  print(k)
  
}




