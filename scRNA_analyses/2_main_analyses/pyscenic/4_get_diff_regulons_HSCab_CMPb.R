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


## Get data from sce object:
seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"

seurat_obj<- subset(
  x = seurat_obj,
  subset = timepoint %in% c("Td0")
)

#Analysis including all donors
auc_mat <- read.csv(file=paste0("scenic_cluster_HSCab_CMPb_auc_mtx.csv"))
df <- seurat_obj@meta.data

celltypes <- c("HSC_a", "HSC_b", "CMP_b")


HSC_obj<- subset(
  x = seurat_obj,
  subset = clust_names %in% c("HSC_a", "HSC_b")
)
HSC_df <- HSC_obj@meta.data

CMP_obj<- subset(
  x = seurat_obj,
  subset = clust_names %in% c("CMP_b")
)
CMP_df <- CMP_obj@meta.data


samples_HSC <- unique(HSC_df$UNIQUE_ID)
for(i in 1:length(samples_HSC))
{
  temp <- row.names(subset(HSC_df, HSC_df$UNIQUE_ID %in% samples_HSC[i]))
  auc_subset <- subset(auc_mat, auc_mat$Cell %in% temp)
  row.names(auc_subset) <- auc_subset$Cell
  auc_subset <- auc_subset[,2:length(colnames(auc_subset))]
  means <- colMeans(auc_subset)
  if(i ==1)
  {
    donor_data_HSC <- means
  }
  else
  {
    donor_data_HSC <- cbind(donor_data_HSC, means)
  }
}
colnames(donor_data_HSC) <- samples_HSC  #donor data is a data frame with mean regulon activity for cells of each sample


samples_CMP <- unique(CMP_df$UNIQUE_ID)
for(i in 1:length(samples_CMP))
{
  temp <- row.names(subset(CMP_df, CMP_df$UNIQUE_ID %in% samples_CMP[i]))
  auc_subset <- subset(auc_mat, auc_mat$Cell %in% temp)
  row.names(auc_subset) <- auc_subset$Cell
  auc_subset <- auc_subset[,2:length(colnames(auc_subset))]
  means <- colMeans(auc_subset)
  if(i ==1)
  {
    donor_data_CMP <- means
  }
  else
  {
    donor_data_CMP <- cbind(donor_data_CMP, means)
  }
}
colnames(donor_data_CMP) <- samples_CMP  #donor data is a data frame with mean regulon activity for cells of each sample


which(row.names(donor_data_HSC) != row.names(donor_data_CMP))


pval_wilcox <- vector(length=length(row.names(donor_data_HSC)))
plot_list <- list()
tf_name <- vector()
plot_count <- 1
direction <- vector(length=length(row.names(donor_data_HSC)))
for(i in 1:length(row.names(donor_data_HSC)))
{
  HSC_vals <- donor_data_HSC[i, ]
  CMP_vals <- donor_data_CMP[i, ]

  direction[i] <- mean(CMP_vals) - mean(HSC_vals)
  
  res_wilcox <- wilcox.test(HSC_vals, CMP_vals)
  pval_wilcox[i] <- res_wilcox$p.value
  
  if(pval_wilcox[i] < 0.00000005)
  {
    df_plot <- data.frame(
      vals <- c(HSC_vals, CMP_vals),
      group <- c(rep("HSC", length(HSC_vals)), rep("CMP", length(CMP_vals)))
    )
    colnames(df_plot) <- c("vals", "group")
    
    p <- ggplot(data=df_plot, aes(x=group, y=vals, fill=group)) + geom_boxplot(aes(fill=group), alpha=0.5)  + geom_jitter(shape=21, size=3, color="black",alpha=1, position=position_jitter(0.2), aes(fill=group))+ theme_light()+
      labs(x="", y="Regulon Activity (HSC vs CMP)", fill="Group", title=paste0(row.names(donor_data_HSC)[i]," Regulon: pval = ", round(pval_wilcox[i], digits=2))) + scale_fill_manual(values=c("#DD4477","#0099C6")) +
      theme(axis.title=element_text(size=10))
    
    
    plot_list[[plot_count]] <- p
    tf_name[plot_count] <- row.names(donor_data_HSC)[i]
    plot_count <- plot_count + 1
    
    
  }
  
}
padj_wilcox <- p.adjust(pval_wilcox, method='BH')
wilcox_data <- data.frame(cbind(pval_wilcox, padj_wilcox, direction))
row.names(wilcox_data) <- row.names(donor_data_HSC)
write.csv(wilcox_data, file=paste0(OUT_DIR, "HSCvsCMP_wilcox_data_no_outliers.csv"))


# for(k in 1:length(plot_list)) 
# {
#   file_name = paste0(OUT_DIR, tf_name[k], "_regulon_HSCvsCMP.tiff")
#   tiff(file_name, units="in", width=3.5, height=3, res=250)
#   print(plot_list[[k]])
#   dev.off()
#   print(k)
#   
# }
