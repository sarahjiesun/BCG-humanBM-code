##correlate alphas and plot

install.packages("devtools") 
devtools::install_github("jakelawlor/PNWColors") 
library(ggplot2)
library(PNWColors)
library(corrplot)
library(RColorBrewer)
library(dplyr)
library(ggrepel)

## Set directory structure
setwd("/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/0.025_elastic_net_model/")
OUT_DIR <- "corr_plots_out/"
dir.create(OUT_DIR)

#my_colors <- pnw_palette(name="Starfish",n=7,type="discrete")
my_colors_Sarah <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1")

## only load selected models
models <- c("IL10_FC_score", "IL1B_FC_score", "IL1RA_FC_score", "IL6_FC_score", "IFNA_FC_score", "IFNG_FC_score", "TNF_FC_score")


alpha=c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

cluster_names <- c("CMP1", "GMP1", "CMP3",  "CMP4", "CMP5",  "GMP3", "GMP2")
plot_list <- list()

for(j in 1:length(cluster_names)){
  
  
  in_dir <- paste0("/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/0.025_elastic_net_model/",cluster_names[j],"_EN_outputs/")
  
  model_dfs <- list()
  data_models_list <-list()
  
  
  ## load metadata
  metadata <- read.table(paste0(cluster_names[j],"_cyt_FC_final.txt"))
  
  for (m in 1:length(models)){
    print(models[m])
    model_m <- models[m]
    
    #get files from folder and add predicted values to a data frame as columns
    all_alphas <- list.files(paste0(in_dir), pattern = paste0(models[m], "_predicted"), full.names = T)
    
    model_dfs[[m]] <- as.data.frame(do.call(cbind,lapply(all_alphas,function(fn)read.table(fn,header=T, sep=",")[,3])))
    alpha_labels <- as.vector(list.files(paste0(in_dir), pattern = paste0(models[m], "_predicted"), full.names = F))    
    
    alpha_labels<- gsub('.{4}$', '', alpha_labels)
    alpha_labels <- sub('.*\\_', '', alpha_labels)
    colnames(model_dfs[[m]])  <- alpha_labels
    
    ## try to improve this code to reduce error
    rownames(model_dfs[[m]]) <- metadata$Donor
    
    #add real values from metadata file
    if (model_m=="IL10_FC_score"){
      model_dfs[[m]]$real_value <- metadata$IL10_FC   
    } else if (model_m=="IL1B_FC_score"){
      model_dfs[[m]]$real_value <- metadata$IL1B_FC
    } else if (model_m=="IL1RA_FC_score"){
      model_dfs[[m]]$real_value <- metadata$IL1RA_FC
    } else if (model_m=="IL6_FC_score"){
      model_dfs[[m]]$real_value <- metadata$IL6_FC
    } else if (model_m=="IFNA_FC_score"){
      model_dfs[[m]]$real_value <- metadata$IFNA_FC
    } else if (model_m=="IFNG_FC_score"){
      model_dfs[[m]]$real_value <- metadata$IFNG_FC
    } else if (model_m=="TNF_FC_score"){
      model_dfs[[m]]$real_value <- metadata$TNF_FC
    }
  }
  names(model_dfs) <- models
  data_models_list <- model_dfs
  
  
  ######################################################################
  ## Correlate actual value with predicted value and select top alpha ##
  ######################################################################
  
  
  for (m in 1:length(models)){
    corr_df <- data_models_list[[m]]
    
    correlation_results <- data.frame(R=numeric(ncol(corr_df)-1), P.value=numeric(ncol(corr_df)-1))
    
    for (i in 1:(ncol(corr_df)-1)){
      test <- cor.test(corr_df[,"real_value"], corr_df[,i], method = c("spearman"))
      correlation_results$R[i] = test$estimate
      correlation_results$P.value[i] = test$p.value
      
      # if(colnames(corr_df)=="a0.3")
      # {
      #   meta_data <- read.csv(file="../meta_data_limma_betas_interaction.csv")
      #   treatment <- vector(length=length(row.names(corr_df)))
      #   donor_name <- vector(length=length(row.names(corr_df)))
      #   for(d in 1:length(row.names(corr_df)))
      #   {
      #     #GET BCG/CTL info for each donor 
      #     treatment[d] <- as.character(unique(meta_data$group[which(meta_data$donor == row.names(corr_df)[d])]))
      #     donor_name[d] <- as.character(unique(meta_data$donor[which(meta_data$donor == row.names(corr_df)[d])]))
      #   }
      #   df_scatter_individual <- data.frame(
      #     real <- corr_df[,"real_value"],
      #     pred <- corr_df[,i],
      #     g <- treatment,
      #     dn <- donor_name
      #   )
      #   
      #   ann_text<-data.frame(
      #     x = 0.9, y = 1.4,
      #     label = paste0("Spearman corr = ",round(correlation_results$R[i],3),"\npvalue = ",round(correlation_results$P.value[i],4))
      #   )
      #   
      #   p1 <- ggplot(df_scatter_individual, aes(x=real, y=pred)) + geom_point(shape=21, size=2, aes(fill=g)) + theme_classic()+
      #     geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1) + labs(x='Real IL1B FC', y="Predicted IL1B FC")+
      #     labs(fill="Group")
      #   p2 <- ggplot(df_scatter_individual, aes(x=real, y=pred)) + geom_point(shape=21, size=3, fill="blue", alpha=0.7) + theme_classic()+
      #     geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1) + labs(x='Real IL1B FC', y="Predicted IL1B FC")+
      #     labs(fill="Group")+ geom_text(data=ann_text, aes(x=0.9, y=1.4, label=label))
      #   
      #   tiff(paste0(OUT_DIR, "ATAC_DApeaks_il1b_GMP3_a0.3_scatter_nocol.tiff"), units="in", width=3.5, height=3, res=250)
      #   p2
      #   dev.off()
      #   
      #   
      # }
    }
    
    correlation_results$alpha <- colnames(corr_df[1:(ncol(corr_df)-1)])
    correlation_results$model <- models[m]
    
    if(m==1){
      model_results <- correlation_results
    } else {
      model_results <- rbind(model_results, correlation_results)
    }
  }
  
  model_results$datatype <- paste0("ATAC_", cluster_names[j])
  
  
  all_results <- model_results
  
  all_results_no_negs <- all_results
  
  ## force all negative numbers to be 0
  all_results_no_negs$R[all_results_no_negs$R<0] <- 0
  ## remove a from alpha labels
  all_results_no_negs$alpha <- gsub("a", "", all_results_no_negs$alpha)
  
  #add lab for p value
  p_group <- rep("none", length(row.names(all_results_no_negs)))
  p_group[which(all_results_no_negs$P.value < 0.05)] <- "sig"
  p_group[which(all_results_no_negs$R == 0)] <- "none"
  all_results_no_negs$p_group <- p_group
  
  lab <- rep("",length(p_group))
  lab[which(p_group=="sig")] <- paste0("p = ", as.character(round(all_results_no_negs$P.value[which(p_group=="sig")],3)))
  all_results_no_negs$l <- lab
  
  ## plot
  p <- ggplot(all_results_no_negs, aes(x=alpha, y = R, group=model, color=p_group, fill=model, label=l))+
    geom_point(alpha = 1, position=position_dodge(width=.3), cex=3.4, shape=21, stroke=1) +   scale_fill_manual(values=my_colors_Sarah) +
    xlab("alpha") + ylab(expression(paste("Spearman's ", italic("ρ"))))+ labs(title=paste0("ATAC ",cluster_names[j]))+
    theme_classic() + scale_color_manual(values=c("white", "black"))+ geom_text_repel(aes(label=l))+
    theme(legend.position = "bottom", text=element_text(size=14), axis.text.x = element_text(size=10, angle=25), 
          legend.title = element_text(size=6), legend.text = element_text(size = 5)) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .2)) + guides(col=FALSE)
  
  plot_list[[j]] <- p
  
  print(j)
}

for(k in 1:length(plot_list)) 
{
  tiff(paste0(OUT_DIR, "ATAC_DApeaks_",cluster_names[k],"_EN.tiff"), units="in", width=4, height=4, res=250)
  print(plot_list[[k]])
  dev.off()
  print(k)
  
}




## plot il1b CMP5_a0.7

corr_df <- data_models_list[[2]]  ### IL1B
correlation_results <- data.frame(R=numeric(ncol(corr_df)-1), P.value=numeric(ncol(corr_df)-1))
test <- cor.test(corr_df[,"real_value"], corr_df[,7], method = c("spearman")) ## a0.7
correlation_results$R[7] = test$estimate ### a0.7
correlation_results$P.value[7] = test$p.value ### a0.7

meta_data <- read.csv(file="../../analysis_FINAL_Clusters_harmony/meta_data_limma_betas_interaction.csv")
treatment <- vector(length=length(row.names(corr_df)))
donor_name <- vector(length=length(row.names(corr_df)))
for(d in 1:length(row.names(corr_df))){
  #GET BCG/CTL info for each donor 
  treatment[d] <- as.character(unique(meta_data$group[which(meta_data$donor == row.names(corr_df)[d])]))
  donor_name[d] <- as.character(unique(meta_data$donor[which(meta_data$donor == row.names(corr_df)[d])]))
}
df_scatter_individual <- data.frame(
  real <- corr_df[,"real_value"],
  pred <- corr_df[,7],
  g <- treatment,
  dn <- donor_name
)

ann_text<-data.frame(
  x = 0.9, y = 1.8,
  label = paste0("Spearman corr = ",round(correlation_results$R[7],3),"\npvalue = ",round(correlation_results$P.value[7],4))
)

p2 <- ggplot(df_scatter_individual, aes(x=real, y=pred)) + geom_point(shape=21, size=4, fill="red", alpha=0.7, color="grey") + theme_classic()+
  geom_smooth(method=lm, color="grey", alpha=0.2, fullrange=TRUE, size=1) + labs(x='Real IL1B FC', y="Predicted IL1B FC")+
  labs(fill="Group")+ geom_text(data=ann_text, aes(x=1, y=1.65, label=label))

tiff(paste0(OUT_DIR, "ATAC_DApeaks_il1b_CMP5_a0.7_scatter_nocol.tiff"), units="in", width=3.5, height=3, res=250)
p2
dev.off()



## plot il6 CMP5_a0

corr_df <- data_models_list[[4]]  ### IL6
correlation_results <- data.frame(R=numeric(ncol(corr_df)-1), P.value=numeric(ncol(corr_df)-1))
test <- cor.test(corr_df[,"real_value"], corr_df[,10], method = c("spearman")) ## a0
correlation_results$R[10] = test$estimate ### a0
correlation_results$P.value[10] = test$p.value ### a0

meta_data <- read.csv(file="../../analysis_FINAL_Clusters_harmony/meta_data_limma_betas_interaction.csv")
treatment <- vector(length=length(row.names(corr_df)))
donor_name <- vector(length=length(row.names(corr_df)))
for(d in 1:length(row.names(corr_df))){
  #GET BCG/CTL info for each donor 
  treatment[d] <- as.character(unique(meta_data$group[which(meta_data$donor == row.names(corr_df)[d])]))
  donor_name[d] <- as.character(unique(meta_data$donor[which(meta_data$donor == row.names(corr_df)[d])]))
}
df_scatter_individual <- data.frame(
  real <- corr_df[,"real_value"],
  pred <- corr_df[,'a0'],
  g <- treatment,
  dn <- donor_name
)

ann_text<-data.frame(
  x = 0.9, y = 1.8,
  label = paste0("Spearman corr = ",round(correlation_results$R[10],3),"\npvalue = ",round(correlation_results$P.value[10],4))
)

p2 <- ggplot(df_scatter_individual, aes(x=real, y=pred)) + geom_point(shape=21, size=4, fill="red", alpha=0.7, color="grey") + theme_classic()+
  geom_smooth(method=lm, color="grey", alpha=0.2, fullrange=TRUE, size=1) + labs(x='Real IL6 FC', y="Predicted IL6 FC")+
  labs(fill="Group")+ geom_text(data=ann_text, aes(x=1, y=1.65, label=label))

tiff(paste0(OUT_DIR, "ATAC_DApeaks_il6_CMP5_a0_scatter_nocol.tiff"), units="in", width=3.5, height=3, res=250)
p2
dev.off()
