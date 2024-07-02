#R version 4.1.0 loaded
#rstudio/1.4
#geos/3.9.1 loaded
#hdf5/1.12.0 loaded 
#cmake loaded

.libPaths("/project/lbarreiro/USERS/sarah/Rlibs_new")

library(devtools)
library(velocyto.R)
library(pagoda2)
library(Seurat)
library(hdf5r)
library(SeuratDisk)
library(SingleCellExperiment)
library(ggplot2)
library(ggrepel)
library(stringr)


setwd("/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis12_label_transfer_emmreml_edited/") 

IN_DIR <- "cellRank/wVelocity_input/"
IN_DIR2 <- "cellRank/terminal_state_probs/"
OUT_DIR <- "cellRank/terminal_state_analysis_output/"
dir.create(OUT_DIR)


donors <- c("S1", "S2", "S3", "S5", "S6", "S7", "S8", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S20", "S21")
cmpa_change <- vector(length=length(donors))
cmpb_change <- vector(length=length(donors))
cmpc_change <- vector(length=length(donors))
gmpa_change <- vector(length=length(donors))
gmpb_change <- vector(length=length(donors))
mlpa_change <- vector(length=length(donors))
mlpb_change <- vector(length=length(donors))
mepa_change <- vector(length=length(donors))
mepb_change <- vector(length=length(donors))
mepc_change <- vector(length=length(donors))
preBNK_change <- vector(length=length(donors))
lymphoid_change <- vector(length=length(donors))
myeloid_change <- vector(length=length(donors))
cmp_change <- vector(length=length(donors))
erythroid_change <- vector(length=length(donors))
cmpa_d90 <- vector(length=length(donors))
cmpb_d90 <- vector(length=length(donors))
cmpc_d90 <- vector(length=length(donors))
gmpa_d90 <- vector(length=length(donors))
gmpb_d90 <- vector(length=length(donors))
mlpa_d90 <- vector(length=length(donors))
mlpb_d90 <- vector(length=length(donors))
mepa_d90 <- vector(length=length(donors))
mepb_d90 <- vector(length=length(donors))
mepc_d90 <- vector(length=length(donors))
cmp_d90 <- vector(length=length(donors))
preBNK_d90 <- vector(length=length(donors))
lymphoid_d90 <- vector(length=length(donors))
myeloid_d90 <- vector(length=length(donors))
erythroid_d90 <- vector(length=length(donors))


cmpa_d0 <- vector(length=length(donors))
cmpb_d0 <- vector(length=length(donors))
cmpc_d0 <- vector(length=length(donors))
gmpa_d0 <- vector(length=length(donors))
gmpb_d0 <- vector(length=length(donors))
mlpa_d0 <- vector(length=length(donors))
mlpb_d0 <- vector(length=length(donors))
mepa_d0 <- vector(length=length(donors))
mepb_d0 <- vector(length=length(donors))
mepc_d0 <- vector(length=length(donors))
preBNK_d0 <- vector(length=length(donors))
erythroid_d0 <- vector(length=length(donors))
lymphoid_d0 <- vector(length=length(donors))
myeloid_d0 <- vector(length=length(donors))

plot_list <- list()
k <- 1
plot_list2 <- list()

for(i in 1:length(donors)){
  time <- c("Td0", "Tm3")
  time_label <- c("D0", "D90")
  
  cmpa_perc <- vector(length=length(time))
  cmpb_perc <- vector(length=length(time))
  cmpc_perc <- vector(length=length(time))
  gmpa_perc <- vector(length=length(time))
  gmpb_perc <- vector(length=length(time))
  mlpa_perc <- vector(length=length(time))
  mlpb_perc <- vector(length=length(time))
  mepa_perc <- vector(length=length(time))
  mepb_perc <- vector(length=length(time))
  mepc_perc <- vector(length=length(time))
  preBNK_perc <- vector(length=length(time))
  lymphoid_perc <- vector(length=length(time))
  myeloid_perc <- vector(length=length(time))
  cmp_perc <- vector(length=length(time))
  erythroid_perc <- vector(length=length(time))
  for(j in 1:length(time))
  {
    coordinates <- read.csv(file=paste0(IN_DIR, "cell_embeddings_", donors[i], "_", time[j],".csv"), header=TRUE)
    id <- read.csv(file=paste0(IN_DIR, "clusters_", donors[i], "_", time[j], ".csv"), header=TRUE)
    terminal_states <- read.csv(file=paste0(IN_DIR2, donors[i],"_",time[j],"_to_terminal_states.csv"), header=TRUE)
    
    id <- id[which(id$X %in% terminal_states$Cell.ID),]
    id <- id[order(match(as.character(id$X), as.character(terminal_states$Cell.ID))),]
    which(id$X != terminal_states$Cell.ID)
    colnames(id) <- c("X", "clust_names")
    
    row.names(coordinates) <- coordinates$X
    coordinates <- coordinates[which(coordinates$X %in% terminal_states$Cell.ID),]
    coordinates <- coordinates[order(match(as.character(coordinates$X), as.character(terminal_states$Cell.ID))),]
    which(coordinates$X != terminal_states$Cell.ID)
    
    terminal_states <- subset(terminal_states, select = -X)
    terminal_states <- subset(terminal_states, select = -Cell.ID)
    most_likely <- colnames(terminal_states)[max.col(terminal_states, ties.method = "first")]
    
    
    df <- data.frame(
      name <- id$X,
      clust <- id$clust_names,
      x <- coordinates$UMAP_1,
      y <- coordinates$UMAP_2,
      col <- most_likely
    )
    
    colnames(df) <- c("name", "clust", "x", "y", "col")
    df_final <- subset(df, df$clust %in% c("HSC_a", "HSC_b"))
    df_final$col <- factor(df_final$col, colnames(terminal_states))
    
    # myColorScale <- c("#FF7F0E","#FFD147","#CC9900", "#33CC00", "#99CC00", "#8C5782", "#990080", "#BA68C8",
    #                   "#46B8DA", "#0288D1", "#1A0099")
    
    
    myColorScale <- c("#006DDB", "#FF7F0E", "#00BB00", "#FF0000", "#BB00BB", "#924900", "#FF6DB6", "#749B58", "#00A9CC", "#AEC7E8", "#FFCC65")
    names(myColorScale) <- colnames(terminal_states)
    
    p <- ggplot(df_final, aes(x=x, y=y, color=col)) + geom_point(aes(fill=col),shape=21, size=2, alpha=1, stroke=0.4, color="black") + theme_void() +
     labs(title=paste0("Donor ", donors[i]," ",time_label[j]," \nPredicted terminal state"), x= "UMAP1", y="UMAP2") + scale_fill_manual(values=myColorScale) + theme(legend.position = "none")   
    
    plot_list[[k]] <- p
    k <- k + 1
    
    
    cmpa_perc[j] <- length(df_final$col[which(df_final$col == "CMP_a")])/length(df_final$col)
    cmpb_perc[j] <- length(df_final$col[which(df_final$col == "CMP_b")])/length(df_final$col)
    cmpc_perc[j] <- length(df_final$col[which(df_final$col == "CMP_c")])/length(df_final$col)
    gmpa_perc[j] <- length(df_final$col[which(df_final$col == "GMP_a")])/length(df_final$col)
    gmpb_perc[j] <- length(df_final$col[which(df_final$col == "GMP_b")])/length(df_final$col)
    mlpa_perc[j] <- length(df_final$col[which(df_final$col == "MLP_a")])/length(df_final$col)
    mlpb_perc[j] <- length(df_final$col[which(df_final$col == "MLP_b")])/length(df_final$col)
    preBNK_perc[j] <- length(df_final$col[which(df_final$col == "PreBNK")])/length(df_final$col)
    mepa_perc[j] <- length(df_final$col[which(df_final$col == "MEP_a")])/length(df_final$col)
    mepb_perc[j] <- length(df_final$col[which(df_final$col == "MEP_b")])/length(df_final$col)
    mepc_perc[j] <- length(df_final$col[which(df_final$col == "MEP_c")])/length(df_final$col)
    myeloid_perc[j] <- cmpa_perc[j] + cmpb_perc[j] + cmpc_perc[j] + gmpa_perc[j] + gmpb_perc[j]
    cmp_perc[j] <- cmpa_perc[j] + cmpb_perc[j] + cmpc_perc[j]
    lymphoid_perc[j] <- mlpa_perc[j] + mlpb_perc[j] + preBNK_perc[j]
    erythroid_perc[j] <- mepa_perc[j] + mepb_perc[j] + mepc_perc[j]
  }
  cmpa_change[i] <- cmpa_perc[2] - cmpa_perc[1]
  cmpb_change[i] <- cmpb_perc[2] - cmpb_perc[1]
  cmpc_change[i] <- cmpc_perc[2] - cmpc_perc[1]
  gmpa_change[i] <- gmpa_perc[2] - gmpa_perc[1]
  gmpb_change[i] <- gmpb_perc[2] - gmpb_perc[1]
  mlpa_change[i] <- mlpa_perc[2] - mlpa_perc[1]
  mlpb_change[i] <- mlpb_perc[2] - mlpb_perc[1]
  preBNK_change[i] <- preBNK_perc[2] - preBNK_perc[1]
  mepa_change[i] <- mepa_perc[2] - mepa_perc[1]
  mepb_change[i] <- mepb_perc[2] - mepb_perc[1]
  mepc_change[i] <- mepc_perc[2] - mepc_perc[1]
  myeloid_change[i] <- myeloid_perc[2] - myeloid_perc[1]
  cmp_change[i] <- cmp_perc[2] - cmp_perc[1]
  lymphoid_change[i] <- lymphoid_perc[2] - lymphoid_perc[1]
  erythroid_change[i] <- erythroid_perc[2] - erythroid_perc[1]
  
  df2 <- data.frame(
    timepoint <- c(rep("D0", length(colnames(terminal_states))), rep("D90", length(colnames(terminal_states))) ),
    Celltype <- rep(paste0(colnames(terminal_states)), 2),
    values <- c(cmpa_perc[1], cmpb_perc[1], cmpc_perc[1], gmpa_perc[1], gmpb_perc[1], mepa_perc[1], mepb_perc[1], mepc_perc[1], mlpa_perc[1], mlpb_perc[1], preBNK_perc[1], cmpa_perc[2], cmpb_perc[2], cmpc_perc[2], gmpa_perc[2], gmpb_perc[2], mepa_perc[2], mepb_perc[2], mepc_perc[2], mlpa_perc[2], mlpb_perc[2], preBNK_perc[2])
  )
  colnames(df2) <- c("timepoint", "Celltype", "values")
  
  myColorScale <- c("#006DDB", "#FF7F0E", "#00BB00", "#FF0000", "#BB00BB", "#924900", "#FF6DB6", "#749B58", "#00A9CC", "#AEC7E8", "#FFCC65")
  names(myColorScale) <- colnames(terminal_states)
  
  # ggplot(data=df2, aes(x=timepoint, y=values, fill=Celltype)) +
  #   geom_bar(stat="identity") + theme_light() + scale_fill_manual(values=myColorScale) + labs(x="", y="Percentage")+ coord_flip()
  
  df2$timepoint <- factor(df2$timepoint, c("D0", "D90"))
  p2 <- ggplot(df2, aes(x=Celltype, y=values, group=timepoint)) +
    geom_line(aes(linetype=timepoint, color=timepoint))+ 
    geom_point(size=3, aes(color=Celltype)) + theme_test() + labs(x="", y="Percentage", linetype="Timepoint", title=paste0("Donor ", donors[i], " HSC \nPredicted terminal state")) + scale_color_manual(values=myColorScale)+guides(color=FALSE)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  plot_list2[[i]] <- p2
  
  #scale_color_manual(values=c("#A50021", "#290AD8")) +  new_scale("color") +
  
  #record day90 percentages
  cmpa_d90[i] <- cmpa_perc[2] 
  cmpb_d90[i] <- cmpb_perc[2] 
  cmpc_d90[i] <- cmpc_perc[2] 
  gmpa_d90[i] <- gmpa_perc[2] 
  gmpb_d90[i] <- gmpb_perc[2] 
  mlpa_d90[i] <- mlpa_perc[2] 
  mlpb_d90[i] <- mlpb_perc[2] 
  mepa_d90[i] <- mepa_perc[2] 
  mepb_d90[i] <- mepb_perc[2] 
  mepc_d90[i] <- mepc_perc[2] 
  preBNK_d90[i] <- preBNK_perc[2] 
  myeloid_d90[i] <- myeloid_perc[2] 
  lymphoid_d90[i] <- lymphoid_perc[2] 
  erythroid_d90[i] <- erythroid_perc[2] 
  cmp_d90[i] <- cmp_perc[2] 
  
  #record day0 percentages
  cmpa_d0[i] <- cmpa_perc[1] 
  cmpb_d0[i] <- cmpb_perc[1] 
  cmpc_d0[i] <- cmpc_perc[1] 
  gmpa_d0[i] <- gmpa_perc[1] 
  gmpb_d0[i] <- gmpb_perc[1] 
  mlpa_d0[i] <- mlpa_perc[1] 
  mlpb_d0[i] <- mlpb_perc[1] 
  preBNK_d0[i] <- preBNK_perc[1] 
  mepa_d0[i] <- mepa_perc[1] 
  mepb_d0[i] <- mepb_perc[1] 
  mepc_d0[i] <- mepc_perc[1] 
  myeloid_d90[i] <- myeloid_perc[1] 
  lymphoid_d90[i] <- lymphoid_perc[1] 
  erythroid_d90[i] <- erythroid_perc[1] 
 
  
  print(i)
}

donors <- c("S1", "S2", "S3", "S5", "S6", "S7", "S8", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S20", "S21")
time <- c("Td0", "Tm3")
k <- 1
for(i in 1:length(donors)) {
  for(j in 1:length(time))
  {
    file_name = paste0(OUT_DIR, donors[i],"_",time[j], "_term_state_umap.tiff")
    tiff(file_name, units="in", width=3, height=3, res=350)
    print(plot_list[[k]])
    dev.off()
    print(k)
    k <- k+1
  }
  
}

for(i in 1:length(donors)) 
{

    file_name = paste0(OUT_DIR, donors[i],"_percentage_lineplot.tiff")
    tiff(file_name, units="in", width=7, height=3, res=350)
    print(plot_list2[[i]])
    dev.off()
    print(i)
}


wilcox_pvals <- vector(length=11)
names(wilcox_pvals) <- c('CMP_a', 'CMP_b', 'CMP_c', 'GMP_a', "GMP_b", 'MEP_a', "MEP_b", "MEP_c", "MLP_a", "MLP_b", "PreBNK")

# CMP_b boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_cmpb <- data.frame(
  val <- cmpb_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
#wilcox.test(df_cmpb$val[which(df_cmpb$group == "BCG")], df_cmpb$val[which(df_cmpb$group == "CTL")])

tiff(paste0(OUT_DIR, "cmp_b_term_state.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpb, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_b Terminal state prediction") + scale_fill_manual(values=c("#FFD320", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black"))
dev.off()
tiff(paste0(OUT_DIR, "cmp_b_term_state_labels.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpb, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_b Terminal state prediction") + scale_fill_manual(values=c("#FFD320", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()


df_cmpb <- subset(df_cmpb, df_cmpb$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_cmpb) <- c("val", "d", "group", "l")

tiff(paste0(OUT_DIR, "cmp_b_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpb, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_b Terminal state prediction") + scale_fill_manual(values=c("#FFD320", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

#values shifted relative to placebo
df_cmpb_subset <- subset(df_cmpb, df_cmpb$group == "CTL")
cmpb_subset_mean <- median(df_cmpb_subset$val)


new_vals <- df_cmpb$val - cmpb_subset_mean
df_cmpb_new <- data.frame(
  val_corrected <- new_vals,
  d <- df_cmpb$d,
  group <- df_cmpb$group,
  l <- df_cmpb$l
)

w_test <- wilcox.test(df_cmpb_new$val_corrected[which(df_cmpb_new$group == "BCG")], df_cmpb_new$val_corrected[which(df_cmpb_new$group == "CTL")])$p.value

tiff(paste0(OUT_DIR, "cmp_b_term_state_labels_2_3_12_removed_relative.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpb_new, aes(x=group, y=val_corrected, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change (Tm3-Td0) relative to placebo", fill="Group", title=paste0("CMP_b Terminal state prediction; \n pval = ", round(w_test,3))) + scale_fill_manual(values=c("#FFD320", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

tiff(paste0(OUT_DIR, "cmp_b_term_state_labels_2_3_12_removed_relative_no_highlight.tiff"), units="in", width=2.4, height=3, res=300)
ggplot(data=df_cmpb_new, aes(x=group, y=val_corrected)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_classic()+labs(x="", y="∆(%D90-%D0) normalized", fill="Group", title=paste0("∆Percentage CMP_b; \n pval = ", round(w_test,3))) + scale_fill_manual(values=c("#FF7F0E", "#004586"))+
   theme(legend.position="none")
dev.off()

#geom_hline(yintercept=0, color="brown") + 


#d90 BCG vs. placebo comparison plots

df_cmpb_d90 <- data.frame(
  vals_d90 <- cmpb_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

df_cmpb_d90 <- subset(df_cmpb_d90, df_cmpb_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_cmpb_d90) <- c("vals_d90", "d", "group", "l")

wilcox_pvals['CMP_b'] <- wilcox_test <- wilcox.test(df_cmpb_d90$vals_d90[which(df_cmpb_d90$group == "BCG")], df_cmpb_d90$vals_d90[which(df_cmpb_d90$group == "CTL")])$p.value 

tiff(paste0(OUT_DIR, "cmp_b_term_state_d90.tiff"), units="in", width=2.4, height=3, res=250)
ggplot(data=df_cmpb_d90, aes(x=group, y=vals_d90)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_classic()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent CMPb biased HSCs", fill="Group", title=paste0("CMP_b bias Day 90; \npval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#FFD320", "#004586"))+
   theme(legend.position = "none")
dev.off()
#version no dots
tiff(paste0(OUT_DIR, "cmp_b_term_state_d90_no_points.tiff"), units="in", width=2.4, height=3, res=250)
ggplot(data=df_cmpb_d90, aes(x=group, y=vals_d90)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + theme_classic()+
  geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent CMPb biased HSCs", fill="Group", title=paste0("CMP_b bias Day 90; \npval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#FFD320", "#004586"))+
  theme(legend.position = "none")
dev.off()



# CMP_a boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_cmpa <- data.frame(
  val <- cmpa_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

tiff(paste0(OUT_DIR, "cmp_a_term_state.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpa, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_a Terminal state prediction") + scale_fill_manual(values=c("#FF7F0E", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black"))
dev.off()
tiff(paste0(OUT_DIR, "cmp_a_term_state_labels.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpa, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_a Terminal state prediction") + scale_fill_manual(values=c("#FF7F0E", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

df_cmpa <- subset(df_cmpa, df_cmpa$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_cmpa) <- c("val", "d", "group", "l")

tiff(paste0(OUT_DIR, "cmp_a_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpa, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_a Terminal state prediction") + scale_fill_manual(values=c("#FF7F0E", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()


#values shifted relative to placebo
df_cmpa_subset <- subset(df_cmpa, df_cmpa$group == "CTL")
cmpa_subset_mean <- median(df_cmpa_subset$val)


new_vals <- df_cmpa$val - cmpa_subset_mean
df_cmpa_new <- data.frame(
  val_corrected <- new_vals,
  d <- df_cmpa$d,
  group <- df_cmpa$group,
  l <- df_cmpa$l
)

w_test <- wilcox.test(df_cmpa_new$val_corrected[which(df_cmpa_new$group == "BCG")], df_cmpa_new$val_corrected[which(df_cmpa_new$group == "CTL")])$p.value

tiff(paste0(OUT_DIR, "cmp_a_term_state_labels_2_3_12_removed_relative.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpa_new, aes(x=group, y=val_corrected, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change (Tm3-Td0) relative to placebo", fill="Group", title=paste0("CMP_a Terminal state prediction; \n pval = ", round(w_test,3))) + scale_fill_manual(values=c("#FF7F0E", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

tiff(paste0(OUT_DIR, "cmp_a_term_state_labels_2_3_12_removed_relative_no_highlight.tiff"), units="in", width=2.4, height=3, res=300)
ggplot(data=df_cmpb_new, aes(x=group, y=val_corrected)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_classic()+labs(x="", y="∆(%D90-%D0) normalized", fill="Group", title=paste0("∆Percentage CMP_a; \n pval = ", round(w_test,3))) + scale_fill_manual(values=c("#006DDB", "#004586"))+
  theme(legend.position="none")
dev.off()

#d90 BCG vs. placebo comparison plots

df_cmpa_d90 <- data.frame(
  vals_d90 <- cmpa_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

df_cmpa_d90 <- subset(df_cmpa_d90, df_cmpa_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_cmpa_d90) <- c("vals_d90", "d", "group", "l")

wilcox_pvals['CMP_a'] <- wilcox.test(df_cmpa_d90$vals_d90[which(df_cmpa_d90$group == "BCG")], df_cmpa_d90$vals_d90[which(df_cmpa_d90$group == "CTL")])$p.value 

tiff(paste0(OUT_DIR, "cmp_a_term_state_d90.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpa_d90, aes(x=group, y=vals_d90, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent cmpa biased HSCs", fill="Group", title=paste0("CMP_a bias Day 90; pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#FF7F0E", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()



# CMP_c boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_cmpc <- data.frame(
  val <- cmpc_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

tiff(paste0(OUT_DIR, "cmp_c_term_state.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpc, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_c Terminal state prediction") + scale_fill_manual(values=c("#CC9900", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black"))
dev.off()
tiff(paste0(OUT_DIR, "cmp_c_term_state_labels.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpc, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_c Terminal state prediction") + scale_fill_manual(values=c("#CC9900", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

df_cmpc <- subset(df_cmpc, df_cmpc$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_cmpc) <- c("val", "d", "group", "l")

tiff(paste0(OUT_DIR, "cmp_c_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpc, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_c Terminal state prediction") + scale_fill_manual(values=c("#CC9900", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

#values shifted relative to placebo
df_cmpc_subset <- subset(df_cmpc, df_cmpc$group == "CTL")
cmpc_subset_mean <- median(df_cmpc_subset$val)


new_vals <- df_cmpc$val - cmpc_subset_mean
df_cmpc_new <- data.frame(
  val_corrected <- new_vals,
  d <- df_cmpc$d,
  group <- df_cmpc$group,
  l <- df_cmpc$l
)

w_test <- wilcox.test(df_cmpc_new$val_corrected[which(df_cmpc_new$group == "BCG")], df_cmpc_new$val_corrected[which(df_cmpc_new$group == "CTL")])$p.value

tiff(paste0(OUT_DIR, "cmp_c_term_state_labels_2_3_12_removed_relative.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpc_new, aes(x=group, y=val_corrected, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change (Tm3-Td0) relative to placebo", fill="Group", title=paste0("CMP_c Terminal state prediction; \n pval = ", round(w_test,3))) + scale_fill_manual(values=c("#CC9900", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

tiff(paste0(OUT_DIR, "cmp_c_term_state_labels_2_3_12_removed_relative_no_highlight.tiff"), units="in", width=2.4, height=3, res=300)
ggplot(data=df_cmpb_new, aes(x=group, y=val_corrected)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_classic()+labs(x="", y="∆(%D90-%D0) normalized", fill="Group", title=paste0("∆Percentage CMP_c; \n pval = ", round(w_test,3))) + scale_fill_manual(values=c("#00BB00", "#004586"))+
   theme(legend.position="none")
dev.off()

#d90 BCG vs. placebo comparison plots

df_cmpc_d90 <- data.frame(
  vals_d90 <- cmpc_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

df_cmpc_d90 <- subset(df_cmpc_d90, df_cmpc_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_cmpc_d90) <- c("vals_d90", "d", "group", "l")

wilcox_pvals['CMP_c'] <- wilcox.test(df_cmpc_d90$vals_d90[which(df_cmpc_d90$group == "BCG")], df_cmpc_d90$vals_d90[which(df_cmpc_d90$group == "CTL")])$p.value 

tiff(paste0(OUT_DIR, "cmp_c_term_state_d90.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpc_d90, aes(x=group, y=vals_d90, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent cmpc biased HSCs", fill="Group", title=paste0("CMP_c bias Day 90; pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#CC9900", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()


# CMP boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_cmp <- data.frame(
  val <- cmp_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

tiff(paste0(OUT_DIR, "cmp_term_state.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmp, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP Terminal state prediction") + scale_fill_manual(values=c("#CC9900", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black"))
dev.off()
tiff(paste0(OUT_DIR, "cmp_term_state_labels.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmp, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP Terminal state prediction") + scale_fill_manual(values=c("#CC9900", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

df_cmp <- subset(df_cmp, df_cmp$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_cmp) <- c("val", "d", "group", "l")
write.csv(df_cmp, file=paste0(OUT_DIR, "df_cmp_raw.csv"))

tiff(paste0(OUT_DIR, "cmp_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmp, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP Terminal state prediction") + scale_fill_manual(values=c("#CC9900", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

#values shifted relative to placebo
df_cmp_subset <- subset(df_cmp, df_cmp$group == "CTL")
cmp_subset_mean <- median(df_cmp_subset$val)


new_vals <- df_cmp$val - cmp_subset_mean
df_cmp_new <- data.frame(
  val_corrected <- new_vals,
  d <- df_cmp$d,
  group <- df_cmp$group,
  l <- df_cmp$l
)
colnames(df_cmp_new) <- c("val_corrected", "d", "group", "l")
write.csv(df_cmp_new, file=paste0(OUT_DIR, "df_cmp_new.csv"))

w_test <- wilcox.test(df_cmp_new$val_corrected[which(df_cmp_new$group == "BCG")], df_cmp_new$val_corrected[which(df_cmp_new$group == "CTL")])$p.value

tiff(paste0(OUT_DIR, "cmp_term_state_labels_2_3_12_removed_relative.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmp_new, aes(x=group, y=val_corrected, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change (Tm3-Td0) relative to placebo", fill="Group", title=paste0("CMP Terminal state prediction; \n pval = ", round(w_test,3))) + scale_fill_manual(values=c("#CC9900", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

tiff(paste0(OUT_DIR, "cmp_term_state_labels_2_3_12_removed_relative_no_highlight.tiff"), units="in", width=2.4, height=3, res=300)
ggplot(data=df_cmp_new, aes(x=group, y=val_corrected)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_classic()+labs(x="", y="∆(%D90-%D0) normalized", fill="Group", title=paste0("∆Percentage CMP; \n pval = ", round(w_test,3))) + scale_fill_manual(values=c("#F72836", "#004586"))+
   theme(legend.position="none")
dev.off()

#d90 BCG vs. placebo comparison plots

df_cmp_d90 <- data.frame(
  vals_d90 <- cmp_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

df_cmp_d90 <- subset(df_cmp_d90, df_cmp_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_cmp_d90) <- c("vals_d90", "d", "group", "l")

wilcox_test <- wilcox.test(df_cmp_d90$vals_d90[which(df_cmp_d90$group == "BCG")], df_cmp_d90$vals_d90[which(df_cmp_d90$group == "CTL")])$p.value 

tiff(paste0(OUT_DIR, "cmp_c_term_state_d90.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmp_d90, aes(x=group, y=vals_d90, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent cmp biased HSCs", fill="Group", title=paste0("CMP_c bias Day 90; pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#CC9900", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()



# Myeloid boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_myeloid <- data.frame(
  val <- myeloid_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
wilcox.test(df_myeloid$val[which(df_myeloid$group == "BCG")], df_myeloid$val[which(df_myeloid$group == "CTL")])

tiff(paste0(OUT_DIR, "myeloid_term_state.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_myeloid, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="Myeloid Terminal state prediction") + scale_fill_manual(values=c("#B82E2E", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black"))
dev.off()
tiff(paste0(OUT_DIR, "myeloid_term_state_labels.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_myeloid, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="Myeloid Terminal state prediction") + scale_fill_manual(values=c("#B82E2E", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

df_myeloid <- subset(df_myeloid, df_myeloid$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_myeloid) <- c("val", "d", "group", "l")

tiff(paste0(OUT_DIR, "myeloid_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_myeloid, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="Myeloid Terminal state prediction") + scale_fill_manual(values=c("#B82E2E", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

#values shifted relative to placebo
df_myeloid_subset <- subset(df_myeloid, df_myeloid$group == "CTL")
myeloid_subset_mean <- median(df_myeloid_subset$val)


new_vals <- df_myeloid$val - myeloid_subset_mean
df_myeloid_new <- data.frame(
  val_corrected <- new_vals,
  d <- df_myeloid$d,
  group <- df_myeloid$group,
  l <- df_myeloid$l
)

wilcox_test <- wilcox.test(df_myeloid_new$val_corrected[which(df_myeloid_new$group == "BCG")], df_myeloid_new$val_corrected[which(df_myeloid_new$group == "CTL")])$p.value

tiff(paste0(OUT_DIR, "myeloid_term_state_labels_2_3_12_removed_relative.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_myeloid_new, aes(x=group, y=val_corrected, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change (Tm3-Td0) relative to placebo", fill="Group", title=paste0("Myeloid Terminal state prediction; \n pval = ", round(wilcox_test,3)) )+ scale_fill_manual(values=c("#B82E2E", "#004586"))+
  scale_color_manual(values=c(  "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

tiff(paste0(OUT_DIR, "myeloid_term_state_labels_2_3_12_removed_relative_no_highlight.tiff"), units="in", width=2.8, height=3, res=300)
ggplot(data=df_cmp_new, aes(x=group, y=val_corrected, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="∆(%D90-%D0) relative to placebo", fill="Group", title=paste0("∆Percentage myeloid; \n pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#B82E2E", "#004586"))+
  geom_text_repel(size=2) + theme(legend.position="none")
dev.off()

#d90 BCG vs. placebo comparison plots

df_myeloid_d90 <- data.frame(
  vals_d90 <- myeloid_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

df_myeloid_d90 <- subset(df_myeloid_d90, df_myeloid_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_myeloid_d90) <- c("vals_d90", "d", "group", "l")

wilcox_test <- wilcox.test(df_myeloid_d90$vals_d90[which(df_myeloid_d90$group == "BCG")], df_myeloid_d90$vals_d90[which(df_myeloid_d90$group == "CTL")])$p.value 

tiff(paste0(OUT_DIR, "myeloid_term_state_d90.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_myeloid_d90, aes(x=group, y=vals_d90, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent myeloid biased HSCs", fill="Group", title=paste0("Myeloid bias Day 90; pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#B82E2E", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()



# lymphoid boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_lymphoid <- data.frame(
  val <- lymphoid_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

tiff(paste0(OUT_DIR, "lymphoid_term_state.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_lymphoid, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="lymphoid Terminal state prediction") + scale_fill_manual(values=c("#0072B2", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black"))
dev.off()
tiff(paste0(OUT_DIR, "lymphoid_term_state_labels.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_lymphoid, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="lymphoid Terminal state prediction") + scale_fill_manual(values=c("#0072B2", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

df_lymphoid <- subset(df_lymphoid, df_lymphoid$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_lymphoid) <- c("val", "d", "group", "l")
write.csv(df_lymphoid, file=paste0(OUT_DIR, "df_lymphoid_raw.csv"))

tiff(paste0(OUT_DIR, "lymphoid_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_lymphoid, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="lymphoid Terminal state prediction") + scale_fill_manual(values=c("#0072B2", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

#values shifted relative to placebo
df_lymphoid_subset <- subset(df_lymphoid, df_lymphoid$group == "CTL")
lymphoid_subset_mean <- median(df_lymphoid_subset$val)


new_vals <- df_lymphoid$val - lymphoid_subset_mean
df_lymphoid_new <- data.frame(
  val_corrected <- new_vals,
  d <- df_lymphoid$d,
  group <- df_lymphoid$group,
  l <- df_lymphoid$l
)
colnames(df_lymphoid_new) <- c("val_corrected", "d", "group", "l")
write.csv(df_lymphoid_new, file=paste0(OUT_DIR, "df_lymphoid_new.csv"))

wilcox_test <- wilcox.test(df_lymphoid_new$val_corrected[which(df_lymphoid_new$group == "BCG")], df_lymphoid_new$val_corrected[which(df_lymphoid_new$group == "CTL")])$p.value

tiff(paste0(OUT_DIR, "lymphoid_term_state_labels_2_3_12_removed_relative.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_lymphoid_new, aes(x=group, y=val_corrected, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change (Tm3-Td0) relative to placebo", fill="Group", title=paste0("Lymphoid Terminal state prediction; \n pval = ", round(wilcox_test,3)) )+ scale_fill_manual(values=c("#0072B2", "#004586"))+
  scale_color_manual(values=c(  "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

tiff(paste0(OUT_DIR, "lymphoid_term_state_labels_2_3_12_removed_relative_no_highlight.tiff"), units="in", width=2.4, height=3, res=300)
ggplot(data=df_lymphoid_new, aes(x=group, y=val_corrected)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_classic()+labs(x="", y="∆(%D90-%D0) normalized", fill="Group", title=paste0("∆Percentage lymphoid; \n pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#00A9CC", "#004586"))+
  theme(legend.position="none")
dev.off()

#d90 BCG vs. placebo comparison plots

df_lymphoid_d90 <- data.frame(
  vals_d90 <- lymphoid_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

df_lymphoid_d90 <- subset(df_lymphoid_d90, df_lymphoid_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_lymphoid_d90) <- c("vals_d90", "d", "group", "l")

wilcox_test <- wilcox.test(df_lymphoid_d90$vals_d90[which(df_lymphoid_d90$group == "BCG")], df_lymphoid_d90$vals_d90[which(df_lymphoid_d90$group == "CTL")])$p.value 

tiff(paste0(OUT_DIR, "lymphoid_term_state_d90.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_lymphoid_d90, aes(x=group, y=vals_d90, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent lymphoid biased HSCs", fill="Group", title=paste0("lymphoid bias Day 90; pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#0072B2", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()


# erythroid boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_erythroid <- data.frame(
  val <- erythroid_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

tiff(paste0(OUT_DIR, "erythroid_term_state.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_erythroid, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="erythroid Terminal state prediction") + scale_fill_manual(values=c("#990099", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black"))
dev.off()
tiff(paste0(OUT_DIR, "erythroid_term_state_labels.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_erythroid, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="erythroid Terminal state prediction") + scale_fill_manual(values=c("#990099", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

df_erythroid <- subset(df_erythroid, df_erythroid$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_erythroid) <- c("val", "d", "group", "l")
write.csv(df_erythroid, file=paste0(OUT_DIR, "df_erythroid_raw.csv"))

tiff(paste0(OUT_DIR, "erythroid_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_erythroid, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="erythroid Terminal state prediction") + scale_fill_manual(values=c("#990099", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

#values shifted relative to placebo
df_erythroid_subset <- subset(df_erythroid, df_erythroid$group == "CTL")
erythroid_subset_mean <- median(df_erythroid_subset$val)


new_vals <- df_erythroid$val - erythroid_subset_mean
df_erythroid_new <- data.frame(
  val_corrected <- new_vals,
  d <- df_erythroid$d,
  group <- df_erythroid$group,
  l <- df_erythroid$l
)
colnames(df_erythroid_new) <- c("val_corrected", "d", "group", "l")
write.csv(df_erythroid_new, file=paste0(OUT_DIR, "df_erythroid_new.csv"))

wilcox_test <- wilcox.test(df_erythroid_new$val_corrected[which(df_erythroid_new$group == "BCG")], df_erythroid_new$val_corrected[which(df_erythroid_new$group == "CTL")])$p.value

tiff(paste0(OUT_DIR, "erythroid_term_state_labels_2_3_12_removed_relative.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_erythroid_new, aes(x=group, y=val_corrected, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change (Tm3-Td0) relative to placebo", fill="Group", title=paste0("Erythroid Terminal state prediction; \n pval = ", round(wilcox_test,3)) )+ scale_fill_manual(values=c("#990099", "#004586"))+
  scale_color_manual(values=c(  "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

tiff(paste0(OUT_DIR, "erythroid_term_state_labels_2_3_12_removed_relative_no_highlight.tiff"), units="in", width=2.4, height=3, res=300)
ggplot(data=df_erythroid_new, aes(x=group, y=val_corrected)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_classic()+ labs(x="", y="∆(%D90-%D0) normalized", fill="Group", title=paste0("∆Percentage erythroid; \n pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c( "#924900", "#004586"))+
  theme(legend.position="none")
dev.off()

#d90 BCG vs. placebo comparison plots

df_erythroid_d90 <- data.frame(
  vals_d90 <- erythroid_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

df_erythroid_d90 <- subset(df_erythroid_d90, df_erythroid_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_erythroid_d90) <- c("vals_d90", "d", "group", "l")

wilcox_test <- wilcox.test(df_erythroid_d90$vals_d90[which(df_erythroid_d90$group == "BCG")], df_erythroid_d90$vals_d90[which(df_erythroid_d90$group == "CTL")])$p.value 

tiff(paste0(OUT_DIR, "erythroid_term_state_d90.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_erythroid_d90, aes(x=group, y=vals_d90, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent erythroid biased HSCs", fill="Group", title=paste0("erythroid bias Day 90; pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#990099", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()



# mepa boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_mepa <- data.frame(
  val <- mepa_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

tiff(paste0(OUT_DIR, "mepa_term_state.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_mepa, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="MEP_a Terminal state prediction") + scale_fill_manual(values=c("#CE6DBD", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black"))
dev.off()
tiff(paste0(OUT_DIR, "mepa_term_state_labels.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_mepa, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="MEP_a Terminal state prediction") + scale_fill_manual(values=c("#CE6DBD", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

df_mepa <- subset(df_mepa, df_mepa$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_mepa) <- c("val", "d", "group", "l")

tiff(paste0(OUT_DIR, "mepa_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_mepa, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="MEP_a Terminal state prediction") + scale_fill_manual(values=c("#CE6DBD", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

#values shifted relative to placebo
df_mepa_subset <- subset(df_mepa, df_mepa$group == "CTL")
mepa_subset_mean <- median(df_mepa_subset$val)


new_vals <- df_mepa$val - mepa_subset_mean
df_mepa_new <- data.frame(
  val_corrected <- new_vals,
  d <- df_mepa$d,
  group <- df_mepa$group,
  l <- df_mepa$l
)

wilcox_test <- wilcox.test(df_mepa_new$val_corrected[which(df_mepa_new$group == "BCG")], df_mepa_new$val_corrected[which(df_mepa_new$group == "CTL")])$p.value

tiff(paste0(OUT_DIR, "mepa_term_state_labels_2_3_12_removed_relative.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_mepa_new, aes(x=group, y=val_corrected, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change (Tm3-Td0) relative to placebo", fill="Group", title=paste0("mepa Terminal state prediction; \n pval = ", round(wilcox_test,3)) )+ scale_fill_manual(values=c("#CE6DBD", "#004586"))+
  scale_color_manual(values=c(  "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

tiff(paste0(OUT_DIR, "mepa_term_state_labels_2_3_12_removed_relative_no_highlight.tiff"), units="in", width=2.8, height=3, res=300)
ggplot(data=df_mepa_new, aes(x=group, y=val_corrected, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="∆(%D90-%D0) relative to placebo", fill="Group", title=paste0("∆Percentage MEPa; \n pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#CE6DBD", "#004586"))+
  geom_text_repel(size=2) + theme(legend.position="none")
dev.off()

#d90 BCG vs. placebo comparison plots

df_mepa_d90 <- data.frame(
  vals_d90 <- mepa_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
wilcox_pvals['MEP_a'] <- wilcox.test(df_mepa_d90$vals_d90[which(df_mepa_d90$group == "BCG")], df_mepa_d90$vals_d90[which(df_mepa_d90$group == "CTL")])$p.value 

ggplot(data=df_mepa_d90, aes(x=group, y=vals_d90, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent mepa biased HSCs", fill="Group", title=paste0("mepa bias Day 90; pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#CE6DBD", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)



# mepb boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_mepb <- data.frame(
  val <- mepb_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
#d90 BCG vs. placebo comparison plots

df_mepb_d90 <- data.frame(
  vals_d90 <- mepb_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
wilcox_pvals['MEP_b'] <- wilcox.test(df_mepb_d90$vals_d90[which(df_mepb_d90$group == "BCG")], df_mepb_d90$vals_d90[which(df_mepb_d90$group == "CTL")])$p.value 


# mepc boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_mepc <- data.frame(
  val <- mepc_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
#d90 BCG vs. placebo comparison plots

df_mepc_d90 <- data.frame(
  vals_d90 <- mepc_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
wilcox_pvals['MEP_c'] <- wilcox.test(df_mepc_d90$vals_d90[which(df_mepc_d90$group == "BCG")], df_mepc_d90$vals_d90[which(df_mepc_d90$group == "CTL")])$p.value 




# GMP_a boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_gmpa <- data.frame(
  val <- gmpa_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

tiff(paste0(OUT_DIR, "gmp_a_term_state.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_gmpa, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="GMP_a Terminal state prediction") + scale_fill_manual(values=c("#33CC00", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black"))
dev.off()
tiff(paste0(OUT_DIR, "gmp_a_term_state_labels.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_gmpa, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="GMP_a Terminal state prediction") + scale_fill_manual(values=c("#33CC00", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

df_gmpa <- subset(df_gmpa, df_gmpa$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_gmpa) <- c("val", "d", "group", "l")

tiff(paste0(OUT_DIR, "gmp_a_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_gmpa, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="GMP_a Terminal state prediction") + scale_fill_manual(values=c("#33CC00", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()


#d90 BCG vs. placebo comparison plots

df_gmpa_d90 <- data.frame(
  vals_d90 <- gmpa_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

df_gmpa_d90 <- subset(df_gmpa_d90, df_gmpa_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_gmpa_d90) <- c("vals_d90", "d", "group", "l")

wilcox_pvals['GMP_a'] <- wilcox.test(df_gmpa_d90$vals_d90[which(df_gmpa_d90$group == "BCG")], df_gmpa_d90$vals_d90[which(df_gmpa_d90$group == "CTL")])$p.value 

tiff(paste0(OUT_DIR, "gmpa_term_state_d90.tiff"), units="in", width=5, height=4, res=250)
ggplot(data=df_gmpa_d90, aes(x=group, y=vals_d90, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent gmpa biased HSCs", fill="Group", title=paste0("GMP_a bias Day 90; pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#33CC00", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()


# gmpb boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_gmpb <- data.frame(
  val <- gmpb_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

tiff(paste0(OUT_DIR, "gmpb_term_state.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_gmpb, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="GMP_b Terminal state prediction") + scale_fill_manual(values=c("#99CC00", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black"))
dev.off()
tiff(paste0(OUT_DIR, "gmpb_term_state_labels.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_gmpb, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="GMP_b Terminal state prediction") + scale_fill_manual(values=c("#99CC00", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

df_gmpb <- subset(df_gmpb, df_gmpb$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_gmpb) <- c("val", "d", "group", "l")

tiff(paste0(OUT_DIR, "gmpb_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_gmpb, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="GMP_b Terminal state prediction") + scale_fill_manual(values=c("#99CC00", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

#d90 BCG vs. placebo comparison plots

df_gmpb_d90 <- data.frame(
  vals_d90 <- gmpb_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

df_gmpb_d90 <- subset(df_gmpb_d90, df_gmpb_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_gmpb_d90) <- c("vals_d90", "d", "group", "l")

wilcox_pvals['GMP_b'] <- wilcox.test(df_gmpb_d90$vals_d90[which(df_gmpb_d90$group == "BCG")], df_gmpb_d90$vals_d90[which(df_gmpb_d90$group == "CTL")])$p.value 

tiff(paste0(OUT_DIR, "gmpb_term_state_d90.tiff"), units="in", width=5, height=4, res=250)
ggplot(data=df_gmpb_d90, aes(x=group, y=vals_d90, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent gmpb biased HSCs", fill="Group", title=paste0("GMP_b bias Day 90; pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#99CC00", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

# mlpb boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_mlpb <- data.frame(
  val <- mlpb_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

tiff(paste0(OUT_DIR, "mlpb_term_state.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_mlpb, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="MLP_b Terminal state prediction") + scale_fill_manual(values=c("#83CAFF", "#004586"))+
  scale_color_manual(values=c( "#DD4477", "black"))
dev.off()
tiff(paste0(OUT_DIR, "mlpb_term_state_labels.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_mlpb, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="MLP_b Terminal state prediction") + scale_fill_manual(values=c("#83CAFF", "#004586"))+
  scale_color_manual(values=c( "#DD4477", "black")) + geom_text_repel(size=2)
dev.off()

df_mlpb <- subset(df_mlpb, df_mlpb$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_mlpb) <- c("val", "d", "group", "l")

tiff(paste0(OUT_DIR, "mlpb_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_mlpb, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="MLP_b Terminal state prediction") + scale_fill_manual(values=c("#83CAFF", "#004586"))+
  scale_color_manual(values=c( "#DD4477", "black")) + geom_text_repel(size=2)
dev.off()

#values shifted relative to placebo
df_mlpb_subset <- subset(df_mlpb, df_mlpb$group == "CTL")
mlpb_subset_mean <- median(df_mlpb_subset$val)


new_vals <- df_mlpb$val - mlpb_subset_mean
df_mlpb_new <- data.frame(
  val_corrected <- new_vals,
  d <- df_mlpb$d,
  group <- df_mlpb$group,
  l <- df_mlpb$l
)

wilcox_test <- wilcox.test(df_mlpb_new$val_corrected[which(df_mlpb_new$group == "BCG")], df_mlpb_new$val_corrected[which(df_mlpb_new$group == "CTL")])$p.value

tiff(paste0(OUT_DIR, "mlpb_term_state_labels_2_3_12_removed_relative.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_mlpb_new, aes(x=group, y=val_corrected, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change (Tm3-Td0) relative to placebo", fill="Group", title=paste0("mlpb Terminal state prediction; \n pval = ", round(wilcox_test,3)) )+ scale_fill_manual(values=c("#83CAFF", "#004586"))+
  scale_color_manual(values=c(  "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

tiff(paste0(OUT_DIR, "mlpb_term_state_labels_2_3_12_removed_relative_no_highlight.tiff"), units="in", width=2.8, height=3, res=300)
ggplot(data=df_mlpb_new, aes(x=group, y=val_corrected, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="∆(%D90-%D0) relative to placebo", fill="Group", title=paste0("∆Percentage MLPb; \n pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#83CAFF", "#004586"))+
  geom_text_repel(size=2) + theme(legend.position="none")
dev.off()




# mlpa boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_mlpa <- data.frame(
  val <- mlpa_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

tiff(paste0(OUT_DIR, "mlpa_term_state.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_mlpa, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="MLP_a Terminal state prediction") + scale_fill_manual(values=c("#83CAFF", "#004586"))+
  scale_color_manual(values=c( "#DD4477", "black"))
dev.off()

#d90 BCG vs. placebo comparison plots

df_mlpa_d90 <- data.frame(
  vals_d90 <- mlpa_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

df_mlpa_d90 <- subset(df_mlpa_d90, df_mlpa_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_mlpa_d90) <- c("vals_d90", "d", "group", "l")

wilcox_pvals['MLP_a'] <- wilcox.test(df_mlpa_d90$vals_d90[which(df_mlpa_d90$group == "BCG")], df_mlpa_d90$vals_d90[which(df_mlpa_d90$group == "CTL")])$p.value 



# mlpb boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_mlpb <- data.frame(
  val <- mlpb_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
#d90 BCG vs. placebo comparison plots

df_mlpb_d90 <- data.frame(
  vals_d90 <- mlpb_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

df_mlpb_d90 <- subset(df_mlpb_d90, df_mlpb_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_mlpb_d90) <- c("vals_d90", "d", "group", "l")

wilcox_pvals['MLP_b'] <- wilcox.test(df_mlpb_d90$vals_d90[which(df_mlpb_d90$group == "BCG")], df_mlpb_d90$vals_d90[which(df_mlpb_d90$group == "CTL")])$p.value 



# pre-BNK boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_preBNK <- data.frame(
  val <- preBNK_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

tiff(paste0(OUT_DIR, "preBNK_term_state.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_preBNK, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="PreBNK Terminal state prediction") + scale_fill_manual(values=c("#83CAFF", "#004586"))+
  scale_color_manual(values=c( "#DD4477", "black"))
dev.off()


#d90 BCG vs. placebo comparison plots

df_preBNK_d90 <- data.frame(
  vals_d90 <- preBNK_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

df_preBNK_d90 <- subset(df_preBNK_d90, df_preBNK_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_preBNK_d90) <- c("vals_d90", "d", "group", "l")

wilcox_pvals['PreBNK'] <- wilcox.test(df_preBNK_d90$vals_d90[which(df_preBNK_d90$group == "BCG")], df_preBNK_d90$vals_d90[which(df_preBNK_d90$group == "CTL")])$p.value 




## D90 plot for all celltypes comparing BCG and placebo


df_combined <- data.frame(
  values <- c(cmpa_d90, cmpb_d90, cmpc_d90, gmpa_d90, gmpb_d90, mlpa_d90, mlpb_d90,
              mepa_d90, mepb_d90, mepc_d90, preBNK_d90),
  donor <- rep(c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),11),
  celltype <- rep(c("CMP_a", "CMP_b", "CMPc", "GMP_a", "GMP_b", "MLP_a", "MLP_b", "MEP_a", "MEP_b",
                  "MEP_c", "PreBNK"), each=18),
  don_id <- rep(donors,11)
)
colnames(df_combined) <- c("values", "donor", "celltype", "don_id")

df_combined_final <- df_combined[which(df_combined$don_id %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21")),]

myColorScale <- c('BCG'='#EF5350','CTL'='#0288D1')

tiff(paste0(OUT_DIR, "test.tiff"), units="in", width=7, height=2.5, res=350)
ggplot(df_combined_final, aes(x=celltype, y=values, fill=donor)) + geom_boxplot(alpha=0.7) + geom_point(shape=21, aes(fill=donor), position=position_dodge(0.7)) + theme_test() +
  scale_fill_manual(values=myColorScale) + labs(x="", y="Percent of HSCs", fill="") +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=10))
dev.off()


#save same plot as RDS
p <- ggplot(df_combined_final, aes(x=celltype, y=values, fill=donor)) + geom_boxplot(alpha=0.7) + geom_point(shape=21, aes(fill=donor), position=position_dodge(0.7)) + theme_test() +
  scale_fill_manual(values=myColorScale) + labs(x="", y="Percent of HSCs", fill="") +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=10))

saveRDS(p, file=paste0(OUT_DIR,"Figure_3C_cell_fates.rds"))

saveRDS(wilcox_pvals, file=paste0(OUT_DIR,"Figure_3C_wilcox_pvalues.rds"))




# 
# df_combined_broad <- data.frame(
#   values <- c(myeloid_d90, lymphoid_d90, erythroid_d90),
#   donors <- rep(c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),3),
#   celltype <- rep(c("myeloid", "lymphoid", "erythroid"), each=18)
# )
# colnames(df_combined_broad) <- c("values", "donors", "celltype")
# 
# tiff(paste0(OUT_DIR, "test_broad.tiff"), units="in", width=7, height=3, res=250)
# ggplot(df_combined_broad, aes(x=celltype, y=values, fill=donors)) + geom_boxplot()
# dev.off()


## D0 plot for all celltypes comparing BCG and placebo



df_combined <- data.frame(
  values <- c(cmpa_d0, cmpb_d0, cmpc_d0, gmpa_d0, gmpb_d0, mlpa_d0, mlpb_d0,
              mepa_d0, mepb_d0, mepc_d0, preBNK_d0),
  donor <- rep(c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),11),
  celltype <- rep(c("CMP_a", "CMP_b", "CMPc", "GMP_a", "GMP_b", "MLP_a", "MLP_b", "MEP_a", "MEP_b",
                    "MEP_c", "PreBNK"), each=18),
  don_id <- rep(donors,11)
)
colnames(df_combined) <- c("values", "donor", "celltype", "don_id")

df_combined_final <- df_combined[which(df_combined$don_id %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21")),]

myColorScale <- c('BCG'='#EF5350','CTL'='#0288D1')

tiff(paste0(OUT_DIR, "test_d0.tiff"), units="in", width=7, height=2.5, res=350)
ggplot(df_combined_final, aes(x=celltype, y=values, fill=donor)) + geom_boxplot(alpha=0.7) + geom_point(shape=21, aes(fill=donor), position=position_dodge(0.7)) + theme_test() +
  scale_fill_manual(values=myColorScale) + labs(x="", y="Percent of HSCs", fill="") +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=10))
dev.off()


# df_combined_broad <- data.frame(
#   values <- c(myeloid_d0, lymphoid_d0, erythroid_d0),
#   donors <- rep(c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),3),
#   celltype <- rep(c("myeloid", "lymphoid", "erythroid"), each=18)
# )
# colnames(df_combined_broad) <- c("values", "donors", "celltype")
# 
# tiff(paste0(OUT_DIR, "test_d0_broad.tiff"), units="in", width=7, height=3, res=250)
# ggplot(df_combined_broad, aes(x=celltype, y=values, fill=donors)) + geom_boxplot()
# dev.off()


# # Composite plot showing results for all possible terminal fates relative to placebo --------
# 
# 
# df_composite <- data.frame(
#   celltype_name <- c(rep("CMP_b", length(row.names(df_cmpb_new))), rep("CMP_a", length(row.names(df_cmpa_new))), rep("CMP_c", length(row.names(df_cmpc_new))), rep("MEP_a", length(row.names(df_mepa_new))), rep("MLP_b", length(row.names(df_mlpb_new)))),
#   vals <- c(df_cmpb_new$val_corrected, df_cmpa_new$val_corrected, df_cmpc_new$val_corrected, df_mepa_new$val_corrected, df_mlpb_new$val_corrected),
#   group_lab <- c(df_cmpb_new$group, df_cmpa_new$group, df_cmpc_new$group, df_mepa_new$group, df_mlpb_new$group)
# )
# 
# unique_lab <- paste0(df_composite$celltype_name, "_", df_composite$group_lab)
# df_composite$unique_lab <- unique_lab
# 
# df_composite_subset <- subset(df_composite, df_composite$group_lab == "BCG")
# colnames(df_composite_subset) <- c("celltype_name", "vals", "group_lab", "unique_lab")
# 
# ggplot(df_composite_subset, aes(x=celltype_name, y=vals)) +  theme_classic() + geom_boxplot(aes(fill=celltype_name), alpha=0.6, color="black") + geom_jitter(shape=21, aes(fill=celltype_name)) +
#   geom_hline(yintercept=0) + labs(x="") + theme (axis.text.x = element_text (angle = 45, vjust=1, hjust=1))
