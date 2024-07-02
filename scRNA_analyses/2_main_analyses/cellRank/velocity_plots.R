library(devtools)
library(velocyto.R)
library(pagoda2)
library(Seurat)
library(hdf5r)
library(SeuratDisk)
library(SingleCellExperiment)
library(ggplot2)
library(ggrepel)


setwd("/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis12_label_transfer_emmreml_edited/") 

IN_DIR <- "cellRank/wVelocity_input/"
IN_DIR2 <- "cellRank/velocity_plots/"
OUT_DIR <- "cellRank/velocity_analysis_output/"
dir.create(OUT_DIR)

donors <- c("S1", "S2", "S3", "S5", "S6", "S7", "S8", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S20", "S21")
cmpa_change <- vector(length=length(donors))
cmpb_change <- vector(length=length(donors))
mlp_change <- vector(length=length(donors))
mep_change <- vector(length=length(donors))
stable_change <- vector(length=length(donors))
myeloid_change <- vector(length=length(donors))

plot_list <- list()
k <- 1
for(i in 1:length(donors))
{
  time <- c("Td0", "Tm3")

    cmpa_perc <- vector(length=length(time))
    cmpb_perc <- vector(length=length(time))
    mlp_perc <- vector(length=length(time))
    mep_perc <- vector(length=length(time))
    stable_perc <- vector(length=length(time))
    myeloid_perc <- vector(length=length(time))
    for(j in 1:length(time))
    {
      coordinates <- read.csv(file=paste0(IN_DIR, "cell_embeddings_", donors[i], "_", time[j],".csv"), header=TRUE)
      id <- read.csv(file=paste0(IN_DIR, "clusters_", donors[i], "_", time[j], ".csv"), header=TRUE)
      velocity_vectors <- read.csv(file=paste0(IN_DIR2, donors[i],"_",time[j],"_velocity_umap_coords.csv"), header=TRUE)
      
      id <- id[which(id$X %in% velocity_vectors$Cell.ID),]
      id <- id[order(match(as.character(id$X), as.character(velocity_vectors$Cell.ID))),]
      which(id$X != velocity_vectors$Cell.ID)
      colnames(id) <- c("X", "clust_names")
      
      row.names(coordinates) <- coordinates$X
      coordinates <- coordinates[which(coordinates$X %in% velocity_vectors$Cell.ID),]
      coordinates <- coordinates[order(match(as.character(coordinates$X), as.character(velocity_vectors$Cell.ID))),]
      which(coordinates$X != velocity_vectors$Cell.ID)
      
      l <- sqrt(velocity_vectors$coord1^2 + velocity_vectors$coord2^2)
      #hist(l)
      col <- rep("no", length(l))
      col[which(l > 0.05)] <- "yes"
      
      x_dir <- vector()
      y_dir <- vector()
      angle <- vector()
      
      x_dir[which(velocity_vectors$coord1 >= 0)] <- "pos"
      x_dir[which(velocity_vectors$coord1 < 0)] <- "neg"
      y_dir[which(velocity_vectors$coord2 >= 0)] <- "pos"
      y_dir[which(velocity_vectors$coord2 < 0)] <- "neg"
      angle <- paste0(x_dir,"_",y_dir)
      angle[which(col=="no")] <- "Stable"
      angle[which(angle == "pos_pos")] <- "MEP"
      angle[which(angle == "neg_pos")] <- "CMP_a"
      angle[which(angle == "pos_neg")] <- "MLP"
      angle[which(angle == "neg_neg")] <- "CMP_b"
      
      df <- data.frame(
        x_vel <- velocity_vectors$coord1,
        y_vel <- velocity_vectors$coord2,
        name <- velocity_vectors$Cell.ID,
        clust <- id$clust_names,
        x <- coordinates$UMAP_1,
        y <- coordinates$UMAP_2,
        c <- angle,
        s <- l
      )
      
      colnames(df) <- c("x_vel", "y_vel", "name", "clust", "x", "y", "Direction", "s")
      df_final <- subset(df, df$clust %in% c("HSC_a", "HSC_b"))
      df_final$Direction <- factor(df_final$Direction, c("Stable", "CMP_a", "MEP", "MLP", "CMP_b"))
      
      # ggplot(df_final, aes(x=x, y=y, fill=c, size=s)) + geom_point(aes(fill=c, size=s), shape=21, color="black",  alpha=0.8) + theme_light() +
      #   geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") + scale_size_continuous(range=c(1,5)) +
      #   scale_fill_manual(values=c("grey","#FF6F00", "#8A4198", "#71D0F5", "#FFCD00"))
      
      myColorScale <- c("grey","#FF6F00", "#8A4198", "#71D0F5", "#FFCD00")
      names(myColorScale) <- c("Stable", "CMP_a", "MEP", "MLP", "CMP_b")


      p <- ggplot(df_final, aes(x=x, y=y, color=Direction)) + geom_point(aes(color=Direction),shape=16, size=2, alpha=0.7, stroke=0.3) + theme_light() +
        geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") +
        scale_color_manual(values=myColorScale) + labs(title=paste0(donors[i]," ",time[j]," Velocity UMAP"), x= "UMAP1", y="UMAP2")
      plot_list[[k]] <- p
      k <- k + 1

    
      cmpa_perc[j] <- length(df_final$Direction[which(df_final$Direction == "CMP_a")])/length(df_final$Direction)
      cmpb_perc[j] <- length(df_final$Direction[which(df_final$Direction == "CMP_b")])/length(df_final$Direction)
      mlp_perc[j] <- length(df_final$Direction[which(df_final$Direction == "MLP")])/length(df_final$Direction)
      mep_perc[j] <- length(df_final$Direction[which(df_final$Direction == "MEP")])/length(df_final$Direction)
      stable_perc[j] <- length(df_final$Direction[which(df_final$Direction == "Stable")])/length(df_final$Direction)  
      myeloid_perc[j] <- cmpa_perc[j] + cmpb_perc[j]
    }
    cmpa_change[i] <- cmpa_perc[2] - cmpa_perc[1]
    cmpb_change[i] <- cmpb_perc[2] - cmpb_perc[1]
    mlp_change[i] <- mlp_perc[2] - mlp_perc[1]
    mep_change[i] <- mep_perc[2] - mep_perc[1]
    stable_change[i] <- stable_perc[2] - stable_perc[1]
    myeloid_change[i] <- myeloid_perc[2] - myeloid_perc[1]
}

donors <- c("S1", "S2", "S3", "S5", "S6", "S7", "S8", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S20", "S21")
time <- c("Td0", "Tm3")
k <- 1
for(i in 1:length(donors)) 
{
  for(j in 1:length(time))
  {
    file_name = paste0(OUT_DIR, donors[i],"_",time[j], "_velocity_umap.tiff")
    tiff(file_name, units="in", width=3.5, height=3, res=250)
    print(plot_list[[k]])
    dev.off()
    print(k)
    k <- k+1
  }

}





#Proportion with CMP_a, CMP_b, or MLP direction
lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_cmpb <- data.frame(
  val <- cmpb_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)

tiff(paste0(OUT_DIR, "cmp_b_velocity_change.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpb, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_b directed velocity vectors") + scale_fill_manual(values=c("#FFD320", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black"))
dev.off()
tiff(paste0(OUT_DIR, "cmp_b_velocity_change_labels.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpb, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_b directed velocity vectors") + scale_fill_manual(values=c("#FFD320", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()



df_cmpa <- data.frame(
  val <- cmpa_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)

tiff(paste0(OUT_DIR, "cmp_a_velocity_change.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpa, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, color="black",alpha=1, position=position_jitter(0.2), aes(fill=group)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_a directed velocity vectors") + scale_fill_manual(values=c("#FF420E", "#004586"))
dev.off()
tiff(paste0(OUT_DIR, "cmp_a_velocity_change_labels.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpa, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, color="black",alpha=1, position=position_jitter(0.2), aes(fill=group)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_a directed velocity vectors") + scale_fill_manual(values=c("#FF420E", "#004586"))+
  geom_text_repel(size=2)
dev.off()




df_mlp <- data.frame(
  val <- mlp_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)
tiff(paste0(OUT_DIR, "mlp_velocity_change.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_mlp, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, color="black",alpha=1, position=position_jitter(0.2), aes(fill=group)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="MLP directed velocity vectors") + scale_fill_manual(values=c("#0084D1", "#004586"))
dev.off()
tiff(paste0(OUT_DIR, "mlp_velocity_change_labels.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_mlp, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, color="black",alpha=1, position=position_jitter(0.2), aes(fill=group)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="MLP directed velocity vectors") + scale_fill_manual(values=c("#0084D1", "#004586")) +
  geom_text_repel(size=2)
dev.off()





lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_mep <- data.frame(
  val <- mep_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
tiff(paste0(OUT_DIR, "mep_velocity_change.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_mep, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3,alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group, color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="MEP directed velocity vectors") + scale_fill_manual(values=c("#4B1F6F", "#004586"))+
  scale_color_manual(values=c( "#7E0021", "black"))
dev.off()
tiff(paste0(OUT_DIR, "mep_velocity_change_label.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_mep, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3,alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group, color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="MEP directed velocity vectors") + scale_fill_manual(values=c("#4B1F6F", "#004586"))+
  scale_color_manual(values=c( "#7E0021", "black")) + geom_text_repel(size=2)
dev.off()





df_stable <- data.frame(
  val <- stable_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)
tiff( paste0(OUT_DIR, "stable_velocity_change.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_stable, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, color="black",alpha=1, position=position_jitter(0.2), aes(fill=group)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="Stable velocity vectors") + scale_fill_manual(values=c("grey", "#004586"))
dev.off()
tiff( paste0(OUT_DIR, "stable_velocity_change_label.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_stable, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, color="black",alpha=1, position=position_jitter(0.2), aes(fill=group)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="Stable velocity vectors") + scale_fill_manual(values=c("grey", "#004586"))+
  geom_text_repel(size=2)
dev.off()






lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_myeloid <- data.frame(
  val <- myeloid_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
tiff( paste0(OUT_DIR, "myeloid_velocity_change.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_stable, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3,alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group, color=lab)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP velocity vectors") + scale_fill_manual(values=c("#B82E2E", "#004586")) +
  scale_color_manual(values=c( "#FFD320", "black"))
dev.off()
tiff( paste0(OUT_DIR, "myeloid_velocity_change_label.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_stable, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3,alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group, color=lab)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP velocity vectors") + scale_fill_manual(values=c("#B82E2E", "#004586")) +
  scale_color_manual(values=c( "#FFD320", "black")) + geom_text_repel(size=2)
dev.off()




#Proportion no S2, S3, S12
lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_cmpb <- data.frame(
  val <- cmpb_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
df_cmpb <- subset(df_cmpb, df_cmpb$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_cmpb) <- c("val", "d", "group", "l")

tiff(paste0(OUT_DIR, "cmp_b_velocity_change_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpb, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_b directed velocity vectors") + scale_fill_manual(values=c("#FFD320", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()












lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_mep <- data.frame(
  val <- mep_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
df_mep <- subset(df_mep, df_mep$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_mep) <- c("val", "d", "group", "l")


tiff(paste0(OUT_DIR, "mep_velocity_change_label_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_mep, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3,alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group, color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="MEP directed velocity vectors") + scale_fill_manual(values=c("#4B1F6F", "#004586"))+
  scale_color_manual(values=c( "#7E0021", "black")) + geom_text_repel(size=2)
dev.off()







lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_myeloid <- data.frame(
  val <- myeloid_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
df_myeloid <- subset(df_myeloid, df_myeloid$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_myeloid) <- c("val", "d", "group", "l")

tiff( paste0(OUT_DIR, "myeloid_velocity_change_label_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_myeloid, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3,alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group, color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP velocity vectors") + scale_fill_manual(values=c("#B82E2E", "#004586")) +
  scale_color_manual(values=c( "#FFD320", "black")) + geom_text_repel(size=2)
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

wilcox.test(df_myeloid_new$val_corrected[which(df_myeloid_new$group == "BCG")], df_myeloid_new$val_corrected[which(df_myeloid_new$group == "CTL")])

ggplot(data=df_myeloid_new, aes(x=group, y=val_corrected, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3,alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group, color=l)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0) relative to placebo", fill="Group", title="CMP velocity vectors") + scale_fill_manual(values=c("#B82E2E", "#004586")) +
  scale_color_manual(values=c( "#FFD320", "black")) + geom_text_repel(size=2)
