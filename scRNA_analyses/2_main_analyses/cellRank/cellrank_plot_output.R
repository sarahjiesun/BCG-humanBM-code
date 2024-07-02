if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pcaMethods")
install_github("velocyto-team/velocyto.R")
install.packages("pagoda2")

library(devtools)
library(velocyto.R)
library(pagoda2)
library(Seurat)
library(hdf5r)
library(SeuratDisk)
library(SingleCellExperiment)
library(ggplot2)


OUT_DIR <- "/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis12_label_transfer_emmreml_edited/cellRank/wVelocity_input/"


trans_prob <- read.csv(file=paste0(OUT_DIR,"S8_Td0_to_terminal_states.csv"), header=TRUE)
coordinates <- read.csv(file=paste0(OUT_DIR, "cell_embeddings_S8_Td0.csv"), header=TRUE)
id <- read.csv(file=paste0(OUT_DIR, "clusters_S8_Td0.csv"), header=TRUE)
velocity_vectors <- read.csv(file=paste0(OUT_DIR,"S8_Td0_velocity_umap_coords.csv"), header=TRUE)

id <- id[which(id$X %in% trans_prob$Cell.ID),]
id <- id[order(match(as.character(id$X), as.character(trans_prob$Cell.ID))),]
which(id$X != trans_prob$Cell.ID)

row.names(coordinates) <- coordinates$X
coordinates <- coordinates[which(coordinates$X %in% trans_prob$Cell.ID),]
coordinates <- coordinates[order(match(as.character(coordinates$X), as.character(trans_prob$Cell.ID))),]
which(coordinates$X != trans_prob$Cell.ID)

velocity_vectors <- velocity_vectors[which(velocity_vectors$Cell.ID %in% trans_prob$Cell.ID),]
velocity_vectors <- velocity_vectors[order(match(as.character(velocity_vectors$Cell.ID), as.character(trans_prob$Cell.ID))),]
which(velocity_vectors$Cell.ID != trans_prob$Cell.ID)

#Plot by percentages for a specific terminal state
df <- data.frame(
  x <- coordinates$UMAP_1,
  y <- coordinates$UMAP_2,
  name <- coordinates$X,
  val <- trans_prob$MEP_a,
  clust <- id$S15_Tm3_obj.meta.data.clust_names
)
colnames(df) <- c("x", "y", "name", "val", "clust")


df_final <- subset(df, df$clust %in% c("HSC_a", "HSC_b"))
hist(df_final$val) 

 
ggplot(df_final, aes(x=x, y=y, fill=val)) + geom_point(aes(fill=val), shape=21, color="black", size=2, alpha=0.8) + theme_light() + 
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                        limits = c(0, 0.3), oob = scales::squish)


#Plot by most likely terminal state

row.names(trans_prob) <- trans_prob$Cell.ID

trans_prob <- subset(trans_prob, select = -X)
trans_prob <- subset(trans_prob, select = -Cell.ID)
most_likely <- colnames(trans_prob)[max.col(trans_prob, ties.method = "first")]

df <- data.frame(
  x <- coordinates$UMAP_1,
  y <- coordinates$UMAP_2,
  name <- coordinates$X,
  term <- most_likely,
  clust <- id$S8_Td0_obj.meta.data.clust_names
)
colnames(df) <- c("x", "y", "name", "term", "clust")

df_final <- subset(df, df$clust %in% c("HSC_a", "HSC_b"))
df_final$term <- factor(df_final$term, unique(df$term))


cols_clust = c('HSC_a','HSC_b','CMP_a','CMP_b','CMP_c','GMP_a', 'GMP_b' , 'MEP_a', 'MEP_b', 'MEP_c','MLP_a', 'MLP_b', 'PreBNK')
cols_val <- c('#EF5350','#CC0C00','#FF7F0E','#FFD147','#CC9900', '#33CC00', '#99CC00', '#8C5782', '#990080', '#BA68C8', '#46B8DA', '#0288D1', '#1A0099')
col <- cols_val[which(cols_clust %in% df_final$term)]

ggplot(df_final, aes(x=x, y=y, fill=term)) + geom_point(aes(fill=term), shape=21, color="black", size=2.5, alpha=0.8) + theme_light() + scale_fill_manual(values=col)

#Plot according to velocity vectors
row.names(trans_prob) <- trans_prob$Cell.ID

trans_prob <- subset(trans_prob, select = -X)
trans_prob <- subset(trans_prob, select = -Cell.ID)
most_likely <- colnames(trans_prob)[max.col(trans_prob, ties.method = "first")]

df <- data.frame(
  x <- velocity_vectors$coord1,
  y <- velocity_vectors$coord2,
  name <- velocity_vectors$Cell.ID,
  term <- most_likely,
  clust <- id$S8_Td0_obj.meta.data.clust_names
)
colnames(df) <- c("x", "y", "name", "term", "clust")


df_final <- subset(df, df$clust %in% c("HSC_a", "HSC_b"))
df_final <- subset(df, df$clust %in% c("MLP_b"))
df_final$term <- factor(df_final$term, unique(df$term))


cols_clust = c('HSC_a','HSC_b','CMP_a','CMP_b','CMP_c','GMP_a', 'GMP_b' , 'MEP_a', 'MEP_b', 'MEP_c','MLP_a', 'MLP_b', 'PreBNK')
cols_val <- c('#EF5350','#CC0C00','#FF7F0E','#FFD147','#CC9900', '#33CC00', '#99CC00', '#8C5782', '#990080', '#BA68C8', '#46B8DA', '#0288D1', '#1A0099')
col <- cols_val[which(cols_clust %in% df_final$term)]

ggplot(df_final, aes(x=x, y=y, fill=term)) + geom_point(aes(fill=term), shape=21, color="black", size=2.5, alpha=0.8) + theme_light() + scale_fill_manual(values=col)+
  geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed")

x_end <- sum(df_final$x)
y_end <- sum(df_final$y)

ggplot(df_final, aes(x=x, y=y, fill=term)) + geom_point(aes(fill=term), shape=21, color="black", size=2.5, alpha=0.8) + theme_light() + scale_fill_manual(values=col)+
  geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") + geom_segment(aes(x = 0,
                                                                                                               y = 0,
                                                                                                               xend = x_end,
                                                                                                               yend = y_end),
                                                                                                           arrow = arrow(length = unit(0.5, "cm")))
