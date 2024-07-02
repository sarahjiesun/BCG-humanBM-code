
#R version 4.1.0 loaded
#rstudio/1.4
#geos/3.9.1 loaded
#hdf5/1.12.0 loaded 
#cmake loaded

.libPaths("/project/lbarreiro/USERS/sarah/Rlibs_new")

library(ggplot2)
library(DESeq2)
library(statmod)
library(RColorBrewer)
library(edgeR)
library(devtools)
library(cluster)
library(ggrepel)
library(EMMREML)
library(ggpointdensity)
library(ashr)
library(mashr)
library(dplyr)
library(plyr)
library(grid)
library(gridExtra)
library(Seurat)


setwd("/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis_Raul")

OUT_DIR <- c("HSC_subcluster_analysis/Hallmark_pathways_dotplot/")
dir.create(OUT_DIR)


# read in gsea results for each cluster 
clusters <- c("c0", "c1", "c2", "c3", "c4", "c5", "c6", 'c7', "c8")
target_pathways_alt <- c("HALLMARK_APOPTOSIS", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                     "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_COMPLEMENT", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_HYPOXIA",
                     "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_CHOLESTEROL_HOMEOSTASIS")


target_pathways <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_APOPTOSIS",
                     "HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_DNA_REPAIR",
                     "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_UV_RESPONSE_UP", "HALLMARK_HEME_METABOLISM", "HALLMARK_KRAS_SIGNALING_DN")

for(i in 1:length(clusters)){
  dat <- read.csv(file=paste0("HSC_subcluster_analysis/9_gsea/NOT_SORTED_gsea_result_",clusters[i],"_no_correlations"))
  pval_values <- dat$pval[which(dat$pathways %in% target_pathways)]
  names(pval_values) <- dat$pathways[which(dat$pathways %in% target_pathways)]
  assign(paste0(clusters[i]), -log10(pval_values))
  NES_values <- dat$NES[which(dat$pathways %in% target_pathways)]
  names(NES_values) <- dat$pathways[which(dat$pathways %in% target_pathways)]
  assign(paste0(clusters[i],"_NES"), NES_values)
}

matrix_pval <- cbind(c0, c1, c2, c3, c4, c5, c6, c7, c8)
matrix_pval_scaled <- data.frame(t(scale(t(matrix_pval))))
matrix_NES <- cbind(c0_NES, c1_NES, c2_NES, c3_NES, c4_NES, c5_NES, c6_NES, c7_NES, c8_NES)
colnames(matrix_NES) <-  c("c0", "c1", "c2", "c3", "c4", "c5", "c6", 'c7', "c8")
matrix_NES_scaled <- scale(t(matrix_NES))


pval <- c(c0, c1, c2, c3, c4, c5, c6, c7, c8)
pval_binary <- rep("not_sig", length(pval))
pval_binary[which(pval>=(-log10(0.05)))] <- "sig"

df <- data.frame(
  pathway <- rep(row.names(matrix_pval),length(clusters)),
  clust_name <- rep(clusters, each=length(row.names(matrix_pval))),
  padj <- c(matrix_pval_scaled$c0, matrix_pval_scaled$c1, matrix_pval_scaled$c2, matrix_pval_scaled$c3, matrix_pval_scaled$c4,
            matrix_pval_scaled$c5, matrix_pval_scaled$c6, matrix_pval_scaled$c7, matrix_pval_scaled$c8),
  NES <- c(c0_NES, c1_NES, c2_NES, c3_NES, c4_NES, c5_NES, c6_NES, c7_NES, c8_NES)
)
colnames(df) <- c("pathway", "clust_name", "pval", "NES")
#df$padj_bin <- factor(df$padj_bin, c("not_sig", "sig"))
df$clust_name <- factor(df$clust_name, c("c0", "c4", "c5", "c7", "c3", "c1", "c6", "c2", "c8"))
## Remove "HALLMARK" from pathway names
df$pathway <- gsub("HALLMARK_", "", df$pathway)
df$pathway <- gsub("_", " ", df$pathway)
df$pathway <- factor(df$pathway, c("OXIDATIVE PHOSPHORYLATION", "TNFA SIGNALING VIA NFKB", "REACTIVE OXYGEN SPECIES PATHWAY", "MYC TARGETS V1", "MYC TARGETS V2", "DNA REPAIR",
                                   "HYPOXIA", "ALLOGRAFT REJECTION", "HEME METABOLISM", "APOPTOSIS", "IL2 STAT5 SIGNALING", 
                                   "EPITHELIAL MESENCHYMAL TRANSITION", "UV RESPONSE UP",  "KRAS SIGNALING DN"))


df$pval[which(df$pval < 1)] <- 0


tiff(paste0(OUT_DIR, "subcluster_hallmark_enrichment_pval.tiff"), units="in", width=4.9, height=3, res=300)
ggplot(df, aes(x=clust_name, y=pathway)) + geom_point(shape=4,stroke=2, aes(color=pval, size=pval, alpha=pval)) + scale_color_gradient(low="white", high="#FF0000")+
  theme_test() + scale_size_continuous(range=c(0.5,3))+ theme(legend.position = "none") + labs(x="") + theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=7))
dev.off()



# ALTERNATIVE KEEP ONLY THE TOP ENRICHMENT

clusters <- c("c0", "c1", "c2", "c3", "c4", "c5", "c6", 'c7', "c8")
target_pathways_alt <- c("HALLMARK_APOPTOSIS", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                         "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_COMPLEMENT", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_HYPOXIA",
                         "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_CHOLESTEROL_HOMEOSTASIS")


target_pathways <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_APOPTOSIS",
                     "HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_DNA_REPAIR",
                     "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_UV_RESPONSE_UP", "HALLMARK_HEME_METABOLISM", "HALLMARK_KRAS_SIGNALING_DN")

for(i in 1:length(clusters)){
  dat <- read.csv(file=paste0("HSC_subcluster_analysis/9_gsea/NOT_SORTED_gsea_result_",clusters[i],"_no_correlations"))
  pval_values <- dat$pval[which(dat$pathways %in% target_pathways)]
  names(pval_values) <- dat$pathways[which(dat$pathways %in% target_pathways)]
  assign(paste0(clusters[i]), -log10(pval_values))
  NES_values <- dat$NES[which(dat$pathways %in% target_pathways)]
  names(NES_values) <- dat$pathways[which(dat$pathways %in% target_pathways)]
  assign(paste0(clusters[i],"_NES"), NES_values)
}

matrix_pval <- cbind(c0, c1, c2, c3, c4, c5, c6, c7, c8)
matrix_pval_scaled <- data.frame(t(scale(t(matrix_pval))))
matrix_NES <- cbind(c0_NES, c1_NES, c2_NES, c3_NES, c4_NES, c5_NES, c6_NES, c7_NES, c8_NES)
colnames(matrix_NES) <-  c("c0", "c1", "c2", "c3", "c4", "c5", "c6", 'c7', "c8")
matrix_NES_scaled <- scale(t(matrix_NES))


pval <- c(c0, c1, c2, c3, c4, c5, c6, c7, c8)
pval_binary <- rep("not_sig", length(pval))
pval_binary[which(pval>=(-log10(0.05)))] <- "sig"

df <- data.frame(
  pathway <- rep(row.names(matrix_pval),length(clusters)),
  clust_name <- rep(clusters, each=length(row.names(matrix_pval))),
  padj <- c(matrix_pval_scaled$c0, matrix_pval_scaled$c1, matrix_pval_scaled$c2, matrix_pval_scaled$c3, matrix_pval_scaled$c4,
            matrix_pval_scaled$c5, matrix_pval_scaled$c6, matrix_pval_scaled$c7, matrix_pval_scaled$c8),
  NES <- c(c0_NES, c1_NES, c2_NES, c3_NES, c4_NES, c5_NES, c6_NES, c7_NES, c8_NES)
)
colnames(df) <- c("pathway", "clust_name", "pval", "NES")
#df$padj_bin <- factor(df$padj_bin, c("not_sig", "sig"))
df$clust_name <- factor(df$clust_name, c("c0", "c4", "c5", "c7", "c3", "c1", "c6", "c2", "c8"))
## Remove "HALLMARK" from pathway names
df$pathway <- gsub("HALLMARK_", "", df$pathway)
df$pathway <- factor(df$pathway, c("OXIDATIVE_PHOSPHORYLATION", "TNFA_SIGNALING_VIA_NFKB", "REACTIVE_OXYGEN_SPECIES_PATHWAY", "MYC_TARGETS_V1", "MYC_TARGETS_V2", "DNA_REPAIR",
                                   "HYPOXIA", "ALLOGRAFT_REJECTION", "HEME_METABOLISM", "APOPTOSIS", "IL2_STAT5_SIGNALING", 
                                   "EPITHELIAL_MESENCHYMAL_TRANSITION", "UV_RESPONSE_UP",  "KRAS_SIGNALING_DN"))


pval_new <- rep(0, length=length(df$pval))
df$p_new <- pval_new

for(i in 1:length(unique(df$pathway))){
  subset <- subset(df, df$pathway == unique(df$pathway)[i])
  max_val <- max(subset$pval)
  df$p_new[which(df$pval == max_val & df$pathway == unique(df$pathway)[i])] <- max_val
}



tiff(paste0(OUT_DIR, "subcluster_hallmark_enrichment_pval_ONE.tiff"), units="in", width=5.2, height=3, res=300)
ggplot(df, aes(x=clust_name, y=pathway)) + geom_point(shape=23,stroke=2, aes(color=p_new, fill=p_new, size=p_new, alpha=p_new)) + scale_color_gradient(low="white", high="#FF0000")+
  scale_fill_gradient(low="white", high="#FF0000")+ theme_test() + scale_size_continuous(range=c(0.5,3))+ theme(legend.position = "none") + labs(x="") + theme(axis.text.x = element_text(size=10))
dev.off()











# ALTERNATIVE NOT SCALED
clusters <- c("c0", "c1", "c2", "c3", "c4", "c5", "c6", 'c7', "c8")
target_pathways_alt <- c("HALLMARK_APOPTOSIS", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                         "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_COMPLEMENT", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_HYPOXIA",
                         "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_CHOLESTEROL_HOMEOSTASIS")


target_pathways <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_APOPTOSIS",
                     "HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_DNA_REPAIR",
                     "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_UV_RESPONSE_UP", "HALLMARK_HEME_METABOLISM", "HALLMARK_KRAS_SIGNALING_DN")

for(i in 1:length(clusters)){
  dat <- read.csv(file=paste0("HSC_subcluster_analysis/9_gsea/NOT_SORTED_gsea_result_",clusters[i],"_no_correlations"))
  pval_values <- dat$pval[which(dat$pathways %in% target_pathways)]
  names(pval_values) <- dat$pathways[which(dat$pathways %in% target_pathways)]
  assign(paste0(clusters[i]), -log10(pval_values))
  NES_values <- dat$NES[which(dat$pathways %in% target_pathways)]
  names(NES_values) <- dat$pathways[which(dat$pathways %in% target_pathways)]
  assign(paste0(clusters[i],"_NES"), NES_values)
}

matrix_pval <- data.frame(cbind(c0, c1, c2, c3, c4, c5, c6, c7, c8))
matrix_NES <- cbind(c0_NES, c1_NES, c2_NES, c3_NES, c4_NES, c5_NES, c6_NES, c7_NES, c8_NES)
colnames(matrix_NES) <-  c("c0", "c1", "c2", "c3", "c4", "c5", "c6", 'c7', "c8")
matrix_NES_scaled <- scale(t(matrix_NES))


pval <- c(c0, c1, c2, c3, c4, c5, c6, c7, c8)
pval_binary <- rep("not_sig", length(pval))
pval_binary[which(pval>=(-log10(0.05)))] <- "sig"

df <- data.frame(
  pathway <- rep(row.names(matrix_pval),length(clusters)),
  clust_name <- rep(clusters, each=length(row.names(matrix_pval))),
  padj <- c(matrix_pval$c0, matrix_pval$c1, matrix_pval$c2, matrix_pval$c3, matrix_pval$c4,
            matrix_pval$c5, matrix_pval$c6, matrix_pval$c7, matrix_pval$c8),
  NES <- c(c0_NES, c1_NES, c2_NES, c3_NES, c4_NES, c5_NES, c6_NES, c7_NES, c8_NES)
)
colnames(df) <- c("pathway", "clust_name", "pval", "NES")
#df$padj_bin <- factor(df$padj_bin, c("not_sig", "sig"))
df$clust_name <- factor(df$clust_name, c("c0", "c4", "c5", "c7", "c3", "c1", "c6", "c2", "c8"))
## Remove "HALLMARK" from pathway names
df$pathway <- gsub("HALLMARK_", "", df$pathway)
df$pathway <- factor(df$pathway, c("OXIDATIVE_PHOSPHORYLATION", "TNFA_SIGNALING_VIA_NFKB", "REACTIVE_OXYGEN_SPECIES_PATHWAY", "MYC_TARGETS_V1", "MYC_TARGETS_V2", "DNA_REPAIR",
                                   "HYPOXIA", "ALLOGRAFT_REJECTION", "HEME_METABOLISM", "APOPTOSIS", "IL2_STAT5_SIGNALING", 
                                   "EPITHELIAL_MESENCHYMAL_TRANSITION", "UV_RESPONSE_UP",  "KRAS_SIGNALING_DN"))



tiff(paste0(OUT_DIR, "subcluster_hallmark_enrichment_pval_NOTSCALED.tiff"), units="in", width=5.2, height=3, res=300)
ggplot(df, aes(x=clust_name, y=pathway)) + geom_point(shape=4,stroke=2, aes(color=pval, size=pval, alpha=pval)) + scale_color_gradient(low="white", high="#FF0000")+
  theme_test() + scale_size_continuous(range=c(0.5,3))+ theme(legend.position = "none") + labs(x="") + theme(axis.text.x = element_text(size=10))
dev.off()
