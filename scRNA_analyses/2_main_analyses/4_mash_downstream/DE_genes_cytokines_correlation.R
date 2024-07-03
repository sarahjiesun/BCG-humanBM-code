#R version 4.1.0
install.packages("VennDiagram")
library(ggplot2)
library(cowplot)
library(reticulate)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(Seurat)
library(sctransform)
library(data.table)
library(VennDiagram)

# setup ---------------------------------------------------------
setwd("/scRNA_analyses/2_main_analyses/4_mash_downstream/")
OUT_DIR <- "DE_genes_cytokines_correlation/"
dir.create(OUT_DIR)


# Read in cytokine data ---------------------------------------------------

cytokine_data <- data.frame(fread("../PBMC_cytokine_data_CHM.txt"))

# Read in the corrected expression matrices -------------------------------

celltypes <- c("HSC_a", "HSC_b", "CMP_a", "CMP_b", "CMP_c", "GMP_b1", "GMP_b2", "MEP_a1", "MEP_a2", "MEP_a3", "mixed_a", "mixed_b", "PreBNK")
data_genes <- vector()
for(i in 1:length(celltypes)){
  c <- celltypes[i]
  
  data <- read.table(file=paste0("emmreml_betas/raw_expression/", c, "_raw_expression.txt"))
  data$X <- row.names(data)
  #Get genes that are significant MASH no correlations with custom matrices
  PM <- read.table(file=paste0("MASH/mash_results/posteriorMeans_wcustom.txt"), header=TRUE)
  lfsr <- read.table(file=paste0("MASH/mash_results/lfsr_wcustom_output.txt"), header=TRUE)
  
  DE_df <- data.frame(cbind(PM[,paste0(celltypes[i])], lfsr[,paste0(celltypes[i])]))
  row.names(DE_df) <- row.names(PM)
  colnames(DE_df) <- c("beta", "lfsr")
  
  DE_subset_genes <- row.names(subset(DE_df, DE_df$lfsr < 0.1))
  
  data_subset <- subset(data, data$X %in% DE_subset_genes)
  
  assign(paste0(c), data_subset)
}

expression_data <- list(HSC_a, HSC_b, CMP_a, CMP_b, CMP_c, GMP_b1, GMP_b2, MEP_a1, MEP_a2, MEP_a3, mixed_a, mixed_b, PreBNK)
expression_data_names <- c("HSC_a", "HSC_b", "CMP_a", "CMP_b", "CMP_c", "GMP_b1", "GMP_b2", "MEP_a1", "MEP_a2", "MEP_a3", "mixed_a", "mixed_b","PreBNK")


meta_data <- read.csv(file="../Analysis12_label_transfer_emmreml_edited/all_samples_meta_data.csv")
meta_data$Sample <- gsub('-','\\.', meta_data$Sample)


Tm3_all <- c("SNG.LB.SS.1S.RS.S1.CD34neg_S1_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S2.CD34neg_S2_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S5.CD34neg_S4_R1_001_Tm3_BCG",
             "SNG.LB.SS.1S.RS.S6.CD34neg_S5_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S8.CD34neg_S7_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S9.CD34neg_S8_R1_001_Tm3_BCG",
             "SNG.LB.SS.1S.RS.S10.CD34neg_S9_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S12.CD34neg_S11_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S14.CD34neg_S13_R1_001_Tm3_BCG",
             "SNG.LB.SS.1S.RS.S15.CD34neg_S14_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S16.CD34neg_S15_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S17.CD34neg_S16_R1_001_Tm3_BCG",
             "SNG.LB.SS.1S.RS.S19.CD34neg_S18_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S20.CD34neg_S19_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S21.CD34neg_S20_R1_001_Tm3_BCG", 
             "SNG.LB.SS.1S.RS.S3.CD34neg_S3_R1_001_Tm3_CTL", "SNG.LB.SS.1S.RS.S7.CD34neg_S6_R1_001_Tm3_CTL", "SNG.LB.SS.1S.RS.S11.CD34neg_S10_R1_001_Tm3_CTL",
             "SNG.LB.SS.1S.RS.S13.CD34neg_S12_R1_001_Tm3_CTL", "SNG.LB.SS.1S.RS.S18.CD34neg_S17_R1_001_Tm3_CTL")


Td0_all <- c("SNG.LB.SS.1S.RS.S1.CD34neg_S1_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S2.CD34neg_S2_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S5.CD34neg_S4_R1_001_Td0_BCG",
             "SNG.LB.SS.1S.RS.S6.CD34neg_S5_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S8.CD34neg_S7_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S9.CD34neg_S8_R1_001_Td0_BCG",
             "SNG.LB.SS.1S.RS.S10.CD34neg_S9_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S12.CD34neg_S11_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S14.CD34neg_S13_R1_001_Td0_BCG",
             "SNG.LB.SS.1S.RS.S15.CD34neg_S14_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S16.CD34neg_S15_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S17.CD34neg_S16_R1_001_Td0_BCG",
             "SNG.LB.SS.1S.RS.S19.CD34neg_S18_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S20.CD34neg_S19_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S21.CD34neg_S20_R1_001_Td0_BCG", 
             "SNG.LB.SS.1S.RS.S3.CD34neg_S3_R1_001_Td0_CTL", "SNG.LB.SS.1S.RS.S7.CD34neg_S6_R1_001_Td0_CTL", "SNG.LB.SS.1S.RS.S11.CD34neg_S10_R1_001_Td0_CTL",
             "SNG.LB.SS.1S.RS.S13.CD34neg_S12_R1_001_Td0_CTL", "SNG.LB.SS.1S.RS.S18.CD34neg_S17_R1_001_Td0_CTL")

logfc_mat_final_list <- list()
for(k in 1:length(expression_data)){
  #split each corrected expression matrix into 2 parts Td0 and Tm3
  celltype <- celltypes[k]
  mat <- expression_data[[k]]
  Td0_cols <- which(grepl("Td0", colnames(mat))==TRUE)
  Td0_subset <- mat[Td0_cols]
  row.names(Td0_subset) <- mat$X
  
  Tm3_cols <- which(grepl("Tm3", colnames(mat))==TRUE)
  Tm3_subset <- mat[Tm3_cols]
  row.names(Tm3_subset) <- mat$X
  
  #for each Tm3 sample get the corresponsing Td0 sample and subtract columns 
  
  logfc_mat <- data.frame()
  count <- 0
  name_short <- vector()
  name_count <- 1
  for(i in 1:length(Tm3_all))
  {
    name <- Tm3_all[i]
    if(name %in% colnames(Tm3_subset))
    {
      d <- name
      Td0_pair_name <- Td0_all[which(Tm3_all == d)]
      Td0_col <- Td0_subset[Td0_pair_name]
      Tm3_col <- Tm3_subset[d]
      log2fc <- Tm3_col - Td0_col
      name_short[name_count] <- subset(meta_data, meta_data$Sample %in% as.character(d))$donor
      name_count <- name_count+1
      if(count==0)
      {
        logfc_mat <- log2fc[1]
      }
      else
      {
        logfc_mat <- cbind(logfc_mat, log2fc[1])
      }
      
      count <- count + 1
    }
    
  }
  colnames(logfc_mat) <- name_short
  logfc_mat_final <- t(logfc_mat)
  logfc_mat_final_list[[k]] <- logfc_mat_final
  names(logfc_mat_final_list)[k] <- expression_data_names[k]
  
  
  il10 <- vector(length=length(name_short))
  il1b <- vector(length=length(name_short))
  il1ra <- vector(length=length(name_short))
  il6 <- vector(length=length(name_short))
  ifna <- vector(length=length(name_short))
  ifng <- vector(length=length(name_short))
  tnfa <- vector(length=length(name_short))
  
  
  for(i in 1:length(name_short))
  {
    cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) == name_short[i])
    cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
    il10[i] <- as.numeric(cyt_subset2$IL10_FC)
    il1b[i] <- as.numeric(cyt_subset2$IL1B_FC)
    il1ra[i] <- as.numeric(cyt_subset2$IL1RA_FC)
    il6[i] <- as.numeric(cyt_subset2$IL6_FC)
    ifna[i] <- as.numeric(cyt_subset2$IFNA_FC)
    ifng[i] <- as.numeric(cyt_subset2$IFNG_FC)
    tnfa[i] <- as.numeric(cyt_subset2$TNF_FC)
    print(i)
  }
  
  cyt_mat <- cbind(il10, il1b, il1ra, il6, ifna, ifng, tnfa)
  row.names(cyt_mat) <- name_short
  
  count <- 1
  gene_name <- vector()
  cyt_name <- vector()
  pval <- vector()
  rho <- vector()
  for(i in 1:length(colnames(logfc_mat_final)))
  {
    gene_vals <- logfc_mat_final[,i]
    for(j in 1:length(colnames(cyt_mat)))
    {
      gene_name[count] <- colnames(logfc_mat_final)[i]
      cyt_vals <- cyt_mat[,j]
      cyt_name[count] <- colnames(cyt_mat)[j]
      val <- cor.test(cyt_vals, gene_vals, method="spearman")
      pval[count] <- val$p.value
      rho[count] <- val$estimate
      
      count <- count + 1
    }
  }
  
  df <- data.frame(cbind(gene_name, cyt_name, pval, rho))
  df_list <- split(df, f = df$cyt_name)
  for(j in 1:length(df_list))
  {
    cyt_name <- names(df_list)[j]
    dat_cyt <- df_list[[j]]
    padj <- p.adjust(dat_cyt[,3], method="BH")
    dat_cyt <- cbind(dat_cyt, padj)
    assign(paste0(celltype, "_",cyt_name, "_result"), dat_cyt)
  }
  print(k)
}
saveRDS(logfc_mat_final_list, file=paste0(OUT_DIR, "logfc_mat_final_list_DR_only.rds"))



ifna_result_list <- list(HSC_a_ifna_result, HSC_b_ifna_result, CMP_a_ifna_result, CMP_b_ifna_result, CMP_c_ifna_result, GMP_b1_ifna_result, GMP_b2_ifna_result, MEP_a1_ifna_result, MEP_a2_ifna_result, MEP_a3_ifna_result, mixed_a_ifna_result, mixed_b_ifna_result, PreBNK_ifna_result)
saveRDS(ifna_result_list, file=paste0(OUT_DIR, "ifna_result_list_DE_genes_only.rds"))
ifng_result_list <- list(HSC_a_ifng_result, HSC_b_ifng_result, CMP_a_ifng_result, CMP_b_ifng_result, CMP_c_ifng_result, GMP_b1_ifng_result, GMP_b2_ifng_result, MEP_a1_ifng_result, MEP_a2_ifng_result, MEP_a3_ifng_result, mixed_a_ifng_result, mixed_b_ifng_result, PreBNK_ifng_result)
saveRDS(ifng_result_list, file=paste0(OUT_DIR, "ifng_result_list_DE_genes_only.rds"))
il10_result_list <- list(HSC_a_il10_result, HSC_b_il10_result, CMP_a_il10_result, CMP_b_il10_result, CMP_c_il10_result, GMP_b1_il10_result, GMP_b2_il10_result, MEP_a1_il10_result, MEP_a2_il10_result, MEP_a3_il10_result, mixed_a_il10_result, mixed_b_il10_result, PreBNK_il10_result)
saveRDS(il10_result_list, file=paste0(OUT_DIR, "il10_result_list_DE_genes_only.rds"))
il1b_result_list <- list(HSC_a_il1b_result, HSC_b_il1b_result, CMP_a_il1b_result, CMP_b_il1b_result, CMP_c_il1b_result, GMP_b1_il1b_result, GMP_b2_il1b_result, MEP_a1_il1b_result, MEP_a2_il1b_result, MEP_a3_il1b_result, mixed_a_il1b_result, mixed_b_il1b_result, PreBNK_il1b_result)
saveRDS(il1b_result_list, file=paste0(OUT_DIR, "il1b_result_list_DE_genes_only.rds"))
il1ra_result_list <- list(HSC_a_il1ra_result, HSC_b_il1ra_result, CMP_a_il1ra_result, CMP_b_il1ra_result, CMP_c_il1ra_result, GMP_b1_il1ra_result, GMP_b2_il1ra_result, MEP_a1_il1ra_result, MEP_a2_il1ra_result, MEP_a3_il1ra_result, mixed_a_il1ra_result, mixed_b_il1ra_result, PreBNK_il1ra_result)
saveRDS(il1ra_result_list, file=paste0(OUT_DIR, "il1ra_result_list_DE_genes_only.rds"))
il6_result_list <- list(HSC_a_il6_result, HSC_b_il6_result, CMP_a_il6_result, CMP_b_il6_result, CMP_c_il6_result, GMP_b1_il6_result, GMP_b2_il6_result, MEP_a1_il6_result, MEP_a2_il6_result, MEP_a3_il6_result, mixed_a_il6_result, mixed_b_il6_result, PreBNK_il6_result)
saveRDS(il6_result_list, file=paste0(OUT_DIR, "il6_result_list_DE_genes_only.rds"))
tnfa_result_list <- list(HSC_a_tnfa_result, HSC_b_tnfa_result, CMP_a_tnfa_result, CMP_b_tnfa_result, CMP_c_tnfa_result, GMP_b1_tnfa_result, GMP_b2_tnfa_result, MEP_a1_tnfa_result, MEP_a2_tnfa_result, MEP_a3_tnfa_result, mixed_a_tnfa_result, mixed_b_tnfa_result, PreBNK_tnfa_result)
saveRDS(tnfa_result_list, file=paste0(OUT_DIR, "tnfa_result_list_DE_genes_only.rds"))



# Bargraph summarizing all results ----------------------------------------

celltypes <- c("HSC_a", "HSC_b", "CMP_a", "CMP_b", "CMP_c", "GMP_b1", "GMP_b2", "MEP_a1", "MEP_a2", "MEP_a3", "mixed_a", "mixed_b", "PreBNK")
cytokines <- c("ifna", "ifng", "il10", "il1b", "il1ra", "il6", "tnfa")


##ifna
num_genes_ifna <- vector(length=length(celltypes))
for(j in 1:length(celltypes))
{
  dat <- data.frame(ifna_result_list[[j]])
  num_genes_ifna[j] <- length(row.names(subset(dat, dat$padj < 0.1)))
}

##ifng
num_genes_ifng <- vector(length=length(celltypes))
for(j in 1:length(celltypes))
{
  dat <- data.frame(ifng_result_list[[j]])
  num_genes_ifng[j] <- length(row.names(subset(dat, dat$padj < 0.1)))
}

##il10
num_genes_il10 <- vector(length=length(celltypes))
for(j in 1:length(celltypes))
{
  dat <- data.frame(il10_result_list[[j]])
  num_genes_il10[j] <- length(row.names(subset(dat, dat$padj < 0.1)))
}

##il1b
num_genes_il1b <- vector(length=length(celltypes))
for(j in 1:length(celltypes))
{
  dat <- data.frame(il1b_result_list[[j]])
  num_genes_il1b[j] <- length(row.names(subset(dat, dat$padj < 0.1)))
}

##il1ra
num_genes_il1ra <- vector(length=length(celltypes))
for(j in 1:length(celltypes))
{
  dat <- data.frame(il1ra_result_list[[j]])
  num_genes_il1ra[j] <- length(row.names(subset(dat, dat$padj < 0.1)))
}

##il6
num_genes_il6 <- vector(length=length(celltypes))
for(j in 1:length(celltypes))
{
  dat <- data.frame(il6_result_list[[j]])
  num_genes_il6[j] <- length(row.names(subset(dat, dat$padj < 0.1)))
}

##tnfa
num_genes_tnfa <- vector(length=length(celltypes))
for(j in 1:length(celltypes))
{
  dat <- data.frame(tnfa_result_list[[j]])
  num_genes_tnfa[j] <- length(row.names(subset(dat, dat$padj < 0.1)))
}

df_summary <- data.frame(
  values <- c(num_genes_ifna, num_genes_ifng, num_genes_il10, num_genes_il1b, num_genes_il1ra, num_genes_il6, num_genes_tnfa),
  cytokine_name <- c(rep("ifna", length(celltypes)), rep("ifng", length(celltypes)), rep("il10", length(celltypes)), rep("il1b", length(celltypes)), rep("il1ra", length(celltypes)), rep("il6", length(celltypes)), rep("tnfa", length(celltypes))),
  celltype_name <- c(rep(c("HSC c1", "HSC c2", "CMP c1", "CMP c2", "CMP c3", "GMP c1", "GMP c2", "MEP c1", "MEP c2", "MEP c3", "MLP c1", "MLP c2", "PreBNK"),7))
)
colnames(df_summary) <- c("values", "cytokine_name", "celltype_name")
df_summary$celltype_name <- factor(df_summary$celltype_name, c("HSC c1", "HSC c2", "CMP c1", "CMP c2", "CMP c3", "GMP c1", "GMP c2", "MEP c1", "MEP c2", "MEP c3", "MLP c1", "MLP c2", "PreBNK"))

#tiff(paste0(OUT_DIR, "summary_df_DE_genes_only.tiff"), units="in", width=4, height=3, res=350)
p <- ggplot(data=df_summary, aes(x=celltype_name, y=values, fill=cytokine_name)) +
  geom_bar(stat="identity", size=0.5)+theme_linedraw()+ylab("Correlated genes (padj < 0.1)")+xlab("")+ labs(fill="")+
  theme(axis.text.x = element_text( angle=45, hjust=1, size=10)) +  theme(axis.title.y = element_text(size=12)) + scale_fill_manual(values=c("#00468B", "#ED0000", "#42B540","#0099B4", "#925E9F", "#FDAF91", "#AD002A")) +
  theme (legend.key.size = unit (0.5, 'cm'))
#dev.off()

final_file <- file.path(OUT_DIR,paste0("summary_df_DE_genes_only.pdf"))
ggsave(filename = final_file, plot = p, height = 3, width = 4)



# Code to make scatterplots for individual gene-cytokine pairs ------------

celltype <- "HSC_a" ###CHOOSE celltype
mat <- expression_data[[1]] ##CHOOSE data
Td0_cols <- which(grepl("Td0", colnames(mat))==TRUE)
Td0_subset <- mat[Td0_cols]
row.names(Td0_subset) <- mat$X

Tm3_cols <- which(grepl("Tm3", colnames(mat))==TRUE)
Tm3_subset <- mat[Tm3_cols]
row.names(Tm3_subset) <- mat$X

#for each Tm3 sample get the corresponsing Td0 sample and subtract columns 

logfc_mat <- data.frame()
count <- 0
name_short <- vector()
name_count <- 1
for(i in 1:length(Tm3_all)){
  name <- Tm3_all[i]
  if(name %in% colnames(Tm3_subset))
  {
    d <- name
    Td0_pair_name <- Td0_all[which(Tm3_all == d)]
    Td0_col <- Td0_subset[Td0_pair_name]
    Tm3_col <- Tm3_subset[d]
    log2fc <- Tm3_col - Td0_col
    name_short[name_count] <- subset(meta_data, meta_data$Sample %in% as.character(d))$donor
    name_count <- name_count+1
    if(count==0)
    {
      logfc_mat <- log2fc[1]
    }
    else
    {
      logfc_mat <- cbind(logfc_mat, log2fc[1])
    }
    
    count <- count + 1
  }
  
}
colnames(logfc_mat) <- name_short
logfc_mat_final <- t(logfc_mat)


gene_vec <- logfc_mat_final[,which(colnames(logfc_mat_final) == "FOSB")] ##CHOOSE GENE NAME

cytokine <- vector(length=length(name_short))
treatment <- vector(length=length(name_short))
meta_data <- read.csv(file="../Analysis12_label_transfer_emmreml_edited/all_samples_meta_data.csv")
meta_data$Sample <- gsub('-','\\.', meta_data$Sample)
for(i in 1:length(name_short)){
  cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) == name_short[i])
  cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
  cytokine[i] <- as.numeric(cyt_subset2$IL1B_FC)  ##CHOOSE CYTOKINE
  
  #GET BCG/CTL info for each donor in name_short
  treatment[i] <- unique(meta_data$vaccination[which(meta_data$donor == name_short[i])])
  print(i)
}

val <- cor.test(cytokine, gene_vec, method="spearman")

df <- data.frame(
  g <- gene_vec,
  c <- cytokine,
  t <- treatment
)
colnames(df) <- c("g", "c", "Group")

#tiff(paste0(OUT_DIR, "HSC_a_KLF6_IL1B.tiff"), units="in", width=4, height=3.2, res=350)
p <- ggplot(df, aes(x=g, y=c)) + geom_point(shape=21, size=4,aes(fill=Group)) + theme_linedraw() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab("Log2FC (Tm3 vs Td0)") +
  ylab("Fold Change IL1B") + ggtitle("FOSB in HSC c1")+
  scale_fill_manual(values=c("#D91B60", "grey")) + ##EDIT GRAPH LABS
  theme(axis.title = element_text(size=14)) + theme(plot.title=element_text(size=15))
#dev.off()


final_file <- file.path(OUT_DIR,paste0("HSC_a_FOSB_IL1B.pdf"))
ggsave(filename = final_file, plot = p, height = 3.2, width = 4)



# correlation of IL1B and IL6 production by PBMCs -------------------------
cytokine_data <- data.frame(fread("PBMC_cytokine_data_CHM.txt"))
meta_data <- read.csv(file="all_samples_meta_data.csv")
donors <- unique(cytokine_data$Donor)
il1b <- vector(length=length(donors))
il6 <- vector(length=length(donors))
treatment <- vector(length=length(donors))

for(i in 1:length(donors))
{
  cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) == donors[i])
  cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
  il1b[i] <- as.numeric(cyt_subset2$IL1B_FC)
  il6[i] <- as.numeric(cyt_subset2$IL6_FC)
  treatment[i] <- unique(meta_data$vaccination[which(meta_data$donor == donors[i])])
  print(i)
}

val <- cor.test(il1b, il6, method="spearman")
pval <- val$p.value

df <- data.frame(
  IL1B <- il1b,
  IL6 <- il6,
  t <- treatment,
  lab <- donors
)
row.names(df) <- donors
colnames(df) <- c("IL1B", "IL6", "Group", "Label")

tiff(paste0(OUT_DIR, "IL6_IL1B_corr.tiff"), units="in", width=5, height=4, res=250)
ggplot(df, aes(x=IL1B, y=IL6)) + geom_point(shape=21, size=4,aes(fill=Group)) + theme_test() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
  geom_vline(xintercept = 1, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab("Fold Change IL1B production of PBMCs") +
  ylab("Fold Change IL6 production of PBMCs") + ggtitle("IL1B vs. IL6")+
  annotate("text", x=0.9, y=2.9, label= paste0("Spearman corr = ",round(val$estimate,2),"\npvalue = ",round(pval,2))) + scale_fill_manual(values=c("#D91B60", "#1E88E5")) +
  geom_text_repel(aes(label=Label))
dev.off()



#gsea on IL1B correlated genes in HSCa and HSCb
library ("DESeq2")
library(statmod)
library(goseq)
library("RColorBrewer")
library(edgeR)
library(devtools)
library(cluster)
library(fgsea)
library(msigdbr)

set.seed(831)

clusters_gsea <- c("HSC_a", "HSC_b")
for(i in 1:length(clusters_gsea))
{
  name <- clusters_gsea[i]
  dat <- il1b_result_list[[i]]
  betas <- as.numeric(dat$rho)
  lfsr_vals <- as.numeric(dat$pval)
  rank_vals <- -log10(lfsr_vals)*betas
  names(rank_vals) <- dat$gene_name
  sorted_ranked_list <- sort(rank_vals, decreasing=FALSE)
  assign(paste0(name,"_table"), sorted_ranked_list)
  print(i)
}


df_list <- list(HSC_a_table, HSC_b_table)


# Gene set enrichment of Hallmark pathways --------------------------------

for(i in 1:length(df_list))
{
  name <- clusters_gsea[i]
  sorted_ranked_list <- df_list[[i]]
  msigdbr_show_species()
  h_df = msigdbr(species = "Homo sapiens", category = "H")  #Hallmark categories n=50
  h_list = h_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  fgseaRes <- fgsea(pathways = h_list, stats= sorted_ranked_list,
                    maxSize  = 500, nperm=100000)
  gsea_result_df<- data.frame(
    pathways <- fgseaRes$pathway,
    pval <- fgseaRes$pval,
    NES <- fgseaRes$NES
  )
  colnames(gsea_result_df) <- c("pathways","pval","NES")
  write.csv(gsea_result_df, file=paste0(OUT_DIR,"NOT_SORTED_gsea_result_",name))
  gsea_result_df <- gsea_result_df[order(gsea_result_df$pval), ]
  assign(paste0(name,"_gsea_result"), gsea_result_df)
  write.csv(gsea_result_df, file=paste0(OUT_DIR,"SORTED_gsea_result_",name))
  print(i)
}

#gsea result as a barchart

hsc_a_gsea <- read.csv(file=paste0("emmreml_downstream/DE_genes_cytokines_correlation/SORTED_gsea_result_HSC_a"))
hsc_b_gsea <- read.csv(file=paste0("emmreml_downstream/DE_genes_cytokines_correlation/SORTED_gsea_result_HSC_b"))

hsc_a_gsea_subset <- subset(hsc_a_gsea, hsc_a_gsea$pval <= 0.05)
hsc_b_gsea_subset <- subset(hsc_b_gsea, hsc_b_gsea$pval <= 0.05)

df_gsea <- data.frame(
  pathway_name <- c(hsc_a_gsea_subset$pathways, hsc_b_gsea_subset$pathways),
  p <- c(-log10(hsc_a_gsea_subset$pval), -log10(hsc_b_gsea_subset$pval)),
  cell <- c(rep("HSC_a", length(hsc_a_gsea_subset$pathways)), rep("HSC_b", length(hsc_b_gsea_subset$pathways)))
)
colnames(df_gsea) <- c("pathway_name", "p", "cell")

ggplot(data=df_gsea, aes(x=pathway_name, y=p, fill=cell)) +
  geom_bar(stat="identity", position=position_dodge()) + theme_test() + coord_flip() + labs(x="Pathway")

immune <- rep("no", length(hsc_a_gsea_subset$pathways))
immune[which(grepl("TNF", hsc_a_gsea_subset$pathways) == TRUE)] <- "yes"
immune[which(grepl("INTERFERON", hsc_a_gsea_subset$pathways) == TRUE)] <- "yes"
immune[which(grepl("TGF", hsc_a_gsea_subset$pathways) == TRUE)] <- "yes"
immune[which(grepl("IL2", hsc_a_gsea_subset$pathways) == TRUE)] <- "yes"
immune[which(grepl("INFLAMMATORY", hsc_a_gsea_subset$pathways) == TRUE)] <- "yes"

df_hsca_gsea <- data.frame(
  pathway_name <- gsub("HALLMARK_", "", hsc_a_gsea_subset$pathways),
  p <- -log10(hsc_a_gsea_subset$pval),
  cell <- rep("HSC_a", length(hsc_a_gsea_subset$pathways)),
  col <- immune
)

colnames(df_hsca_gsea) <- c("pathway_name", "p", "cell", "col")

tiff(paste0(OUT_DIR, "HSC_a_gsea_bargraph.tiff"), units="in", width=6, height=3.5, res=300)
ggplot(data=df_hsca_gsea, aes(x=pathway_name, y=p, fill=col)) +
  geom_bar(stat="identity", aes(fill=col)) + theme_test() + coord_flip() + labs(x="Pathway", y="-log10(pval)") + theme(legend.position = "none") + geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  scale_fill_manual(values=c("#E7BA52", "#D6616B"))
dev.off()




# calculate HSCa gene score and correlate with IL6 and IL1B -----------------------------------------------

##Get names of all HSCa significantly correlated genes

il1b_result_list <- readRDS(file=paste0(OUT_DIR, "il1b_result_list_DE_genes_only.rds"))
names(il1b_result_list) <- c("HSC_a", "HSC_b", "CMP_a", "CMP_b", "CMP_c", "GMP_b1", "GMP_b2", "MEP_a1", "MEP_a2", "MEP_a3", "MLP_a", "MLP_b", "PreBNK")
hsca <- il1b_result_list[["HSC_a"]]

hsca_sig <- subset(hsca, hsca$padj <= 0.1)  
hsca_sig_up <- subset(hsca_sig, hsca_sig$rho > 0)
hsca_sig_down <- subset(hsca_sig, hsca_sig$rho < 0)

#get z-scores for each gene
celltype <- "HSC_a" ###CHOOSE celltype
mat <- expression_data[[1]] ##CHOOSE data
Td0_cols <- which(grepl("Td0", colnames(mat))==TRUE)
Td0_subset <- mat[Td0_cols]
row.names(Td0_subset) <- mat$X

Tm3_cols <- which(grepl("Tm3", colnames(mat))==TRUE)
Tm3_subset <- mat[Tm3_cols]
row.names(Tm3_subset) <- mat$X

logfc_mat <- data.frame()
count <- 0
name_short <- vector()
name_count <- 1
for(i in 1:length(Tm3_all))
{
  name <- Tm3_all[i]
  if(name %in% colnames(Tm3_subset))
  {
    d <- name
    Td0_pair_name <- Td0_all[which(Tm3_all == d)]
    Td0_col <- Td0_subset[Td0_pair_name]
    Tm3_col <- Tm3_subset[d]
    log2fc <- Tm3_col - Td0_col
    name_short[name_count] <- subset(meta_data, meta_data$Sample %in% as.character(d))$donor
    name_count <- name_count+1
    if(count==0)
    {
      logfc_mat <- log2fc[1]
    }
    else
    {
      logfc_mat <- cbind(logfc_mat, log2fc[1])
    }
    
    count <- count + 1
  }
  
}
colnames(logfc_mat) <- name_short
logfc_mat_final <- t(logfc_mat)


gene_vec_up <- logfc_mat_final[,which(colnames(logfc_mat_final) %in% hsca_sig_up$gene_name)]

scaled_gene_vec_up <- scale(gene_vec_up)
hsca_score_up <- rowMeans(scaled_gene_vec_up)

gene_vec_down <- logfc_mat_final[,which(colnames(logfc_mat_final) %in% hsca_sig_down$gene_name)]

scaled_gene_vec_down <- scale(gene_vec_down)
hsca_score_down <- rowMeans(scaled_gene_vec_down)


#plot these average UP scores versus IL1B
cytokine_data <- data.frame(fread("PBMC_cytokine_data_CHM.txt"))
meta_data <- read.csv(file="all_samples_meta_data.csv")
il1b <- vector(length=length(hsca_score_up))
treatment <- vector(length=length(hsca_score_up))

for(i in 1:length(hsca_score_up))
{
  cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) == names(hsca_score_up)[i])
  cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
  il1b[i] <- as.numeric(cyt_subset2$IL1B_FC)
  treatment[i] <- unique(meta_data$vaccination[which(meta_data$donor == names(hsca_score_up)[i])])
  print(i)
}

val <- cor.test(hsca_score_up, il1b, method="spearman")
pval_up_il1b <- val$p.value
rho_up_il1b <- val$estimate

df_hsca_up <- data.frame(
  up_score <- hsca_score_up,
  c <- il1b,
  t <- treatment
)
colnames(df_hsca_up) <- c("up", "c", "Group")

tiff(paste0(OUT_DIR, "HSCa_up_score_IL1B_corr.tiff"), units="in", width=5, height=4, res=300)
ggplot(df_hsca_up, aes(x=up, y=c)) + geom_point(shape=21, size=4,aes(fill=Group)) + theme_test() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab("HSCa gene score") +
  ylab("Fold Change IL1B production of PBMCs") + ggtitle("Positive gene score vs. IL1B")+
  annotate("text", x=-0.6, y=2, label= paste0("Spearman corr = ",round(val$estimate,2),"\npvalue = ",round(pval_up_il1b,6))) + scale_fill_manual(values=c("#8AE67E", "#004586"))
dev.off()


#plot these average UP scores versus IL6
cytokine_data <- data.frame(fread("PBMC_cytokine_data_CHM.txt"))
meta_data <- read.csv(file="all_samples_meta_data.csv")
il6 <- vector(length=length(hsca_score_up))
treatment <- vector(length=length(hsca_score_up))

for(i in 1:length(hsca_score_up))
{
  cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) == names(hsca_score_up)[i])
  cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
  il6[i] <- as.numeric(cyt_subset2$IL6_FC)
  treatment[i] <- unique(meta_data$vaccination[which(meta_data$donor == names(hsca_score_up)[i])])
  print(i)
}

val <- cor.test(hsca_score_up, il6, method="spearman")
pval_up_il6 <- val$p.value
rho_up_il6 <- val$estimate

df_hsca_up <- data.frame(
  up_score <- hsca_score_up,
  c <- il6,
  t <- treatment
)
colnames(df_hsca_up) <- c("up", "c", "Group")

tiff(paste0(OUT_DIR, "HSCa_up_score_IL6_corr.tiff"), units="in", width=5, height=4, res=300)
ggplot(df_hsca_up, aes(x=up, y=c)) + geom_point(shape=21, size=4,aes(fill=Group)) + theme_test() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab("HSCa gene score") +
  ylab("Fold Change IL6 production of PBMCs") + ggtitle("Positive gene score vs. IL6")+
  annotate("text", x=-0.55, y=2, label= paste0("Spearman corr = ",round(val$estimate,2),"\npvalue = ",round(pval_up_il6,6))) + scale_fill_manual(values=c("#F72836", "#004586"))
dev.off()


#plot these average DOWN scores versus IL1B
cytokine_data <- data.frame(fread("PBMC_cytokine_data_CHM.txt"))
meta_data <- read.csv(file="all_samples_meta_data.csv")
il1b <- vector(length=length(hsca_score_down))
treatment <- vector(length=length(hsca_score_down))

for(i in 1:length(hsca_score_down))
{
  cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) == names(hsca_score_down)[i])
  cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
  il1b[i] <- as.numeric(cyt_subset2$IL1B_FC)
  treatment[i] <- unique(meta_data$vaccination[which(meta_data$donor == names(hsca_score_down)[i])])
  print(i)
}

val <- cor.test(hsca_score_down, il1b, method="spearman")
pval_down_il1b <- val$p.value
rho_down_il1b <- val$estimate

df_hsca_down <- data.frame(
  down_score <- hsca_score_down,
  c <- il1b,
  t <- treatment
)
colnames(df_hsca_down) <- c("down", "c", "Group")

tiff(paste0(OUT_DIR, "HSCa_down_score_IL1B_corr.tiff"), units="in", width=5, height=4, res=300)
ggplot(df_hsca_down, aes(x=down, y=c)) + geom_point(shape=21, size=4,aes(fill=Group)) + theme_test() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab("HSCa gene score") +
  ylab("Fold Change IL1B production of PBMCs") + ggtitle("Repressive gene score vs. IL1B")+
  annotate("text", x=0.7, y=2, label= paste0("Spearman corr = ",round(val$estimate,2),"\npvalue = ",round(pval_down_il1b,6))) + scale_fill_manual(values=c("#8AE67E", "#004586"))
dev.off()

#plot these average DOWN scores versus IL6
cytokine_data <- data.frame(fread("PBMC_cytokine_data_CHM.txt"))
meta_data <- read.csv(file="all_samples_meta_data.csv")
il6 <- vector(length=length(hsca_score_down))
treatment <- vector(length=length(hsca_score_down))

for(i in 1:length(hsca_score_down))
{
  cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) == names(hsca_score_down)[i])
  cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
  il6[i] <- as.numeric(cyt_subset2$IL6_FC)
  treatment[i] <- unique(meta_data$vaccination[which(meta_data$donor == names(hsca_score_down)[i])])
  print(i)
}

val <- cor.test(hsca_score_down, il6, method="spearman")
pval_down_il6 <- val$p.value
rho_down_il6 <- val$estimate

df_hsca_down <- data.frame(
  down_score <- hsca_score_down,
  c <- il6,
  t <- treatment
)
colnames(df_hsca_down) <- c("down", "c", "Group")

tiff(paste0(OUT_DIR, "HSCa_down_score_IL6_corr.tiff"), units="in", width=5, height=4, res=300)
ggplot(df_hsca_down, aes(x=down, y=c)) + geom_point(shape=21, size=4, aes(fill=Group)) + theme_test() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab("HSCa gene score") +
  ylab("Fold Change IL6 production of PBMCs") + ggtitle("Repressive gene score vs. IL6")+
  annotate("text", x=0.7, y=2, label= paste0("Spearman corr = ",round(val$estimate,2),"\npvalue = ",round(pval_down_il6,6))) + scale_fill_manual(values=c("#F72836", "#004586"))
dev.off()

# bargraph to compare IL1B and IL6 rho values -----------------------------

df_bargraph <- data.frame(
  rho <- c(rho_up_il6, rho_up_il1b, rho_down_il6, rho_down_il1b),
  cyt <- c("IL6", "IL1B", "IL6", "IL1B"),
  p <- c(pval_up_il6, pval_up_il1b, pval_down_il6, pval_down_il1b),
  group <- c("+Gene score IL6", "+Gene score IL1B", "-Gene score IL6", "-Gene score IL1B")
)  
colnames(df_bargraph) <- c("rho", "cyt", "pval", "group")
df_bargraph$cyt <- factor(df_bargraph$cyt, c("IL6", "IL1B"))
df_bargraph$group <- factor(df_bargraph$group, c("+Gene score IL6", "+Gene score IL1B", "-Gene score IL6", "-Gene score IL1B"))

tiff(paste0(OUT_DIR, "rho_comparison_bargraph.tiff"), units="in", width=6, height=4, res=300)
ggplot(df_bargraph, aes(x=group, y=rho, fill=cyt)) + geom_bar(stat="identity", aes(fill=cyt), alpha=0.7) + theme_test() + 
  labs(x="", y="Rho", fill="Cytokine") + scale_fill_manual(values=c("#F72836", "#8AE67E")) + geom_text(aes(label=round(pval,6)),vjust=ifelse(sign(rho)>0, 1.7, -1.7),color="black", size=4)
dev.off()


#BCG samples ONLY
gene_vec_up <- logfc_mat_final[,which(colnames(logfc_mat_final) %in% hsca_sig_up$gene_name)]
BCG_donors <- unique(meta_data$donor[which(meta_data$vaccination == "BCG")])
scaled_gene_vec_up <- scale(subset(gene_vec_up, row.names(gene_vec_up) %in% BCG_donors))
hsca_score_up_BCG <- rowMeans(scaled_gene_vec_up)

gene_vec_down <- logfc_mat_final[,which(colnames(logfc_mat_final) %in% hsca_sig_down$gene_name)]
scaled_gene_vec_down <- scale(subset(gene_vec_down, row.names(gene_vec_down) %in% BCG_donors))
hsca_score_down_BCG <- rowMeans(scaled_gene_vec_down)

## UP scores versus IL6
cytokine_data <- data.frame(fread("PBMC_cytokine_data_CHM.txt"))
meta_data <- read.csv(file="all_samples_meta_data.csv")
il6 <- vector(length=length(hsca_score_up_BCG))
treatment <- vector(length=length(hsca_score_up_BCG))

for(i in 1:length(hsca_score_up_BCG))
{
  cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) == names(hsca_score_up_BCG)[i])
  cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
  il6[i] <- as.numeric(cyt_subset2$IL6_FC)
  treatment[i] <- unique(meta_data$vaccination[which(meta_data$donor == names(hsca_score_up_BCG)[i])])
  print(i)
}

val <- cor.test(hsca_score_up_BCG, il6, method="spearman")
pval_up_il6 <- val$p.value
rho_up_il6 <- val$estimate

df_hsca_up <- data.frame(
  up_score <- hsca_score_up_BCG,
  c <- il6,
  t <- treatment
)
colnames(df_hsca_up) <- c("up", "c", "Group")

tiff(paste0(OUT_DIR, "HSCa_up_score_IL6_corr_BCG.tiff"), units="in", width=4.5, height=4, res=300)
ggplot(df_hsca_up, aes(x=up, y=c)) + geom_point(shape=21, size=4, fill="#F72836") + theme_test() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab("HSCa gene score") +
  ylab("Fold Change IL6 production of PBMCs") + ggtitle("Positive gene score vs. IL6")+
  annotate("text", x=0.8, y=2, label= paste0("Spearman corr = ",round(val$estimate,2),"\npvalue = ",round(pval_up_il6,4))) 
dev.off()


## UP scores versus IL1B
cytokine_data <- data.frame(fread("PBMC_cytokine_data_CHM.txt"))
meta_data <- read.csv(file="all_samples_meta_data.csv")
il1b <- vector(length=length(hsca_score_up_BCG))
treatment <- vector(length=length(hsca_score_up_BCG))

for(i in 1:length(hsca_score_up_BCG))
{
  cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) == names(hsca_score_up_BCG)[i])
  cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
  il1b[i] <- as.numeric(cyt_subset2$IL1B_FC)
  treatment[i] <- unique(meta_data$vaccination[which(meta_data$donor == names(hsca_score_up_BCG)[i])])
  print(i)
}

val <- cor.test(hsca_score_up_BCG, il1b, method="spearman")
pval_up_il1b <- val$p.value
rho_up_il1b <- val$estimate

df_hsca_up <- data.frame(
  up_score <- hsca_score_up_BCG,
  c <- il1b,
  t <- treatment
)
colnames(df_hsca_up) <- c("up", "c", "Group")

tiff(paste0(OUT_DIR, "HSCa_up_score_IL1B_corr_BCG.tiff"), units="in", width=4.5, height=4, res=300)
ggplot(df_hsca_up, aes(x=up, y=c)) + geom_point(shape=21, size=4, fill="#8AE67E") + theme_test() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab("HSCa gene score") +
  ylab("Fold Change IL1B production of PBMCs") + ggtitle("Positive gene score vs. IL1B")+
  annotate("text", x=-0.1, y=2.1, label= paste0("Spearman corr = ",round(val$estimate,2),"\npvalue = ",round(pval_up_il1b,4))) 
dev.off()


## DOWN scores versus IL6
cytokine_data <- data.frame(fread("PBMC_cytokine_data_CHM.txt"))
meta_data <- read.csv(file="all_samples_meta_data.csv")
il6 <- vector(length=length(hsca_score_down_BCG))
treatment <- vector(length=length(hsca_score_down_BCG))

for(i in 1:length(hsca_score_down_BCG))
{
  cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) == names(hsca_score_down_BCG)[i])
  cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
  il6[i] <- as.numeric(cyt_subset2$IL6_FC)
  treatment[i] <- unique(meta_data$vaccination[which(meta_data$donor == names(hsca_score_down_BCG)[i])])
  print(i)
}

val <- cor.test(hsca_score_down_BCG, il6, method="spearman")
pval_down_il6 <- val$p.value
rho_down_il6 <- val$estimate

df_hsca_down <- data.frame(
  down_score <- hsca_score_down_BCG,
  c <- il6,
  t <- treatment
)
colnames(df_hsca_down) <- c("down", "c", "Group")

tiff(paste0(OUT_DIR, "HSCa_down_score_IL6_corr_BCG.tiff"), units="in", width=4.5, height=4, res=300)
ggplot(df_hsca_down, aes(x=down, y=c)) + geom_point(shape=21, size=4, fill="#F72836") + theme_test() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab("HSCa gene score") +
  ylab("Fold Change IL6 production of PBMCs") + ggtitle("Repressive gene score vs. IL6")+
  annotate("text", x=-0.6, y=2, label= paste0("Spearman corr = ",round(val$estimate,2),"\npvalue = ",round(pval_down_il6,4))) 
dev.off()


## DOWN scores versus IL1B
cytokine_data <- data.frame(fread("PBMC_cytokine_data_CHM.txt"))
meta_data <- read.csv(file="all_samples_meta_data.csv")
il1b <- vector(length=length(hsca_score_down_BCG))
treatment <- vector(length=length(hsca_score_down_BCG))

for(i in 1:length(hsca_score_down_BCG))
{
  cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) == names(hsca_score_down_BCG)[i])
  cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
  il1b[i] <- as.numeric(cyt_subset2$IL1B_FC)
  treatment[i] <- unique(meta_data$vaccination[which(meta_data$donor == names(hsca_score_down_BCG)[i])])
  print(i)
}

val <- cor.test(hsca_score_down_BCG, il1b, method="spearman")
pval_down_il1b <- val$p.value
rho_down_il1b <- val$estimate

df_hsca_down <- data.frame(
  down_score <- hsca_score_down_BCG,
  c <- il1b,
  t <- treatment
)
colnames(df_hsca_down) <- c("down", "c", "Group")

tiff(paste0(OUT_DIR, "HSCa_down_score_IL1B_corr_BCG.tiff"), units="in", width=4.5, height=4, res=300)
ggplot(df_hsca_down, aes(x=down, y=c)) + geom_point(shape=21, size=4, fill="#8AE67E") + theme_test() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab("HSCa gene score") +
  ylab("Fold Change IL1B production of PBMCs") + ggtitle("Repressive gene score vs. IL1B")+
  annotate("text", x=0.2, y=2, label= paste0("Spearman corr = ",round(val$estimate,2),"\npvalue = ",round(pval_down_il1b,4))) 
dev.off()



# bargraph to compare IL1B and IL6 rho values -----------------------------

df_bargraph <- data.frame(
  rho <- c(rho_up_il6, rho_up_il1b, rho_down_il6, rho_down_il1b),
  cyt <- c("IL6", "IL1B", "IL6", "IL1B"),
  p <- c(pval_up_il6, pval_up_il1b, pval_down_il6, pval_down_il1b),
  group <- c("+Gene score IL6", "+Gene score IL1B", "-Gene score IL6", "-Gene score IL1B")
)  
colnames(df_bargraph) <- c("rho", "cyt", "pval", "group")
df_bargraph$cyt <- factor(df_bargraph$cyt, c("IL6", "IL1B"))
df_bargraph$group <- factor(df_bargraph$group, c("+Gene score IL6", "+Gene score IL1B", "-Gene score IL6", "-Gene score IL1B"))

tiff(paste0(OUT_DIR, "rho_comparison_bargraph_BCG.tiff"), units="in", width=6, height=3.5, res=300)
ggplot(df_bargraph, aes(x=group, y=rho, fill=cyt)) + geom_bar(stat="identity", aes(fill=cyt), alpha=0.7) + theme_test() + 
  labs(x="", y="Rho", fill="Cytokine") + scale_fill_manual(values=c("#F72836", "#8AE67E")) + geom_text(aes(label=paste0("p = ",round(pval,4))), vjust=ifelse(sign(rho)>0, 1.3, -0.9),color="black", size=4)
dev.off()
