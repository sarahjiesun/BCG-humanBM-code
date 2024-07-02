.libPaths("/project/lbarreiro/USERS/sarah/Rlibs_new")
library(ggplot2)
library (DESeq2)
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

setwd("/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis_Raul/") 

OUT_DIR <- "HSC_subcluster_analysis/5_mash/"
dir.create(OUT_DIR)

# read in results for each cluster ----------------------------------------

#exclude cluster c9 due to not enough samples 
clusters <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")

genes_common <- vector()
for(i in 1:length(clusters)){
  name <- clusters[i]
  dat <- read.csv(file=paste0("HSC_subcluster_analysis/4_emmreml_betas/", name, "_res_full.csv"))
  if(i == 1)
  {
    genes_common <- dat$X
  }
  else
  {
    genes_common <- intersect(genes_common, dat$X)
  }
  assign(paste0("c",name,"_res"), dat)
}

results <- list(c0_res, c1_res, c2_res, c3_res, c4_res, c5_res, c6_res, c7_res, c8_res)


# make table of betas and standard errors -----------------------------------------------------

for(i in 1:length(results)){
  name <- clusters[i]
  dat <- results[[i]]
  dat_subset <- subset(dat, dat$X %in% genes_common)
  print(which(dat_subset$X != genes_common))
  assign(paste0("c",name, "_subset"), dat_subset)
}


beta_table <- cbind(c0_subset$beta_groupBCG.timepoint, c1_subset$beta_groupBCG.timepoint, c2_subset$beta_groupBCG.timepoint, c3_subset$beta_groupBCG.timepoint,
                    c4_subset$beta_groupBCG.timepoint, c5_subset$beta_groupBCG.timepoint, c6_subset$beta_groupBCG.timepoint , c7_subset$beta_groupBCG.timepoint ,
                    c8_subset$beta_groupBCG.timepoint)
colnames(beta_table) <- c("c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8")
row.names(beta_table) <- genes_common


SE_table <- cbind(sqrt(c0_subset$sdev_groupBCG.timepoint), sqrt(c1_subset$sdev_groupBCG.timepoint), sqrt(c2_subset$sdev_groupBCG.timepoint), sqrt(c3_subset$sdev_groupBCG.timepoint),
                  sqrt(c4_subset$sdev_groupBCG.timepoint),  sqrt(c5_subset$sdev_groupBCG.timepoint),  sqrt(c6_subset$sdev_groupBCG.timepoint),  sqrt(c7_subset$sdev_groupBCG.timepoint),
                  sqrt(c8_subset$sdev_groupBCG.timepoint))
colnames(SE_table) <- c("c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8")
row.names(SE_table) <- genes_common



# Run MASH with correlations ----------------------------------------------

data = mash_set_data(beta_table, SE_table)

# estimate null correlations and add to data (decreases number of false positive from my test)
V = estimate_null_correlation_simple(data)
data_V = mash_update_data(data, V = V)

#set up data driven covariance matrices
m.1by1 = mash_1by1(data_V)
strong = get_significant_results(m.1by1,0.05)
U.pca = cov_pca(data_V,5,subset=strong)
print(names(U.pca))
U.ed = cov_ed(data_V, U.pca, subset=strong)

#canonical covariance matrix
U.c = cov_canonical(data_V)
print(names(U.c))

#mash
m_corr   = mash(data_V, c(U.c,U.ed))


get_significant_results(m_corr)
print(length(get_significant_results(m_corr)))


## write posterior outs
write.table(get_lfsr(m_corr), paste0(OUT_DIR,"mash_results/lfsr_wcorrelations_output.txt"), quote = FALSE)
write.table(get_pm(m_corr), paste0(OUT_DIR, "mash_results/posteriorMeans_wcorrelations.txt"), quote = FALSE)
write.table(get_psd(m_corr), paste0(OUT_DIR, "mash_results/posteriorStandardDevs_wcorrelations.txt"), quote = FALSE)



# Run MASH without correlations ----------------------------------------------------------------
data = mash_set_data(beta_table, SE_table)

#set up data driven covariance matrices
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.05)
U.pca = cov_pca(data,5,subset=strong)
print(names(U.pca))
U.ed = cov_ed(data, U.pca, subset=strong)

#canonical covariance matrix
U.c = cov_canonical(data)
print(names(U.c))


#mash
m   = mash(data, c(U.c,U.ed))


get_significant_results(m)
print(length(get_significant_results(m)))


## write posterior outs
write.table(get_lfsr(m), paste0(OUT_DIR,"mash_results/lfsr_output.txt"), quote = FALSE)
write.table(get_pm(m), paste0(OUT_DIR, "mash_results/posteriorMeans.txt"), quote = FALSE)
write.table(get_psd(m), paste0(OUT_DIR, "mash_results/posteriorStandardDevs.txt"), quote = FALSE)