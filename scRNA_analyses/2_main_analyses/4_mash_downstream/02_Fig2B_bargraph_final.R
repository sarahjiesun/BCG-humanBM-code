
#R version 4.1.0 loaded
#rstudio/1.4
#geos/3.9.1 loaded
#hdf5/1.12.0 loaded 
#cmake loaded

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
library(cowplot)

setwd("/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis_Raul")

OUT_DIR <- c("MASH_emmreml_downstream/DE_genes_bargraphs/")

source("bcf_csf_utils.R")

theme_set(theme_minimal_grid() +
            theme(
              #base_size=2,
              legend.position = "bottom" , 
              text = element_text(family = "Helvetica"), 
              legend.text = element_text(size=8),
              legend.title =  element_text(size=8),
              axis.title = element_text(size=8),
              axis.text.y = element_text(size=8),
              strip.text.x = element_text(size=8, angle=0, vjust = 0.5, hjust = 0.5),
              strip.text.y = element_text(size=8, angle=0, vjust = 0.5, hjust = 0.5),
              axis.text.x = element_text(size=8))
)




barplot_file <- "/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis_Raul/MASH_emmreml_downstream/DE_genes_bargraphs/Figure_2A_bargraph.rds"
barplot <- readRDS(barplot_file)

# myColorScale <- c('HSC_a'='#EF5350','HSC_b'='#CC0C00','CMP_a'='#FF7F0E','CMP_b'='#FFD147','CMP_c'='#CC9900',
#                   'GMP_a'='#33CC00', 'GMP_b' ='#99CC00', 'MEP_a'='#8C5782', 'MEP_b'='#990080', 'MEP_c'='#BA68C8',
#                   'MLP_a'='#46B8DA', 'MLP_b'='#0288D1', 'PreBNK'='#1A0099')
# 
# 
# ct_info <- data.frame(ct=names(myColorScale)) %>% 
#   mutate(type= gsub("_.*","",ct)) 
# ct_info$type[ct_info$type == "PreBNK"] <- "MLP"
# 

barplot$data$cluster <- factor(barplot$data$cluster, c("HSC_a", "HSC_b", "CMP_a", "CMP_b", "CMP_c", "GMP_a", "GMP_b", "MEP_a", "MEP_b", "MEP_c", "MLP_a", "MLP_b", "PreBNK"))

barplot_new <- barplot$data %>% 
  mutate(fill_col= ifelse(group == "up","BCG","Placebo")) %>% 
  mutate(cluster_type = ct_info$type[match(cluster,ct_info$ct)]) %>% 
  mutate(cluster_type = factor(cluster_type, levels= names(ct_type_scale))) %>% 
  ggplot(aes(x = cluster, y = values)) + 
  geom_bar(aes(fill= fill_col), stat="identity", alpha=0.9) + 
  facet_grid(~cluster_type, scales = "free",space = "free") + 
  #scale_fill_met_d(name = "Benedictus") + 
  scale_fill_manual(values = vx_colors) + 
  ylab("Number of DR genes (lfsr < 0.01)") + 
  xlab("Cell type") + 
  theme(text = element_text(size=8, family = "Helvetica"),   
        strip.text = element_text(size=8, family = "Helvetica",angle=0, hjust = 0.5),
        #axis.text = element_text(size=9, family = "Helvetica"),
        axis.text.x = element_text(size=7,angle=45, hjust = 1), 
        axis.text.y = element_text(size=8), 
        strip.background = element_blank(),
        legend.position="none") 


#transposed 
t_barplot_new <- barplot$data %>% 
  mutate(fill_col= ifelse(group == "up","BCG","Placebo")) %>% 
  mutate(cluster_type = ct_info$type[match(cluster,ct_info$ct)]) %>% 
  mutate(cluster_type = factor(cluster_type, levels= names(ct_type_scale))) %>% 
  ggplot(aes(y = cluster, x = values)) + 
  geom_bar(aes(fill= fill_col), stat="identity", alpha=0.9) + 
  facet_grid(cluster_type~., scales = "free",space = "free") + 
  #scale_fill_met_d(name = "Benedictus") + 
  scale_fill_manual(values = vx_colors) + 
  xlab("N DRG (lfsr < 0.01)") + 
  ylab("Cell type") + 
  theme(text = element_text(size=8, family = "Helvetica"),   
        strip.text = element_text(size=8, family = "Helvetica",angle=0, hjust = 0.5),
        #axis.text = element_text(size=9, family = "Helvetica"),
        axis.text.x = element_text(size=7,angle=45, hjust = 1), 
        axis.text.y = element_text(size=8), 
        strip.background = element_blank(),
        legend.position="none") + scale_y_discrete(limits=rev)

t_barplot_file <- file.path(OUT_DIR, 
                            paste0("Fig2_t_barplot",".pdf"))
ggsave(filename = t_barplot_file, plot = t_barplot_new, height = 2.5 , width = 1.8)
