library(ArchR)
addArchRThreads(threads = 16) 
addArchRGenome("hg38")
library(readr)
library(harmony)
library(Seurat)

# Run noninteractively ----------------------------------------------------

#See callpeaks_all_combined.sbatch

setwd("/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis1_clusters_harmony/")


# Add Clusters harmony labels to undited object ---------------------------

#DO NOT CHANGE THIS OBJECT
projBM_labels <- loadArchRProject(path = "../ArchR_initial_processing/3_Save5_ProjBM5_scRNA_integrated_clusters_remapped", force = FALSE, showLogo = TRUE)

#Unedited object
projBM <- loadArchRProject(path = "0_Save1-ProjBM1", force = FALSE, showLogo = TRUE)


df_labels <- getCellColData(projBM_labels)
df <- getCellColData(projBM)


#Cells are in the same order!
which(row.names(df) != row.names(df_labels))

projBM$Clusters_harmony <- projBM_labels$Clusters_harmony

df_updated <- getCellColData(projBM)

projBM <- addGroupCoverages(ArchRProj = projBM, groupBy = "Clusters_harmony")

# Call peaks --------------------------------------------------------------

pathToMacs2 <- findMacs2()

projBM <- addReproduciblePeakSet(
  ArchRProj = projBM, 
  groupBy = "Clusters_harmony", 
  pathToMacs2 = pathToMacs2
)

getPeakSet(projBM)

projBM <- addPeakMatrix(projBM)

getAvailableMatrices(projBM)

saveArchRProject(ArchRProj = projBM, outputDirectory = "Unedited_with_Clusters_harmony_peaks_combined", load = FALSE)
