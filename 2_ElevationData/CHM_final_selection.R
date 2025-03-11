rm(list=ls())
gc()

library(readxl)
library(raster)
library(fs)

base_path <- "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/2_ElevationData"
excel_path <- file.path(base_path, "selection_final_CHM.xlsx")
dest_folder <- file.path(base_path, "CHM_final")

dir.create(dest_folder, showWarnings = FALSE)

df <- read_excel(excel_path) %>% filter(keep %in% c(1,2))

for(i in 1:nrow(df)) {
  source_folder <- file.path(base_path, df$origin_CHM[i])
  source_file <- file.path(source_folder, paste0(df$plot_ref[i], ".tif"))
  dest_file <- file.path(dest_folder, paste0(df$plot_ref[i], ".tif"))
  
  if(file.exists(source_file)) {
    file.copy(source_file, dest_file, overwrite = TRUE)
    cat("Copié:", df$plot_ref[i], "\n")
  } else {
    cat("ERREUR - Fichier non trouvé:", source_file, "\n")
  }
}

cat("\nTerminé! Les fichiers ont été copiés dans:", dest_folder)
