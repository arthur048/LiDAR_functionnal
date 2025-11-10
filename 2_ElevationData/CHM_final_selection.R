rm(list=ls())
gc()

library(readxl)
library(raster)
library(fs)
library(here)
library(dplyr)

excel_path <- here("2_ElevationData", "selection_final_CHM.xlsx")
dest_folder <- here("2_ElevationData", "CHM_final")

dir.create(dest_folder, showWarnings = FALSE, recursive = TRUE)

df <- read_excel(excel_path) %>% filter(keep %in% c(1,2))

for(i in 1:nrow(df)) {
  source_folder <- here("2_ElevationData", df$origin_CHM[i])
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
