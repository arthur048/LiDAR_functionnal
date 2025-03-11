rm(list=ls())
gc()

# Chargement des packages nécessaires
library(tidyverse)
library(sf)
library(terra)

# Lecture des données
plots_info <- read_csv2("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/final/plots_info.csv") %>%
  filter(source == "AfriSAR")

# Chemins des données
plots_path <- "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/final/plots_unique"
chm_dir <- "E:/Arthur/OneDrive2/R/DoctoratGIS/OriginalDataFiles/AfriSAR_TropiSAR_CHM_FieldInventories/CHM"

# Liste des CHM
chm_files <- list.files(chm_dir, pattern = "\\.tif$", full.names = TRUE)

# Création du dossier de sortie avec la date du jour
output_dir <- file.path("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/2_ElevationData", 
                        paste0("CHM_", format(Sys.Date(), "%Y%m%d"), "_", sprintf("%03d", sample(0:999, 1))))
dir.create(output_dir, recursive = TRUE)

# Pour chaque plot
for(i in 1:nrow(plots_info)) {
  plot_ref <- plots_info$plot_ref[i]
  plot_path <- file.path(plots_path, paste0(plot_ref, ".gpkg"))
  
  plot_polygon <- st_read(plot_path, quiet = TRUE) %>%
    st_transform(plots_info$EPSG[i]) %>%
    st_buffer(200)
  
  for(chm_file in chm_files) {
    chm <- terra::rast(chm_file)
    chm_epsg <- as.numeric(gsub("EPSG:", "", terra::crs(chm, describe=TRUE)$code))
    
    if(chm_epsg == st_crs(plot_polygon)$epsg) {
      chm_bbox <- st_bbox(terra::ext(chm)) %>% 
        st_as_sfc() %>%
        st_set_crs(chm_epsg)  # On assigne explicitement le CRS au bbox
      
      if(st_intersects(plot_polygon, chm_bbox, sparse = FALSE)[1]) {
        # Crop et mask
        chm_cropped <- terra::mask(
          terra::crop(chm, vect(plot_polygon)),
          vect(plot_polygon)
        )
        
        # Sauvegarde
        terra::writeRaster(chm_cropped, 
                           file.path(output_dir, paste0(plot_ref, ".tif")),
                           overwrite = TRUE)
        
        cat("Plot", plot_ref, "traité avec", basename(chm_file), "\n")
        break
      }
    }
  }
}

cat("\nTraitement terminé! Les résultats sont dans:", output_dir, "\n")
