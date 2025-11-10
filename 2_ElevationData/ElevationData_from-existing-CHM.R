rm(list=ls())
gc()

# Chargement des packages
library(tidyverse)
library(sf)
library(terra)
library(foreach)
library(doParallel)
library(here)

# Fonction pour obtenir l'EPSG à partir de la zone UTM
get_epsg_from_utm <- function(utm_zone) {
  epsg <- utm_epsg[utm_zone]
  if(is.na(epsg)) warning(paste("EPSG not found for", utm_zone))
  return(as.integer(epsg))
}

# Fonction pour obtenir le chemin du CHM
get_chm_path <- function(source, utm, folder_name) {
  # EXTERNAL DATA PATH: Pre-processed LiDAR rasters from NASA 2015 campaign
  # Adjust to your local setup if needed
  raster_base <- "E:/Arthur/Doctorat_DataAnnexe/LiDAR_RDC_2015/001_Lidar_Rasters/Lidar_2meter"
  if(source == "transects") {
    return(file.path(raster_base, "CHM", paste0(utm, "_", folder_name, "_CHM.tif")))
  } else {
    return(file.path(raster_base, "Ferry/CHM", paste0(utm, "_", folder_name, "_CHM.tif")))
  }
}

# Configuration des chemins
input_base <- here("1_LiDAR_extraction", "output")
output_base <- here("2_ElevationData")

# Création du dossier CHM_NASA2015 dans ElevationData
output_chm <- file.path(output_base, "CHM_NASA2015")
dir.create(output_chm, showWarnings = FALSE, recursive = TRUE)

# Définition des EPSG UTM
utm_epsg <- c(
  "UTM33S" = 32733,
  "UTM34N" = 32634,
  "UTM34S" = 32734,
  "UTM35N" = 32635,
  "UTM35S" = 32735,
  "UTM36N" = 32636
)

# Lecture du fichier plots_info
plots_info <- read_csv2(file.path(input_base, "plots_info.csv"))

# Configuration de la parallélisation
n_cores <- parallel::detectCores() / 2
cl <- makeCluster(floor(n_cores))
registerDoParallel(cl)

# Définition des packages à exporter vers les workers
packages_to_export <- c("terra", "sf", "tidyverse")

# Traitement parallèle
results <- foreach(
  i = 1:nrow(plots_info),
  .packages = packages_to_export,
  .errorhandling = "pass"
) %dopar% {
  tryCatch({
    plot_info <- plots_info[i,]
    plot_ref <- plot_info$plot_ref
    target_epsg <- get_epsg_from_utm(plot_info$utm)
    
    # Obtention des chemins
    gpkg_path <- here("0_Inventories_plot_preparation", "final", "plots_unique", paste0(plot_ref, ".gpkg"))
    chm_path <- get_chm_path(plot_info$source, plot_info$utm, plot_info$folder_name)
    
    # Vérification de l'existence des fichiers
    if(!file.exists(chm_path)) {
      return(list(
        status = "error",
        plot_ref = plot_ref,
        message = paste("CHM file not found:", chm_path)
      ))
    }
    if(is.na(gpkg_path)) {
      return(list(
        status = "error",
        plot_ref = plot_ref,
        message = "Shapefile not found"
      ))
    }
    
    # Lecture des données
    chm <- terra::rast(chm_path)
    plot_boundary <- st_read(gpkg_path, quiet = TRUE) %>% st_buffer(dist = 200)
    
    # Vérification et assignation du CRS pour le CHM
    if(is.na(terra::crs(chm))) {
      terra::crs(chm) <- st_crs(target_epsg)$wkt
    }
    
    # Transformation du shapefile dans le CRS du CHM
    plot_boundary <- st_transform(plot_boundary, terra::crs(chm))
    
    # Découpage du CHM selon les limites du plot
    chm_masked <- terra::mask(chm, vect(plot_boundary))
    chm_cropped <- terra::crop(chm_masked, vect(plot_boundary))
    
    # Export du CHM découpé
    output_file <- file.path(output_chm, paste0(plot_ref, "_CHM.tif"))
    terra::writeRaster(chm_cropped, output_file, overwrite = TRUE)
    
    return(list(
      status = "success",
      plot_ref = plot_ref,
      original_crs = terra::crs(chm)
    ))
    
  }, error = function(e) {
    return(list(
      status = "error",
      plot_ref = plot_ref,
      message = as.character(e)
    ))
  })
}

# Arrêt du cluster
stopCluster(cl)

# Analyse des résultats
results_df <- bind_rows(lapply(results, function(x) {
  if(is.list(x)) {
    as_tibble(x)
  } else {
    tibble(status = "error", message = "Unknown error")
  }
}))

# Affichage des résultats
cat("\nTraitement terminé!\n")
cat("\nRésumé:\n")
cat("Succès:", sum(results_df$status == "success"), "\n")
cat("Erreurs:", sum(results_df$status == "error"), "\n")

# Affichage des erreurs s'il y en a
if(any(results_df$status == "error")) {
  cat("\nDétails des erreurs:\n")
  print(results_df %>% filter(status == "error") %>% select(plot_ref, message))
}

# Sauvegarde des résultats
write_csv(results_df, file.path(output_chm, "chm_extraction_results.csv"))