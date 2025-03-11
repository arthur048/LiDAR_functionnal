rm(list=ls())
gc()

# Chargement des packages ----
library(tidyverse)
library(sf)
library(lidR)
library(terra)
library(foreach)
library(doParallel)

# Fonctions utilitaires ----
# Fonction pour filtrer le bruit
filter_noise = function(las, sensitivity) {
  if (is(las, "LAS")) {
    p95 <- pixel_metrics(las, ~quantile(Z, probs = 0.95), 10)
    las <- merge_spatial(las, p95, "p95")
    las <- filter_poi(las, Z < p95*sensitivity)
    las$p95 <- NULL
    return(las)
  }
  
  if (is(las, "LAScatalog")) {
    options <- list(
      need_output_file = TRUE,
      need_buffer = TRUE)
    res <- catalog_map(las, filter_noise, sensitivity = sensitivity, .options = options)
    return(res)
  }
}

# Fonction pour obtenir l'EPSG à partir de la zone UTM
get_epsg_from_utm <- function(utm_zone) {
  epsg <- utm_epsg[utm_zone]
  if(is.na(epsg)) warning(paste("EPSG not found for", utm_zone))
  return(as.integer(epsg))
}

# Fonction pour interpoler les valeurs manquantes dans un raster
interpolate_na = function(rast) {
  # Création d'une matrice de voisinage 3x3
  w <- matrix(1, 3, 3)
  
  # Tant qu'il y a des NA, continuer l'interpolation
  while(any(is.na(values(rast)))) {
    focal_mean <- focal(rast, w = w, fun = mean, na.rm = TRUE)
    na_cells <- is.na(values(rast))
    values(rast)[na_cells] <- values(focal_mean)[na_cells]
  }
  
  return(rast)
}

# Fonction pour vérifier et corriger le CRS
ensure_correct_crs <- function(obj, target_epsg) {
  if(is(obj, "LAS")) {
    if(epsg(obj) != target_epsg) {
      st_crs(obj) <- target_epsg
    }
    return(obj)
  } else if(is(obj, "SpatRaster")) {
    if(terra::crs(obj, proj=TRUE) != st_crs(target_epsg)$wkt) {
      terra::crs(obj) <- st_crs(target_epsg)$wkt
    }
    return(obj)
  }
  return(obj)
}

# Fonction pour obtenir le masque de référence
get_reference_mask <- function(source, utm, name, target_epsg) {
  base_path <- "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/shapefile/lidar"
  
  ref_path <- if(source == "transects") {
    file.path(base_path, "DRC01_lidar_transects.shp")
  } else {
    file.path(base_path, "DRC02_ferry_lines.shp")
  }
  
  ref_shapefile <- st_read(ref_path, quiet = TRUE) %>%
    filter(Name == paste0(utm, "_", name)) %>%
    st_transform(target_epsg)
  
  if(nrow(ref_shapefile) == 0) {
    stop(paste("No polygon found for:", paste0(utm, "_", name)))
  }
  
  return(ref_shapefile)
}

# Fonction pour traiter les rasters
process_raster <- function(raster, ref_mask, plot_mask) {
  # Interpolation des valeurs manquantes si nécessaire
  if(any(is.na(values(raster)))) {
    raster <- interpolate_na(raster)
  }
  
  # Application séquentielle des masques
  raster_masked <- crop(raster, vect(ref_mask), mask = TRUE)
  raster_final <- crop(raster_masked, vect(plot_mask), mask = TRUE)
  
  return(raster_final)
}

# Configuration ----
# Chemins
path_base <- "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles"
input_base <- file.path(path_base, "LiDAR_functionnal/1_LiDAR_extraction/output")
output_base <- file.path(path_base, "LiDAR_functionnal/2_ElevationData")

# Définition des EPSG UTM
utm_epsg <- c(
  "UTM33S" = 32733,
  "UTM34N" = 32634,
  "UTM34S" = 32734,
  "UTM35N" = 32635,
  "UTM35S" = 32735,
  "UTM36N" = 32636
)

# Création des dossiers de sortie
paths <- list(
  chm = file.path(output_base, "CHM"),
  dtm = file.path(output_base, "DTM"),
  nlas = file.path(output_base, "LAS_normalized")
)

# Création des dossiers s'ils n'existent pas
walk(paths, ~dir.create(., showWarnings = FALSE, recursive = TRUE))

# Lecture du fichier plots_info
plots_info <- read_csv2("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/final/plots_info_NASA2015.csv")

# Configuration de la parallélisation
n_cores <- parallel::detectCores() / 2
cl <- makeCluster(floor(n_cores))
registerDoParallel(cl)

# Traitement parallèle ----
results <- foreach(
  i = 1:nrow(plots_info),
  .packages = c("lidR", "terra", "tidyverse", "sf"),
  .errorhandling = "pass"
) %dopar% {
  tryCatch({
    # Initialisation des données pour le plot
    plot_info <- plots_info[i,]
    plot_ref <- plot_info$plot_ref
    target_epsg <- get_epsg_from_utm(plot_info$utm)
    plot_dir <- file.path(input_base, plot_ref)
    
    # Lecture du fichier LAS
    las_path <- list.files(plot_dir, pattern = "\\.las$", full.names = TRUE)[1]
    if(is.na(las_path)) {
      return(list(status = "error", 
                  plot_ref = plot_ref,
                  message = "LAS file not found"))
    }
    
    # Lecture du shapefile du plot
    shp_path <- list.files(plot_dir, pattern = "\\.shp$", full.names = TRUE)[1]
    if(is.na(shp_path)) {
      return(list(status = "error", 
                  plot_ref = plot_ref,
                  message = "Shapefile not found"))
    }
    
    # Lecture des données
    las <- readLAS(las_path, select = "xyzcrn")
    plot_boundary <- st_read(shp_path, quiet = TRUE)
    
    # Assignation et transformation des CRS
    st_crs(las) <- target_epsg
    plot_boundary <- st_transform(plot_boundary, target_epsg)
    
    # Obtention du masque de référence
    ref_mask <- get_reference_mask(
      plot_info$source,
      plot_info$utm,
      plot_info$transect_name,
      target_epsg
    )
    
    # Création et traitement du DTM
    dtm <- rasterize_terrain(las, algorithm = tin())
    dtm <- ensure_correct_crs(dtm, target_epsg)
    
    # Normalisation et filtrage du nuage de points
    nlas <- normalize_height(las, tin())
    nlas <- ensure_correct_crs(nlas, target_epsg)
    nlas_filtered <- filter_noise(nlas, 1.2)
    nlas_filtered <- ensure_correct_crs(nlas_filtered, target_epsg)
    
    # Création du CHM
    chm <- rasterize_canopy(
      nlas_filtered,
      res = 1.0,
      algorithm = pitfree(
        thresholds = c(0,2,5,10,20,30,40,50,60),
        max_edge = c(0,3),
        subcircle = 0.3
      )
    )
    chm <- ensure_correct_crs(chm, target_epsg)
    
    # Application des masques aux rasters
    chm_final <- process_raster(chm, ref_mask, plot_boundary)
    dtm_final <- process_raster(dtm, ref_mask, plot_boundary)
    
    # Export des résultats
    terra::writeRaster(
      chm_final,
      file.path(paths$chm, paste0(plot_ref, ".tif")),
      overwrite = TRUE
    )
    
    terra::writeRaster(
      dtm_final,
      file.path(paths$dtm, paste0(plot_ref, ".tif")),
      overwrite = TRUE
    )
    
    lidR::writeLAS(
      nlas_filtered,
      file.path(paths$nlas, paste0(plot_ref, ".las")),
      index = TRUE
    )
    
    return(list(
      status = "success",
      plot_ref = plot_ref
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
write_csv(results_df, file.path(output_base, "processing_results.csv"))