rm(list=ls())
gc()

library(tidyverse)
library(sf)
library(lidR)

# Fonction pour extraire l'UTM et le nom du transect
extract_utm_plot <- function(name) {
  utm <- str_extract(name, "UTM\\d{2}[A-Z]")
  plot_name <- str_extract(name, "Plot.*$")
  return(list(utm = utm, plot_name = plot_name))
}

# Fonction pour formater le nom du dossier ferry
format_ferry_folder <- function(plot_name) {
  if(str_detect(plot_name, "to")) {
    nums <- str_extract_all(plot_name, "\\d+")[[1]]
    return(paste0("Plot", nums[1], " to Plot", nums[2]))
  }
  return(plot_name)
}

# Fonction pour obtenir le chemin LiDAR
get_lidar_path <- function(source, utm, folder_name) {
  if(source == "transects") {
    return(file.path(lidar_base, "Plots", utm, folder_name, "Las"))
  } else {
    return(file.path(lidar_base, "Ferry Data", utm, folder_name))
  }
}

# Fonction pour obtenir l'EPSG à partir de la zone UTM
get_epsg_from_utm <- function(utm_zone) {
  epsg <- utm_epsg[utm_zone]
  if(is.na(epsg)) warning(paste("EPSG not found for", utm_zone))
  return(as.integer(epsg))
}

utm_epsg <- c(
  "UTM33S" = 32733,
  "UTM34N" = 32634,
  "UTM34S" = 32734,
  "UTM35N" = 32635,
  "UTM35S" = 32735,
  "UTM36N" = 32636
)

# 1. Lecture et préparation des données ----
path_base <- "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles"
lidar_base <- "E:/Arthur/Doctorat_DataAnnexe/LiDAR_RDC_2015/00-raw-data"

# Lecture des shapefiles avec correction de géométrie
plots <- st_read(paste0(path_base, "/LiDAR_functionnal/0_Inventories_plot_preparation/final/plots_limites_NASA2015.gpkg"))
transects <- st_read(paste0(path_base, "/LiDAR_functionnal/0_Inventories_plot_preparation/shapefile/lidar/DRC01_lidar_transects.shp")) %>%
  st_make_valid()

ferry_lines <- st_read(paste0(path_base, "/LiDAR_functionnal/0_Inventories_plot_preparation/shapefile/lidar/DRC02_ferry_lines.shp"), 
                       options = "OGR_GEOMETRY_ACCEPT_UNCLOSED_RING=NO") %>%
  st_buffer(0) %>%  # Force la fermeture des polygones
  st_make_valid()

# Harmonisation des CRS
transects <- st_transform(transects, st_crs(plots))
ferry_lines <- st_transform(ferry_lines, st_crs(plots))

# 2. Création du tibble d'information ----
# Intersection avec les transects et ferry lines
plots_transects <- st_intersection(plots, transects) %>%
  mutate(source = "transects")
plots_ferry <- st_intersection(plots, ferry_lines) %>%
  mutate(source = "ferry")

plots_info <- bind_rows(
  plots_transects %>% st_drop_geometry(),
  plots_ferry %>% st_drop_geometry()
) %>%
  mutate(
    utm_info = map(Name, extract_utm_plot),
    utm = map_chr(utm_info, "utm"),
    transect_name = map_chr(utm_info, "plot_name"),
    folder_name = if_else(source == "ferry",
                          map_chr(transect_name, format_ferry_folder),
                          transect_name)
  ) %>%
  select(plot_ref, source, utm, transect_name, folder_name) %>%
  distinct()

# 3. Extraction des données LiDAR ----
# Création du dossier de sortie
output_base <- paste0(path_base, "/LiDAR_functionnal/1_LiDAR_extraction/output")
dir.create(output_base, showWarnings = FALSE, recursive = TRUE)

# Paramètres de traitement
nb.coeurs <- parallel::detectCores()
set_lidr_threads(nb.coeurs / 2)

# Boucle d'extraction pour chaque plot
for(i in 1:nrow(plots_info)) {
  plot_info <- plots_info[i,]
  cat("Processing plot:", plot_info$plot_ref, "\n")
  
  # Obtention du chemin LiDAR
  lidar_path <- get_lidar_path(plot_info$source, plot_info$utm, plot_info$folder_name)
  
  # Vérification de l'existence du dossier
  if(!dir.exists(lidar_path)) {
    warning(paste("Directory not found:", lidar_path))
    next
  }
  
  epsg <- get_epsg_from_utm(plot_info$utm)
  ctg <- readLAScatalog(lidar_path)
  st_crs(ctg) <- epsg
  
  las_check(ctg)
  
  # Obtention du plot correspondant et reprojection
  plot_shp <- plots %>%
    filter(plot_ref == plot_info$plot_ref) %>%
    st_transform(st_crs(ctg))
  
  # Création du buffer
  plot_buffer <- st_buffer(plot_shp, 200)
  
  # Création du dossier de sortie pour ce plot
  plot_output <- file.path(output_base, plot_info$plot_ref)
  dir.create(plot_output, showWarnings = FALSE)
  
  # Export du shapefile du plot
  st_write(plot_buffer, 
           paste0(plot_output, "/plot_", plot_info$plot_ref, ".shp"), 
           append = FALSE)
  
  # Configuration du catalogue et extraction
  opt_output_files(ctg) <- file.path(plot_output, paste0("plot_", plot_info$plot_ref))
  opt_progress(ctg) <- TRUE
  opt_chunk_buffer(ctg) <- 200
  
  # Extraction des données LiDAR
  clip_roi(ctg, plot_buffer, append = FALSE)
}

output_info = paste0(path_base, "/LiDAR_functionnal/0_Inventories_plot_preparation/final")
write_csv2(plots_info, paste0(output_info, "/plots_info_NASA2015.csv"))
