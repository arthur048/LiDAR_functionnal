rm(list=ls())
gc()

library(tidyverse)
library(sf)
library(readxl)

# Définition des EPSG UTM complète
utm_epsg <- c(
  "UTM32N" = 32632,
  "UTM32S" = 32732,
  "UTM33S" = 32733,
  "UTM34N" = 32634,
  "UTM34S" = 32734,
  "UTM35N" = 32635,
  "UTM35S" = 32735,
  "UTM36N" = 32636
)

# 1. Lecture des données ----
# Lecture des shapefiles principaux et extraction des IDs

rabi_sf <- st_read("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/shapefile/Rabi_AfriSAR_1ha.shp")
mondah_sf <- st_read("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/shapefile/Mondah_AfriSAR_1ha.shp")
mabounie_sf <- st_read("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/shapefile/Mabounie_AfriSAR_1ha.shp")
lope_sf <- st_read("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/shapefile/Lope_AfriSAR_1ha.shp")

plots_sf <- st_read("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/shapefile/plots_under_lidar_without_ituri_final.shp")
plots_sf_id <- plots_sf %>% pull(ID)
plots_sf_id2 <- plots_sf %>% filter(is.na(ID)) %>% pull(ID2)

# Lecture du fichier plots_info existant
plots_info <- rio::import("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/final/plots_info_NASA2015.csv", sep = ";", dec = ",") %>% as_tibble()

# Lecture des données avec vérification des séparateurs
data_excel <- read_excel("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/field_data/NASA_JPL/db_Africa_NASA.xlsx") %>%
  filter(ID %in% plots_sf_id2, dbh >= 10) %>%
  select(ID, site, species, dbh, h) %>%
  rename(plot_ref = ID, sp = species)

data_malebo <- read_delim("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/field_data/RDC_others/MaleboFieldPlots.csv", 
                          delim = ";", locale = locale(decimal_mark = ",")) %>%
  filter(plot_ref %in% plots_sf_id, dbh >= 10) %>%
  select(plot_ref, site, sp, dbh, h)

data_yangambi_nestor <- read_delim("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/field_data/RDC_others/Yangambi_Nestor_FieldPlots.csv",
                                   delim = ";", locale = locale(decimal_mark = ",")) %>%
  filter(plot_ref %in% plots_sf_id, dbh >= 10) %>%
  select(plot_ref, site, sp, dbh, h)

data_yangambi <- read_excel("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/field_data/RDC_others/YangambiFieldPlots.xlsx") %>%
  filter(plot_ref %in% plots_sf_id, dbh >= 10) %>%
  select(plot_ref, site, sp, dbh, h)

# 2. Combinaison et filtrage des données d'inventaires de terrain ----
data_combined <- bind_rows(
  data_excel,
  data_yangambi,
  data_malebo,
  data_yangambi_nestor
) %>%
  filter(dbh >= 10, (h <= 65 | is.na(h)))

# 3. informations sur les utm et projection ----

# Création des données pour chaque site avec leur UTM correspondant
rabi_info <- tibble(
  plot_ref = rabi_sf$areacode,
  utm = "UTM32S"
)

mondah_info <- tibble(
  plot_ref = mondah_sf$areacode,
  utm = "UTM32N"
)

mabounie_info <- tibble(
  plot_ref = mabounie_sf$areacode,
  utm = "UTM32S"
)

lope_info <- tibble(
  plot_ref = lope_sf$areacode,
  utm = "UTM32S"
)

# Combinaison de toutes les nouvelles informations
new_plots_info <- bind_rows(
  rabi_info,
  mondah_info,
  mabounie_info,
  lope_info
) %>% mutate(source = "AfriSAR")

# Mise à jour du fichier plots_info
plots_info_updated <- plots_info %>%
  # Supprime les anciennes entrées si elles existent
  anti_join(new_plots_info, by = "plot_ref") %>%
  # Ajoute les nouvelles entrées
  bind_rows(new_plots_info) %>%
  # Ajoute la colonne EPSG
  mutate(EPSG = utm_epsg[utm] %>% as.integer()) %>%
  # Supprime les colonnes inutiles
  select(-transect_name, -folder_name)


# 4. traitement spatial et uniformisation des couches ----

# Récupération de l'EPSG de plots_sf pour harmonisation
target_epsg <- st_crs(plots_sf)$epsg

# Reprojection des couches AfriSAR
rabi_sf_reproj <- st_transform(rabi_sf, target_epsg)
mondah_sf_reproj <- st_transform(mondah_sf, target_epsg)
mabounie_sf_reproj <- st_transform(mabounie_sf, target_epsg)
lope_sf_reproj <- st_transform(lope_sf, target_epsg)

# Préparation de plots_sf
plots_sf <- plots_sf %>%
  mutate(plot_ref = if_else(!is.na(ID), ID, ID2))

# Préparation des couches AfriSAR avec plot_ref
afrisar_combined <- bind_rows(
  rabi_sf_reproj %>% 
    mutate(plot_ref = areacode) %>%
    select(plot_ref),
  mondah_sf_reproj %>% 
    mutate(plot_ref = areacode) %>%
    select(plot_ref),
  mabounie_sf_reproj %>% 
    mutate(plot_ref = areacode) %>%
    select(plot_ref),
  lope_sf_reproj %>% 
    mutate(plot_ref = areacode) %>%
    select(plot_ref)
)

# Combinaison de toutes les couches
all_plots_sf <- bind_rows(
  plots_sf %>% select(plot_ref),
  afrisar_combined
) %>%
  # Jointure avec plots_info_updated pour utm et epsg
  left_join(plots_info_updated %>% select(plot_ref, utm, EPSG), by = "plot_ref") %>%
  # Calcul des centroïdes et leurs coordonnées
  mutate(
    centroid = st_centroid(geometry),
    x_coord_EPSG3857 = st_coordinates(centroid)[,1],
    y_coord_EPSG3857 = st_coordinates(centroid)[,2]
  ) %>%
  # Suppression de la colonne centroid temporaire
  select(-centroid) %>%
  # filtre les NA
  filter(!is.na(utm))

# 5. Export des données ----
# Création des dossiers
export_path <- "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/final/"
dir.create(export_path, showWarnings = FALSE)
dir.create(paste0(export_path, "plots_unique"), showWarnings = FALSE)

# Export en GeoPackage
st_write(all_plots_sf, 
         paste0(export_path, "plots_limites.gpkg"), 
         append = FALSE)

# Export des plots uniques avec reprojection dans leur EPSG respectif
plots_list <- unique(all_plots_sf$plot_ref)

walk(plots_list[!is.na(plots_list)], function(plot_id) {
  # Sélection du plot
  single_plot <- all_plots_sf %>% 
    filter(plot_ref == plot_id)
  
  # Récupération de l'EPSG correspondant
  target_epsg <- single_plot$EPSG
  
  # Reprojection dans le bon EPSG
  single_plot_reproj <- single_plot %>%
    st_transform(target_epsg)
  
  # Export en GeoPackage
  st_write(single_plot_reproj, 
           paste0(export_path, "plots_unique/", plot_id, ".gpkg"),
           append = FALSE)
})

# 6. Traitement des données Wood Density ----

# Lecture des données
afrisar_agb <- rio::import("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/field_data/AfriSAR/AfriSAR_plotbased_AGB.csv",
                          dec = ".", sep = ";") %>%
  as_tibble() %>%
  select(plot_ref = Area_code, WD)

# Mise à jour de plots_info avec WD
plots_AfriSAR_WD <- plots_info_updated %>%
  left_join(afrisar_agb, by = "plot_ref") %>%
  filter(!is.na(WD)) %>%
  select(plot_ref, WD)

# Export des données mises à jour
rio::export(plots_info_updated, 
            "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/final/plots_info.csv",
            sep = ";", dec = ",", append = FALSE)

rio::export(plots_AfriSAR_WD, 
            "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/final/plots_AfriSAR_WD.csv",
            sep = ";", dec = ",", append = FALSE)


