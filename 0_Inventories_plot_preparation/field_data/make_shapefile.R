rm(list=ls())
gc()

#load("")

#installr::updateR()
#update.packages()

# Open the file in RStudio to edit it : ---> usethis::edit_file("~/AppData/Roaming/RStudio/templates/default.R") <---

# Initialisation ----------------------------------------------------------

pkgs = c("future","sysfonts", "showtext","rio","foreach", "doParallel","lidR", "terra", "sf", "tidyverse", "here")

to_install = !pkgs %in% installed.packages() ; if(any(to_install)) {install.packages(pkgs[to_install])} ; inst = lapply(pkgs, library, character.only = TRUE) # load them

# EXTERNAL PATH: Font directory - adjust to your local setup
font_paths(new = "E:/Arthur/OneDrive2/R"); font_add("LM12", regular = "lmroman12-regular.otf", bold = "lmroman12-bold.otf") ; showtext_auto()

# EXTERNAL DATA PATH: Field data from NASA 2015 campaign - adjust to your local setup
path0 = "E:/Arthur/Doctorat_DataAnnexe/LiDAR_RDC_2015/006_Field_Data/" ; setwd(path0)


# -------------------------------------------------------------------------
if (!require(readxl)) install.packages("readxl")

col_types <- c("guess", "guess", "text", "guess", "guess", "guess", "guess", "guess", "guess", "guess", "guess")
db_africa <- read_excel("E:/Arthur/Doctorat_DataAnnexe/LiDAR_RDC_2015/006_Field_Data/db_Africa_NASA.xlsx", col_types = col_types) %>%
  as_tibble()

db_plot_africa = rio::import("E:/Arthur/Doctorat_DataAnnexe/LiDAR_RDC_2015/006_Field_Data/db_plot_Africa_NASA.csv", dec = ".") %>% as_tibble() 
db_meta_africa = rio::import("E:/Arthur/Doctorat_DataAnnexe/LiDAR_RDC_2015/006_Field_Data/meta_Africa_NASA.xlsx", dec = ",") %>% as_tibble()

unique(db_africa$local_ID)

db_plot_africa <- db_plot_africa %>% 
  filter(!is.na(Long) & !is.na(Lat))%>%
  mutate(Long = as.numeric(Long), Lat = as.numeric(Lat))

sf_plot_africa <- st_as_sf(db_plot_africa, coords = c("Long", "Lat"), crs = 4326)
st_write(sf_plot_africa, paste0(path0, "db_Africa_NASA.shp"), delete_layer = TRUE)

# Lire les shapefiles
sf_lidar_transects <- st_read(here("0_Inventories_plot_preparation", "shapefile", "lidar", "DRC01_lidar_transects.shp"))
sf_ferry_lines <- st_read(here("0_Inventories_plot_preparation", "shapefile", "lidar", "DRC02_ferry_lines.shp"))
sf_ferry_lines <- st_make_valid(sf_ferry_lines)

# Effectuer les intersections et ne garder que les colonnes de sf_plot_africa
sf_intersect_transects <- st_intersection(sf_plot_africa, sf_lidar_transects) %>% select(names(sf_plot_africa))
sf_intersect_ferry <- st_intersection(sf_plot_africa, sf_ferry_lines) %>% select(names(sf_plot_africa))

# Ajouter les lignes où la colonne ID est "Ituri Edoro1" ou "Ituri Lenda2"
additional_rows <- sf_plot_africa %>% filter(ID %in% c("Ituri Edoro1", "Ituri Lenda2"))

# Combiner les deux intersections et les lignes additionnelles
sf_intersect_combined <- rbind(sf_intersect_transects, sf_intersect_ferry, additional_rows)

st_write(sf_intersect_combined, here("0_Inventories_plot_preparation", "plots_under_lidar.shp"), delete_layer = TRUE)
