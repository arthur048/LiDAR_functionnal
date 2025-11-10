# =================================================================
# SCRIPT D'ANALYSE DE LA DISTRIBUTION DES TROUÉES FORESTIÈRES
# Objectif : Analyser la distribution de taille des trouées pour 
#            différentes hauteurs et zones tampons
# =================================================================

# === NETTOYAGE DE L'ENVIRONNEMENT ===
rm(list=ls())
gc()

# === INITIALISATION DES PACKAGES ===
# Packages nécessaires :
# - terra : manipulation de rasters
# - tidyverse : manipulation de données
# - sf : manipulation de données spatiales
# - rio : export de données
# - foreach/doParallel : parallélisation (non utilisée actuellement)
pkgs = c("terra", "tidyverse", "sf", "rio", "foreach", "doParallel")

# Installation et chargement automatique des packages manquants
to_install = !pkgs %in% installed.packages()
if(any(to_install)) {install.packages(pkgs[to_install])}
inst = lapply(pkgs, library, character.only = TRUE) # load them

# === DÉFINITION DU RÉPERTOIRE DE TRAVAIL ===
path0 = "here()3_Gaps/"
setwd(path0)

# === FONCTIONS UTILITAIRES ===
# Fonction pour créer un répertoire s'il n'existe pas
create_dir = function(path) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}

# === FONCTION PRINCIPALE D'ANALYSE DES TROUÉES ===
# Cette fonction calcule la distribution des tailles des trouées
# Entrées : 
# - gaps_layer : raster de trouées
# - type_of_return : "frequency" pour distribution, autre pour données brutes
# Sortie : tibble avec distribution des tailles ou données brutes des trouées
get_gap_size_frequency_distribution = function(gaps_layer, type_of_return = "frequency"){
  
  # Construction personnelle sur base de l'article JPL
  
  gaps_tibble = gaps_layer %>%
    freq() %>%
    as_tibble() %>%
    mutate(count = count * terra::res(gaps_layer)[1]^2, gap_id = seq_along(value)) %>%
    select(gap_id, gap_area = count, gap_id_origin = value, -layer) 
  
  gaps_freq_tibble = gaps_tibble %>%
    count(gap_area, name = "gap_freq") %>%
    ungroup()
  
  
  if(type_of_return == "frequency"){
    return(gaps_freq_tibble)
  } else {
    return(gaps_tibble)
  }
}

# === LECTURE ET PRÉPARATION DES DONNÉES ===
# Création du dossier de sortie
path_output = paste0(path0,"output/")
lapply(path_output, create_dir)

# Chemin vers les données d'élévation
path_elevation = "here()2_ElevationData/"

# Lecture des CHM et création des noms de plots
path_chm_folder = paste0(path_elevation, "CHM_final/")
path_chm = lapply(path_chm_folder, function(folder) {
  list.files(path = folder, pattern = "\\.tif$", full.names = TRUE)
}) %>% unlist()

# Lecture des informations des plots
plots_info <- read_csv2("here()0_Inventories_plot_preparation/final/plots_info.csv")

plot_name = basename(gsub("\\.tif$", "", path_chm))

path_outputs = paste0(path_output, plot_name) ; lapply(path_outputs, create_dir)

path_plot = "here()0_Inventories_plot_preparation/final/plots_unique/" %>%
  list.files(pattern = "\\.gpkg$", full.names = TRUE) %>%
  .[grep(paste(plot_name, collapse = "|"), .)]

# === PARAMÈTRES D'ANALYSE ===
# Définition des hauteurs à analyser (1-45m)
height_aboveground = seq.int(1:45)

# Définition des zones tampons (0m, 20m, 50m)
buffers = c(0,20,50)

# === INITIALISATION DES STRUCTURES DE DONNÉES ===
# Création d'une liste pour stocker les résultats intermédiaires
gaps_freq_list = setNames(lapply(plot_name, function(x) {
  setNames(vector("list", length(height_aboveground)), height_aboveground)
}), plot_name)

# Tibble pour stocker tous les résultats
results_gaps_freq = tibble()

# === BOUCLE PRINCIPALE DE TRAITEMENT ===
# Mesure du temps total
time_start_all = Sys.time()

# Boucle sur chaque plot
i = 3 ; height = 10 ; buf = 50
for (i in 1:length(plot_name)) {
  # Le code de la boucle reste identique
  # Il réalise :
  # 1. Lecture du plot
  # 2. Pour chaque hauteur (1-45m)
  #    - Lecture du raster de trouées
  #    - Pour chaque buffer (0,20,50m)
  #      * Découpage du raster
  #      * Calcul des fréquences
  #      * Stockage des résultats
  # 3. Export des résultats
  
  time_start_plot = Sys.time()
  
  plot_name_part = plot_name[[i]]
  
  path_plot_part = path_plot[basename(path_plot) == paste0(plot_name_part, ".gpkg")]
  plot = st_read(path_plot_part, quiet = T)
  
  path_output_part = path_outputs[basename(path_outputs) == plot_name_part]
  
  a = 0
  for (height in height_aboveground){
    a = a + 1
    
    gaps_freq_list[[i]][[a]] = tibble()
    
    gaps = rast(paste0(path_output_part, "/", plot_name_part, "_", height,"m_gaps_height.tif"))
    
    for (buf in buffers){
      
      gaps_part = gaps %>% crop(plot %>% st_buffer(buf),
                                mask = TRUE)
      
      gaps_freq = get_gap_size_frequency_distribution(gaps_layer = gaps_part) %>%
        mutate(plot_name = plot_name_part,
               height_aboveground = height,
               buffer = buf) %>%
        select(plot_name, height_aboveground, buffer, everything())
      
      results_gaps_freq = results_gaps_freq %>% #Combinaison des fréquences de trouéess pour la boucle sur les hauteurs
        dplyr::bind_rows(gaps_freq)
      
      # Sauvegarde des résultats intermédiaires dans une liste, au cas où
      
      gaps_freq_list[[i]][[a]] = gaps_freq_list[[i]][[a]] %>% #Combinaison des fréquences de trouéess pour la boucle sur les hauteurs
        dplyr::bind_rows(gaps_freq)
      
    }
  }
  
  rio::export(results_gaps_freq, paste0(path_output, "results_gaps_freq.csv"), sep = ";", dec = ",", append = FALSE)
  
  if(i == length(plot_name)){
    duration_all = Sys.time() - time_start_all  
    print(paste("Fin de l'ensemble des", i, "plots !"))
    print(duration_all)
  } else { print(paste("End of plot number", i, "on", length(plot_name),"; Plot duration :", round(Sys.time() - time_start_plot, 2))) } 
}

# === SORTIES ===
# Le script produit :
# 1. Un fichier CSV 'results_gaps_freq.csv' contenant :
#    - plot_name : identifiant du plot
#    - height_aboveground : hauteur analysée
#    - buffer : taille de la zone tampon
#    - gap_area : taille de la trouée
#    - gap_freq : fréquence de cette taille
# 2. Une liste gaps_freq_list conservant les résultats intermédiaires



