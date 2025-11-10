# =================================================================
# SCRIPT DE CALCUL DE MÉTRIQUE RELATIVES AUX TROUÉES
# Objectif : Calculer la proportion de pixels considérés comme des trouées 
# en fonction du buffer et du seuil de hauteur
# =================================================================

#### 1️⃣ INITIALISATION ####
# Nettoyage de l'environnement et de la mémoire
rm(list = ls())
gc()

# Définition des packages nécessaires
pkgs <- c("terra", "tidyverse", "sf", "rio", "foreach", "doParallel", "here", "glue", "testthat", "progress")

# Installation et chargement automatique des packages manquants
to_install <- !pkgs %in% installed.packages()
if(any(to_install)) {install.packages(pkgs[to_install])}
inst <- lapply(pkgs, library, character.only = TRUE)

# Définition du répertoire de travail avec here
path0 <- here()
setwd(path0)

#### 2️⃣ FONCTIONS ####
#' Crée un répertoire s'il n'existe pas
#' 
#' @param path Chemin du répertoire à créer
#' @return NULL, crée le répertoire si nécessaire
create_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

#' Calcule la proportion de trouées dans un CHM pour un seuil de hauteur donné
#' 
#' @param chm SpatRaster du Canopy Height Model
#' @param height_threshold Seuil de hauteur en dessous duquel un pixel est considéré comme une trouée
#' @return Proportion de pixels considérés comme des trouées (entre 0 et 1)
calculate_gap_proportion <- function(chm, height_threshold) {
  # Classification binaire : TRUE pour les trouées, FALSE pour la canopée
  gap_mask <- chm < height_threshold
  
  # Calcul des fréquences (nombre de pixels TRUE et FALSE)
  freq_table <- terra::freq(gap_mask) %>%
    as_tibble()
  
  # Cas spécial : aucun pixel n'est une trouée (tous au-dessus du seuil)
  if(!any(freq_table$value == 1)) {
    return(0.0)
  }
  
  # Cas spécial : tous les pixels sont des trouées (tous en dessous du seuil)
  if(!any(freq_table$value == 0)) {
    return(1.0)
  }
  
  # Si aucun pixel n'est NA, on peut calculer directement la proportion
  if(nrow(freq_table) == 2) {
    # Récupération du nombre de pixels considérés comme des trouées (TRUE)
    gap_pixels <- freq_table %>%
      filter(value == 1) %>%
      pull(count)
    
    # Nombre total de pixels
    total_pixels <- sum(freq_table$count)
    
    # Calcul de la proportion
    proportion <- gap_pixels / total_pixels
  } else {
    # Si présence de NA, gestion spéciale
    gap_pixels <- freq_table %>%
      filter(value == 1) %>%
      pull(count)
    
    valid_pixels <- freq_table %>%
      filter(!is.na(value)) %>%
      summarise(total = sum(count)) %>%
      pull(total)
    
    # Calcul de la proportion sur les pixels valides uniquement
    proportion <- gap_pixels / valid_pixels
  }
  
  return(proportion)
}

#### 3️⃣ CHARGEMENT DES DONNÉES ####
# Création du répertoire de sortie
path_output <- file.path(path0, "output")
create_dir(path_output)

# Chemin vers les données d'élévation
path_elevation = "here()2_ElevationData/"

# Chemin vers les CHM
path_chm_folder <- file.path(path_elevation, "CHM_final")
path_chm <- list.files(path = path_chm_folder, 
                       pattern = "\\.tif$", 
                       full.names = TRUE)

# Extraction des noms des plots
plot_name_chm <- basename(gsub("\\.tif$", "", path_chm))

# Lecture des informations des plots
plots_info <- read_csv2("here()3_Gaps/plots_infos.csv") %>% filter(plot_ref %in% plot_name_chm)

plot_name = plots_info$plot_ref

# Chemin vers les polygones des plots
path_plot = "here()0_Inventories_plot_preparation/final/plots_unique/" %>%
  list.files(pattern = "\\.gpkg$", full.names = TRUE) %>%
  .[grep(paste(plot_name, collapse = "|"), .)]

#### 4️⃣ PRÉTRAITEMENT ####
# Définition des paramètres
height_aboveground <- seq.int(1, 60)  # Seuils de hauteur à tester
buffers <- c(0, 20, 50, 100, 200)               # Tailles de buffer à appliquer

# Préparation du dataframe pour stocker les résultats
results_gaps_metrics <- tibble(
  plot_name = character(),
  height_aboveground = numeric(),
  buffer = numeric(),
  proportion = numeric()
)

#### 5️⃣ ANALYSE ####
# Initialisation du temps de calcul
time_start_all <- Sys.time()

# Configuration de la parallélisation
num_cores <- parallel::detectCores()
num_workers <- max(1, floor(num_cores * 0.75))  # Utilisation de 75% des cœurs disponibles
cat(glue::glue("Démarrage du traitement parallèle avec {num_workers} cœurs ({round(num_workers/num_cores*100)}% des {num_cores} cœurs disponibles)\n"))

# Initialisation du cluster
cl <- parallel::makeCluster(num_workers)
doParallel::registerDoParallel(cl)

# Affichage de l'information avant le traitement
cat(glue::glue("Traitement parallèle de {length(plot_name)} plots...\n"))
cat(glue::glue("Hauteurs testées : {min(height_aboveground)}-{max(height_aboveground)} m ({length(height_aboveground)} valeurs)\n"))
cat(glue::glue("Buffers testés : {paste(buffers, collapse = ', ')} m\n"))
cat(glue::glue("Nombre total de combinaisons : {length(plot_name) * length(height_aboveground) * length(buffers)}\n\n"))

# Boucle parallèle sur les plots avec foreach
results_list <- foreach::foreach(
  i = seq_along(plot_name),
  .packages = c("terra", "tidyverse", "sf", "glue"),
  .export = c("calculate_gap_proportion", "height_aboveground", "buffers")
) %dopar% {
  time_start_plot <- Sys.time()
  
  # Extraction du nom du plot courant
  plot_name_part <- plot_name[[i]]
  
  # Préparation du dataframe pour les résultats de ce plot
  plot_results <- tibble(
    plot_name = character(),
    height_aboveground = numeric(),
    buffer = numeric(),
    proportion = numeric()
  )
  
  # Chargement du CHM et du polygone du plot
  path_chm_part <- path_chm[basename(path_chm) == paste0(plot_name_part, ".tif")]
  path_plot_part <- path_plot[basename(path_plot) == paste0(plot_name_part, ".gpkg")]
  
  chm_part <- terra::rast(path_chm_part)
  plot_part <- sf::st_read(path_plot_part, quiet = TRUE) %>%
    sf::st_transform(crs = sf::st_crs(chm_part))
  
  # Boucle sur les hauteurs seuils
  for (height in height_aboveground) {
    # Boucle sur les tailles de buffer
    for (buf in buffers) {
      # Découpage du CHM selon le buffer
      chm <- terra::crop(chm_part, plot_part %>% sf::st_buffer(buf), mask = TRUE)
      
      # Calcul de la proportion de trouées
      gap_proportion <- calculate_gap_proportion(chm, height)
      
      # Ajout des résultats au dataframe
      plot_results <- plot_results %>%
        bind_rows(tibble(
          plot_name = plot_name_part,
          height_aboveground = height,
          buffer = buf,
          proportion = gap_proportion
        ))
    }
  }
  
  # Calcul de la durée pour ce plot
  plot_duration <- as.numeric(difftime(Sys.time(), time_start_plot, units = "secs"))
  
  # Retourner les résultats et les informations sur ce plot
  list(
    results = plot_results,
    plot_name = plot_name_part,
    duration = plot_duration
  )
}

# Fermeture du cluster
parallel::stopCluster(cl)

# Combinaison des résultats de tous les workers
cat("Combinaison des résultats de tous les workers...\n")
for (result in results_list) {
  # Ajout des résultats de ce plot au dataframe global
  results_gaps_metrics <- results_gaps_metrics %>%
    bind_rows(result$results)
  
  # Affichage des informations sur ce plot
  cat(glue::glue("Plot '{result$plot_name}' traité en {round(result$duration, 2)} secondes\n"))
}

# Affichage de la progression finale
duration_all <- difftime(Sys.time(), time_start_all, units = "mins")
cat("\n")
cat(glue::glue("========================================================\n"))
cat(glue::glue("Fin de l'ensemble des {length(plot_name)} plots !\n"))
cat(glue::glue("Durée totale d'exécution: {round(as.numeric(duration_all), 2)} minutes\n"))
cat(glue::glue("Gain de parallélisation estimé: {round(num_workers)} x plus rapide\n"))
cat(glue::glue("========================================================\n"))

#### 6️⃣ EXPORT ####
# Export des résultats au format CSV
rio::export(
  results_gaps_metrics,
  file.path(path_output, "results_gaps_metrics.csv"),
  sep = ";",
  dec = ",",
  append = FALSE
)

# Affichage d'un message de fin
cat(glue::glue("Calcul des proportions de trouées terminé.\n"))
cat(glue::glue("Résultats détaillés sauvegardés dans {file.path(path_output, 'results_gaps_metrics.csv')}\n"))

rio::export(
  tibble(plot_name),
  file.path("here()final_plot_name.csv"),
  sep = ";",
  dec = ",",
  append = FALSE
)
