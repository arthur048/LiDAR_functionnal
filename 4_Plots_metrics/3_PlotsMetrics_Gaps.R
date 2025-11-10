# =================================================================
# CALCUL DES MÉTRIQUES DE TROUÉES POUR L'ANALYSE DE LA STRUCTURE FORESTIÈRE
# Objectif : Calculer une suite simplifiée de métriques caractérisant
#            la structure verticale forestière en relation avec la densité du bois
# =================================================================

#### 1️⃣ INITIALISATION ####
# Nettoyage de l'environnement et de la mémoire
rm(list = ls())
gc()

# Définition des packages nécessaires
pkgs <- c("terra", "tidyverse", "sf", "rio", "foreach", "doParallel",
          "mgcv", "minpack.lm", "segmented", "splines", "waveslim",
          "fields", "zoo", "progress", "here")

# Installation et chargement automatique des packages manquants
to_install <- !pkgs %in% installed.packages()
if(any(to_install)) {install.packages(pkgs[to_install])}
inst <- lapply(pkgs, library, character.only = TRUE)

# Définition des répertoires de travail
path_base <- here()
path_gaps <- here("3_Gaps")
path_metrics <- here("4_Plots_metrics")
dir.create(path_metrics, recursive = TRUE, showWarnings = FALSE)

setwd(path_metrics)

#### 2️⃣ FONCTIONS ####
#' Crée un répertoire s'il n'existe pas
#'
#' @param path Chemin du répertoire à créer
#' @return NULL
create_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

#' Calcule la hauteur à laquelle la proportion de trouées atteint un seuil donné
#' Cette fonction est essentielle pour déterminer les points caractéristiques
#' de la courbe de transition (h10, h25, h50, h75, h90) qui décrivent la 
#' structure verticale de la canopée forestière.
#' 
#' @param data Données d'un plot filtré 
#' @param threshold Seuil de proportion (entre 0 et 1)
#' @return Hauteur à laquelle la proportion atteint le seuil
gap_height_at_threshold <- function(data, threshold) {
  # Vérifier que les données ont au moins 2 lignes
  if (nrow(data) < 2) return(NA)
  
  # Ordonner les données par hauteur
  data <- data %>% arrange(height_aboveground)
  
  # Rechercher où la proportion passe au-dessus du seuil
  for (i in 2:nrow(data)) {
    if (data$proportion[i-1] < threshold && data$proportion[i] >= threshold) {
      # Interpolation linéaire
      height1 <- data$height_aboveground[i-1]
      height2 <- data$height_aboveground[i]
      prop1 <- data$proportion[i-1]
      prop2 <- data$proportion[i]
      
      # Formule d'interpolation
      height_at_threshold <- height1 + (threshold - prop1) * (height2 - height1) / (prop2 - prop1)
      return(height_at_threshold)
    }
  }
  
  # Si tous les points sont en dessous du seuil (jamais atteint)
  if (all(data$proportion < threshold)) return(NA)
  
  # Si tous les points sont au-dessus du seuil (atteint dès le début)
  if (all(data$proportion >= threshold)) return(min(data$height_aboveground))
  
  return(NA)
}

#' Calcule l'aire sous la courbe de proportion vs hauteur
#' L'AUC représente l'intégrale de la proportion de trouées sur une plage de hauteur,
#' fournissant une mesure de l'ouverture de la canopée dans la structure verticale spécifiée.
#' 
#' @param data Données d'un plot filtré
#' @param min_prop Proportion minimale à considérer (défaut: 0)
#' @param max_prop Proportion maximale à considérer (défaut: 1)
#' @param min_height Hauteur minimale à considérer (si NULL, déterminée par min_prop)
#' @param max_height Hauteur maximale à considérer (si NULL, déterminée par max_prop)
#' @param use_proportion Booléen indiquant si on utilise min_prop/max_prop (TRUE) ou min_height/max_height (FALSE)
#' @param normalize Booléen indiquant si l'AUC doit être normalisée (défaut: TRUE)
#' @return Aire sous la courbe, normalisée ou non selon le paramètre
calculate_auc <- function(data, min_prop = 0, max_prop = 1, min_height = NULL, max_height = NULL, 
                          use_proportion = TRUE, normalize = TRUE) {
  # Ordonner les données par hauteur
  data <- data %>% arrange(height_aboveground)
  
  # Si on utilise les proportions pour définir les bornes
  if (use_proportion) {
    # Convertir les proportions en hauteurs
    min_height <- gap_height_at_threshold(data, min_prop)
    max_height <- gap_height_at_threshold(data, max_prop)
    
    # Si les hauteurs ne peuvent pas être déterminées, utiliser les extrêmes
    if (is.na(min_height)) min_height <- min(data$height_aboveground)
    if (is.na(max_height)) max_height <- max(data$height_aboveground)
  } else {
    # Utiliser les hauteurs spécifiées, ou les extrêmes si non spécifiées
    if (is.null(min_height)) min_height <- min(data$height_aboveground)
    if (is.null(max_height)) max_height <- max(data$height_aboveground)
  }
  
  # S'assurer que min_height et max_height sont dans les limites des données
  min_height <- max(min_height, min(data$height_aboveground))
  max_height <- min(max_height, max(data$height_aboveground))
  
  # Filtrer les données pour ne garder que celles dans la plage spécifiée
  filtered_data <- data %>% 
    filter(height_aboveground >= min_height, 
           height_aboveground <= max_height)
  
  # Si aucune donnée après filtrage, retourner 0
  if (nrow(filtered_data) <= 1) {
    return(0)
  }
  
  # Calculer l'AUC avec la méthode des trapèzes
  auc <- 0
  for (i in 2:nrow(filtered_data)) {
    h1 <- filtered_data$height_aboveground[i-1]
    h2 <- filtered_data$height_aboveground[i]
    p1 <- filtered_data$proportion[i-1]
    p2 <- filtered_data$proportion[i]
    
    # Aire du trapèze
    auc <- auc + (h2 - h1) * (p1 + p2) / 2
  }
  
  # Normaliser par la plage de hauteur si demandé
  if (normalize) {
    height_range <- max_height - min_height
    auc <- auc / height_range
  }
  
  return(auc)
}

#' Décompose une courbe de proportion de trouées en segments fonctionnels
#' Simplifiée pour ne garder que l'épaisseur et la pente de chaque zone
#'
#' @param height Vecteur des hauteurs
#' @param prop Vecteur des proportions de trouées
#' @param threshold_low Seuil inférieur pour la zone basale (défaut: 0.25)
#' @param threshold_high Seuil supérieur pour la zone de transition (défaut: 0.75)
#' @return Liste contenant les segments et métriques simplifiées par zone
decompose_curve_segments <- function(height, prop, threshold_low = 0.25, threshold_high = 0.75) {
  # S'assurer que les données sont triées par hauteur
  data <- tibble(height = height, prop = prop) %>%
    arrange(height)
  
  # Ajuster une spline pour l'interpolation
  spline_fit <- smooth.spline(data$height, data$prop)
  
  # Prédire sur une grille régulière pour faciliter l'analyse
  height_grid <- seq(min(data$height), max(data$height), length.out = 200)
  prop_grid <- predict(spline_fit, height_grid)$y
  
  # Trouver les hauteurs correspondant aux seuils
  h_low <- approx(prop_grid, height_grid, xout = threshold_low)$y
  h_high <- approx(prop_grid, height_grid, xout = threshold_high)$y
  
  # Diviser la courbe en segments
  segment_basal <- tibble(
    height = height_grid[height_grid <= h_low],
    prop = prop_grid[height_grid <= h_low],
    segment = "basal"
  )
  
  segment_transition <- tibble(
    height = height_grid[height_grid > h_low & height_grid <= h_high],
    prop = prop_grid[height_grid > h_low & height_grid <= h_high],
    segment = "transition"
  )
  
  segment_upper <- tibble(
    height = height_grid[height_grid > h_high],
    prop = prop_grid[height_grid > h_high],
    segment = "upper"
  )
  
  # Combiner les segments
  segments <- bind_rows(segment_basal, segment_transition, segment_upper)
  
  # Calculer des métriques simplifiées pour chaque segment
  metrics <- tibble(
    zone = c("basal", "transition", "upper"),
    thickness = c(
      max(segment_basal$height) - min(segment_basal$height),
      max(segment_transition$height) - min(segment_transition$height),
      max(segment_upper$height) - min(segment_upper$height)
    )
  )
  
  # Convertir en format pour la fonction calculate_auc
  data_for_auc <- tibble(
    height_aboveground = data$height,
    proportion = data$prop
  )
  
  # Calculer l'AUC pour chaque segment, à la fois normalisée et non normalisée
  basal_auc_norm <- calculate_auc(data_for_auc, 
                                  min_prop = 0, 
                                  max_prop = threshold_low, 
                                  use_proportion = TRUE, 
                                  normalize = TRUE)
  
  basal_auc <- calculate_auc(data_for_auc, 
                             min_prop = 0, 
                             max_prop = threshold_low, 
                             use_proportion = TRUE, 
                             normalize = FALSE)
  
  transition_auc_norm <- calculate_auc(data_for_auc, 
                                       min_prop = threshold_low, 
                                       max_prop = threshold_high, 
                                       use_proportion = TRUE, 
                                       normalize = TRUE)
  
  transition_auc <- calculate_auc(data_for_auc, 
                                  min_prop = threshold_low, 
                                  max_prop = threshold_high, 
                                  use_proportion = TRUE, 
                                  normalize = FALSE)
  
  upper_auc_norm <- calculate_auc(data_for_auc, 
                                  min_prop = threshold_high, 
                                  max_prop = 1, 
                                  use_proportion = TRUE, 
                                  normalize = TRUE)
  
  upper_auc <- calculate_auc(data_for_auc, 
                             min_prop = threshold_high, 
                             max_prop = 1, 
                             use_proportion = TRUE, 
                             normalize = FALSE)
  
  # Ajouter l'AUC à la table des métriques
  metrics <- metrics %>%
    mutate(
      auc_norm = c(basal_auc_norm, transition_auc_norm, upper_auc_norm),
      auc = c(basal_auc, transition_auc, upper_auc)
    )
  
  # Calculer les pentes pour chaque zone de façon simple
  if (nrow(segment_basal) >= 2) {
    basal_slope <- lm(prop ~ height, data = segment_basal)
    metrics$slope_basal <- coef(basal_slope)[2]
  } else {
    metrics$slope_basal <- NA
  }
  
  if (nrow(segment_transition) >= 2) {
    transition_slope <- lm(prop ~ height, data = segment_transition)
    metrics$slope_transition <- coef(transition_slope)[2]
  } else {
    metrics$slope_transition <- NA
  }
  
  if (nrow(segment_upper) >= 2) {
    upper_slope <- lm(prop ~ height, data = segment_upper)
    metrics$slope_upper <- coef(upper_slope)[2]
  } else {
    metrics$slope_upper <- NA
  }
  
  # Retourner les segments et les métriques
  return(list(segments = segments, metrics = metrics))
}

#' Calcule les dérivées d'une courbe de proportion de trouées
#' Simplifiée pour ne garder que height max d1/d2 et max d1/d2
#'
#' @param height Vecteur des hauteurs
#' @param prop Vecteur des proportions de trouées
#' @param smooth_span Paramètre de lissage pour la spline (défaut: 0.4)
#' @return Dataframe avec métriques simplifiées
calculate_derivatives <- function(height, prop, smooth_span = 0.4) {
  # S'assurer que les données sont triées par hauteur
  data <- tibble(height = height, prop = prop) %>%
    arrange(height)
  
  # Ajuster une spline pour le lissage
  spline_fit <- smooth.spline(data$height, data$prop, spar = smooth_span)
  
  # Prédire sur une grille régulière
  height_grid <- seq(min(data$height), max(data$height), length.out = 200)
  prop_smooth <- predict(spline_fit, height_grid)$y
  
  # Calculer la dérivée première
  d1 <- diff(prop_smooth) / diff(height_grid)
  height_d1 <- (height_grid[-1] + height_grid[-length(height_grid)]) / 2
  
  # Calculer la dérivée seconde
  d2 <- diff(d1) / diff(height_d1)
  height_d2 <- (height_d1[-1] + height_d1[-length(height_d1)]) / 2
  
  # Trouver les points caractéristiques
  max_d1_idx <- which.max(d1)
  max_d2_idx <- which.max(d2)
  min_d2_idx <- which.min(d2)
  
  height_max_d1 <- height_d1[max_d1_idx]
  height_max_d2 <- height_d2[max_d2_idx]
  height_min_d2 <- height_d2[min_d2_idx]
  
  # Calculer des métriques simplifiées sur les dérivées
  deriv_metrics <- tibble(
    height_max_d1 = height_max_d1,
    max_d1 = d1[max_d1_idx],
    height_max_d2 = height_max_d2,
    max_d2 = d2[max_d2_idx],
    height_min_d2 = height_min_d2,
    min_d2 = d2[min_d2_idx]
  )
  
  # Retourner les métriques
  return(deriv_metrics)
}

#### 3️⃣ CHARGEMENT DES DONNÉES ####
# Création du répertoire de sortie
create_dir(path_metrics)

# Charger les métriques de trouées
gaps_metrics <- rio::import(file.path(path_gaps, "output", "results_gaps_metrics.csv"), 
                            sep = ";", dec = ",") %>%
  as_tibble()

# Charger les informations des plots
plots_info <- rio::import(file.path(path_metrics, "PlotsMetrics_FieldDataInventories.csv"), 
                          sep = ";", dec = ",") %>%
  as_tibble() %>%
  rename_with(~ ifelse(tolower(.) == "plot_ref", "plot_name", .))

# Récupérer les noms des plots
plot_names <- read.csv2(here("final_plot_name.csv")) %>%
  pull(plot_name)

# Filtrer les informations des plots pour ne garder que ceux de notre étude
plots_info <- plots_info %>%
  filter(plot_name %in% plot_names)

#### 4️⃣ PRÉTRAITEMENT ####
# Définir le buffer à utiliser pour l'analyse principale
buffer_main <- 0  # Buffer de 0m pour l'analyse principale

# Filtrer les données pour le buffer principal
buffer_data <- gaps_metrics %>%
  filter(buffer == buffer_main)

# Vérifier que tous les plots ont des données pour ce buffer
missing_plots <- setdiff(plot_names, unique(buffer_data$plot_name))
if (length(missing_plots) > 0) {
  warning(paste("Les plots suivants n'ont pas de données pour le buffer", buffer_main, "m:", 
                paste(missing_plots, collapse = ", ")))
}

# Joindre avec les infos de plots pour intégrer la densité du bois (WD)
if ("WD_BA" %in% colnames(plots_info)) {
  combined_data <- buffer_data %>%
    left_join(plots_info %>% dplyr::select(plot_name, WD_BA), by = "plot_name")
} else {
  warning("La colonne WD_BA n'est pas disponible dans le fichier plots_info.")
  combined_data <- buffer_data
}

# Création du dataframe unique pour stocker les résultats
plot_metrics <- tibble(
  plot_name = character(),
  WD_BA = numeric(),
  
  # Métriques de base
  prop_at_5m = numeric(),
  prop_at_10m = numeric(),
  prop_at_15m = numeric(),
  prop_at_20m = numeric(),
  prop_at_25m = numeric(),
  prop_at_30m = numeric(),
  prop_at_35m = numeric(),
  prop_at_40m = numeric(),
  h10 = numeric(),
  h25 = numeric(),
  h50 = numeric(),
  h75 = numeric(),
  h90 = numeric(),
  h100 = numeric(),
  
  # Aires sous la courbe (normalisées et non normalisées)
  auc_norm = numeric(),      # AUC normalisée pour la courbe entière
  auc = numeric(),           # AUC non normalisée pour la courbe entière
  auc_0_25_norm = numeric(), # AUC normalisée entre 0% et 25% de trouées
  auc_0_25 = numeric(),      # AUC non normalisée entre 0% et 25% de trouées
  auc_25_75_norm = numeric(), # AUC normalisée entre 25% et 75% de trouées
  auc_25_75 = numeric(),      # AUC non normalisée entre 25% et 75% de trouées
  auc_75_100_norm = numeric(), # AUC normalisée entre 75% et 100% de trouées
  auc_75_100 = numeric(),      # AUC non normalisée entre 75% et 100% de trouées
  
  # Métriques zonales simplifiées
  basal_thickness = numeric(),
  transition_thickness = numeric(),
  upper_thickness = numeric(),
  basal_auc_norm = numeric(),
  basal_auc = numeric(),
  transition_auc_norm = numeric(),
  transition_auc = numeric(),
  upper_auc_norm = numeric(),
  upper_auc = numeric(),
  slope_basal = numeric(),
  slope_transition = numeric(),
  slope_upper = numeric(),
  
  # Métriques différentielles simplifiées
  height_max_d1 = numeric(),
  max_d1 = numeric(),
  height_max_d2 = numeric(),
  max_d2 = numeric(),
  height_min_d2 = numeric(),
  min_d2 = numeric()
)

#### 5️⃣ ANALYSE ####
# Initialisation de la progression
cat(paste0("Calcul des métriques pour ", length(plot_names), " plots...\n"))
pb <- progress_bar$new(
  format = "  traitement [:bar] :percent | Plots traités: :current/:total | Temps restant: :eta",
  total = length(plot_names),
  clear = FALSE,
  width = 100
)

# Mettre à jour la boucle principale sur les plots
for (plot_i in plot_names) {
  # Filtrer les données pour ce plot
  plot_data <- combined_data %>%
    filter(plot_name == plot_i)
  
  # Récupérer WD_BA
  WD_BA_val <- NA
  if ("WD_BA" %in% names(plot_data)) {
    WD_BA_val <- plot_data$WD_BA[1]
  }
  
  #### Calcul des métriques de base ####
  
  # 1. Proportion à hauteurs spécifiques
  prop_5m <- plot_data %>% 
    filter(height_aboveground == 5) %>% 
    pull(proportion)
  
  prop_10m <- plot_data %>% 
    filter(height_aboveground == 10) %>% 
    pull(proportion)
  
  prop_15m <- plot_data %>% 
    filter(height_aboveground == 15) %>% 
    pull(proportion)
  
  prop_20m <- plot_data %>% 
    filter(height_aboveground == 20) %>% 
    pull(proportion)
  
  prop_25m <- plot_data %>% 
    filter(height_aboveground == 25) %>% 
    pull(proportion)
  
  prop_30m <- plot_data %>% 
    filter(height_aboveground == 30) %>% 
    pull(proportion)
  
  prop_35m <- plot_data %>% 
    filter(height_aboveground == 35) %>% 
    pull(proportion)
  
  prop_40m <- plot_data %>% 
    filter(height_aboveground == 40) %>% 
    pull(proportion)
  
  # 2. Hauteurs à seuils spécifiques (métriques fondamentales)
  h10 <- gap_height_at_threshold(plot_data, 0.10)
  h25 <- gap_height_at_threshold(plot_data, 0.25)
  h50 <- gap_height_at_threshold(plot_data, 0.50)
  h75 <- gap_height_at_threshold(plot_data, 0.75)
  h90 <- gap_height_at_threshold(plot_data, 0.90)
  h100 <- gap_height_at_threshold(plot_data, 1)
  
  # 3. Aires sous la courbe (normalisées et non normalisées)
  # Pour la courbe entière
  auc_norm <- calculate_auc(plot_data, min_prop = 0, max_prop = 1, use_proportion = TRUE, normalize = TRUE)
  auc <- calculate_auc(plot_data, min_prop = 0, max_prop = 1, use_proportion = TRUE, normalize = FALSE)
  
  # Pour les segments 0-25%, 25-75%, 75-100%
  auc_0_25_norm <- calculate_auc(plot_data, min_prop = 0, max_prop = 0.25, use_proportion = TRUE, normalize = TRUE)
  auc_0_25 <- calculate_auc(plot_data, min_prop = 0, max_prop = 0.25, use_proportion = TRUE, normalize = FALSE)
  
  auc_25_75_norm <- calculate_auc(plot_data, min_prop = 0.25, max_prop = 0.75, use_proportion = TRUE, normalize = TRUE)
  auc_25_75 <- calculate_auc(plot_data, min_prop = 0.25, max_prop = 0.75, use_proportion = TRUE, normalize = FALSE)
  
  auc_75_100_norm <- calculate_auc(plot_data, min_prop = 0.75, max_prop = 1, use_proportion = TRUE, normalize = TRUE)
  auc_75_100 <- calculate_auc(plot_data, min_prop = 0.75, max_prop = 1, use_proportion = TRUE, normalize = FALSE)
  
  # 4. Décomposition fonctionnelle
  segments_result <- decompose_curve_segments(plot_data$height_aboveground, plot_data$proportion)
  
  # 5. Analyse différentielle
  derivatives_result <- calculate_derivatives(plot_data$height_aboveground, plot_data$proportion)
  
  # Ajouter au dataframe de métriques
  plot_metrics <- plot_metrics %>%
    bind_rows(tibble(
      plot_name = plot_i,
      WD_BA = WD_BA_val,
      
      # Métriques de base
      prop_at_5m = prop_5m,
      prop_at_10m = prop_10m,
      prop_at_15m = prop_15m,
      prop_at_20m = prop_20m,
      prop_at_25m = prop_25m,
      prop_at_30m = prop_30m,
      prop_at_35m = prop_35m,
      prop_at_40m = prop_40m,
      h10 = h10,
      h25 = h25,
      h50 = h50,
      h75 = h75,
      h90 = h90,
      h100 = h100,
      
      # Aires sous la courbe
      auc_norm = auc_norm,
      auc = auc,
      auc_0_25_norm = auc_0_25_norm,
      auc_0_25 = auc_0_25,
      auc_25_75_norm = auc_25_75_norm,
      auc_25_75 = auc_25_75,
      auc_75_100_norm = auc_75_100_norm,
      auc_75_100 = auc_75_100,
      
      # Métriques zonales simplifiées
      basal_thickness = segments_result$metrics %>% filter(zone == "basal") %>% pull(thickness),
      transition_thickness = segments_result$metrics %>% filter(zone == "transition") %>% pull(thickness),
      upper_thickness = segments_result$metrics %>% filter(zone == "upper") %>% pull(thickness),
      basal_auc_norm = segments_result$metrics %>% filter(zone == "basal") %>% pull(auc_norm),
      basal_auc = segments_result$metrics %>% filter(zone == "basal") %>% pull(auc),
      transition_auc_norm = segments_result$metrics %>% filter(zone == "transition") %>% pull(auc_norm),
      transition_auc = segments_result$metrics %>% filter(zone == "transition") %>% pull(auc),
      upper_auc_norm = segments_result$metrics %>% filter(zone == "upper") %>% pull(auc_norm),
      upper_auc = segments_result$metrics %>% filter(zone == "upper") %>% pull(auc),
      slope_basal = segments_result$metrics$slope_basal[1],
      slope_transition = segments_result$metrics$slope_transition[2],
      slope_upper = segments_result$metrics$slope_upper[3],
      
      # Métriques différentielles simplifiées
      height_max_d1 = derivatives_result$height_max_d1,
      max_d1 = derivatives_result$max_d1,
      height_max_d2 = derivatives_result$height_max_d2,
      max_d2 = derivatives_result$max_d2,
      height_min_d2 = derivatives_result$height_min_d2,
      min_d2 = derivatives_result$min_d2
    ))
  
  # Mettre à jour la barre de progression
  pb$tick()
}

#### 6️⃣ EXPORT ET VALIDATION ####
# Exporter les métriques calculées
rio::export(
  plot_metrics,
  file.path(path_metrics, "plot_metrics_gaps.csv"),
  sep = ";",
  dec = ",",
  append = FALSE
)

# Vérification rapide des corrélations avec WD
if (any(!is.na(plot_metrics$WD_BA))) {
  # Calculer les corrélations avec WD pour toutes les métriques numériques
  correlation_vars <- names(plot_metrics)[sapply(plot_metrics, is.numeric)]
  correlation_vars <- setdiff(correlation_vars, c("WD_BA"))
  
  correlations <- data.frame(
    metric = character(),
    correlation = numeric(),
    p_value = numeric()
  )
  
  for (var in correlation_vars) {
    # Données valides
    valid_data <- !is.na(plot_metrics[[var]]) & !is.na(plot_metrics$WD_BA)
    
    # Vérifier suffisamment de données
    if (sum(valid_data) >= 3) {
      # Pearson
      pearson_test <- tryCatch({
        cor.test(plot_metrics[[var]][valid_data], plot_metrics$WD_BA[valid_data], method = "pearson")
      }, error = function(e) {
        list(estimate = NA, p.value = NA)
      })
      
      # Ajouter les résultats
      correlations <- correlations %>%
        bind_rows(tibble(
          metric = var,
          correlation = ifelse(is.na(pearson_test$estimate), NA, pearson_test$estimate),
          p_value = ifelse(is.na(pearson_test$p.value), NA, pearson_test$p.value)
        ))
    }
  }
  
  # Afficher les 10 métriques les plus corrélées
  cat("\nTop 10 métriques les plus corrélées avec WD_BA:\n")
  print(correlations %>% 
          filter(!is.na(correlation)) %>%
          arrange(desc(abs(correlation))) %>% 
          head(10))
} else {
  warning("Aucune valeur de WD_BA disponible pour calculer les corrélations.")
}
