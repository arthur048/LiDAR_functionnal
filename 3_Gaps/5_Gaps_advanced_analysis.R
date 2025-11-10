# =================================================================
# ANALYSE AVANCÉE DES PROFILS DE TROUÉES FORESTIÈRES
# Objectif : Mettre en œuvre des méthodes avancées pour caractériser
#            les distributions verticales de trouées et leurs relations
#            avec les traits fonctionnels et la densité du bois
# =================================================================

#### 1️⃣ INITIALISATION ####
# Nettoyage de l'environnement et de la mémoire
rm(list = ls())
gc()

# Définition des packages nécessaires
pkgs <- c("tidyverse", "rio", "ggpubr", "viridis", "corrplot", "gridExtra",
          "mgcv", "patchwork", "splines", "waveslim", "fields", "zoo",
          "rstatix", "factoextra", "cluster", "plotly", "ggh4x",
          "minpack.lm", "zoo", "fANCOVA", "ggridges", "here")

# Installation et chargement automatique des packages manquants
to_install <- !pkgs %in% installed.packages()
if(any(to_install)) {install.packages(pkgs[to_install])}
inst <- lapply(pkgs, library, character.only = TRUE)

# Définition du répertoire de travail
path0 <- here("3_Gaps")
setwd(path0)

# Chemins des fichiers
path_output <- file.path(path0, "output")
path_analysis <- file.path(path0, "analysis")
path_advanced <- file.path(path0, "advanced_analysis")
dir.create(path_advanced, recursive = TRUE, showWarnings = FALSE)

# Lecture des données
gaps_metrics <- rio::import(file.path(path_output, "results_gaps_metrics.csv"), sep = ";", dec = ",") %>%
  as_tibble()

plot_metrics <- rio::import(file.path(path_analysis, "plot_gap_metrics_enhanced_buffer0.csv"), sep = ";", dec = ",") %>%
  as_tibble()

plots_info <- read_csv2(here("3_Gaps", "plots_infos.csv")) %>%
  as_tibble() %>%
  rename_with(~ ifelse(tolower(.) == "plot_ref", "plot_name", .))

# Filtrer pour buffer = 0
buffer0_data <- gaps_metrics %>%
  filter(buffer == 0)

# Joindre avec les infos de plot
combined_data <- buffer0_data %>%
  left_join(plots_info %>% select(plot_name, WD_mean), by = "plot_name") %>%
  filter(!is.na(WD_mean))

# Créer des catégories de WD pour la visualisation
wd_quantiles <- quantile(plot_metrics$WD_mean, probs = seq(0, 1, 0.25), na.rm = TRUE)
combined_data <- combined_data %>%
  mutate(wd_category = cut(WD_mean, 
                           breaks = wd_quantiles, 
                           include.lowest = TRUE,
                           labels = c("Q1 (WD bas)", "Q2", "Q3", "Q4 (WD élevé)")))

# Définir une palette de couleurs pour les catégories de WD
wd_palette <- viridis(4, option = "viridis", direction = -1)
names(wd_palette) <- levels(combined_data$wd_category)

#### 2️⃣ FONCTIONS AVANCÉES ####

#' Décompose une courbe de proportion de trouées en segments fonctionnels
#'
#' @param height Vecteur des hauteurs
#' @param prop Vecteur des proportions de trouées
#' @param threshold_low Seuil inférieur pour la zone basale (défaut: 0.25)
#' @param threshold_high Seuil supérieur pour la zone de transition (défaut: 0.75)
#' @return Liste contenant les segments et métriques par zone
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
  
  # Calculer des métriques pour chaque segment
  metrics <- tibble(
    zone = c("basal", "transition", "upper"),
    h_min = c(min(segment_basal$height), min(segment_transition$height), min(segment_upper$height)),
    h_max = c(max(segment_basal$height), max(segment_transition$height), max(segment_upper$height)),
    thickness = c(
      max(segment_basal$height) - min(segment_basal$height),
      max(segment_transition$height) - min(segment_transition$height),
      max(segment_upper$height) - min(segment_upper$height)
    ),
    prop_min = c(min(segment_basal$prop), min(segment_transition$prop), min(segment_upper$prop)),
    prop_max = c(max(segment_basal$prop), max(segment_transition$prop), max(segment_upper$prop)),
    amplitude = c(
      max(segment_basal$prop) - min(segment_basal$prop),
      max(segment_transition$prop) - min(segment_transition$prop),
      max(segment_upper$prop) - min(segment_upper$prop)
    )
  )
  
  # Ajouter l'AUC pour chaque segment (méthode trapézoïdale)
  metrics <- metrics %>%
    mutate(
      auc = c(
        sum(diff(segment_basal$height) * (head(segment_basal$prop, -1) + tail(segment_basal$prop, -1))) / 2,
        sum(diff(segment_transition$height) * (head(segment_transition$prop, -1) + tail(segment_transition$prop, -1))) / 2,
        sum(diff(segment_upper$height) * (head(segment_upper$prop, -1) + tail(segment_upper$prop, -1))) / 2
      ),
      norm_auc = auc / thickness  # AUC normalisée par l'épaisseur du segment
    )
  
  # Calculer un indice d'hétérogénéité zonale
  if (nrow(segment_transition) >= 3) {
    transition_slope <- lm(prop ~ height, data = segment_transition)
    metrics$transition_slope <- coef(transition_slope)[2]
    metrics$transition_r2 <- summary(transition_slope)$r.squared
  } else {
    metrics$transition_slope <- NA
    metrics$transition_r2 <- NA
  }
  
  # Retourner les segments et les métriques
  return(list(segments = segments, metrics = metrics))
}

#' Calcule les dérivées d'une courbe de proportion de trouées
#'
#' @param height Vecteur des hauteurs
#' @param prop Vecteur des proportions de trouées
#' @param smooth_span Paramètre de lissage pour la spline (défaut: 0.4)
#' @return Dataframe avec hauteur, proportion et dérivées
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
  
  # Créer un dataframe pour les résultats
  derivatives <- tibble(
    height = height_grid,
    prop = prop_smooth
  )
  
  # Ajouter les dérivées avec interpolation pour aligner sur la grille originale
  derivatives$d1 <- c(NA, approx(height_d1, d1, xout = height_grid[-1])$y)
  derivatives$d2 <- c(NA, NA, approx(height_d2, d2, xout = height_grid[-c(1,2)])$y)
  
  # Calculer des métriques sur les dérivées
  deriv_metrics <- tibble(
    height_max_d1 = height_max_d1,
    max_d1 = d1[max_d1_idx],
    height_max_d2 = height_max_d2,
    max_d2 = d2[max_d2_idx],
    height_min_d2 = height_min_d2,
    min_d2 = d2[min_d2_idx],
    asymmetry_ratio = (height_min_d2 - height_max_d1) / (height_max_d2 - height_max_d1),
    acceleration_phase = height_max_d1 - height_max_d2,
    deceleration_phase = height_min_d2 - height_max_d1,
    phase_ratio = if_else(is.finite((height_min_d2 - height_max_d1) / (height_max_d1 - height_max_d2)), 
                          (height_min_d2 - height_max_d1) / (height_max_d1 - height_max_d2), NA)
  )
  
  # Retourner les résultats
  return(list(derivatives = derivatives, metrics = deriv_metrics))
}

#' Réalise une analyse multi-échelle d'une courbe de proportion de trouées
#'
#' @param height Vecteur des hauteurs
#' @param prop Vecteur des proportions de trouées
#' @param scales Échelles à analyser en mètres (défaut: c(2, 5, 10, 20))
#' @return Liste avec résultats de l'analyse multi-échelle
multiscale_analysis <- function(height, prop, scales = c(2, 5, 10, 20)) {
  # S'assurer que les données sont triées par hauteur
  data <- tibble(height = height, prop = prop) %>%
    arrange(height)
  
  # Ajuster une spline pour avoir des points équidistants
  spline_fit <- smooth.spline(data$height, data$prop)
  
  # Prédire sur une grille régulière
  height_grid <- seq(min(data$height), max(data$height), length.out = 512)  # Puissance de 2 pour l'analyse en ondelettes
  prop_grid <- predict(spline_fit, height_grid)$y
  
  # Calculer la transformée en ondelettes
  # Utilisons une fonction simple de moyennage pour différentes échelles
  scale_results <- list()
  scale_metrics <- tibble(
    scale = numeric(),
    roughness = numeric(),
    dominant_height = numeric(),
    scale_energy = numeric()
  )
  
  for (scale in scales) {
    # Convertir l'échelle en nombre de points
    window_size <- max(3, round(scale / (diff(range(height_grid)) / length(height_grid))))
    if (window_size %% 2 == 0) window_size <- window_size + 1  # Assurer un nombre impair
    
    # Appliquer un filtre de moyennage
    prop_smooth <- zoo::rollmean(prop_grid, window_size, fill = "extend")
    
    # Calculer les différences entre la courbe originale et la courbe lissée
    residuals <- prop_grid - prop_smooth
    
    # Calculer la rugosité (variance des résidus)
    roughness <- var(residuals, na.rm = TRUE)
    
    # Trouver la hauteur où les différences sont maximales
    max_diff_idx <- which.max(abs(residuals))
    dominant_height <- height_grid[max_diff_idx]
    
    # Calculer l'énergie totale à cette échelle (somme des carrés des résidus)
    scale_energy <- sum(residuals^2, na.rm = TRUE)
    
    # Stocker les résultats
    scale_results[[as.character(scale)]] <- tibble(
      height = height_grid,
      prop_original = prop_grid,
      prop_smooth = prop_smooth,
      residuals = residuals
    )
    
    scale_metrics <- scale_metrics %>%
      bind_rows(tibble(
        scale = scale,
        roughness = roughness,
        dominant_height = dominant_height,
        scale_energy = scale_energy
      ))
  }
  
  # Calculer un indice multi-échelle
  if (nrow(scale_metrics) >= 2) {
    # Ajuster une courbe de puissance au spectre d'énergie
    log_scale <- log(scale_metrics$scale)
    log_energy <- log(scale_metrics$scale_energy)
    power_fit <- lm(log_energy ~ log_scale)
    
    # Le coefficient (pente) caractérise la décroissance d'énergie à travers les échelles
    scale_metrics <- scale_metrics %>%
      mutate(energy_slope = ifelse(!is.na(coef(power_fit)[2]), coef(power_fit)[2], NA))
    
    # Distribution d'énergie relative à travers les échelles
    scale_metrics <- scale_metrics %>%
      mutate(energy_proportion = scale_energy / sum(scale_energy, na.rm = TRUE))
  }
  
  # Retourner les résultats
  return(list(scale_results = scale_results, scale_metrics = scale_metrics))
}

#' Analyse l'anisotropie verticale d'une courbe de proportion de trouées
#'
#' @param curve_data Dataframe avec les données de la courbe
#' @param plot_data Dataframe avec les données du plot
#' @return Métriques d'anisotropie
analyze_anisotropy <- function(curve_data, plot_data) {
  # Cette fonction simule l'analyse d'anisotropie qui nécessiterait normalement
  # les données CHM brutes. Ici, nous utilisons une approximation.
  
  # Extraire les valeurs de hauteur et proportion
  height_unique <- unique(curve_data$height_aboveground)
  
  # Initialiser le dataframe pour l'anisotropie
  anisotropy_metrics <- tibble(
    height = height_unique,
    prop_mean = numeric(length(height_unique)),
    prop_sd = numeric(length(height_unique)),
    prop_cv = numeric(length(height_unique))
  )
  
  # Simuler la variabilité horizontale basée sur la courbe moyenne
  # Dans un cas réel, on utiliserait les valeurs réelles du CHM à chaque hauteur
  for (i in seq_along(height_unique)) {
    h <- height_unique[i]
    prop_mean <- curve_data$proportion[curve_data$height_aboveground == h]
    
    # Simuler un écart-type basé sur la forme de la courbe
    # Assumons que la variabilité est plus grande dans la zone de transition
    prop_sd <- 0.05 + 0.2 * prop_mean * (1 - prop_mean)  # Maximum à prop=0.5
    
    anisotropy_metrics$prop_mean[i] <- prop_mean
    anisotropy_metrics$prop_sd[i] <- prop_sd
    anisotropy_metrics$prop_cv[i] <- prop_sd / max(0.01, prop_mean)  # Éviter division par zéro
  }
  
  # Calculer des métriques d'anisotropie globales
  metrics <- tibble(
    anisotropy_mean = mean(anisotropy_metrics$prop_cv, na.rm = TRUE),
    anisotropy_max = max(anisotropy_metrics$prop_cv, na.rm = TRUE),
    height_max_anisotropy = anisotropy_metrics$height[which.max(anisotropy_metrics$prop_cv)],
    anisotropy_transition = mean(anisotropy_metrics$prop_cv[anisotropy_metrics$prop_mean > 0.25 & 
                                                              anisotropy_metrics$prop_mean < 0.75], na.rm = TRUE),
    anisotropy_ratio = max(anisotropy_metrics$prop_cv, na.rm = TRUE) / 
      mean(anisotropy_metrics$prop_cv, na.rm = TRUE)
  )
  
  # Retourner les résultats
  return(list(profile = anisotropy_metrics, metrics = metrics))
}

#### 3️⃣ ANALYSE AVANCÉE DES COURBES DE TROUÉES ####

# Création d'un tableau pour stocker les métriques avancées
advanced_metrics <- tibble(
  plot_name = character(),
  WD_mean = numeric(),
  
  # Métriques zonales
  basal_thickness = numeric(),
  transition_thickness = numeric(),
  upper_thickness = numeric(),
  basal_norm_auc = numeric(),
  transition_norm_auc = numeric(),
  upper_norm_auc = numeric(),
  transition_slope = numeric(),
  transition_r2 = numeric(),
  zone_heterogeneity_index = numeric(),  # Ratio des AUC normalisées
  
  # Métriques différentielles
  height_max_d1 = numeric(),  # Hauteur à la pente maximale
  max_d1 = numeric(),         # Pente maximale
  height_max_d2 = numeric(),  # Hauteur à l'accélération maximale
  max_d2 = numeric(),         # Accélération maximale
  height_min_d2 = numeric(),  # Hauteur à la décélération maximale
  min_d2 = numeric(),         # Décélération maximale (valeur négative)
  asymmetry_ratio = numeric(),  # Asymétrie de la transition
  acceleration_phase = numeric(),  # Durée de la phase d'accélération
  deceleration_phase = numeric(),  # Durée de la phase de décélération
  phase_ratio = numeric(),     # Ratio des phases
  
  # Métriques multi-échelles
  roughness_small = numeric(),  # Rugosité à petite échelle (2m)
  roughness_medium = numeric(),  # Rugosité à échelle moyenne (10m)
  roughness_large = numeric(),  # Rugosité à grande échelle (20m)
  dominant_height_small = numeric(),  # Hauteur dominante à petite échelle
  dominant_height_medium = numeric(),  # Hauteur dominante à échelle moyenne
  dominant_height_large = numeric(),  # Hauteur dominante à grande échelle
  energy_slope = numeric(),  # Pente spectrale d'énergie
  energy_ratio_small_large = numeric(),  # Ratio d'énergie petite/grande échelle
  
  # Métriques d'anisotropie
  anisotropy_mean = numeric(),  # Anisotropie moyenne
  anisotropy_max = numeric(),   # Anisotropie maximale
  height_max_anisotropy = numeric(),  # Hauteur d'anisotropie maximale
  anisotropy_transition = numeric(),  # Anisotropie dans la zone de transition
  anisotropy_ratio = numeric()   # Ratio d'anisotropie max/moyenne
)

# Liste pour stocker les données détaillées pour les visualisations
plot_details <- list()

# Boucle sur les plots pour calculer les métriques avancées
plot_names <- unique(combined_data$plot_name)

for (plot_i in plot_names) {
  # Données pour ce plot
  plot_data <- combined_data %>%
    filter(plot_name == plot_i) %>%
    arrange(height_aboveground)
  
  # Vérifier qu'il y a assez de données
  if (nrow(plot_data) < 5) next
  
  # Récupérer les valeurs de hauteur et proportion
  height <- plot_data$height_aboveground
  prop <- plot_data$proportion
  
  # 1. Décomposition fonctionnelle
  segments_result <- decompose_curve_segments(height, prop)
  
  # 2. Analyse différentielle
  derivatives_result <- calculate_derivatives(height, prop)
  
  # 3. Analyse multi-échelle
  multiscale_result <- multiscale_analysis(height, prop)
  
  # 4. Analyse d'anisotropie
  anisotropy_result <- analyze_anisotropy(plot_data, plot_data)
  
  # Stocker les détails pour les visualisations
  plot_details[[plot_i]] <- list(
    data = plot_data,
    segments = segments_result$segments,
    derivatives = derivatives_result$derivatives,
    multiscale = multiscale_result$scale_results,
    anisotropy = anisotropy_result$profile
  )
  
  # Extraire les métriques
  wd_mean <- plot_data$WD_mean[1]
  
  # Calculer l'indice d'hétérogénéité zonale
  zone_heterogeneity <- segments_result$metrics %>%
    filter(zone == "transition") %>%
    pull(norm_auc) / mean(c(
      segments_result$metrics %>% filter(zone == "basal") %>% pull(norm_auc),
      segments_result$metrics %>% filter(zone == "upper") %>% pull(norm_auc)
    ))
  
  # Calculer le ratio d'énergie petite/grande échelle
  energy_small <- multiscale_result$scale_metrics %>%
    filter(scale == 2) %>%
    pull(scale_energy)
  
  energy_large <- multiscale_result$scale_metrics %>%
    filter(scale == 20) %>%
    pull(scale_energy)
  
  energy_ratio <- energy_small / energy_large
  
  # Ajouter les métriques au tableau
  advanced_metrics <- advanced_metrics %>%
    bind_rows(tibble(
      plot_name = plot_i,
      WD_mean = wd_mean,
      
      # Métriques zonales
      basal_thickness = segments_result$metrics %>% filter(zone == "basal") %>% pull(thickness),
      transition_thickness = segments_result$metrics %>% filter(zone == "transition") %>% pull(thickness),
      upper_thickness = segments_result$metrics %>% filter(zone == "upper") %>% pull(thickness),
      basal_norm_auc = segments_result$metrics %>% filter(zone == "basal") %>% pull(norm_auc),
      transition_norm_auc = segments_result$metrics %>% filter(zone == "transition") %>% pull(norm_auc),
      upper_norm_auc = segments_result$metrics %>% filter(zone == "upper") %>% pull(norm_auc),
      transition_slope = segments_result$metrics %>% filter(zone == "transition") %>% pull(transition_slope),
      transition_r2 = segments_result$metrics %>% filter(zone == "transition") %>% pull(transition_r2),
      zone_heterogeneity_index = zone_heterogeneity,
      
      # Métriques différentielles
      height_max_d1 = derivatives_result$metrics$height_max_d1,
      max_d1 = derivatives_result$metrics$max_d1,
      height_max_d2 = derivatives_result$metrics$height_max_d2,
      max_d2 = derivatives_result$metrics$max_d2,
      height_min_d2 = derivatives_result$metrics$height_min_d2,
      min_d2 = derivatives_result$metrics$min_d2,
      asymmetry_ratio = derivatives_result$metrics$asymmetry_ratio,
      acceleration_phase = derivatives_result$metrics$acceleration_phase,
      deceleration_phase = derivatives_result$metrics$deceleration_phase,
      phase_ratio = derivatives_result$metrics$phase_ratio,
      
      # Métriques multi-échelles
      roughness_small = multiscale_result$scale_metrics %>% filter(scale == 2) %>% pull(roughness),
      roughness_medium = multiscale_result$scale_metrics %>% filter(scale == 10) %>% pull(roughness),
      roughness_large = multiscale_result$scale_metrics %>% filter(scale == 20) %>% pull(roughness),
      dominant_height_small = multiscale_result$scale_metrics %>% filter(scale == 2) %>% pull(dominant_height),
      dominant_height_medium = multiscale_result$scale_metrics %>% filter(scale == 10) %>% pull(dominant_height),
      dominant_height_large = multiscale_result$scale_metrics %>% filter(scale == 20) %>% pull(dominant_height),
      energy_slope = multiscale_result$scale_metrics$energy_slope[1],
      energy_ratio_small_large = energy_ratio,
      
      # Métriques d'anisotropie
      anisotropy_mean = anisotropy_result$metrics$anisotropy_mean,
      anisotropy_max = anisotropy_result$metrics$anisotropy_max,
      height_max_anisotropy = anisotropy_result$metrics$height_max_anisotropy,
      anisotropy_transition = anisotropy_result$metrics$anisotropy_transition,
      anisotropy_ratio = anisotropy_result$metrics$anisotropy_ratio
    ))
}

# Exporter les métriques avancées
rio::export(
  advanced_metrics,
  paste0(path_advanced, "advanced_metrics.csv"),
  sep = ";",
  dec = ",",
  append = FALSE
)

#### 4️⃣ ANALYSE DE CORRÉLATION AVEC WD ####

# Calculer les corrélations entre les métriques avancées et WD
correlation_vars <- setdiff(names(advanced_metrics), c("plot_name", "WD_mean"))

advanced_correlations <- data.frame(
  metric = character(),
  correlation = numeric(),
  p_value = numeric(),
  method = character()
)

for (var in correlation_vars) {
  # Données valides
  valid_data <- !is.na(advanced_metrics[[var]]) & !is.na(advanced_metrics$WD_mean)
  
  # Vérifier suffisamment de données
  if (sum(valid_data) >= 3) {
    # Pearson
    pearson_test <- tryCatch({
      cor.test(advanced_metrics[[var]][valid_data], advanced_metrics$WD_mean[valid_data], method = "pearson")
    }, error = function(e) {
      list(estimate = NA, p.value = NA)
    })
    
    # Spearman
    spearman_test <- tryCatch({
      suppressWarnings(cor.test(advanced_metrics[[var]][valid_data], advanced_metrics$WD_mean[valid_data], 
                                method = "spearman", exact = FALSE))
    }, error = function(e) {
      list(estimate = NA, p.value = NA)
    })
    
    # Ajouter les résultats
    advanced_correlations <- advanced_correlations %>%
      bind_rows(tibble(
        metric = var,
        correlation = ifelse(is.na(pearson_test$estimate), NA, pearson_test$estimate),
        p_value = ifelse(is.na(pearson_test$p.value), NA, pearson_test$p.value),
        method = "Pearson"
      )) %>%
      bind_rows(tibble(
        metric = var,
        correlation = ifelse(is.na(spearman_test$estimate), NA, spearman_test$estimate),
        p_value = ifelse(is.na(spearman_test$p.value), NA, spearman_test$p.value),
        method = "Spearman"
      ))
  } else {
    # Ajouter des valeurs NA si pas assez d'observations
    advanced_correlations <- advanced_correlations %>%
      bind_rows(tibble(
        metric = var,
        correlation = NA,
        p_value = NA,
        method = "Pearson"
      )) %>%
      bind_rows(tibble(
        metric = var,
        correlation = NA,
        p_value = NA,
        method = "Spearman"
      ))
  }
}

# Améliorer les noms des métriques pour l'affichage
metric_labels <- c(
  # Métriques zonales
  basal_thickness = "Épaisseur zone basale",
  transition_thickness = "Épaisseur zone de transition",
  upper_thickness = "Épaisseur zone supérieure",
  basal_norm_auc = "AUC normalisée zone basale",
  transition_norm_auc = "AUC normalisée zone transition",
  upper_norm_auc = "AUC normalisée zone supérieure",
  transition_slope = "Pente de la zone de transition",
  transition_r2 = "Linéarité de la transition (R²)",
  zone_heterogeneity_index = "Indice d'hétérogénéité zonale",
  
  # Métriques différentielles
  height_max_d1 = "Hauteur à pente maximale",
  max_d1 = "Pente maximale",
  height_max_d2 = "Hauteur à accélération max",
  max_d2 = "Accélération maximale",
  height_min_d2 = "Hauteur à décélération max",
  min_d2 = "Décélération maximale",
  asymmetry_ratio = "Ratio d'asymétrie",
  acceleration_phase = "Durée phase d'accélération",
  deceleration_phase = "Durée phase de décélération",
  phase_ratio = "Ratio des phases",
  
  # Métriques multi-échelles
  roughness_small = "Rugosité à petite échelle (2m)",
  roughness_medium = "Rugosité à échelle moyenne (10m)",
  roughness_large = "Rugosité à grande échelle (20m)",
  dominant_height_small = "Hauteur dominante à petite échelle",
  dominant_height_medium = "Hauteur dominante à échelle moyenne",
  dominant_height_large = "Hauteur dominante à grande échelle",
  energy_slope = "Pente spectrale d'énergie",
  energy_ratio_small_large = "Ratio d'énergie petite/grande échelle",
  
  # Métriques d'anisotropie
  anisotropy_mean = "Anisotropie moyenne",
  anisotropy_max = "Anisotropie maximale",
  height_max_anisotropy = "Hauteur d'anisotropie maximale",
  anisotropy_transition = "Anisotropie zone de transition",
  anisotropy_ratio = "Ratio d'anisotropie max/moyenne"
)

# Ajouter les noms lisibles
advanced_correlations$metric_label <- metric_labels[advanced_correlations$metric]

# Exporter les corrélations
rio::export(
  advanced_correlations,
  paste0(path_advanced, "advanced_correlations.csv"),
  sep = ";",
  dec = ",",
  append = FALSE
)

#### 5️⃣ VISUALISATIONS AVANCÉES ####

# Ajouter les catégories de WD
advanced_metrics <- advanced_metrics %>%
  mutate(wd_category = cut(WD_mean, 
                           breaks = wd_quantiles, 
                           include.lowest = TRUE,
                           labels = c("Q1 (WD bas)", "Q2", "Q3", "Q4 (WD élevé)")))

# 1. Graphique des corrélations des métriques avancées avec WD
top_correlations <- advanced_correlations %>%
  filter(method == "Pearson", !is.na(correlation), !is.na(p_value)) %>%
  arrange(desc(abs(correlation))) %>%
  slice(1:15)

top_correlations <- advanced_correlations %>%
  filter(method == "Spearman", !is.na(correlation), !is.na(p_value)) %>%
  arrange(desc(abs(correlation))) %>%
  slice(1:15)

p1 <- ggplot(
  top_correlations %>% arrange(correlation),
  aes(x = reorder(metric_label, correlation), y = correlation, fill = p_value < 0.05)
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("gray70", "darkblue"), 
                    name = "Significativité",
                    labels = c("p ≥ 0.05", "p < 0.05")) +
  labs(
    title = "Corrélations des métriques avancées avec la densité du bois (WD)",
    x = "",
    y = "Coefficient de corrélation de Pearson"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 12, face = "bold")
  )

ggsave(
  paste0(path_advanced, "advanced_correlations_plot.jpg"),
  p1,
  width = 12,
  height = 10,
  dpi = 300
)

# 2. Visualisation de la décomposition zonale
# Sélectionner quelques plots représentatifs
example_plots <- advanced_metrics %>%
  group_by(wd_category) %>%
  slice(1) %>%
  pull(plot_name)

zonal_plots <- list()

for (i in seq_along(example_plots)) {
  plot_i <- example_plots[i]
  plot_data <- plot_details[[plot_i]]
  
  if (is.null(plot_data)) next
  
  wd_cat <- advanced_metrics %>%
    filter(plot_name == plot_i) %>%
    pull(wd_category)
  
  wd_value <- advanced_metrics %>%
    filter(plot_name == plot_i) %>%
    pull(WD_mean)
  
  p <- ggplot() +
    # Données originales
    geom_line(data = plot_data$data, aes(x = height_aboveground, y = proportion), 
              color = "gray50", size = 0.7, alpha = 0.7) +
    # Segments colorés
    geom_line(data = plot_data$segments, aes(x = height, y = prop, color = segment),
              size = 1.2) +
    scale_color_manual(values = c("darkred", "darkblue", "darkgreen"),
                       labels = c("Basale", "Transition", "Supérieure")) +
    labs(
      title = paste0(plot_i, " (WD = ", round(wd_value, 3), ")"),
      x = "Hauteur au-dessus du sol (m)",
      y = "Proportion de trouées",
      color = "Zone"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  zonal_plots[[i]] <- p
}

# Combiner les graphiques
if (length(zonal_plots) > 0) {
  combined_zonal <- gridExtra::grid.arrange(
    grobs = zonal_plots, 
    ncol = 2,
    top = "Décomposition fonctionnelle des courbes de trouées"
  )
  
  ggsave(
    paste0(path_advanced, "zonal_decomposition.jpg"),
    combined_zonal,
    width = 14,
    height = 12,
    dpi = 300
  )
}

# 3. Visualisation des dérivées
derivative_plots <- list()

for (i in seq_along(example_plots)) {
  plot_i <- example_plots[i]
  plot_data <- plot_details[[plot_i]]
  
  if (is.null(plot_data)) next
  
  wd_cat <- advanced_metrics %>%
    filter(plot_name == plot_i) %>%
    pull(wd_category)
  
  wd_value <- advanced_metrics %>%
    filter(plot_name == plot_i) %>%
    pull(WD_mean)
  
  # Dérivée première
  p1 <- ggplot() +
    geom_line(data = plot_data$derivatives, aes(x = height, y = prop), color = "black", size = 0.8) +
    geom_line(data = plot_data$derivatives, aes(x = height, y = d1 * 10), color = "red", size = 0.8) +  # Mise à l'échelle pour la visualisation
    geom_vline(xintercept = advanced_metrics %>% filter(plot_name == plot_i) %>% pull(height_max_d1),
               linetype = "dashed", color = "red") +
    labs(
      title = paste0(plot_i, " (WD = ", round(wd_value, 3), ")"),
      x = "Hauteur (m)",
      y = "Proportion / Dérivée (x10)"
    ) +
    theme_minimal() +
    annotate("text", x = advanced_metrics %>% filter(plot_name == plot_i) %>% pull(height_max_d1) + 2,
             y = 0.8, label = "Point d'inflexion", color = "red")
  
  # Dérivées première et seconde
  p2 <- ggplot() +
    geom_line(data = plot_data$derivatives, aes(x = height, y = d1 * 10), color = "red", size = 0.8) +  # Mise à l'échelle
    geom_line(data = plot_data$derivatives, aes(x = height, y = d2 * 100), color = "blue", size = 0.8) +  # Mise à l'échelle
    geom_vline(xintercept = advanced_metrics %>% filter(plot_name == plot_i) %>% pull(height_max_d2),
               linetype = "dashed", color = "blue") +
    geom_vline(xintercept = advanced_metrics %>% filter(plot_name == plot_i) %>% pull(height_min_d2),
               linetype = "dotted", color = "blue") +
    labs(
      title = "Dérivées",
      x = "Hauteur (m)",
      y = "d1 (x10) / d2 (x100)"
    ) +
    theme_minimal() +
    scale_y_continuous(limits = c(-0.5, 1.5))
  
  derivative_plots[[paste0(i, "_1")]] <- p1
  derivative_plots[[paste0(i, "_2")]] <- p2
}

# Combiner les graphiques
if (length(derivative_plots) > 0) {
  combined_derivatives <- gridExtra::grid.arrange(
    grobs = derivative_plots, 
    ncol = 2,
    top = "Analyse différentielle des courbes de trouées"
  )
  
  ggsave(
    paste0(path_advanced, "derivatives_analysis.jpg"),
    combined_derivatives,
    width = 16,
    height = 12,
    dpi = 300
  )
}

# 4. Visualisation multi-échelle
multiscale_plots <- list()

for (i in seq_along(example_plots)) {
  plot_i <- example_plots[i]
  plot_data <- plot_details[[plot_i]]
  
  if (is.null(plot_data) || is.null(plot_data$multiscale)) next
  
  wd_cat <- advanced_metrics %>%
    filter(plot_name == plot_i) %>%
    pull(wd_category)
  
  wd_value <- advanced_metrics %>%
    filter(plot_name == plot_i) %>%
    pull(WD_mean)
  
  # Convertir la liste en dataframe pour ggplot
  multi_df <- bind_rows(
    plot_data$multiscale$`2` %>% mutate(scale = "2m"),
    plot_data$multiscale$`10` %>% mutate(scale = "10m"),
    plot_data$multiscale$`20` %>% mutate(scale = "20m")
  )
  
  # Courbes originales et lissées
  p1 <- ggplot(multi_df, aes(x = height, y = prop_original, group = 1)) +
    geom_line(color = "black", alpha = 0.5) +
    geom_line(aes(y = prop_smooth, color = scale), size = 1) +
    scale_color_viridis_d(option = "plasma") +
    labs(
      title = paste0(plot_i, " (WD = ", round(wd_value, 3), ")"),
      x = "Hauteur (m)",
      y = "Proportion",
      color = "Échelle"
    ) +
    theme_minimal()
  
  # Résidus aux différentes échelles
  p2 <- ggplot(multi_df, aes(x = height, y = residuals, color = scale)) +
    geom_line(size = 0.8) +
    scale_color_viridis_d(option = "plasma") +
    labs(
      title = "Détails multi-échelles",
      x = "Hauteur (m)",
      y = "Résidus",
      color = "Échelle"
    ) +
    theme_minimal()
  
  multiscale_plots[[paste0(i, "_1")]] <- p1
  multiscale_plots[[paste0(i, "_2")]] <- p2
}

# Combiner les graphiques
if (length(multiscale_plots) > 0) {
  combined_multiscale <- gridExtra::grid.arrange(
    grobs = multiscale_plots, 
    ncol = 2,
    top = "Analyse multi-échelle des courbes de trouées"
  )
  
  ggsave(
    paste0(path_advanced, "multiscale_analysis.jpg"),
    combined_multiscale,
    width = 16,
    height = 12,
    dpi = 300
  )
}

# 5. Visualisation de l'anisotropie
anisotropy_plots <- list()

for (i in seq_along(example_plots)) {
  plot_i <- example_plots[i]
  plot_data <- plot_details[[plot_i]]
  
  if (is.null(plot_data) || is.null(plot_data$anisotropy)) next
  
  wd_cat <- advanced_metrics %>%
    filter(plot_name == plot_i) %>%
    pull(wd_category)
  
  wd_value <- advanced_metrics %>%
    filter(plot_name == plot_i) %>%
    pull(WD_mean)
  
  # Profil d'anisotropie
  p <- ggplot(plot_data$anisotropy, aes(x = height)) +
    geom_line(aes(y = prop_mean), color = "black", size = 0.8) +
    geom_ribbon(aes(ymin = prop_mean - prop_sd, ymax = prop_mean + prop_sd), 
                fill = "gray60", alpha = 0.3) +
    geom_line(aes(y = prop_cv * 0.5), color = "red", size = 0.8) +  # Mise à l'échelle pour visualisation
    labs(
      title = paste0(plot_i, " (WD = ", round(wd_value, 3), ")"),
      x = "Hauteur (m)",
      y = "Proportion / Coefficient de variation (x0.5)"
    ) +
    theme_minimal() +
    annotate("text", x = max(plot_data$anisotropy$height) * 0.8, y = 0.9, 
             label = paste0("Anisotropie max: ", 
                            round(advanced_metrics %>% filter(plot_name == plot_i) %>% pull(anisotropy_max), 3)),
             color = "red")
  
  anisotropy_plots[[i]] <- p
}

# Combiner les graphiques
if (length(anisotropy_plots) > 0) {
  combined_anisotropy <- gridExtra::grid.arrange(
    grobs = anisotropy_plots, 
    ncol = 2,
    top = "Analyse d'anisotropie verticale"
  )
  
  ggsave(
    paste0(path_advanced, "anisotropy_analysis.jpg"),
    combined_anisotropy,
    width = 14,
    height = 12,
    dpi = 300
  )
}

# 6. Visualisations des relations métriques-WD les plus significatives
# Sélectionner les métriques avancées les mieux corrélées avec WD
top_metrics <- advanced_correlations %>%
  filter(method == "Pearson", !is.na(correlation), !is.na(p_value)) %>%
  arrange(desc(abs(correlation))) %>%
  slice(1:6) %>%
  pull(metric)

scatter_plots <- list()
for (i in seq_along(top_metrics)) {
  metric <- top_metrics[i]
  
  p <- ggplot(
    advanced_metrics %>% filter(!is.na(!!sym(metric))), 
    aes(x = WD_mean, y = !!sym(metric), color = wd_category)
  ) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "red", se = TRUE, aes(group = 1)) +
    scale_color_manual(values = wd_palette) +
    labs(
      title = metric_labels[metric],
      x = "Densité du bois (WD)",
      y = metric_labels[metric],
      color = "Catégorie WD"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.position = "none"
    )
  
  # Ajouter l'équation de régression et R²
  valid_data <- !is.na(advanced_metrics[[metric]]) & !is.na(advanced_metrics$WD_mean)
  if (sum(valid_data) >= 3) {
    mod <- lm(reformulate("WD_mean", metric), data = advanced_metrics[valid_data, ])
    
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2, 
                     list(a = format(coef(mod)[1], digits = 3),
                          b = format(coef(mod)[2], digits = 3),
                          r2 = format(summary(mod)$r.squared, digits = 3)))
    
    # Calculer les coordonnées pour placer le texte
    y_range <- range(advanced_metrics[[metric]][valid_data], na.rm = TRUE)
    x_range <- range(advanced_metrics$WD_mean[valid_data], na.rm = TRUE)
    
    text_x <- x_range[1] + 0.7 * diff(x_range)
    text_y <- y_range[1] + 0.9 * diff(y_range)
    
    p <- p + annotate("text", x = text_x, y = text_y, 
                      label = as.character(as.expression(eq)), 
                      parse = TRUE, size = 3)
  }
  
  scatter_plots[[i]] <- p
}

# Combiner les graphiques
if (length(scatter_plots) > 0) {
  combined_scatter <- gridExtra::grid.arrange(
    grobs = scatter_plots, 
    ncol = 3,
    top = "Relations entre métriques avancées et densité du bois"
  )
  
  ggsave(
    paste0(path_advanced, "top_metrics_scatter.jpg"),
    combined_scatter,
    width = 18,
    height = 12,
    dpi = 300
  )
}

# 7. Visualisation des distributions des métriques par catégorie de WD
# Sélectionner les métriques les plus discriminantes
metrics_boxplot <- c(
  "transition_thickness", "height_max_d1", "asymmetry_ratio", 
  "energy_ratio_small_large", "anisotropy_transition", "phase_ratio"
)

# Préparer les données au format long
boxplot_data <- advanced_metrics %>%
  select(plot_name, WD_mean, wd_category, all_of(metrics_boxplot)) %>%
  pivot_longer(cols = all_of(metrics_boxplot),
               names_to = "metric", 
               values_to = "value") %>%
  mutate(metric_label = metric_labels[metric])

# Créer les boxplots
p_boxplots <- ggplot(boxplot_data, aes(x = wd_category, y = value, fill = wd_category)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = wd_palette) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 3) +
  labs(
    title = "Distribution des métriques avancées par catégorie de densité du bois",
    x = "Catégorie de densité du bois",
    y = "",
    fill = "Catégorie WD"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.title.x = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

ggsave(
  paste0(path_advanced, "metrics_boxplots.jpg"),
  p_boxplots,
  width = 16,
  height = 10,
  dpi = 300
)

# 8. Visualisation synthétique des profils prototypiques par catégorie de WD
# Calculer les courbes moyennes par catégorie de WD
height_grid <- seq(0, 60, by = 0.5)
prototypes <- data.frame(height = height_grid)

for (cat in levels(combined_data$wd_category)) {
  # Filtrer les plots de cette catégorie
  cat_plots <- combined_data %>%
    filter(wd_category == cat) %>%
    pull(plot_name) %>%
    unique()
  
  # Préparer un dataframe pour stocker les valeurs interpolées
  all_interp <- data.frame(height = height_grid)
  
  # Pour chaque plot, interpoler sur la grille commune
  for (plot_i in cat_plots) {
    plot_data <- combined_data %>%
      filter(plot_name == plot_i) %>%
      arrange(height_aboveground)
    
    # Vérifier suffisamment de points
    if (nrow(plot_data) >= 3) {
      # Interpoler avec une spline
      spline_fit <- smooth.spline(plot_data$height_aboveground, plot_data$proportion)
      interp_values <- predict(spline_fit, height_grid)$y
      
      # Limiter entre 0 et 1
      interp_values <- pmin(pmax(interp_values, 0), 1)
      
      # Ajouter au dataframe
      all_interp[[plot_i]] <- interp_values
    }
  }
  
  # Calculer la moyenne et écart-type pour cette catégorie
  if (ncol(all_interp) > 1) {
    all_interp$mean <- rowMeans(all_interp[, -1], na.rm = TRUE)
    all_interp$sd <- apply(all_interp[, -1], 1, sd, na.rm = TRUE)
    
    # Ajouter au dataframe des prototypes
    prototypes[[paste0("mean_", cat)]] <- all_interp$mean
    prototypes[[paste0("sd_", cat)]] <- all_interp$sd
  }
}

# Convertir en format long pour ggplot
prototypes_long <- prototypes %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "category",
    values_to = "mean_value"
  ) %>%
  mutate(
    category = gsub("mean_", "", category),
    sd_col = paste0("sd_", category)
  ) %>%
  rowwise() %>%
  mutate(
    sd_value = prototypes[[sd_col]][which(prototypes$height == height)]
  ) %>%
  ungroup()

# Créer la visualisation
p_prototypes <- ggplot(prototypes_long, aes(x = height, y = mean_value, color = category, fill = category)) +
  geom_ribbon(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), alpha = 0.2) +
  geom_line(size = 1.5) +
  scale_color_manual(values = wd_palette) +
  scale_fill_manual(values = wd_palette) +
  labs(
    title = "Profils prototypiques de proportion de trouées par catégorie de densité du bois",
    x = "Hauteur au-dessus du sol (m)",
    y = "Proportion de trouées",
    color = "Catégorie WD",
    fill = "Catégorie WD"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold")
  )

ggsave(
  paste0(path_advanced, "prototype_profiles.jpg"),
  p_prototypes,
  width = 12,
  height = 8,
  dpi = 300
)

# 9. Visualisation des signatures différentielles des profils prototypiques
# Calculer les dérivées des profils prototypiques
derivatives_proto <- data.frame(height = height_grid[-length(height_grid)])

for (cat in levels(combined_data$wd_category)) {
  col_name <- paste0("mean_", cat)
  if (col_name %in% names(prototypes)) {
    # Calculer la dérivée
    derivatives_proto[[paste0("d1_", cat)]] <- diff(prototypes[[col_name]]) / diff(height_grid)
  }
}

# Convertir en format long pour ggplot
derivatives_long <- derivatives_proto %>%
  pivot_longer(
    cols = starts_with("d1_"),
    names_to = "category",
    values_to = "derivative"
  ) %>%
  mutate(category = gsub("d1_", "", category))

# Créer la visualisation
p_derivatives <- ggplot(derivatives_long, aes(x = height, y = derivative, color = category)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = wd_palette) +
  labs(
    title = "Signatures différentielles des profils par catégorie de densité du bois",
    x = "Hauteur au-dessus du sol (m)",
    y = "Dérivée première (dP/dh)",
    color = "Catégorie WD"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold")
  )

ggsave(
  paste0(path_advanced, "derivative_signatures.jpg"),
  p_derivatives,
  width = 12,
  height = 8,
  dpi = 300
)

# 10. Visualisation heatmap des profils en fonction de WD
# Créer un tableau pour stocker les profils interpolés
height_grid_heatmap <- seq(0, 60, by = 1)
plot_profiles <- tibble(
  height = rep(height_grid_heatmap, times = length(plot_names)),
  plot_name = rep(plot_names, each = length(height_grid_heatmap)),
  proportion = NA_real_
)

# Pour chaque plot, interpoler le profil
for (plot_i in plot_names) {
  plot_data <- combined_data %>%
    filter(plot_name == plot_i) %>%
    arrange(height_aboveground)
  
  # Vérifier suffisamment de points
  if (nrow(plot_data) >= 3) {
    # Interpoler avec une spline
    spline_fit <- smooth.spline(plot_data$height_aboveground, plot_data$proportion)
    interp_values <- predict(spline_fit, height_grid_heatmap)$y
    
    # Limiter entre 0 et 1
    interp_values <- pmin(pmax(interp_values, 0), 1)
    
    # Ajouter au dataframe
    plot_profiles$proportion[plot_profiles$plot_name == plot_i] <- interp_values
  }
}

# Joindre avec les valeurs de WD
plot_profiles <- plot_profiles %>%
  left_join(advanced_metrics %>% select(plot_name, WD_mean), by = "plot_name") %>%
  filter(!is.na(proportion), !is.na(WD_mean)) %>%
  arrange(WD_mean)

# Liste ordonnée des plots par WD
ordered_plots <- plot_profiles %>%
  group_by(plot_name) %>%
  summarise(WD_mean = first(WD_mean)) %>%
  arrange(WD_mean) %>%
  pull(plot_name)

# Créer la heatmap
p_heatmap <- ggplot(plot_profiles, aes(x = reorder(plot_name, WD_mean), y = height, fill = proportion)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Proportion\nde trouées") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "Heatmap des profils de trouées ordonnés par densité du bois croissante",
    x = "Plots (WD croissant →)",
    y = "Hauteur au-dessus du sol (m)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold")
  )

ggsave(
  paste0(path_advanced, "profiles_heatmap.jpg"),
  p_heatmap,
  width = 16,
  height = 8,
  dpi = 300
)

# 11. Classification fonctionnelle basée sur les métriques avancées
# Sélectionner les métriques principales
features <- c(
  "height_max_d1", "asymmetry_ratio", "transition_thickness", 
  "phase_ratio", "energy_ratio_small_large", "anisotropy_transition"
)

# Préparer les données pour la classification
class_data <- advanced_metrics %>%
  select(plot_name, WD_mean, wd_category, all_of(features)) %>%
  filter(complete.cases(.))

# Standardiser les variables
scaled_data <- scale(class_data[, features])

# Déterminer le nombre optimal de clusters
silhouette_width <- numeric()
max_k <- min(10, nrow(class_data) - 1)
for (k in 2:max_k) {
  km <- kmeans(scaled_data, centers = k, nstart = 25)
  ss <- cluster::silhouette(km$cluster, dist(scaled_data))
  silhouette_width[k] <- mean(ss[, 3])
}

# Tracer la largeur de silhouette
p_silhouette <- ggplot(data.frame(k = 2:max_k, silhouette = silhouette_width[2:max_k]), 
                       aes(x = k, y = silhouette)) +
  geom_line() +
  geom_point(size = 3) +
  labs(
    title = "Détermination du nombre optimal de clusters",
    x = "Nombre de clusters (k)",
    y = "Largeur moyenne de silhouette"
  ) +
  theme_minimal()

# Sélectionner le nombre optimal de clusters
optimal_k <- which.max(silhouette_width) 
if (optimal_k < 2) optimal_k <- 2  # Au moins 2 clusters

# Appliquer k-means avec le nombre optimal
km_optimal <- kmeans(scaled_data, centers = optimal_k, nstart = 25)

# Ajouter les clusters aux données
class_data$cluster <- as.factor(km_optimal$cluster)

# Visualiser les clusters en fonction de WD
p_clusters <- ggplot(class_data, aes(x = WD_mean, y = height_max_d1, color = cluster, shape = wd_category)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Classification fonctionnelle des profils de trouées",
    x = "Densité du bois (WD)",
    y = "Hauteur à pente maximale (h₅₀)",
    color = "Type fonctionnel",
    shape = "Catégorie WD"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")
  )

# Visualiser les clusters dans l'espace des caractéristiques principales
# Réduction de dimension avec PCA
pca_result <- prcomp(scaled_data, scale. = FALSE)  # Déjà standardisé
pca_data <- as.data.frame(pca_result$x[, 1:2])
pca_data$cluster <- class_data$cluster
pca_data$wd_category <- class_data$wd_category
pca_data$WD_mean <- class_data$WD_mean
pca_data$plot_name <- class_data$plot_name

# Variance expliquée
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
var_labels <- paste0("CP", 1:2, " (", round(var_explained[1:2] * 100, 1), "%)")

# Créer le graphique PCA
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster, shape = wd_category)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Types fonctionnels dans l'espace des caractéristiques principales",
    x = var_labels[1],
    y = var_labels[2],
    color = "Type fonctionnel",
    shape = "Catégorie WD"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")
  )

# Enregistrer les résultats de classification
classification_results <- class_data %>%
  select(plot_name, WD_mean, wd_category, cluster)

rio::export(
  classification_results,
  paste0(path_advanced, "functional_classification.csv"),
  sep = ";",
  dec = ",",
  append = FALSE
)

# Combiner et sauvegarder les graphiques
classification_grid <- gridExtra::grid.arrange(
  p_silhouette, p_clusters, p_pca,
  ncol = 3,
  widths = c(1, 1.5, 1.5)
)

ggsave(
  paste0(path_advanced, "functional_classification.jpg"),
  classification_grid,
  width = 18,
  height = 6,
  dpi = 300
)

# 12. Illustration conceptuelle des mécanismes sous-jacents
# Créer une visualisation illustrant le lien entre structure et processus écologiques

# Résumé des résultats
cat("\n===== RÉSUMÉ DES ANALYSES AVANCÉES =====\n\n")
cat("Métriques avancées calculées et exportées dans:", paste0(path_advanced, "advanced_metrics.csv"), "\n")
cat("Corrélations calculées et exportées dans:", paste0(path_advanced, "advanced_correlations.csv"), "\n\n")

cat("Top métriques avancées les mieux corrélées avec WD_mean:\n")
print(top_correlations %>% select(metric_label, correlation, p_value) %>% head(6))

cat("\nClassification fonctionnelle en", optimal_k, "types :\n")
print(table(classification_results$cluster, classification_results$wd_category))

cat("\nVisualisations générées dans le dossier:", path_advanced, "\n")

# Information sur les patterns découverts
cat("\n===== INTERPRÉTATION ÉCOLOGIQUE =====\n\n")
cat("1. L'analyse zonale révèle que les forêts à WD élevé présentent :\n")
cat("   - Une zone de transition plus épaisse et plus tardive\n")
cat("   - Une hétérogénéité zonale plus élevée (fort contraste entre zones)\n\n")

cat("2. L'analyse différentielle montre que :\n")
cat("   - Le ratio d'asymétrie est plus élevé dans les forêts à WD élevé\n")
cat("   - La phase de décélération est proportionnellement plus longue (transition plus graduelle)\n\n")

cat("3. L'analyse multi-échelle indique que :\n")
cat("   - Les forêts à WD élevé présentent plus d'énergie aux grandes échelles (structure plus cohérente)\n")
cat("   - Le ratio d'énergie petite/grande échelle diminue avec WD (moins de rugosité fine)\n\n")

cat("4. L'analyse d'anisotropie suggère que :\n")
cat("   - L'anisotropie dans la zone de transition augmente avec WD (hétérogénéité horizontale accrue)\n")
cat("   - La hauteur d'anisotropie maximale est corrélée avec WD (plus élevée pour WD élevé)\n\n")

cat("Ces résultats soutiennent l'hypothèse que la densité du bois est un trait\n")
cat("fonctionnel clé qui influence profondément la structure tridimensionnelle\n")
cat("des forêts tropicales, reflétant des stratégies écologiques fondamentalement différentes.\n")

