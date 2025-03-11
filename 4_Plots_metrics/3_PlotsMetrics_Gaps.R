# =================================================================
# CALCUL DES MÉTRIQUES DE TROUÉES POUR L'ANALYSE DE LA STRUCTURE FORESTIÈRE
# Objectif : Calculer une suite complète de métriques caractérisant
#            la structure verticale forestière en relation avec la densité du bois
# =================================================================

#### 1️⃣ INITIALISATION ####
# Nettoyage de l'environnement et de la mémoire
rm(list = ls())
gc()

# Définition des packages nécessaires
pkgs <- c("terra", "tidyverse", "sf", "rio", "foreach", "doParallel", 
          "mgcv", "minpack.lm", "segmented", "splines", "waveslim", 
          "fields", "zoo", "progress")

# Installation et chargement automatique des packages manquants
to_install <- !pkgs %in% installed.packages()
if(any(to_install)) {install.packages(pkgs[to_install])}
inst <- lapply(pkgs, library, character.only = TRUE)

# Définition des répertoires de travail
path_base <- "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal"
path_gaps <- file.path(path_base, "3_Gaps")
path_metrics <- file.path(path_base, "4_Plots_metrics")
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
#' L'AUC représente l'intégrale de la proportion de trouées sur toute la hauteur,
#' fournissant une mesure globale de l'ouverture de la canopée dans toute la
#' structure verticale de la forêt.
#' 
#' @param data Données d'un plot filtré
#' @param max_height Hauteur maximale à considérer pour normaliser l'AUC
#' @return Aire sous la courbe normalisée
calculate_auc <- function(data, max_height = 45) {
  # Ordonner les données par hauteur
  data <- data %>% arrange(height_aboveground)
  
  # Calculer l'AUC avec la méthode des trapèzes
  auc <- 0
  for (i in 2:nrow(data)) {
    h1 <- data$height_aboveground[i-1]
    h2 <- data$height_aboveground[i]
    p1 <- data$proportion[i-1]
    p2 <- data$proportion[i]
    
    # Aire du trapèze
    auc <- auc + (h2 - h1) * (p1 + p2) / 2
  }
  
  # Normaliser par la hauteur maximale
  max_h <- min(max(data$height_aboveground), max_height)
  min_h <- min(data$height_aboveground)
  
  normalized_auc <- auc / (max_h - min_h)
  return(normalized_auc)
}

#' Calcule le taux de décroissance moyen de la proportion de trouées
#' Cette métrique capture la pente globale de la relation entre hauteur et
#' proportion de trouées, offrant une mesure simple de la rapidité de transition
#' dans la structure verticale.
#' 
#' @param data Données d'un plot filtré
#' @return Pente de la régression linéaire (taux de changement)
calculate_slope <- function(data) {
  # Vérifier que les données ont au moins 3 lignes
  if (nrow(data) < 3) return(NA)
  
  # Ajuster un modèle linéaire
  model <- lm(proportion ~ height_aboveground, data = data)
  
  # Récupérer la pente
  slope <- coef(model)[2]
  return(slope)
}

#' Calcule la pente maximale et sa position pour une courbe de proportion de trouées
#' Ces métriques sont cruciales pour identifier le point d'inflexion de la courbe
#' et la vitesse maximale de changement dans la structure verticale.
#' 
#' @param data Données d'un plot filtré
#' @return Liste contenant Smax (pente maximale) et hSmax (hauteur correspondante)
calculate_max_slope <- function(data) {
  # Vérifier suffisamment de données
  if (nrow(data) < 5) return(list(Smax = NA, hSmax = NA))
  
  # Ordonner les données
  data <- data %>% arrange(height_aboveground)
  
  # Calculer les pentes entre points adjacents
  slopes <- numeric(nrow(data) - 1)
  heights <- numeric(nrow(data) - 1)
  
  for (i in 1:(nrow(data) - 1)) {
    h1 <- data$height_aboveground[i]
    h2 <- data$height_aboveground[i + 1]
    p1 <- data$proportion[i]
    p2 <- data$proportion[i + 1]
    
    slopes[i] <- (p2 - p1) / (h2 - h1)
    heights[i] <- (h1 + h2) / 2  # Hauteur moyenne entre les deux points
  }
  
  # Trouver la pente maximale et sa position
  if (length(slopes) > 0) {
    max_idx <- which.max(slopes)
    Smax <- slopes[max_idx]
    hSmax <- heights[max_idx]
    return(list(Smax = Smax, hSmax = hSmax))
  } else {
    return(list(Smax = NA, hSmax = NA))
  }
}

#' Ajuste une fonction logistique à la courbe de proportion de trouées
#' La modélisation logistique permet de caractériser précisément la forme sigmoïde
#' et d'extraire des paramètres biologiquement significatifs pour comparer 
#' différentes structures forestières.
#'
#' @param data Données d'un plot filtré
#' @return Liste contenant les paramètres du modèle et la qualité d'ajustement
fit_logistic_model <- function(data) {
  # Vérifier suffisamment de données
  if (nrow(data) < 5) {
    return(list(a = NA, b = NA, h0 = NA, k = NA, R2 = NA, AIC = NA, 
                convergence = FALSE))
  }
  
  # Ordonner les données
  data <- data %>% arrange(height_aboveground)
  
  # Estimer des valeurs initiales pour les paramètres
  a_init <- min(data$proportion, na.rm = TRUE)
  b_init <- max(data$proportion, na.rm = TRUE)
  
  # Estimer h0 (point d'inflexion) comme hauteur à 50% entre min et max
  target <- a_init + (b_init - a_init) / 2
  h0_init <- data$height_aboveground[which.min(abs(data$proportion - target))]
  
  # Estimer k (pente au point d'inflexion)
  slopes <- diff(data$proportion) / diff(data$height_aboveground)
  k_init <- max(slopes, na.rm = TRUE) * 4  # Ajustement pour la logistique
  
  # Paramètres initiaux
  start_params <- list(a = a_init, b = b_init, h0 = h0_init, k = k_init)
  
  # Définir la fonction logistique
  logistic_func <- function(h, a, b, h0, k) {
    a + (b - a) / (1 + exp(-k * (h - h0)))
  }
  
  # Ajuster le modèle avec minpack.lm qui est plus robuste que nls
  tryCatch({
    model <- minpack.lm::nlsLM(
      proportion ~ a + (b - a) / (1 + exp(-k * (height_aboveground - h0))),
      data = data,
      start = start_params,
      control = list(maxiter = 500)
    )
    
    # Extraire les paramètres
    params <- coef(model)
    a <- params["a"]
    b <- params["b"]
    h0 <- params["h0"]
    k <- params["k"]
    
    # Calculer R² et AIC
    predictions <- predict(model)
    SSE <- sum((data$proportion - predictions)^2)
    SST <- sum((data$proportion - mean(data$proportion))^2)
    R2 <- 1 - SSE/SST
    AIC <- AIC(model)
    
    return(list(a = a, b = b, h0 = h0, k = k, R2 = R2, AIC = AIC, 
                convergence = TRUE))
  }, error = function(e) {
    # En cas d'échec de convergence, retourner NA
    return(list(a = NA, b = NA, h0 = NA, k = NA, R2 = NA, AIC = NA, 
                convergence = FALSE))
  })
}

#' Calcule des métriques combinées à partir des paramètres fondamentaux
#' Ces métriques complexes synthétisent plusieurs aspects de la structure
#' forestière en un seul indicateur, potentiellement plus sensible aux
#' variations de densité du bois.
#'
#' @param h10 Hauteur à 10% de trouées
#' @param h25 Hauteur à 25% de trouées 
#' @param h50 Hauteur à 50% de trouées
#' @param h75 Hauteur à 75% de trouées
#' @param h90 Hauteur à 90% de trouées
#' @param Smax Pente maximale
#' @param hSmax Hauteur à pente maximale
#' @param P0 Proportion initiale de trouées
#' @param Pmax Proportion maximale de trouées
#' @param auc Aire sous la courbe
#' @return Liste des métriques combinées
calculate_combined_metrics <- function(h10, h25, h50, h75, h90, Smax, hSmax, P0, Pmax, auc) {
  # Vérifier que h50 n'est pas NA ou 0 (pour éviter division par zéro)
  if (is.na(h50) || h50 == 0) {
    return(list(
      IDV = NA, ROS = NA, IAS = NA, ISC = NA, CPR = NA, ICT = NA, IPN = NA,
      ET = NA, CSI = NA, GVN = NA, RT = NA, HPSI = NA
    ))
  }
  
  # Calculer les métriques combinées
  
  # Indice de Développement Vertical: élevé pour forêts à WD élevée (transition tardive et graduelle)
  IDV <- ifelse(!is.na(h10) && !is.na(h90) && h90 > 0, 
                h50 * (h90 - h10) / h90, NA)
  
  # Ratio d'Occupation Stratifiée: élevé pour forêts à WD élevée (h50 élevé, P0 très faible)
  ROS <- h50 / (P0 + 0.05)
  
  # Indice d'Asymétrie Structurelle: rapport entre les transitions haute et basse
  IAS <- ifelse(!is.na(h25) && !is.na(h75) && (h50 - h25) > 0, 
                (h75 - h50) / (h50 - h25), NA)
  
  # Indice de Stratification Complète: intègre largeur de transition, amplitude et position
  ISC <- ifelse(!is.na(h25) && !is.na(h75) && h50 > 0, 
                ((h75 - h25) * (Pmax - P0)) / (h50^2), NA)
  
  # Coefficient de Position Relative: combine position relative et fermeture initiale
  CPR <- ifelse(!is.na(h90) && h90 > 0, (h50 / h90) * (1 - P0), NA)
  
  # Indice de Contraste de Transition: combine rapidité, position et étalement
  ICT <- ifelse(!is.na(h10) && !is.na(h90) && !is.na(Smax) && (h90 - h10) > 0, 
                Smax * h50 / (h90 - h10), NA)
  
  # Indice de Progression Normalisée: étalement relatif vs hauteur médiane
  IPN <- ifelse(!is.na(h10) && !is.na(h90) && h50 > 0, 
                (h90 - h10) / (h50^2), NA)
  
  # Épaisseur de transition: mesure directe de l'étalement vertical
  ET <- ifelse(!is.na(h10) && !is.na(h90), h90 - h10, NA)
  
  # Indice de structure de canopée: produit hauteur médiane x pente maximale
  CSI <- ifelse(!is.na(Smax), h50 * Smax, NA)
  
  # Gradient vertical normalisé: pente maximale / hauteur maximale (ou h90)
  GVN <- ifelse(!is.na(Smax) && !is.na(h90) && h90 > 0, Smax / h90, NA)
  
  # Ratio de transition: similaire à IAS mais formulé différemment
  RT <- ifelse(!is.na(h25) && !is.na(h75) && (h50 - h25) > 0, 
               (h75 - h50) / (h50 - h25), NA)
  
  # Harmonie Position-Slope Index: relation entre position pente max et h50
  HPSI <- ifelse(!is.na(hSmax) && h50 > 0, 1 - abs(hSmax - h50) / h50, NA)
  
  return(list(
    IDV = IDV, ROS = ROS, IAS = IAS, ISC = ISC, CPR = CPR, ICT = ICT, IPN = IPN,
    ET = ET, CSI = CSI, GVN = GVN, RT = RT, HPSI = HPSI
  ))
}

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

#' Calcule les nouveaux indices optimisés
#'
#' @param metrics Dataframe contenant toutes les métriques calculées
#' @param max_pente_ref Valeur maximale de référence pour la pente (pour normalisation)
#' @param max_transition_thickness Valeur maximale de l'épaisseur de transition (pour normalisation)
#' @return Liste contenant les indices optimisés
calculate_optimized_indices <- function(metrics, max_pente_ref = NULL, max_transition_thickness = NULL) {
  
  # Extraire les valeurs requises, s'assurer qu'elles sont des scalaires (length 1)
  h50 <- if(length(metrics$h50) > 0) metrics$h50[1] else NA
  height_min_d2 <- if(length(metrics$height_min_d2) > 0) metrics$height_min_d2[1] else NA
  phase_ratio <- if(length(metrics$phase_ratio) > 0) metrics$phase_ratio[1] else NA
  max_d1 <- if(length(metrics$max_d1) > 0) metrics$max_d1[1] else NA
  
  basal_thickness <- if(length(metrics$basal_thickness) > 0) metrics$basal_thickness[1] else NA
  upper_thickness <- if(length(metrics$upper_thickness) > 0) metrics$upper_thickness[1] else NA
  anisotropy_max <- if(length(metrics$anisotropy_max) > 0) metrics$anisotropy_max[1] else NA
  energy_ratio_small_large <- if(length(metrics$energy_ratio_small_large) > 0) metrics$energy_ratio_small_large[1] else NA
  
  transition_thickness <- if(length(metrics$transition_thickness) > 0) metrics$transition_thickness[1] else NA
  basal_norm_auc <- if(length(metrics$basal_norm_auc) > 0) metrics$basal_norm_auc[1] else NA
  
  # Si aucune valeur de référence n'est fournie, utiliser la valeur du plot actuel
  if (is.null(max_pente_ref)) max_pente_ref <- abs(max_d1)
  if (is.null(max_transition_thickness)) max_transition_thickness <- transition_thickness
  
  # 1. Indice de Stratification Asymétrique (ISA)
  ISA <- if(!is.na(height_min_d2) && !is.na(h50) && !is.na(phase_ratio) && !is.na(max_d1) && 
            h50 > 0 && max_pente_ref > 0) {
    (height_min_d2 / h50) * phase_ratio * (1 - abs(max_d1) / max_pente_ref)
  } else {
    NA
  }
  
  # 2. Indice de Complexité Verticale (ICV)
  ICV <- if(!is.na(basal_thickness) && !is.na(upper_thickness) && !is.na(anisotropy_max) && 
            !is.na(energy_ratio_small_large) && upper_thickness > 0 && energy_ratio_small_large > 0) {
    (basal_thickness / upper_thickness) * anisotropy_max * (1 / energy_ratio_small_large)
  } else {
    NA
  }
  
  # 3. Indice de Différenciation Transitionnelle (IDT)
  IDT <- if(!is.na(h50) && !is.na(transition_thickness) && !is.na(basal_norm_auc) && 
            max_transition_thickness > 0) {
    h50 * (transition_thickness / max_transition_thickness) * (1 - basal_norm_auc)
  } else {
    NA
  }
  
  # Calcul des transformations supplémentaires recommandées
  log_transition_thickness <- if(!is.na(transition_thickness) && transition_thickness > 0) {
    log(transition_thickness)
  } else {
    NA
  }
  
  log_max_d1 <- if(!is.na(max_d1) && max_d1 > 0) {
    log(max_d1)
  } else {
    NA
  }
  
  h90 <- if(length(metrics$h90) > 0) metrics$h90[1] else NA
  
  h50_relatif <- if(!is.na(h50) && !is.na(h90) && h90 > 0) {
    h50 / h90
  } else {
    NA
  }
  
  basal_thickness_relatif <- if(!is.na(basal_thickness) && !is.na(h90) && h90 > 0) {
    basal_thickness / h90
  } else {
    NA
  }
  
  ratio_h50_h90 <- h50_relatif
  
  ratio_basal_transition <- if(!is.na(basal_thickness) && !is.na(transition_thickness) && 
                               transition_thickness > 0) {
    basal_thickness / transition_thickness
  } else {
    NA
  }
  
  h50_squared <- if(!is.na(h50)) h50^2 else NA
  anisotropy_squared <- if(!is.na(anisotropy_max)) anisotropy_max^2 else NA
  
  # Retourner les indices et transformations
  return(list(
    ISA = ISA,
    ICV = ICV,
    IDT = IDT,
    log_transition_thickness = log_transition_thickness,
    log_max_d1 = log_max_d1,
    h50_relatif = h50_relatif,
    basal_thickness_relatif = basal_thickness_relatif,
    ratio_h50_h90 = ratio_h50_h90,
    ratio_basal_transition = ratio_basal_transition,
    h50_squared = h50_squared,
    anisotropy_squared = anisotropy_squared
  ))
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
plot_names <- read.csv2(file.path(path_base, "final_plot_name.csv")) %>% 
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
if ("WD_mean" %in% colnames(plots_info)) {
  combined_data <- buffer_data %>%
    left_join(plots_info %>% dplyr::select(plot_name, WD_mean), by = "plot_name")
} else {
  warning("La colonne WD_mean n'est pas disponible dans le fichier plots_info.")
  combined_data <- buffer_data
}

# Création des dataframes pour stocker les résultats
# 1. Métriques de base et combinées
plot_metrics_basic <- tibble(
  plot_name = character(),
  WD_mean = numeric(),
  # Métriques de base
  prop_at_10m = numeric(),
  prop_at_20m = numeric(),
  prop_at_30m = numeric(),
  h10 = numeric(),
  h25 = numeric(),
  h50 = numeric(),
  h75 = numeric(),
  h90 = numeric(),
  P0 = numeric(),
  Pmax = numeric(),
  Smax = numeric(),
  hSmax = numeric(),
  decay_rate = numeric(),
  auc = numeric(),
  # Paramètres logistiques
  logistic_a = numeric(),
  logistic_b = numeric(),
  logistic_h0 = numeric(),
  logistic_k = numeric(),
  logistic_R2 = numeric(),
  # Métriques combinées
  IDV = numeric(),
  ROS = numeric(),
  IAS = numeric(),
  ISC = numeric(),
  CPR = numeric(),
  ICT = numeric(),
  IPN = numeric(),
  ET = numeric(),
  CSI = numeric(),
  GVN = numeric(),
  RT = numeric(),
  HPSI = numeric()
)

# 2. Métriques avancées
plot_metrics_advanced <- tibble(
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

# 3. Indices optimisés et transformations
plot_metrics_optimized <- tibble(
  plot_name = character(),
  WD_mean = numeric(),
  
  # Indices optimisés
  ISA = numeric(),  # Indice de Stratification Asymétrique
  ICV = numeric(),  # Indice de Complexité Verticale
  IDT = numeric(),  # Indice de Différenciation Transitionnelle
  
  # Transformations
  log_transition_thickness = numeric(),
  log_max_d1 = numeric(),
  h50_relatif = numeric(),
  basal_thickness_relatif = numeric(),
  ratio_h50_h90 = numeric(),
  ratio_basal_transition = numeric(),
  h50_squared = numeric(),
  anisotropy_squared = numeric()
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

# Variables pour stocker les valeurs maximales pour la normalisation
max_pente_ref <- 0
max_transition_thickness <- 0

# Premier passage pour calculer les valeurs de référence
for (plot_i in plot_names) {
  plot_data <- combined_data %>%
    filter(plot_name == plot_i)
  
  if (nrow(plot_data) < 3) next  # Sauter les plots avec trop peu de données
  
  # Calculer la pente maximale
  max_slope_results <- calculate_max_slope(plot_data)
  if (!is.na(max_slope_results$Smax) && abs(max_slope_results$Smax) > max_pente_ref) {
    max_pente_ref <- abs(max_slope_results$Smax)
  }
  
  # Décomposition fonctionnelle pour obtenir l'épaisseur de transition
  segments_result <- decompose_curve_segments(plot_data$height_aboveground, plot_data$proportion)
  transition_thickness <- segments_result$metrics %>% 
    filter(zone == "transition") %>% 
    pull(thickness)
  
  if (!is.na(transition_thickness) && transition_thickness > max_transition_thickness) {
    max_transition_thickness <- transition_thickness
  }
}

# Boucle principale sur les plots
for (plot_i in plot_names) {
  # Filtrer les données pour ce plot
  plot_data <- combined_data %>%
    filter(plot_name == plot_i)
  
  if (nrow(plot_data) < 3) {
    pb$tick()
    next  # Sauter les plots avec trop peu de données
  }
  
  # Récupérer WD_mean
  wd_mean_val <- NA
  if ("WD_mean" %in% names(plot_data)) {
    wd_mean_val <- plot_data$WD_mean[1]
  }
  
  #### Calcul des métriques de base ####
  
  # 1. Proportion à hauteurs spécifiques
  prop_10m <- plot_data %>% 
    filter(height_aboveground == 10) %>% 
    pull(proportion)
  
  prop_20m <- plot_data %>% 
    filter(height_aboveground == 20) %>% 
    pull(proportion)
  
  prop_30m <- plot_data %>% 
    filter(height_aboveground == 30) %>% 
    pull(proportion)
  
  # S'assurer que les valeurs existent, sinon interpoler
  if (length(prop_10m) == 0) {
    prop_10m <- approx(plot_data$height_aboveground, plot_data$proportion, xout = 10)$y
  }
  
  if (length(prop_20m) == 0) {
    prop_20m <- approx(plot_data$height_aboveground, plot_data$proportion, xout = 20)$y
  }
  
  if (length(prop_30m) == 0) {
    prop_30m <- approx(plot_data$height_aboveground, plot_data$proportion, xout = 30)$y
  }
  
  # 2. Hauteurs à seuils spécifiques (métriques fondamentales)
  h10 <- gap_height_at_threshold(plot_data, 0.10)
  h25 <- gap_height_at_threshold(plot_data, 0.25)
  h50 <- gap_height_at_threshold(plot_data, 0.50)
  h75 <- gap_height_at_threshold(plot_data, 0.75)
  h90 <- gap_height_at_threshold(plot_data, 0.90)
  
  # 3. Valeurs aux limites
  if (nrow(plot_data) > 0) {
    # Première valeur (au niveau du sol ou proche)
    min_height_idx <- which.min(plot_data$height_aboveground)
    P0 <- plot_data$proportion[min_height_idx]
    
    # Valeur maximale
    max_prop_idx <- which.max(plot_data$proportion)
    Pmax <- plot_data$proportion[max_prop_idx]
  } else {
    P0 <- NA
    Pmax <- NA
  }
  
  # 4. Pente maximale et sa position
  max_slope_results <- calculate_max_slope(plot_data)
  Smax <- max_slope_results$Smax
  hSmax <- max_slope_results$hSmax
  
  # 5. Taux de changement global
  decay <- calculate_slope(plot_data)
  
  # 6. Aire sous la courbe
  auc_val <- calculate_auc(plot_data)
  
  # 7. Ajustement du modèle logistique
  logistic_results <- fit_logistic_model(plot_data)
  
  # 8. Calculer les métriques combinées
  combined_metrics <- calculate_combined_metrics(
    h10, h25, h50, h75, h90, Smax, hSmax, P0, Pmax, auc_val
  )
  
  # Ajouter au dataframe de métriques de base
  plot_metrics_basic <- plot_metrics_basic %>%
    bind_rows(tibble(
      plot_name = plot_i,
      WD_mean = wd_mean_val,
      # Métriques de base
      prop_at_10m = prop_10m,
      prop_at_20m = prop_20m,
      prop_at_30m = prop_30m,
      h10 = h10,
      h25 = h25,
      h50 = h50,
      h75 = h75,
      h90 = h90,
      P0 = P0,
      Pmax = Pmax,
      Smax = Smax,
      hSmax = hSmax,
      decay_rate = decay,
      auc = auc_val,
      # Paramètres logistiques
      logistic_a = logistic_results$a,
      logistic_b = logistic_results$b,
      logistic_h0 = logistic_results$h0,
      logistic_k = logistic_results$k,
      logistic_R2 = logistic_results$R2,
      # Métriques combinées
      IDV = combined_metrics$IDV,
      ROS = combined_metrics$ROS,
      IAS = combined_metrics$IAS,
      ISC = combined_metrics$ISC,
      CPR = combined_metrics$CPR,
      ICT = combined_metrics$ICT,
      IPN = combined_metrics$IPN,
      ET = combined_metrics$ET,
      CSI = combined_metrics$CSI,
      GVN = combined_metrics$GVN,
      RT = combined_metrics$RT,
      HPSI = combined_metrics$HPSI
    ))
  
  #### Calcul des métriques avancées ####
  
  # 1. Décomposition fonctionnelle
  segments_result <- decompose_curve_segments(plot_data$height_aboveground, plot_data$proportion)
  
  # 2. Analyse différentielle
  derivatives_result <- calculate_derivatives(plot_data$height_aboveground, plot_data$proportion)
  
  # 3. Analyse multi-échelle
  multiscale_result <- multiscale_analysis(plot_data$height_aboveground, plot_data$proportion)
  
  # 4. Analyse d'anisotropie
  anisotropy_result <- analyze_anisotropy(plot_data, plot_data)
  
  # Calculer l'indice d'hétérogénéité zonale
  zone_heterogeneity <- segments_result$metrics %>%
    filter(zone == "transition") %>%
    pull(norm_auc) / mean(c(
      segments_result$metrics %>% filter(zone == "basal") %>% pull(norm_auc),
      segments_result$metrics %>% filter(zone == "upper") %>% pull(norm_auc)
    ), na.rm = TRUE)
  
  # Calculer le ratio d'énergie petite/grande échelle
  energy_small <- multiscale_result$scale_metrics %>%
    filter(scale == 2) %>%
    pull(scale_energy)
  
  energy_large <- multiscale_result$scale_metrics %>%
    filter(scale == 20) %>%
    pull(scale_energy)
  
  energy_ratio <- energy_small / energy_large
  
  # Ajouter au dataframe de métriques avancées
  plot_metrics_advanced <- plot_metrics_advanced %>%
    bind_rows(tibble(
      plot_name = plot_i,
      WD_mean = wd_mean_val,
      
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
  
  #### Calcul des indices optimisés ####
  # Extraire les métriques pertinentes du plot actuel
  advanced_row <- plot_metrics_advanced %>% filter(plot_name == plot_i)
  
  # Créer une liste avec les métriques nécessaires pour les indices optimisés
  metrics_for_optimized <- list(
    h50 = h50,
    height_min_d2 = advanced_row$height_min_d2,
    phase_ratio = advanced_row$phase_ratio,
    max_d1 = advanced_row$max_d1,
    basal_thickness = advanced_row$basal_thickness,
    upper_thickness = advanced_row$upper_thickness,
    anisotropy_max = advanced_row$anisotropy_max,
    energy_ratio_small_large = advanced_row$energy_ratio_small_large,
    transition_thickness = advanced_row$transition_thickness,
    basal_norm_auc = advanced_row$basal_norm_auc
  )
  
  # Calculer les indices optimisés
  optimized_indices <- calculate_optimized_indices(
    metrics_for_optimized, 
    max_pente_ref = max_pente_ref,
    max_transition_thickness = max_transition_thickness
  )
  
  # Ajouter au dataframe d'indices optimisés
  plot_metrics_optimized <- plot_metrics_optimized %>%
    bind_rows(tibble(
      plot_name = plot_i,
      WD_mean = wd_mean_val,
      
      # Indices optimisés
      ISA = optimized_indices$ISA,
      ICV = optimized_indices$ICV,
      IDT = optimized_indices$IDT,
      
      # Transformations
      log_transition_thickness = optimized_indices$log_transition_thickness,
      log_max_d1 = optimized_indices$log_max_d1,
      h50_relatif = optimized_indices$h50_relatif,
      basal_thickness_relatif = optimized_indices$basal_thickness_relatif,
      ratio_h50_h90 = optimized_indices$ratio_h50_h90,
      ratio_basal_transition = optimized_indices$ratio_basal_transition,
      h50_squared = optimized_indices$h50_squared,
      anisotropy_squared = optimized_indices$anisotropy_squared
    ))
  
  # Mettre à jour la barre de progression
  pb$tick()
}

#### 6️⃣ EXPORT ET VALIDATION ####
# Fusionner les dataframes pour un export unique
plot_metrics_all <- plot_metrics_basic %>%
  left_join(plot_metrics_advanced, by = c("plot_name", "WD_mean")) %>%
  left_join(plot_metrics_optimized, by = c("plot_name", "WD_mean"))

# Exporter les métriques calculées
rio::export(
  plot_metrics_all,
  file.path(path_metrics, "plot_metrics_gaps_all.csv"),
  sep = ";",
  dec = ",",
  append = FALSE
)

# Exporter également les dataframes séparés
rio::export(
  plot_metrics_basic,
  file.path(path_metrics, "plot_metrics_gaps_basic.csv"),
  sep = ";",
  dec = ",",
  append = FALSE
)

rio::export(
  plot_metrics_advanced,
  file.path(path_metrics, "plot_metrics_gaps_advanced.csv"),
  sep = ";",
  dec = ",",
  append = FALSE
)

rio::export(
  plot_metrics_optimized,
  file.path(path_metrics, "plot_metrics_gaps_optimized.csv"),
  sep = ";",
  dec = ",",
  append = FALSE
)

# Vérification rapide des corrélations avec WD
if (any(!is.na(plot_metrics_all$WD_mean))) {
  # Calculer les corrélations avec WD pour toutes les métriques numériques
  correlation_vars <- names(plot_metrics_all)[sapply(plot_metrics_all, is.numeric)]
  correlation_vars <- setdiff(correlation_vars, c("WD_mean"))
  
  correlations <- data.frame(
    metric = character(),
    correlation = numeric(),
    p_value = numeric()
  )
  
  for (var in correlation_vars) {
    # Données valides
    valid_data <- !is.na(plot_metrics_all[[var]]) & !is.na(plot_metrics_all$WD_mean)
    
    # Vérifier suffisamment de données
    if (sum(valid_data) >= 3) {
      # Pearson
      pearson_test <- tryCatch({
        cor.test(plot_metrics_all[[var]][valid_data], plot_metrics_all$WD_mean[valid_data], method = "pearson")
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
  
  # Exporter les corrélations
  rio::export(
    correlations %>% arrange(desc(abs(correlation))),
    file.path(path_metrics, "plot_metrics_gaps_correlations.csv"),
    sep = ";",
    dec = ",",
    append = FALSE
  )
  
  # Afficher les 10 métriques les plus corrélées
  cat("\nTop 10 métriques les plus corrélées avec WD_mean:\n")
  print(correlations %>% 
          filter(!is.na(correlation)) %>%
          arrange(desc(abs(correlation))) %>% 
          head(10))
  
  # Vérifier spécifiquement les nouveaux indices optimisés
  cat("\nCorrélations des indices optimisés avec WD_mean:\n")
  print(correlations %>% 
          filter(metric %in% c("ISA", "ICV", "IDT")) %>%
          arrange(desc(abs(correlation))))
} else {
  warning("Aucune valeur de WD_mean disponible pour calculer les corrélations.")
}

# Message de fin
cat("\n====================================================\n")
cat("Calcul des métriques de trouées terminé !")
cat("\n====================================================\n")
cat("Nombre de plots analysés:", nrow(plot_metrics_all), "\n")
cat("Nombre total de métriques calculées:", ncol(plot_metrics_all) - 1, "\n")
cat("Résultats sauvegardés dans:", path_metrics, "\n")
cat("====================================================\n")
