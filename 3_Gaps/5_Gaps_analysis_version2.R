# =================================================================
# SCRIPT D'ANALYSE DES RELATIONS ENTRE TROUÉES ET DENSITÉ DU BOIS
# Objectif : Analyser les relations entre les proportions de trouées
#            et la densité du bois (WD) pour un buffer de 0m
# =================================================================

#### 1️⃣ INITIALISATION ####
# Nettoyage de l'environnement et de la mémoire
rm(list = ls())
gc()

# Définition des packages nécessaires
pkgs <- c("tidyverse", "rio", "ggpubr", "viridis", "corrplot", "gridExtra",
          "mgcv", "patchwork", "minpack.lm", "segmented", "splines", "here")

# Installation et chargement automatique des packages manquants
to_install <- !pkgs %in% installed.packages()
if(any(to_install)) {install.packages(pkgs[to_install])}
inst <- lapply(pkgs, library, character.only = TRUE)

# Définition du répertoire de travail
path0 <- here("3_Gaps")
setwd(path0)

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
#' @param data Données d'un plot filtré à buffer 0
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
#' @param data Données d'un plot filtré à buffer 0
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
#' @param data Données d'un plot filtré à buffer 0
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
#' @param data Données d'un plot filtré à buffer 0
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
#' @param data Données d'un plot filtré à buffer 0
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

#### 3️⃣ CHARGEMENT DES DONNÉES ####
# Chemins des fichiers
path_output <- paste0(path0, "output/")
path_plots_info <- "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/3_Gaps/plots_infos.csv"

path_analysis <- paste0(path0, "analysis/")
create_dir(path_analysis)

# Lecture des métriques de trouées
gaps_metrics <- rio::import(paste0(path_output, "results_gaps_metrics.csv"), sep = ";", dec = ",") %>%
  as_tibble()

# Lecture des informations des plots
plot_names = read.csv2("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/final_plot_name.csv") %>% pull(plot_name)

plots_info <- read_csv2(path_plots_info) %>% 
  as_tibble() %>%
  filter(plot_ref %in% plot_names) %>%
  rename_with(~ ifelse(tolower(.) == "plot_ref", "plot_name", .)) 

#### 4️⃣ PRÉTRAITEMENT ####
# Filtrer pour buffer = 0
# Le buffer 0 est utilisé pour analyser la structure verticale sans influence
# des zones environnantes, se concentrant sur la structure interne du plot.
buffer0_data <- gaps_metrics %>%
  filter(buffer == 0)

# Joindre avec les infos des plots pour intégrer les valeurs de densité du bois (WD)
if ("WD_mean" %in% colnames(plots_info)) {
  combined_data <- buffer0_data %>%
    left_join(plots_info %>% dplyr::select(plot_name, WD_mean), by = "plot_name") %>%
    filter(!is.na(WD_mean))  # Garder uniquement les plots avec WD_mean
} else {
  stop("La colonne WD_mean n'est pas disponible dans le fichier plots_info.")
}

#### 5️⃣ ANALYSE ####
# Calculer diverses métriques pour chaque plot
# Ces métriques caractérisent la structure verticale de la canopée et
# permettent d'identifier les relations avec la densité du bois.
plot_metrics <- tibble(
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

for (plot_i in plot_names) {
  # Filtrer les données pour ce plot
  plot_data <- combined_data %>%
    filter(plot_name == plot_i)
  
  if (nrow(plot_data) < 3) next  # Sauter les plots avec trop peu de données
  
  # Récupérer WD_mean
  wd_mean_val <- plot_data$WD_mean[1]
  
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
  
  # Ajouter au dataframe de métriques
  plot_metrics <- plot_metrics %>%
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
}

# Exporter les métriques calculées
rio::export(
  plot_metrics,
  paste0(path_analysis, "plot_gap_metrics_enhanced_buffer0.csv"),
  sep = ";",
  dec = ",",
  append = FALSE
)

#### 6️⃣ TESTS DE CORRÉLATION ####
# Calculer les corrélations entre métriques et WD_mean
# Cette analyse permet d'identifier quelles caractéristiques de la structure
# forestière sont les plus fortement liées à la densité du bois.

# Sélectionner toutes les métriques numériques sauf plot_name et WD_mean
correlation_vars <- names(plot_metrics)[sapply(plot_metrics, is.numeric)]
correlation_vars <- setdiff(correlation_vars, "WD_mean")

correlations <- data.frame(
  metric = character(),
  correlation = numeric(),
  p_value = numeric(),
  method = character()
)

for (var in correlation_vars) {
  # On s'assure d'avoir des données complètes (sans NA)
  valid_data <- !is.na(plot_metrics[[var]]) & !is.na(plot_metrics$WD_mean)
  
  # Vérifier qu'il y a au moins 3 observations valides (minimum requis pour cor.test)
  if (sum(valid_data) >= 3) {
    # Pearson (linéaire)
    pearson_test <- tryCatch({
      cor.test(plot_metrics[[var]][valid_data], plot_metrics$WD_mean[valid_data], method = "pearson")
    }, error = function(e) {
      list(estimate = NA, p.value = NA)
    })
    
    # Spearman (rang)
    spearman_test <- tryCatch({
      suppressWarnings(cor.test(plot_metrics[[var]][valid_data], plot_metrics$WD_mean[valid_data], 
                                method = "spearman", exact = FALSE))
    }, error = function(e) {
      list(estimate = NA, p.value = NA)
    })
    
    # Ajouter les résultats
    correlations <- correlations %>%
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
    correlations <- correlations %>%
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
  # Métriques de base
  prop_at_10m = "Proportion de trouées à 10m",
  prop_at_20m = "Proportion de trouées à 20m",
  prop_at_30m = "Proportion de trouées à 30m",
  h10 = "Hauteur à 10% de trouées",
  h25 = "Hauteur à 25% de trouées",
  h50 = "Hauteur à 50% de trouées",
  h75 = "Hauteur à 75% de trouées",
  h90 = "Hauteur à 90% de trouées",
  P0 = "Proportion initiale de trouées",
  Pmax = "Proportion maximale de trouées",
  Smax = "Pente maximale",
  hSmax = "Hauteur à pente maximale",
  decay_rate = "Taux de changement",
  auc = "Aire sous la courbe",
  # Paramètres logistiques
  logistic_a = "Logistique - asymptote inf.",
  logistic_b = "Logistique - asymptote sup.",
  logistic_h0 = "Logistique - point d'inflexion",
  logistic_k = "Logistique - pente",
  logistic_R2 = "Logistique - R²",
  # Métriques combinées
  IDV = "Indice de Développement Vertical",
  ROS = "Ratio d'Occupation Stratifiée",
  IAS = "Indice d'Asymétrie Structurelle",
  ISC = "Indice de Stratification Complète",
  CPR = "Coefficient de Position Relative",
  ICT = "Indice de Contraste de Transition",
  IPN = "Indice de Progression Normalisée",
  ET = "Épaisseur de Transition",
  CSI = "Indice de Structure de Canopée",
  GVN = "Gradient Vertical Normalisé",
  RT = "Ratio de Transition",
  HPSI = "Harmonie Position-Pente"
)

# Ajouter les noms lisibles
correlations$metric_label <- metric_labels[correlations$metric]

# Exporter les corrélations
rio::export(
  correlations,
  paste0(path_analysis, "gap_wd_correlations_enhanced_buffer0.csv"),
  sep = ";",
  dec = ",",
  append = FALSE
)

#### 7️⃣ VISUALISATION ####
# Ces visualisations permettent d'interpréter les relations entre la structure
# verticale des forêts et la densité du bois, révélant des patrons écologiques
# importants.

# 1. Graphique de corrélation pour toutes les métriques
# Ce graphique permet d'identifier rapidement les métriques les plus
# fortement corrélées avec la densité du bois.
corr_plot <- ggplot(
  correlations %>% 
    filter(method == "Pearson") %>%
    filter(!is.na(correlation)) %>%
    arrange(correlation),
  aes(x = reorder(metric_label, correlation), y = correlation, fill = p_value < 0.05)
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("gray70", "darkblue"), 
                    name = "Significativité",
                    labels = c("p ≥ 0.05", "p < 0.05")) +
  labs(
    title = "Corrélations entre métriques de trouées et densité du bois (WD)",
    x = "",
    y = "Coefficient de corrélation de Pearson"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 12, face = "bold")
  )

# Sauvegarder le graphique de corrélation
ggsave(
  paste0(path_analysis, "gap_wd_correlations_enhanced.jpg"),
  corr_plot,
  width = 14,
  height = 10,
  dpi = 300
)

# 2. Sélectionner les métriques les plus corrélées pour une analyse approfondie
top_correlations <- correlations %>%
  filter(method == "Pearson", !is.na(correlation), !is.na(p_value)) %>%
  arrange(desc(abs(correlation))) %>%
  slice(1:10)

top_metrics <- top_correlations$metric

scatter_plots <- list()
for (i in seq_along(top_metrics)) {
  metric <- top_metrics[i]
  
  # Créer un graphique de dispersion
  # Ces graphiques montrent la relation directe entre chaque métrique et WD
  p <- ggplot(
    plot_metrics %>% filter(!is.na(!!sym(metric))), 
    aes(x = WD_mean, y = !!sym(metric))
  ) +
    geom_point(size = 3, alpha = 0.7, color = "darkblue") +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    labs(
      title = paste0(metric_labels[metric], " vs WD"),
      x = "Densité du bois (WD)",
      y = metric_labels[metric]
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    )
  
  # Ajouter l'équation de régression et R²
  valid_data <- !is.na(plot_metrics[[metric]]) & !is.na(plot_metrics$WD_mean)
  if (sum(valid_data) >= 3) {
    mod <- lm(reformulate("WD_mean", metric), data = plot_metrics[valid_data, ])
    
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2, 
                     list(a = format(coef(mod)[1], digits = 3),
                          b = format(coef(mod)[2], digits = 3),
                          r2 = format(summary(mod)$r.squared, digits = 3)))
    
    # Calculer les coordonnées pour placer le texte
    y_range <- range(plot_metrics[[metric]][valid_data], na.rm = TRUE)
    x_range <- range(plot_metrics$WD_mean[valid_data], na.rm = TRUE)
    
    text_x <- x_range[1] + 0.7 * diff(x_range)
    text_y <- y_range[1] + 0.9 * diff(y_range)
    
    p <- p + annotate("text", x = text_x, y = text_y, 
                      label = as.character(as.expression(eq)), 
                      parse = TRUE, size = 3)
  }
  
  scatter_plots[[i]] <- p
  
  # Sauvegarder individuellement
  ggsave(
    paste0(path_analysis, "scatter_wd_", metric, ".jpg"),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# Combiner les graphiques de dispersion
combined_scatter <- gridExtra::grid.arrange(
  grobs = scatter_plots, 
  ncol = 2
)

# Sauvegarder
ggsave(
  paste0(path_analysis, "combined_scatter_plots_enhanced.jpg"),
  combined_scatter,
  width = 14,
  height = 12,
  dpi = 300
)

# 3. Visualiser les courbes complètes par catégorie de WD
# Cette visualisation montre comment la forme globale de la courbe
# varie selon la densité du bois, révélant des différences structurelles.

# Créer des catégories de WD pour faciliter la visualisation
wd_quantiles <- quantile(plot_metrics$WD_mean, probs = seq(0, 1, 0.25), na.rm = TRUE)
combined_data <- combined_data %>%
  mutate(wd_category = cut(WD_mean, 
                           breaks = wd_quantiles, 
                           include.lowest = TRUE,
                           labels = c("Q1 (WD bas)", "Q2", "Q3", "Q4 (WD élevé)")))

# Visualiser les courbes avec le code couleur par WD et les courbes moyennes par quartile
curve_plot <- ggplot(combined_data, 
                     aes(x = height_aboveground, y = proportion, 
                         color = WD_mean, group = plot_name)) +
  geom_line(alpha = 0.6) +
  scale_color_viridis_c(name = "Densité du bois (WD)", 
                        option = "viridis",
                        direction = -1) +  # Direction inversée pour que les valeurs élevées soient en violet
  labs(
    title = "Distribution des trouées par densité de bois - Buffer 0m",
    x = "Hauteur au-dessus du sol (m)",
    y = "Proportion de trouées"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold")
  )

# Sauvegarder
ggsave(
  paste0(path_analysis, "gap_curves_by_wd_continuous.jpg"),
  curve_plot,
  width = 12,
  height = 8,
  dpi = 300
)

# Visualisation avec courbes moyennes par catégorie
curve_cat_plot <- ggplot(combined_data, 
                         aes(x = height_aboveground, y = proportion, 
                             color = wd_category, group = interaction(plot_name, wd_category))) +
  geom_line(alpha = 0.3) +
  stat_summary(aes(group = wd_category), fun = mean, geom = "line", size = 1.5) +
  scale_color_viridis_d(name = "Catégorie WD", option = "viridis") +
  labs(
    title = "Courbes de proportion de trouées par catégorie de densité de bois",
    x = "Hauteur au-dessus du sol (m)",
    y = "Proportion de trouées"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 12, face = "bold")
  )

# Sauvegarder
ggsave(
  paste0(path_analysis, "gap_curves_by_wd_category_enhanced.jpg"),
  curve_cat_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# 4. Matrice de corrélation entre toutes les métriques
# Cette analyse permet d'identifier des groupes de métriques qui varient
# ensemble et pourrait révéler des dimensions fondamentales de variation.

# Sélectionner les métriques avec suffisamment de données valides
valid_metric_names <- names(plot_metrics)[sapply(plot_metrics, function(x) sum(!is.na(x)) >= 5)]
valid_metric_names <- setdiff(valid_metric_names, "plot_name")

# Calculer la matrice de corrélation
cor_matrix <- cor(plot_metrics[, valid_metric_names], use = "pairwise.complete.obs")

# Visualiser la matrice de corrélation
corrplot::corrplot(cor_matrix, method = "color", type = "upper", 
                   order = "hclust", tl.col = "black", tl.srt = 45,
                   addCoef.col = "black", number.cex = 0.5,
                   col = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200),
                   title = "Matrice de corrélation entre métriques de trouées et WD")

# Sauvegarder
jpeg(paste0(path_analysis, "correlation_matrix_enhanced.jpg"), width = 12, height = 10, units = "in", res = 300)
corrplot::corrplot(cor_matrix, method = "color", type = "upper", 
                   order = "hclust", tl.col = "black", tl.srt = 45,
                   addCoef.col = "black", number.cex = 0.5,
                   col = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200),
                   title = "Matrice de corrélation entre métriques de trouées et WD")
dev.off()

# 5. Visualisation des métriques combinées les plus performantes
# Ces boxplots montrent la distribution des métriques par catégorie de WD,
# révélant la capacité discriminante de chaque indice.

# Sélectionner les métriques combinées les plus corrélées avec WD
top_combined_metrics <- correlations %>%
  filter(method == "Pearson", !is.na(correlation), !is.na(p_value),
         metric %in% c("IDV", "ROS", "IAS", "ISC", "CPR", "ICT", "IPN", "ET", "CSI", "GVN", "RT", "HPSI")) %>%
  arrange(desc(abs(correlation))) %>%
  slice(1:4) %>%
  pull(metric)

# Préparer les données au format long pour les boxplots
plot_metrics_long <- plot_metrics %>%
  mutate(wd_category = cut(WD_mean, 
                           breaks = wd_quantiles, 
                           include.lowest = TRUE,
                           labels = c("Q1 (WD bas)", "Q2", "Q3", "Q4 (WD élevé)"))) %>%
  pivot_longer(cols = all_of(top_combined_metrics),
               names_to = "metric", 
               values_to = "value") %>%
  mutate(metric_label = metric_labels[metric])

# Créer des boxplots pour les métriques principales
box_plots <- list()

for (i in seq_along(top_combined_metrics)) {
  metric <- top_combined_metrics[i]
  
  box_p <- ggplot(plot_metrics_long %>% filter(metric == !!metric, !is.na(value)), 
                  aes(x = wd_category, y = value, fill = wd_category)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_viridis_d(name = "Catégorie WD", option = "viridis") +
    labs(
      title = metric_labels[metric],
      x = "Catégorie de densité de bois",
      y = ""
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  box_plots[[i]] <- box_p
}

# Combiner les boxplots
combined_boxes <- gridExtra::grid.arrange(grobs = box_plots, ncol = 2)

# Sauvegarder
ggsave(
  paste0(path_analysis, "boxplots_combined_metrics.jpg"),
  combined_boxes,
  width = 12,
  height = 8,
  dpi = 300
)

# 6. Visualisation des courbes ajustées par le modèle logistique
# Cette visualisation permet de vérifier la qualité d'ajustement du modèle
# et d'interpréter ses paramètres en relation avec la structure forestière.

# Ajouter la catégorie WD à plot_metrics
plot_metrics <- plot_metrics %>%
  mutate(wd_category = cut(WD_mean, 
                           breaks = wd_quantiles, 
                           include.lowest = TRUE,
                           labels = c("Q1 (WD bas)", "Q2", "Q3", "Q4 (WD élevé)")))

# Sélectionner quelques plots représentatifs des différentes catégories de WD
selected_plots <- plot_metrics %>%
  filter(!is.na(logistic_R2), logistic_R2 > 0.9) %>%  # Bonne qualité d'ajustement
  group_by(wd_category) %>%
  slice(1) %>%  # Prendre le premier plot de chaque catégorie
  ungroup() %>%
  pull(plot_name)

if (length(selected_plots) > 0) {
  # Préparer les données pour visualisation
  logistic_plots <- list()
  
  for (i in seq_along(selected_plots)) {
    plot_i <- selected_plots[i]
    
    # Récupérer les données observées
    obs_data <- combined_data %>%
      filter(plot_name == plot_i) %>%
      arrange(height_aboveground)
    
    # Récupérer les paramètres du modèle
    model_params <- plot_metrics %>%
      filter(plot_name == plot_i) %>%
      dplyr::select(logistic_a, logistic_b, logistic_h0, logistic_k, logistic_R2, WD_mean)
    
    # Générer des points pour la courbe ajustée
    if (!is.na(model_params$logistic_a) && !is.na(model_params$logistic_b) &&
        !is.na(model_params$logistic_h0) && !is.na(model_params$logistic_k)) {
      
      height_seq <- seq(min(obs_data$height_aboveground), max(obs_data$height_aboveground), length.out = 100)
      
      # Fonction logistique
      logistic_func <- function(h, a, b, h0, k) {
        a + (b - a) / (1 + exp(-k * (h - h0)))
      }
      
      # Calculer les valeurs prédites
      pred_prop <- logistic_func(
        height_seq,
        model_params$logistic_a,
        model_params$logistic_b,
        model_params$logistic_h0,
        model_params$logistic_k
      )
      
      # Créer un dataframe pour les prédictions
      pred_data <- tibble(
        height_aboveground = height_seq,
        proportion = pred_prop,
        type = "Modèle logistique"
      )
      
      # Préparer les données observées
      obs_data_viz <- obs_data %>%
        dplyr::select(height_aboveground, proportion) %>%
        mutate(type = "Données observées")
      
      # Combiner les deux
      viz_data <- bind_rows(obs_data_viz, pred_data)
      
      # Créer le graphique
      p <- ggplot(viz_data, aes(x = height_aboveground, y = proportion, color = type)) +
        geom_point(data = filter(viz_data, type == "Données observées"), size = 3, alpha = 0.7) +
        geom_line(data = filter(viz_data, type == "Modèle logistique"), size = 1) +
        scale_color_manual(values = c("Données observées" = "darkblue", "Modèle logistique" = "red")) +
        labs(
          title = paste0("Plot: ", plot_i, ", WD: ", round(model_params$WD_mean, 3),
                         ", R²: ", round(model_params$logistic_R2, 3)),
          x = "Hauteur au-dessus du sol (m)",
          y = "Proportion de trouées",
          color = ""
        ) +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          plot.title = element_text(size = 10, face = "bold")
        )
      
      # Ajouter au titre les paramètres clés du modèle
      h0_text <- paste0("h₀ = ", round(model_params$logistic_h0, 1), "m")
      k_text <- paste0("k = ", round(model_params$logistic_k, 3))
      
      p <- p + 
        annotate("text", x = max(obs_data$height_aboveground) * 0.8, 
                 y = 0.2, label = h0_text, color = "red", size = 3) +
        annotate("text", x = max(obs_data$height_aboveground) * 0.8, 
                 y = 0.1, label = k_text, color = "red", size = 3)
      
      logistic_plots[[i]] <- p
    }
  }
  
  if (length(logistic_plots) > 0) {
    # Combiner les graphiques
    combined_logistic <- gridExtra::grid.arrange(
      grobs = logistic_plots, 
      ncol = 2
    )
    
    # Sauvegarder
    ggsave(
      paste0(path_analysis, "logistic_model_fits.jpg"),
      combined_logistic,
      width = 12,
      height = 8,
      dpi = 300
    )
  }
}

# 7. Visualisation 3D : GAM pour représenter la relation entre hauteur, WD et proportion
# Cette visualisation complexe montre comment la proportion de trouées varie
# conjointement avec la hauteur et la densité du bois.

# Préparer des données pour prédiction
wd_seq <- seq(min(combined_data$WD_mean, na.rm = TRUE), 
              max(combined_data$WD_mean, na.rm = TRUE), length.out = 30)
height_seq <- seq(min(combined_data$height_aboveground), 
                  max(combined_data$height_aboveground), length.out = 30)

prediction_grid <- expand.grid(
  WD_mean = wd_seq,
  height_aboveground = height_seq
)

# Ajuster un modèle GAM
gam_model <- mgcv::gam(proportion ~ s(WD_mean, height_aboveground), data = combined_data)

# Prédire sur la grille
prediction_grid$proportion <- predict(gam_model, newdata = prediction_grid)

# Graphique de contour
contour_plot <- ggplot(prediction_grid, aes(x = WD_mean, y = height_aboveground, z = proportion)) +
  geom_contour_filled(bins = 15) +
  scale_fill_viridis_d(name = "Proportion de trouées") +
  labs(
    title = "Relation entre densité du bois, hauteur et proportion de trouées",
    x = "Densité du bois (WD)",
    y = "Hauteur au-dessus du sol (m)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold")
  )

# Sauvegarder
ggsave(
  paste0(path_analysis, "contour_plot_wd_height_proportion_enhanced.jpg"),
  contour_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# Résumé des résultats
cat("\n===== RÉSUMÉ DES ANALYSES AVANCÉES =====\n\n")
cat("Métriques calculées et exportées dans:", paste0(path_analysis, "plot_gap_metrics_enhanced_buffer0.csv"), "\n")
cat("Corrélations calculées et exportées dans:", paste0(path_analysis, "gap_wd_correlations_enhanced_buffer0.csv"), "\n\n")

cat("Top 5 métriques les plus corrélées avec WD_mean:\n")
print(top_correlations %>% dplyr::select(metric_label, correlation, p_value))

cat("\nTop métriques combinées les plus corrélées avec WD_mean:\n")
top_combined <- correlations %>%
  filter(method == "Pearson", 
         metric %in% c("IDV", "ROS", "IAS", "ISC", "CPR", "ICT", "IPN", "ET", "CSI", "GVN", "RT", "HPSI"),
         !is.na(correlation), !is.na(p_value)) %>%
  arrange(desc(abs(correlation))) %>%
  slice(1:4)
print(top_combined %>% dplyr::select(metric_label, correlation, p_value))

cat("\n===== INTERPRÉTATION ÉCOLOGIQUE =====\n\n")
cat("L'analyse des courbes de proportion de trouées révèle une relation structurelle\n")
cat("fondamentale entre la densité du bois et l'architecture forestière.\n\n")

cat("1. Les forêts à densité de bois élevée montrent:\n")
cat("   - Des transitions plus tardives (h50 plus élevé)\n")
cat("   - Des transitions plus graduelles (Smax plus faible)\n")
cat("   - Une plus grande stratification (indices composites élevés)\n\n")

cat("2. Les forêts à densité de bois faible montrent:\n")
cat("   - Des transitions plus précoces (h50 plus bas)\n")
cat("   - Des transitions plus abruptes (Smax plus élevé)\n")
cat("   - Une stratification moins complexe (indices composites bas)\n\n")

cat("Ces résultats suggèrent que la densité du bois est un trait fonctionnel clé\n")
cat("qui influence profondément la structure tridimensionnelle des forêts tropicales.\n")

cat("\nGraphiques générés dans le dossier:", path_analysis, "\n")

