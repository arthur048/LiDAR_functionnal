# =================================================================
# ANALYSE APPROFONDIE DES CORRÉLATIONS AVEC LA DENSITÉ DU BOIS
# Objectif : Explorer en détail les relations entre les métriques de structure
#            verticale et différentes variables cibles
# =================================================================

#### 1️⃣ INITIALISATION ####
# Nettoyage de l'environnement et de la mémoire
rm(list = ls())
gc()

# Définition des packages nécessaires
pkgs <- c("tidyverse", "ggpubr", "viridis", "corrplot", "gridExtra",
          "patchwork", "ggridges", "zoo", "fields", "ggrepel", "cowplot", "here")

# Installation et chargement automatique des packages manquants
to_install <- !pkgs %in% installed.packages()
if(any(to_install)) {install.packages(pkgs[to_install])}
inst <- lapply(pkgs, library, character.only = TRUE)

# Configuration de l'environnement
options(scipen = 999)        # Éviter la notation scientifique
set.seed(42)                 # Fixer la graine pour reproductibilité

#### 2️⃣ DÉFINITION DES CHEMINS ET CHARGEMENT DES DONNÉES ####
# Définition du répertoire de base
project_dir <- here()

# Création d'une fonction pour les chemins relatifs au projet
project_path <- function(...) {
  here(...)
}

# Définition des chemins essentiels
path_main <- project_path("5_Main_analysis")
path_input <- file.path(path_main, "input")
path_output <- file.path(path_main, "output")
path_data <- file.path(path_output, "data")
path_figures <- file.path(path_output, "figures")
path_tables <- file.path(path_output, "tables")

# Sous-répertoires pour les analyses de corrélation
path_data_corr <- file.path(path_data, "correlation")
path_figures_corr <- file.path(path_figures, "correlation")
path_tables_corr <- file.path(path_tables, "correlation")

# S'assurer que les répertoires de sortie existent
dirs_to_check <- c(path_data, path_figures, path_tables,
                   path_data_corr, path_figures_corr, path_tables_corr)
for (dir in dirs_to_check) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("Répertoire créé:", dir, "\n")
  }
}

# Chargement des données
plot_metrics <- read_csv2(file.path(path_input, "plot_metrics.csv"))
gaps_metrics <- read_csv2(file.path(path_input, "gaps_metrics.csv"))
combined_data <- read_csv2(file.path(path_input, "combined_data.csv"))
plots_location <- read_csv2(file.path(project_dir, "plots_localisation.csv"))  %>%
  rename(plot_name = plot_ref)

# Charger les catégories de WD et la palette de couleurs
wd_quantiles <- quantile(plot_metrics$WD_BA, probs = seq(0, 1, 0.25), na.rm = TRUE)
wd_palette <- viridis(4, option = "viridis", direction = -1)
names(wd_palette) <- levels(combined_data$wd_category)

# Charger les fonctions utilitaires
source(file.path(path_main, "functions", "metric_labels.R"))
source(file.path(path_main, "functions", "graphics_correlation.R"))

# Obtenir les labels des métriques
metric_labels <- get_metric_labels()

# Paramètres de palette de couleurs
color_palette_discrete <- "turbo"  # Options: "magma", "inferno", "plasma", "viridis", "cividis", "turbo"
color_palette_continuous <- "turbo"  # Options: "magma", "inferno", "plasma", "viridis", "cividis", "turbo"

# Paramètres globaux pour l'analyse
# Ces paramètres contrôlent l'ensemble de l'analyse et peuvent être modifiés facilement
target_vars <- c("WD_BA", "WD_mean", "prop_g_helio", "prop_g_npld", 
                 "prop_g_shade", "prop_ind_helio", "prop_ind_npld", "prop_ind_shade")  # Variables cibles à analyser
comparison_pairs <- list(
  c("WD_BA", "WD_mean"),
  c("WD_BA", "prop_ind_helio"),
  c("WD_BA", "prop_ind_shade")
)  # Paires à comparer
corr_methods <- c("Pearson", "Spearman")  # Méthodes de corrélation à utiliser
top_n_value <- 20  # Nombre de métriques à retenir pour les visualisations
dpi_value <- 600  # Résolution des graphiques exportés

# Variables à exclure des variables explicatives
exclude_vars <- c("plot_name", target_vars)

#### 3️⃣ CRÉER DES SOUS-RÉPERTOIRES PAR VARIABLE CIBLE ####
# Cette section crée les sous-répertoires pour chaque variable cible

# Créer un sous-répertoire pour chaque variable cible
var_dirs <- c()
for (target_var in target_vars) {
  # Créer le répertoire
  path_target_var <- file.path(path_figures_corr, target_var)
  if (!dir.exists(path_target_var)) {
    dir.create(path_target_var, recursive = TRUE)
    cat("Répertoire créé pour", target_var, ":", path_target_var, "\n")
  }
  var_dirs <- c(var_dirs, path_target_var)
}

#### 4️⃣ ANALYSE DE CORRÉLATION AVEC LES VARIABLES CIBLES ####
# Cette section examine les relations entre les différentes métriques 
# structurelles et les variables cibles

# Fonction pour calculer les corrélations pour une variable cible
calculate_correlations <- function(target_var, exclude_vars) {
  correlation_vars <- setdiff(names(plot_metrics), exclude_vars)
  
  correlations <- data.frame()
  
  for (var in correlation_vars) {
    for (method in c("Pearson", "Spearman")) {
      # Données valides
      valid_data <- !is.na(plot_metrics[[var]]) & !is.na(plot_metrics[[target_var]])
      
      # Vérifier suffisamment de données
      if (sum(valid_data) >= 3) {
        if (method == "Pearson") {
          test_result <- tryCatch({
            cor.test(plot_metrics[[var]][valid_data], plot_metrics[[target_var]][valid_data], method = "pearson")
          }, error = function(e) {
            list(estimate = NA, p.value = NA)
          })
        } else {
          test_result <- tryCatch({
            suppressWarnings(cor.test(plot_metrics[[var]][valid_data], plot_metrics[[target_var]][valid_data], 
                                      method = "spearman", exact = FALSE))
          }, error = function(e) {
            list(estimate = NA, p.value = NA)
          })
        }
        
        correlations <- rbind(correlations, data.frame(
          metric = var,
          correlation = ifelse(is.na(test_result$estimate), NA, test_result$estimate),
          p_value = ifelse(is.na(test_result$p.value), NA, test_result$p.value),
          method = method,
          target_variable = target_var
        ))
      } else {
        correlations <- rbind(correlations, data.frame(
          metric = var,
          correlation = NA,
          p_value = NA,
          method = method,
          target_variable = target_var
        ))
      }
    }
  }
  
  return(correlations)
}

# Calculer les corrélations pour chaque variable cible
all_correlations <- data.frame()
for (target_var in target_vars) {
  correlations <- calculate_correlations(target_var, exclude_vars)
  all_correlations <- rbind(all_correlations, correlations)
}

# Convertir all_correlations en tibble pour faciliter les manipulations
all_correlations <- as_tibble(all_correlations)

# Ajouter les noms lisibles aux corrélations
all_correlations$metric_label <- sapply(all_correlations$metric, function(m) {
  if (m %in% names(metric_labels)) {
    return(metric_labels[[m]])
  } else {
    return(m)
  }
})

# Exporter les corrélations complètes uniquement
write_csv2(
  all_correlations,
  file.path(path_tables_corr, "all_metrics_correlations.csv")
)

#### 5️⃣ VISUALISATIONS DES CORRÉLATIONS ####
# Cette section crée différentes visualisations des corrélations

# Créer les graphiques en barres pour chaque variable cible et méthode
for (target_var in target_vars) {
  # Obtenir le chemin du répertoire pour cette variable
  path_target_var <- file.path(path_figures_corr, target_var)
  
  for (corr_method in corr_methods) {
    # Graphique avec noms complets
    p_full <- plot_correlation_bars(
      all_correlations, 
      target_var, 
      corr_method, 
      top_n = top_n_value, 
      use_full_names = TRUE
    )
    
    # Calculer la hauteur en fonction du nombre de métriques pour éviter les problèmes de taille
    plot_height <- min(20, 5 + (top_n_value * 0.3))
    
    # Sauvegarder dans le dossier de la variable cible
    ggsave(
      file.path(path_target_var, paste0("top", top_n_value, "_correlations_", corr_method, ".jpg")),
      p_full,
      width = 14,
      height = plot_height,
      dpi = dpi_value,
      limitsize = FALSE
    )
  }
}

# Générer les corrplots pour les variables cibles
for (target_var in target_vars) {
  # Obtenir le chemin du répertoire pour cette variable
  path_target_var <- file.path(path_figures_corr, target_var)
  
  for (corr_method in corr_methods) {
    # Créer corrplot avec les top métriques
    create_target_corrplot(
      plot_metrics,
      target_var,
      exclude_vars,
      metric_labels,
      path_target_var,
      method = corr_method,
      top_n = top_n_value
    )
  }
}

if(F){
  # Générer les corrplots pour les variables cibles
  for (target_var in target_vars) {
    # Obtenir le chemin du répertoire pour cette variable
    path_target_var <- file.path(path_figures_corr, target_var)
    
    for (corr_method in corr_methods) {
      # Créer corrplot complet
      create_target_corrplot(
        plot_metrics,
        target_var,
        exclude_vars,
        metric_labels,
        path_target_var,
        method = corr_method,
        top_n = NULL  # Toutes les métriques
      )
    }
  }
}


#### 6️⃣ ANALYSE DES RELATIONS INDIVIDUELLES ####
# Cette section génère des scatterplots pour les métriques les plus corrélées

# Générer les scatterplots pour chaque combinaison de variable cible et méthode
for (target_var in target_vars) {
  # Obtenir le chemin du répertoire pour cette variable
  path_target_var <- file.path(path_figures_corr, target_var)
  
  for (corr_method in corr_methods) {
    # Générer les visualisations améliorées
    lm_stats <- generate_enhanced_scatterplots(
      plot_metrics,
      all_correlations,
      target_var,
      metric_labels,
      plots_location,
      path_target_var,
      corr_method = corr_method,
      top_n = top_n_value,
      dpi = dpi_value
    )
  }
}

#### 7️⃣ ANALYSE COMPARATIVE DES VARIABLES CIBLES ####
# Cette section compare les corrélations entre paires de variables cibles

# Pour chaque paire de variables à comparer
for (pair in comparison_pairs) {
  var1 <- pair[1]
  var2 <- pair[2]
  
  # Pour chaque méthode de corrélation
  for (corr_method in corr_methods) {
    # Comparer les corrélations avec le nombre de métriques défini
    comparison_result <- compare_correlations(
      all_correlations,
      var1,
      var2,
      metric_labels,
      corr_method = corr_method,
      top_n = top_n_value,
      use_full_names = TRUE
    )
    
    # Si la comparaison a réussi
    if (!is.null(comparison_result)) {
      # Sauvegarder le graphique dans le dossier de corrélation principal
      ggsave(
        file.path(path_figures_corr, paste0("comparison_", var1, "_vs_", var2, "_top", top_n_value, "_", corr_method, ".jpg")),
        comparison_result$plot,
        width = 10,
        height = 8,
        dpi = dpi_value
      )
    }
  }
}

#### 8️⃣ SYNTHÈSE DES RÉSULTATS ####
# Cette section crée un résumé des principales métriques corrélées

# Générer les résumés pour chaque variable cible et méthode
for (target_var in target_vars) {
  # Obtenir le chemin du répertoire pour cette variable
  path_target_var <- file.path(path_figures_corr, target_var)
  
  for (corr_method in corr_methods) {
    # Obtenir le résumé des corrélations
    top_metrics <- summarize_top_correlations(
      all_correlations,
      target_var,
      corr_method,
      top_n = top_n_value
    )
    
    # Sauvegarder la version spécifique dans le dossier de la variable cible
    write_csv2(
      top_metrics,
      file.path(path_target_var, paste0("top", top_n_value, "_metrics_", corr_method, ".csv"))
    )
  }
}