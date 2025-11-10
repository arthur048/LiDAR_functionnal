# =================================================================
# VISUALISATION DES PROFILS DE TROUÉES EN RELATION AVEC LA DENSITÉ DU BOIS
# Objectif : Générer des visualisations avancées pour caractériser les 
#            relations entre structure verticale des forêts et densité du bois
# =================================================================

#### 1️⃣ INITIALISATION ####
# Nettoyage de l'environnement et de la mémoire
rm(list = ls())
gc()

# Définition des packages nécessaires
pkgs <- c("tidyverse", "ggpubr", "viridis", "gridExtra",
          "patchwork", "ggridges", "zoo", "fields", "ggrepel", "cowplot", "foreach", "doParallel", "here")

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

# Sous-répertoires pour les visualisations de profils
path_figures_profiles <- file.path(path_figures, "profiles")
path_data_profiles <- file.path(path_data, "profiles")

# S'assurer que les répertoires de sortie existent
dirs_to_check <- c(path_data_profiles, path_figures_profiles)
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

# Charger les fonctions de visualisation
source(file.path(path_main, "functions", "graphics_profiles.R"))

# Liste des plots
plot_names <- unique(combined_data$plot_name)

# Créer des sous-dossiers pour chaque plot
plot_folders <- list()
for (plot_name in plot_names) {
  plot_folder <- file.path(path_figures_profiles, plot_name)
  if (!dir.exists(plot_folder)) {
    dir.create(plot_folder, recursive = TRUE)
    cat("Répertoire créé:", plot_folder, "\n")
  }
  plot_folders[[plot_name]] <- plot_folder
}

# Créer les palettes de couleurs pour les deux types de catégorisation en les triant selon WD_BA
# Pour les quartiles
quartile_categories <- combined_data %>%
  group_by(wd_category_quartile) %>%
  summarise(mean_wd = mean(WD_BA, na.rm=TRUE)) %>%
  arrange(mean_wd) %>%
  pull(wd_category_quartile)

wd_palette_quartile <- viridis(length(quartile_categories), option = "viridis", direction = -1)
names(wd_palette_quartile) <- quartile_categories

# Pour les groupes uniformes
uniform_categories <- combined_data %>%
  group_by(wd_category_uniform) %>%
  summarise(mean_wd = mean(WD_BA, na.rm=TRUE)) %>%
  arrange(mean_wd) %>%
  pull(wd_category_uniform)

wd_palette_uniform <- viridis(length(uniform_categories), option = "viridis", direction = -1)
names(wd_palette_uniform) <- uniform_categories

# Paramètres graphiques
dpi_value <- 600  # Résolution des graphiques exportés
color_palette <- "viridis"  # Palette de couleurs pour les visualisations

#### 3️⃣ GÉNÉRATION DES VISUALISATIONS DE PROFILS ####

### 3A. COURBES DE PROPORTION DE TROUÉES ###
# 1. Courbes par quartile de WD
gap_curves_quartile <- plot_gap_curves_by_wd(
  combined_data %>% mutate(wd_category = wd_category_quartile),
  wd_palette_quartile,
  title = "Courbes de proportion de trouées par quartile de densité de bois"
)

# Sauvegarde du graphique
ggsave(
  file.path(path_figures_profiles, "gap_curves_by_wd_quartile.jpg"),
  gap_curves_quartile,
  width = 12,
  height = 8,
  dpi = dpi_value
)

# 2. Courbes par groupe uniforme de WD
gap_curves_uniform <- plot_gap_curves_by_wd(
  combined_data %>% mutate(wd_category = wd_category_uniform),
  wd_palette_uniform,
  title = "Courbes de proportion de trouées par groupe uniforme de densité de bois"
)

# Sauvegarde du graphique
ggsave(
  file.path(path_figures_profiles, "gap_curves_by_wd_uniform.jpg"),
  gap_curves_uniform,
  width = 12,
  height = 8,
  dpi = dpi_value
)

### 3B. HEATMAP DES PROFILS ###
# Heatmap ordonnée par WD
heatmap_plot <- create_profiles_heatmap(
  combined_data,
  plot_names = unique(combined_data$plot_name),
  color_palette = "plasma",
  title = "Heatmap des profils de trouées ordonnés par densité du bois croissante"
)

# Sauvegarde du graphique
ggsave(
  file.path(path_figures_profiles, "profiles_heatmap.jpg"),
  heatmap_plot,
  width = 16,
  height = 8,
  dpi = dpi_value
)

### 3C. PROFILS TYPES ###
# 1. Profils types par quartile - avec écart-type
typical_profiles_quartile_with_sd <- create_typical_profiles(
  combined_data,
  wd_palette_quartile,
  category_column = "wd_category_quartile",
  show_deviation = TRUE,
  title = "Profils types par quartile de densité du bois (avec écart-type)"
)

# Sauvegarde du graphique
ggsave(
  file.path(path_figures_profiles, "typical_profiles_quartile_with_sd.jpg"),
  typical_profiles_quartile_with_sd,
  width = 12,
  height = 8,
  dpi = dpi_value
)

# 2. Profils types par quartile - sans écart-type
typical_profiles_quartile_no_sd <- create_typical_profiles(
  combined_data,
  wd_palette_quartile,
  category_column = "wd_category_quartile",
  show_deviation = FALSE,
  title = "Profils types par quartile de densité du bois"
)

# Sauvegarde du graphique
ggsave(
  file.path(path_figures_profiles, "typical_profiles_quartile_no_sd.jpg"),
  typical_profiles_quartile_no_sd,
  width = 12,
  height = 8,
  dpi = dpi_value
)

# 3. Profils types par groupe uniforme - avec écart-type
typical_profiles_uniform_with_sd <- create_typical_profiles(
  combined_data,
  wd_palette_uniform,
  category_column = "wd_category_uniform",
  show_deviation = TRUE,
  title = "Profils types par groupe uniforme de densité du bois (avec écart-type)"
)

# Sauvegarde du graphique
ggsave(
  file.path(path_figures_profiles, "typical_profiles_uniform_with_sd.jpg"),
  typical_profiles_uniform_with_sd,
  width = 12,
  height = 8,
  dpi = dpi_value
)

# 4. Profils types par groupe uniforme - sans écart-type
typical_profiles_uniform_no_sd <- create_typical_profiles(
  combined_data,
  wd_palette_uniform,
  category_column = "wd_category_uniform",
  show_deviation = FALSE,
  title = "Profils types par groupe uniforme de densité du bois"
)

# Sauvegarde du graphique
ggsave(
  file.path(path_figures_profiles, "typical_profiles_uniform_no_sd.jpg"),
  typical_profiles_uniform_no_sd,
  width = 12,
  height = 8,
  dpi = dpi_value
)

### 3D. SIGNATURES DIFFÉRENTIELLES ###
# 1. Signatures différentielles par quartile avec écart-type
differential_signatures_quartile_with_sd <- create_differential_signatures(
  combined_data %>% mutate(wd_category = wd_category_quartile),
  wd_palette_quartile,
  title = "Signatures différentielles par quartile de densité du bois (avec écart-type)",
  confidence_interval = "sd",
  show_deviation = TRUE
)

ggsave(
  file.path(path_figures_profiles, "differential_signatures_quartile_with_sd.jpg"),
  differential_signatures_quartile_with_sd,
  width = 12,
  height = 8,
  dpi = dpi_value
)

# 2. Signatures différentielles par quartile sans écart-type
differential_signatures_quartile_no_sd <- create_differential_signatures(
  combined_data %>% mutate(wd_category = wd_category_quartile),
  wd_palette_quartile,
  title = "Signatures différentielles par quartile de densité du bois (sans écart-type)",
  confidence_interval = "sd",
  show_deviation = FALSE
)

ggsave(
  file.path(path_figures_profiles, "differential_signatures_quartile_no_sd.jpg"),
  differential_signatures_quartile_no_sd,
  width = 12,
  height = 8,
  dpi = dpi_value
)

# 3. Signatures différentielles par groupe uniforme avec écart-type
differential_signatures_uniform_with_sd <- create_differential_signatures(
  combined_data %>% mutate(wd_category = wd_category_uniform),
  wd_palette_uniform,
  title = "Signatures différentielles par groupe uniforme de densité du bois (avec écart-type)",
  confidence_interval = "sd",
  show_deviation = TRUE
)

ggsave(
  file.path(path_figures_profiles, "differential_signatures_uniform_with_sd.jpg"),
  differential_signatures_uniform_with_sd,
  width = 12,
  height = 8,
  dpi = dpi_value
)

# 4. Signatures différentielles par groupe uniforme sans écart-type
differential_signatures_uniform_no_sd <- create_differential_signatures(
  combined_data %>% mutate(wd_category = wd_category_uniform),
  wd_palette_uniform,
  title = "Signatures différentielles par groupe uniforme de densité du bois (sans écart-type)",
  confidence_interval = "sd",
  show_deviation = FALSE
)

ggsave(
  file.path(path_figures_profiles, "differential_signatures_uniform_no_sd.jpg"),
  differential_signatures_uniform_no_sd,
  width = 12,
  height = 8,
  dpi = dpi_value
)


#### 4️⃣ VISUALISATIONS COMBINÉES ####
# Combinaison d'un même type de graphique avec différentes catégorisations

# A. Profils types comparés (quartiles vs groupes uniformes)
profiles_comparison <- typical_profiles_quartile_no_sd + typical_profiles_uniform_no_sd +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Comparaison des profils types: quartiles vs groupes uniformes",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

# Sauvegarde de la visualisation combinée
ggsave(
  file.path(path_figures_profiles, "profiles_comparison.jpg"),
  profiles_comparison,
  width = 16,
  height = 8,
  dpi = dpi_value
)

# B. Signatures différentielles comparées (quartiles vs groupes uniformes)
signatures_comparison <- differential_signatures_quartile_no_sd + differential_signatures_uniform_no_sd +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Comparaison des signatures différentielles: quartiles vs groupes uniformes",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

# Sauvegarde de la visualisation combinée
ggsave(
  file.path(path_figures_profiles, "signatures_comparison.jpg"),
  signatures_comparison,
  width = 16,
  height = 8,
  dpi = dpi_value
)

#### 5️⃣ TABLEAUX DE BORD COMPLETS ####
# Regroupement par type de catégorisation pour des tableaux de bord complets

# A. Tableau de bord pour les quartiles
dashboard_quartiles <- (gap_curves_quartile + typical_profiles_quartile_no_sd) /
  (heatmap_plot + differential_signatures_quartile_no_sd) +
  plot_annotation(
    title = "Analyse complète des profils de trouées par quartile de densité du bois",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

# Sauvegarde du tableau de bord
ggsave(
  file.path(path_figures_profiles, "dashboard_quartiles.jpg"),
  dashboard_quartiles,
  width = 20,
  height = 16,
  dpi = dpi_value
)

# B. Tableau de bord pour les groupes uniformes
dashboard_uniform <- (gap_curves_uniform + typical_profiles_uniform_no_sd) /
  (heatmap_plot + differential_signatures_uniform_no_sd) +
  plot_annotation(
    title = "Analyse complète des profils de trouées par groupe uniforme de densité du bois",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

# Sauvegarde du tableau de bord
ggsave(
  file.path(path_figures_profiles, "dashboard_uniform.jpg"),
  dashboard_uniform,
  width = 20,
  height = 16,
  dpi = dpi_value
)

#### 6️⃣ ANALYSES AVANCÉES PAR PLOT ####
# Création de visualisations avancées pour chaque plot individuel

for(plot_i in plot_names) {
  # Données pour ce plot
  plot_data <- combined_data %>%
    filter(plot_name == plot_i) %>%
    arrange(height_aboveground)
  
  # Vérifier qu'il y a assez de données
  if (nrow(plot_data) < 5) {
    return(NULL)  # Passer au plot suivant si pas assez de données
  }
  
  # Récupérer les valeurs de hauteur, proportion et WD
  height <- plot_data$height_aboveground
  prop <- plot_data$proportion
  wd_value <- plot_data$WD_BA[1]
  
  # Dossier de sortie pour ce plot
  output_folder <- plot_folders[[plot_i]]
  
  # 1. Décomposition fonctionnelle
  zonal_plot <- plot_decompose_curve_segments(height, prop, plot_i, wd_value)
  
  # Sauvegarde du graphique
  ggsave(
    file.path(output_folder, "zonal_decomposition.jpg"),
    zonal_plot,
    width = 10,
    height = 6,
    dpi = dpi_value
  )
  
  # 2. Analyse différentielle
  deriv_plots <- plot_derivatives(height, prop, plot_i, wd_value)
  
  # Sauvegarder les graphiques
  ggsave(
    file.path(output_folder, "derivatives_curve.jpg"),
    deriv_plots$p1,
    width = 10,
    height = 6,
    dpi = dpi_value
  )
  
  ggsave(
    file.path(output_folder, "derivatives_only.jpg"),
    deriv_plots$p2,
    width = 10,
    height = 6,
    dpi = dpi_value
  )
  
  # 3. Analyse multi-échelle
  multiscale_plots <- plot_multiscale_analysis(height, prop, plot_i, wd_value, 
                                              scales = c(5, 10, 15, 20))
  
  # Sauvegarder les graphiques
  ggsave(
    file.path(output_folder, "multiscale_curves.jpg"),
    multiscale_plots$p1,
    width = 10,
    height = 6,
    dpi = dpi_value
  )
  
  ggsave(
    file.path(output_folder, "multiscale_residuals.jpg"),
    multiscale_plots$p2,
    width = 10,
    height = 6,
    dpi = dpi_value
  )
  
  # 4. Création du dashboard avec titres améliorés
  dashboard <- (zonal_plot + deriv_plots$p1 + deriv_plots$p2) /
               (multiscale_plots$p1 + multiscale_plots$p2) +
    plot_annotation(
      title = paste0("Analyse avancée des profils de trouées: ", plot_i, " (WD = ", round(wd_value, 3), ")"),
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  # Sauvegarder le dashboard
  ggsave(
    file.path(output_folder, "dashboard_advanced.jpg"),
    dashboard,
    width = 20,
    height = 12,
    dpi = dpi_value
  )
}
 
