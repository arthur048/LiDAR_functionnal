# =================================================================
# SCRIPT DE PRÉPARATION DES DONNÉES POUR L'ANALYSE DES RELATIONS
# ENTRE STRUCTURE VERTICALE FORESTIÈRE ET DENSITÉ DU BOIS
# Objectif : Charger, organiser et préparer les données pour les analyses
# =================================================================

#### 1️⃣ INITIALISATION ####
# Nettoyage de l'environnement et de la mémoire
rm(list = ls())
gc()

# Définition des packages nécessaires
pkgs <- c("tidyverse")

# Installation et chargement automatique des packages manquants
to_install <- !pkgs %in% installed.packages()
if(any(to_install)) {install.packages(pkgs[to_install])}
inst <- lapply(pkgs, library, character.only = TRUE)

# Configuration de l'environnement
options(scipen = 999)        # Éviter la notation scientifique
set.seed(42)                 # Fixer la graine pour reproductibilité

#### 2️⃣ DÉFINITION DES CHEMINS ET CRÉATION DE LA STRUCTURE ####
# Définition du répertoire de base
project_dir <- "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal"

# Création d'une fonction pour les chemins relatifs au projet
project_path <- function(...) {
  file.path(project_dir, ...)
}

# Définition des chemins essentiels
path_gaps <- project_path("3_Gaps")
path_metrics <- project_path("4_Plots_metrics")
path_main <- project_path("5_Main_analysis")

# Création des sous-répertoires pour les analyses
dirs_to_create <- c(
  file.path(path_main, "functions"),
  file.path(path_main, "input"), 
  file.path(path_main, "output", "data"),
  file.path(path_main, "output", "figures"),
  file.path(path_main, "output", "tables")
)

# Création des répertoires s'ils n'existent pas déjà
for (dir in dirs_to_create) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("Répertoire créé:", dir, "\n")
  }
}

#### 3️⃣ CHARGEMENT DES DONNÉES ####
# Charger les métriques de trouées calculées précédemment
gaps_metrics <- read_csv2(file.path(path_gaps, "output", "results_gaps_metrics.csv"))

# Charger les métriques de plots déjà calculées
plot_metrics_gaps <- read_csv2(file.path(path_metrics, "plot_metrics_gaps.csv"))

# Charger les noms des plots
plot_names <- read_csv2(file.path(project_dir, "final_plot_name.csv")) %>% 
  pull(plot_name)

# Charger les métriques d'inventaires de terrain, dont le tempérament des espèces
plot_metrics_inv <- read_csv2(file.path(path_metrics, "PlotsMetrics_FieldDataInventories.csv")) %>%
  rename_with(~ ifelse(tolower(.) == "plot_ref", "plot_name", .)) %>%
  filter(plot_name %in% plot_names)

# Charger les données de tempérament
plot_metrics_temperament <- read_csv2(file.path(path_metrics, "PlotsMetrics_Temperament_all.csv")) %>%
  filter(Type_subdivision == 3) %>%
  filter(plot_ref %in% plot_names)

# Charger les inventaires individuels avec tempérament
inventory_temperament <- read_csv2(file.path(path_metrics, "field_data_inventories_with_temperament.csv"))

#### 4️⃣ PRÉPARATION DES DONNÉES ####

# Extraction des proportions de tempérament pour la surface terrière (G.rel)
temperament_g <- plot_metrics_temperament %>%
  filter(Type_Proportion == "G.rel") %>%
  mutate(temp_type = case_when(
    Temperament == "héliophile" ~ "prop_g_helio",
    Temperament == "héliophile non pionnière" ~ "prop_g_npld",
    Temperament == "tolérante à l'ombrage" ~ "prop_g_shade",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(temp_type)) %>%
  dplyr::select(plot_ref, temp_type, Proportion) %>%
  pivot_wider(names_from = temp_type, values_from = Proportion)

# Extraction des proportions de tempérament pour le nombre d'individus (ind.rel)
temperament_ind <- plot_metrics_temperament %>%
  filter(Type_Proportion == "ind.rel") %>%
  mutate(temp_type = case_when(
    Temperament == "héliophile" ~ "prop_ind_helio",
    Temperament == "héliophile non pionnière" ~ "prop_ind_npld",
    Temperament == "tolérante à l'ombrage" ~ "prop_ind_shade",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(temp_type)) %>%
  dplyr::select(plot_ref, temp_type, Proportion) %>%
  pivot_wider(names_from = temp_type, values_from = Proportion)

# Création du dataframe plot_metrics complet
plot_metrics <- plot_metrics_gaps %>%
  left_join(temperament_g %>% rename(plot_name = plot_ref), by = "plot_name") %>%
  left_join(temperament_ind %>% rename(plot_name = plot_ref), by = "plot_name")

# Vérification des données manquantes
missing_data_summary <- plot_metrics %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "missing_count") %>%
  filter(missing_count > 0) %>%
  arrange(desc(missing_count))

# Afficher un résumé des données manquantes s'il y en a
if(nrow(missing_data_summary) > 0) {
  cat("Variables avec données manquantes:\n")
  print(missing_data_summary)
}

# Filtrer pour buffer = 0 (structure forestière interne)
buffer0_data <- gaps_metrics %>%
  filter(buffer == 0)

# Joindre avec les données de WD pour les analyses de courbes
combined_data <- buffer0_data %>%
  left_join(plot_metrics %>% dplyr::select(plot_name, WD_BA), by = "plot_name") %>%
  filter(!is.na(WD_BA))

# 1. Catégorisation par quartiles (4 groupes avec même nombre d'observations)
wd_quartiles <- quantile(plot_metrics$WD_BA, probs = seq(0, 1, 0.25), na.rm = TRUE)
combined_data <- combined_data %>%
  mutate(wd_category_quartile = cut(WD_BA, 
                                    breaks = wd_quartiles, 
                                    include.lowest = TRUE,
                                    labels = c("Q1", "Q2", "Q3", "Q4")))

# 2. Catégorisation par groupes uniformes (5 groupes de même amplitude)
wd_min <- min(plot_metrics$WD_BA, na.rm = TRUE)
wd_max <- max(plot_metrics$WD_BA, na.rm = TRUE)
wd_uniform_breaks <- seq(wd_min, wd_max, length.out = 6)
combined_data <- combined_data %>%
  mutate(wd_category_uniform = cut(WD_BA,
                                   breaks = wd_uniform_breaks,
                                   include.lowest = TRUE,
                                   labels = c("G1 (0.3-0.4)", "G2 (0.4-0.5)", "G3 (0.5-0.6)", 
                                              "G4 (0.6-0.7)", "G5 (0.7-0.8)")))

# Définir des palettes de couleurs pour les deux types de catégories
wd_palette_quartile <- viridis(4, option = "turbo", direction = -1)
names(wd_palette_quartile) <- levels(combined_data$wd_category_quartile)

wd_palette_uniform <- viridis(5, option = "turbo", direction = -1)
names(wd_palette_uniform) <- levels(combined_data$wd_category_uniform)

#### 5️⃣ EXPORT DES DONNÉES VERS LE DOSSIER INPUT ####
path_input = file.path(path_main, "input")

# Export des données vers le dossier input
write_csv2(gaps_metrics, file.path(path_input, "gaps_metrics.csv"))
write_csv2(plot_metrics_gaps, file.path(path_input, "plot_metrics_gaps.csv"))
write_csv2(plot_metrics_inv, file.path(path_input, "plot_metrics_field_inventories.csv"))
write_csv2(plot_metrics_temperament, file.path(path_input, "plot_metrics_temperament.csv"))
write_csv2(inventory_temperament, file.path(path_input, "inventory_with_temperament.csv"))
write_csv2(plot_metrics, file.path(path_input, "plot_metrics.csv"))
write_csv2(combined_data, file.path(path_input, "combined_data.csv"))

