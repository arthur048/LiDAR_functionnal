#### 1️⃣ INITIALIZATION ####
# Clean environment and memory
rm(list = ls())
gc()

# Définition des packages nécessaires
pkgs <- c("tidyverse", "FactoMineR", "factoextra", "ggrepel", "viridis", 
          "corrplot", "patchwork", "ade4", "vegan", "broom", "here")

# Automatic installation and loading of missing packages
to_install <- !pkgs %in% installed.packages()
if(any(to_install)) {install.packages(pkgs[to_install])}
inst <- lapply(pkgs, library, character.only = TRUE)

#### 2️⃣ PATH DEFINITION ####
# Helper function to build paths (adapt to your environment)
project_path <- function(...) {
  here(...)
}

# Base directory definition
project_dir <- here()

# Essential paths
path_main <- project_path("5_Main_analysis")
path_input <- file.path(path_main, "input")
path_output <- file.path(path_main, "output")

# Paths for this specific analysis
path_data <- file.path(path_output, "data/multivariate_suite")
path_figures <- file.path(path_output, "figures/multivariate_suite")
path_table <- file.path(path_output, "tables/multivariate_suite")

# Folder creation
dir.create(path_data, recursive = TRUE, showWarnings = FALSE)
dir.create(path_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(path_table, recursive = TRUE, showWarnings = FALSE)

#### 3️⃣ DATA LOADING ####
# Load main data
plot_metrics <- read_csv2(file.path(path_input, "plot_metrics.csv"))
inventory_temperament <- read_csv2(file.path(path_input, "inventory_with_temperament.csv"))

site_score_nsca_species <- read_csv(file.path(path_main, "output/data/multivariate/nsca_species_site_coordinates.csv"))
site_score_nsca_species_dbh <- read_csv(file.path(path_main, "output/data/multivariate/nsca_species_dbh_site_coordinates.csv"))
site_score_ca_species <- read_csv(file.path(path_main, "output/data/multivariate/ca_species_site_coordinates.csv"))
site_score_nsca_genus <- read_csv(file.path(path_main, "output/data/multivariate/nsca_genus_site_coordinates.csv"))
site_score_nsca_genus_dbh <- read_csv(file.path(path_main, "output/data/multivariate/nsca_genus_dbh_site_coordinates.csv"))

#### 4️⃣ DATA PREPARATION ####
# Filter plots with temperament data
plots_with_temperament <- inventory_temperament %>%
  pull(plot_ref) %>%
  unique()

# Filter plot_metrics to keep only these plots
plot_data_complete <- plot_metrics %>%
  filter(plot_name %in% plots_with_temperament)

# Define columns for analyses
metadata_cols <- c("plot_name", "WD_BA", "prop_g_helio", "prop_g_npld", "prop_g_shade", 
                   "prop_ind_helio", "prop_ind_npld", "prop_ind_shade")

# Canopy structural metrics
gap_metrics_cols <- names(plot_data_complete)[!names(plot_data_complete) %in% metadata_cols]

# Metadata with wood density and succession
metadata <- plot_data_complete %>%
  select(all_of(metadata_cols))

plots_name <- metadata %>%
  pull(plot_name)

#### 5️⃣ CONSTRUCTION DE LA TABLE FINALE ####

# Préparer les données des analyses multivariées (2 premiers axes)
nsca_species_data <- site_score_nsca_species %>%
  filter(plot_ref %in% plots_with_temperament) %>%
  select(plot_ref, NSCA1:NSCA2) %>%
  rename(plot_name = plot_ref)

nsca_species_dbh_data <- site_score_nsca_species_dbh %>%
  filter(plot_ref %in% plots_with_temperament) %>%
  select(plot_ref, NSCA_DBH1:NSCA_DBH2) %>%
  rename(plot_name = plot_ref)

ca_species_data <- site_score_ca_species %>%
  filter(plot_ref %in% plots_with_temperament) %>%
  select(plot_ref, CA1:CA2) %>%
  rename(plot_name = plot_ref)

# nsca_genus_data <- site_score_nsca_genus %>%
#   filter(plot_ref %in% plots_with_temperament) %>%
#   select(plot_ref, NSCA_Genus1:NSCA_Genus2) %>%
#   rename(plot_name = plot_ref)

nsca_genus_dbh_data <- site_score_nsca_genus_dbh %>%
  filter(plot_ref %in% plots_with_temperament) %>%
  select(plot_ref, NSCA_Genus_DBH1:NSCA_Genus_DBH2) %>%
  rename(plot_name = plot_ref)

# Sélectionner les colonnes gap_metrics et les variables cibles de plot_data_complete
structural_and_target_data <- plot_data_complete %>%
  select(plot_name, all_of(gap_metrics_cols), WD_BA, prop_g_helio, prop_g_npld, prop_g_shade) %>%
  # Remplacer la valeur NA dans h100 par 60
  mutate(h100 = replace_na(h100, 60))  # Méthode tidyverse recommandée
# Ou alternative : mutate(h100 = ifelse(is.na(h100), 60, h100))

# Combiner toutes les données
pca_data_final <- structural_and_target_data %>%
  left_join(nsca_species_data, by = "plot_name") %>%
  left_join(nsca_species_dbh_data, by = "plot_name") %>%
  left_join(ca_species_data, by = "plot_name") %>%
  # left_join(nsca_genus_data, by = "plot_name") %>%
  left_join(nsca_genus_dbh_data, by = "plot_name") %>%
  column_to_rownames("plot_name")

pca_data_clean <- pca_data_final %>%
  drop_na()

#### 6️⃣ DIAGNOSTIC DE COLINÉARITÉ ET FILTRAGE DES VARIABLES ####

# Extraire uniquement les variables structurelles pour le diagnostic
structure_data_raw <- pca_data_clean %>%
  select(all_of(intersect(gap_metrics_cols, names(pca_data_clean))))

# =============================================================================
# A.1. DIAGNOSTIC DE COLINÉARITÉ ET FILTRAGE DES VARIABLES
# =============================================================================

# Diagnostic initial de colinéarité
cat("🔍 Diagnostic de colinéarité des variables de structure:\n")

# Matrice de corrélation
cor_matrix_raw <- cor(structure_data_raw, use = "complete.obs")

# Identifier les paires hautement corrélées (|r| > 0.95)
high_cor_pairs <- which(abs(cor_matrix_raw) > 0.95 & abs(cor_matrix_raw) < 1, arr.ind = TRUE)
if(nrow(high_cor_pairs) > 0) {
  cat("⚠️  Variables hautement corrélées (|r| > 0.95):\n")
  for(i in 1:nrow(high_cor_pairs)) {
    var1 <- rownames(cor_matrix_raw)[high_cor_pairs[i,1]]
    var2 <- colnames(cor_matrix_raw)[high_cor_pairs[i,2]]
    cor_val <- cor_matrix_raw[high_cor_pairs[i,1], high_cor_pairs[i,2]]
    cat("  ", var1, "↔", var2, ":", round(cor_val, 3), "\n")
  }
} else {
  cat("✅ Aucune corrélation > 0.95 détectée\n")
}

# Identifier et supprimer les versions non normalisées
# Prioriser les versions normalisées (_norm) et supprimer les versions brutes
vars_all <- colnames(structure_data_raw)
vars_norm <- vars_all[grepl("_norm$", vars_all)]
vars_base <- gsub("_norm$", "", vars_norm)

# Variables à supprimer (versions non normalisées quand version normalisée existe)
vars_to_remove <- vars_base[vars_base %in% vars_all]

# Variables supplémentaires identifiées comme redondantes
additional_removes <- c("upper_thickness", "basal_thickness", "transition_auc_norm", 
                        "basal_auc_norm", "upper_auc_norm")  # Ajuster selon vos données

# Liste finale des variables à supprimer
all_vars_to_remove <- unique(c(vars_to_remove, additional_removes))
all_vars_to_remove <- all_vars_to_remove[all_vars_to_remove %in% vars_all]

cat("\n🗑️  Variables supprimées (redondantes):\n")
if(length(all_vars_to_remove) > 0) {
  cat("  ", paste(all_vars_to_remove, collapse = ", "), "\n")
} else {
  cat("  Aucune variable supprimée\n")
}

# Filtrer les données structurelles
structure_data_filtered <- structure_data_raw[, !colnames(structure_data_raw) %in% all_vars_to_remove]
cat("📊 Variables restantes:", ncol(structure_data_filtered), "sur", ncol(structure_data_raw), "\n")

# Vérification finale de la colinéarité
cor_matrix_filtered <- cor(structure_data_filtered, use = "complete.obs")
high_cor_final <- which(abs(cor_matrix_filtered) > 0.95 & abs(cor_matrix_filtered) < 1, arr.ind = TRUE)

if(nrow(high_cor_final) > 0) {
  cat("⚠️  Corrélations restantes > 0.95:\n")
  for(i in 1:nrow(high_cor_final)) {
    var1 <- rownames(cor_matrix_filtered)[high_cor_final[i,1]]
    var2 <- colnames(cor_matrix_filtered)[high_cor_final[i,2]]
    cor_val <- cor_matrix_filtered[high_cor_final[i,1], high_cor_final[i,2]]
    cat("  ", var1, "↔", var2, ":", round(cor_val, 3), "\n")
  }
} else {
  cat("✅ Aucune colinéarité résiduelle > 0.95\n")
}

# CHOIX INTERACTIF DES VARIABLES À SUPPRIMER SUPPLÉMENTAIRES
cat("\n🎯 Variables structurelles finales disponibles:\n")
cat("  ", paste(sort(colnames(structure_data_filtered)), collapse = ", "), "\n")

# Variables à supprimer manuellement (modifiez cette liste selon vos besoins)
manual_removes <- c("height_max_d2", "height_min_d2", "basal_auc_norm", "upper_auc_norm") 

# Si vous voulez supprimer des variables supplémentaires, modifiez la ligne ci-dessus
if(length(manual_removes) > 0) {
  manual_removes <- manual_removes[manual_removes %in% colnames(structure_data_filtered)]
  structure_data_filtered <- structure_data_filtered[, !colnames(structure_data_filtered) %in% manual_removes]
  cat("🔧 Variables supprimées manuellement:", paste(manual_removes, collapse = ", "), "\n")
}

# Mettre à jour la liste des variables structurelles filtrées
gap_metrics_filtered <- colnames(structure_data_filtered)

# Reconstruire le dataset final avec les variables filtrées
pca_data_final_filtered <- pca_data_clean %>%
  select(all_of(gap_metrics_filtered), WD_BA, prop_g_helio, prop_g_npld, prop_g_shade,
         NSCA_Genus_DBH1, NSCA_Genus_DBH2, CA1, CA2)

cat("📋 Dataset final: ", ncol(pca_data_final_filtered), " variables pour ", nrow(pca_data_final_filtered), " plots\n")

#### 7️⃣ DÉFINITION DES VARIABLES ACTIVES ET SUPPLÉMENTAIRES ####

# Variables structurelles filtrées - VARIABLES ACTIVES pour la PCA
structural_vars_filtered <- gap_metrics_filtered

# Variables supplémentaires (spécifiées par l'utilisateur)
supplementary_vars <- c("WD_BA", "prop_g_helio", "prop_g_npld", "prop_g_shade",
                        "NSCA_Genus_DBH1", "NSCA_Genus_DBH2",
                        "CA1", "CA2")

# Vérifier quelles variables sont effectivement présentes
supplementary_vars_present <- intersect(supplementary_vars, names(pca_data_final_filtered))

# Créer un vecteur d'indices pour les variables supplémentaires
supp_indices <- which(names(pca_data_final_filtered) %in% supplementary_vars_present)

#### 8️⃣ ANALYSE EN COMPOSANTES PRINCIPALES ####

# Réaliser la PCA avec variables supplémentaires
pca_result <- PCA(pca_data_final_filtered, 
                  scale.unit = TRUE, 
                  quanti.sup = supp_indices,  # Variables supplémentaires
                  graph = FALSE)

# Résumé de la PCA
print(summary(pca_result))

# Extraire les résultats
pca_individuals <- get_pca_ind(pca_result)
pca_variables <- get_pca_var(pca_result)

# Extraire les variables supplémentaires si elles existent
if(length(supp_indices) > 0) {
  pca_supp_quanti <- pca_result$quanti.sup
} else {
  pca_supp_quanti <- NULL
}

#### 9️⃣ VISUALISATIONS ####

# Préparer les données pour WD_BA (coloration des individus)
ind_coords_df <- as.data.frame(pca_individuals$coord[, 1:2]) %>%
  rownames_to_column("plot_name") %>%
  left_join(
    pca_data_final_filtered %>% 
      rownames_to_column("plot_name") %>% 
      select(plot_name, WD_BA), 
    by = "plot_name"
  )

# Créer un dataframe pour les variables ACTIVES
var_contrib_df <- as.data.frame(pca_variables$contrib[, 1:2]) %>%
  rownames_to_column("variable") %>%
  mutate(type = "Active") %>%
  rename(Dim1_contrib = Dim.1, Dim2_contrib = Dim.2)

# Ajouter les coordonnées des variables actives
var_coord_df <- as.data.frame(pca_variables$coord[, 1:2]) %>%
  rownames_to_column("variable") %>%
  rename(Dim1_coord = Dim.1, Dim2_coord = Dim.2)

# Combiner variables actives
plot_vars_df <- var_contrib_df %>%
  left_join(var_coord_df, by = "variable")

# Ajouter les variables supplémentaires si elles existent
if(!is.null(pca_supp_quanti)) {
  supp_coord_df <- as.data.frame(pca_supp_quanti$coord[, 1:2]) %>%
    rownames_to_column("variable") %>%
    mutate(
      type = "Supplementary",
      Dim1_contrib = NA,  # Les variables supplémentaires ne contribuent pas
      Dim2_contrib = NA
    ) %>%
    rename(Dim1_coord = Dim.1, Dim2_coord = Dim.2)
  
  # Combiner avec les variables actives
  plot_vars_df <- bind_rows(plot_vars_df, supp_coord_df)
}

# 1. Graphique des individus coloré par WD_BA
p_individuals <- ind_coords_df %>%
  ggplot(aes(x = Dim.1, y = Dim.2, color = WD_BA)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = plot_name), size = 3, max.overlaps = 15, 
                  color = "black", alpha = 0.7) +
  scale_color_viridis_c(name = "WD_BA", option = "viridis") +
  labs(
    title = "PCA - Individuals colored by Wood Density (WD_BA)",
    subtitle = paste0("Based on canopy structural metrics only - Dim 1 (", 
                      round(pca_result$eig[1, 2], 1), "%) vs Dim 2 (", 
                      round(pca_result$eig[2, 2], 1), "%)"),
    x = paste0("Dim 1 (", round(pca_result$eig[1, 2], 1), "%)"),
    y = paste0("Dim 2 (", round(pca_result$eig[2, 2], 1), "%)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  ) +
  coord_equal()

# 2. Graphique des variables (actives et supplémentaires)
p_variables <- plot_vars_df %>%
  ggplot(aes(x = Dim1_coord, y = Dim2_coord, color = type)) +
  geom_segment(aes(xend = 0, yend = 0), 
               arrow = arrow(length = unit(0.3, "cm")), alpha = 0.7) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = variable), size = 3, max.overlaps = 25) +
  scale_color_manual(
    values = c("Active" = "#2E8B57", "Supplementary" = "#CD853F"),
    name = "Variable Type",
    labels = c("Active (Structural)", "Supplementary")
  ) +
  labs(
    title = "PCA - Variables Factor Map",
    subtitle = "Active variables (structural) vs Supplementary variables",
    x = paste0("Dim 1 (", round(pca_result$eig[1, 2], 1), "%)"),
    y = paste0("Dim 2 (", round(pca_result$eig[2, 2], 1), "%)")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  coord_equal()

# 3. Contribution des variables ACTIVES aux axes
active_vars_df <- plot_vars_df %>% filter(type == "Active")

p_contrib_dim1 <- active_vars_df %>%
  arrange(desc(Dim1_contrib)) %>%
  slice_head(n = 15) %>%
  ggplot(aes(x = reorder(variable, Dim1_contrib), y = Dim1_contrib)) +
  geom_col(fill = "#2E8B57", alpha = 0.8) +
  coord_flip() +
  labs(
    title = "Contribution to Dimension 1 (Active Variables Only)",
    x = "Structural Variables",
    y = "Contribution (%)"
  ) +
  theme_minimal()

p_contrib_dim2 <- active_vars_df %>%
  arrange(desc(Dim2_contrib)) %>%
  slice_head(n = 15) %>%
  ggplot(aes(x = reorder(variable, Dim2_contrib), y = Dim2_contrib)) +
  geom_col(fill = "#2E8B57", alpha = 0.8) +
  coord_flip() +
  labs(
    title = "Contribution to Dimension 2 (Active Variables Only)",
    x = "Structural Variables",
    y = "Contribution (%)"
  ) +
  theme_minimal()

# 4. Graphique spécial pour les variables supplémentaires
if(!is.null(pca_supp_quanti)) {
  supp_vars_df <- plot_vars_df %>% filter(type == "Supplementary")
  
  p_supplementary <- supp_vars_df %>%
    ggplot(aes(x = Dim1_coord, y = Dim2_coord)) +
    geom_segment(aes(xend = 0, yend = 0), 
                 arrow = arrow(length = unit(0.3, "cm")), 
                 color = "#CD853F", alpha = 0.8, size = 1) +
    geom_point(size = 3, color = "#CD853F") +
    geom_text_repel(aes(label = variable), size = 3.5, color = "#8B4513") +
    labs(
      title = "Supplementary Variables in PCA Space",
      subtitle = "Projection of WD_BA, succession groups, and multivariate axes",
      x = paste0("Dim 1 (", round(pca_result$eig[1, 2], 1), "%)"),
      y = paste0("Dim 2 (", round(pca_result$eig[2, 2], 1), "%)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    coord_equal()
} else {
  p_supplementary <- ggplot() + 
    geom_text(aes(x = 0, y = 0, label = "No supplementary variables found")) +
    theme_void()
}

# 4. Scree plot
scree_data <- as.data.frame(pca_result$eig) %>%
  rownames_to_column("dimension") %>%
  mutate(dimension = factor(dimension, levels = dimension))

p_scree <- scree_data %>%
  ggplot(aes(x = dimension, y = `percentage of variance`)) +
  geom_col(fill = "#4472C4", alpha = 0.7) +
  geom_line(aes(group = 1), color = "#E74C3C", size = 1) +
  geom_point(color = "#E74C3C", size = 2) +
  labs(
    title = "Scree Plot",
    x = "Dimensions",
    y = "Percentage of Variance (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Combiner les graphiques principaux
if(!is.null(pca_supp_quanti)) {
  combined_plot <- (p_individuals | p_variables) / (p_contrib_dim1 | p_contrib_dim2)
  supplementary_plot <- p_supplementary
} else {
  combined_plot <- (p_individuals | p_variables) / (p_contrib_dim1 | p_contrib_dim2)
  supplementary_plot <- NULL
}

# Afficher les graphiques
print(combined_plot)
print(p_scree)
if(!is.null(supplementary_plot)) {
  print(supplementary_plot)
}

#### 🔟 SAUVEGARDE DES RÉSULTATS ####

# Sauvegarder les données filtrées
write_csv(pca_data_final_filtered %>% rownames_to_column("plot_name"), 
          file.path(path_data, "pca_final_dataset_filtered.csv"))

# Sauvegarder les résultats du diagnostic de colinéarité
write_csv(data.frame(
  variable_removed = all_vars_to_remove,
  reason = "Colinearity or redundancy"
), file.path(path_data, "variables_removed_colinearity.csv"))

# Sauvegarder les résultats PCA
write_csv(ind_coords_df, file.path(path_data, "pca_individuals_coordinates.csv"))
write_csv(plot_vars_df, file.path(path_data, "pca_variables_results.csv"))
write_csv(scree_data, file.path(path_data, "pca_eigenvalues.csv"))

# Sauvegarder informations sur les variables
variable_info <- data.frame(
  variable = names(pca_data_final_filtered),
  type = ifelse(names(pca_data_final_filtered) %in% supplementary_vars_present, 
                "Supplementary", "Active"),
  category = case_when(
    names(pca_data_final_filtered) %in% structural_vars_filtered ~ "Structural_metrics",
    names(pca_data_final_filtered) %in% c("WD_BA", "prop_g_helio", "prop_g_npld", "prop_g_shade") ~ "Wood_density_succession",
    str_detect(names(pca_data_final_filtered), "NSCA|CA") ~ "Multivariate_axes",
    TRUE ~ "Other"
  )
)
write_csv(variable_info, file.path(path_data, "pca_variable_classification.csv"))

# Sauvegarder les graphiques
ggsave(file.path(path_figures, "pca_combined_plot.png"), 
       combined_plot, width = 16, height = 12, dpi = 300)
ggsave(file.path(path_figures, "pca_scree_plot.png"), 
       p_scree, width = 10, height = 6, dpi = 300)
ggsave(file.path(path_figures, "pca_individuals_WD_BA.png"), 
       p_individuals, width = 12, height = 10, dpi = 300)

if(!is.null(supplementary_plot)) {
  ggsave(file.path(path_figures, "pca_supplementary_variables.png"), 
         supplementary_plot, width = 12, height = 10, dpi = 300)
}

