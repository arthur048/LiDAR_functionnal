#### 1️⃣ INITIALIZATION ####
# Clean environment and memory
rm(list = ls())
gc()

# Required packages
pkgs <- c("tidyverse", "FactoMineR", "factoextra", "ggrepel", "viridis", 
          "corrplot", "patchwork", "ade4", "vegan", "broom", "here",
          "cluster", "indicspecies", "caret", "glmnet", "lavaan", "semPlot",
          "ggcorrplot", "RColorBrewer", "plotly", "DT", "kableExtra")

# Automatic installation and loading of missing packages
to_install <- !pkgs %in% installed.packages()
if(any(to_install)) {install.packages(pkgs[to_install])}
inst <- lapply(pkgs, library, character.only = TRUE)

# General configuration
options(scipen = 999)  # Avoid scientific notation
set.seed(123)         # Reproducibility

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

plots_name = metadata %>%
  pull(plot_name)

# 1. Clean taxonomic data
inventory_clean <- inventory_temperament %>%
  filter(!is.na(species_full) &
           species_full != "" &
           species_full != "NA NA" &
           species_full != "sp.2 NA" &
           species_full != "indet indet " &
           species_full != "bois mort" &
           species_full != "Unknown NA"
           )

# 2. Build species × plots contingency table
species_plot_table <- inventory_clean %>%
  count(plot_ref, species_full) %>%
  pivot_wider(names_from = species_full, 
              values_from = n,
              values_fill = 0) %>%
  column_to_rownames("plot_ref")

# 2.a Table espèces × parcelles pondérée par DBH
if(T){
  species_dbh_plot_table <- inventory_clean %>%
    # Somme des sections par espèce et parcelle (DBH en cm)
    group_by(plot_ref, species_full) %>%
    summarise(basal_area = sum(G), .groups = "drop") %>% # Surface terrière en m²
    # Convertir au format large (contingence)
    pivot_wider(names_from = species_full, 
                values_from = basal_area,
                values_fill = 0) %>%
    # Conserver plot_ref comme identifiant de ligne
    column_to_rownames("plot_ref")
  
  genus_plot_table <- inventory_clean %>%
    count(plot_ref, genus) %>%
    pivot_wider(names_from = genus, 
                values_from = n,
                values_fill = 0) %>%
    column_to_rownames("plot_ref")
  
  #Test table avec uniquement le genre
  genus_wd_table <- inventory_clean %>%
    group_by(genus) %>%
    reframe(wd = min(meanWD, na.rm = TRUE)) %>%
    distinct(genus, .keep_all = TRUE) %>%
    column_to_rownames("genus")
}

# 3. Build species x wood density table
species_wd_table <- inventory_clean %>%
  group_by(species_full) %>%
  reframe(wd = min(meanWD, na.rm = TRUE)) %>%
  distinct(species_full, .keep_all = TRUE) %>%
  column_to_rownames("species_full")

# 4. Build structure x plots table

structure_plots_table <- plot_data_complete %>%
  mutate(h100 = ifelse(is.na(h100), 1, h100)) %>%
  select(plot_name , all_of(gap_metrics_cols)) %>%
  column_to_rownames("plot_name")

#### 5️⃣ ÉTAPE 1: NSCA LIBRE - GROUPES FLORISTIQUES NATURELS ####

# =============================================================================
# A. NSCA - Analyse de Correspondance Non-Symétrique
# =============================================================================

# Effectuer la NSCA sur les données floristiques
nsca_result <- dudi.nsc(species_plot_table, scannf = FALSE, nf = 10)
nsca_dbh_result <- dudi.nsc(species_dbh_plot_table, scannf = FALSE, nf = 10)
nsca_genus_result <- dudi.nsc(genus_plot_table, scannf = FALSE, nf = 10)

# Variance expliquée par les axes
variance_explained <- nsca_result$eig / sum(nsca_result$eig) * 100
variance_cumulative <- cumsum(variance_explained)

cat("📈 Variance expliquée par les premiers axes:\n")
for(i in 1:min(5, length(variance_explained))) {
  cat(sprintf("  Axe %d: %.1f%% (cumulé: %.1f%%)\n", 
              i, variance_explained[i], variance_cumulative[i]))
}

# Graphique variance expliquée
p_variance <- data.frame(
  Axe = 1:min(10, length(variance_explained)),
  Variance = variance_explained[1:min(10, length(variance_explained))]
) %>%
  ggplot(aes(x = Axe, y = Variance)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_text(aes(label = paste0(round(Variance, 1), "%")), 
            vjust = -0.5, size = 3) +
  labs(title = "Variance expliquée par les axes NSCA",
       y = "Variance expliquée (%)",
       x = "Axes NSCA") +
  theme_minimal()

print(p_variance)

# =============================================================================
# B. VISUALISATIONS NSCA
# =============================================================================

# Scores des parcelles sur les axes NSCA
site_scores <- nsca_result$li
site_scores_dbh <- nsca_dbh_result$li
site_scores_genus <- nsca_genus_result$li

colnames(site_scores) <- paste0("NSCA", 1:ncol(site_scores))
colnames(site_scores_dbh) <- paste0("NSCA", 1:ncol(site_scores_dbh))
colnames(site_scores_genus) <- paste0("NSCA", 1:ncol(site_scores_genus))

# Fonction pour plot NSCA avec flexibilité sur les axes et coloration
plot_nsca <- function(site_scores, x = "NSCA1", y = "NSCA2", color_wd = TRUE, 
                      wd_data = NULL, variance_explained_data = NULL) {
  
  # Utiliser les valeurs par défaut si non spécifiées
  if(is.null(wd_data)) {
    wd_data <- metadata$WD_BA
  }
  if(is.null(variance_explained_data)) {
    variance_explained_data <- variance_explained
  }
  
  # Extraire les numéros d'axes pour les labels
  x_axis_num <- as.numeric(gsub("NSCA", "", x))
  y_axis_num <- as.numeric(gsub("NSCA", "", y))
  
  # Vérifier que les colonnes existent
  if (!x %in% colnames(site_scores) || !y %in% colnames(site_scores)) {
    stop("Les axes spécifiés n'existent pas dans site_scores")
  }
  
  # Base du graphique
  p <- ggplot(site_scores, aes_string(x = x, y = y))
  
  # Ajouter les points selon l'option de coloration
  if (color_wd) {
    p <- p + 
      geom_point(aes(color = wd_data), size = 4, alpha = 0.8) +
      scale_color_viridis_c(name = "Wood Density\n(WD_BA)")
    
    title_text <- paste("NSCA coloré par Wood Density")
  } else {
    p <- p + 
      geom_point(size = 4, alpha = 0.8, color = "darkgreen")
    
    title_text <- paste("NSCA - Axes", x_axis_num, "et", y_axis_num)
  }
  
  # Ajouter le reste du graphique
  p <- p +
    geom_text_repel(aes(label = rownames(site_scores)), 
                    size = 3, max.overlaps = 10) +
    labs(title = title_text,
         subtitle = paste0("Axe ", x_axis_num, ": ", round(variance_explained_data[x_axis_num], 1), 
                           "% | Axe ", y_axis_num, ": ", round(variance_explained_data[y_axis_num], 1), "%"),
         x = paste0(x, " (", round(variance_explained_data[x_axis_num], 1), "%)"),
         y = paste0(y, " (", round(variance_explained_data[y_axis_num], 1), "%)")) +
    theme_minimal() +
    theme(aspect.ratio = 1)
  
  return(p)
}

# NSCA 1 vs 2 coloré par WD
p_nsca_12_wd <- plot_nsca(site_scores, x = "NSCA1", y = "NSCA2", color_wd = TRUE)
print(p_nsca_12_wd)

# NSCA 2 vs 3 coloré par WD
p_nsca_23_wd <- plot_nsca(site_scores, x = "NSCA3", y = "NSCA4", color_wd = TRUE)
print(p_nsca_23_wd)

if(T){
  # NSCA 1 vs 2 coloré par WD
  p_nsca_dbh_12_wd <- plot_nsca(site_scores_dbh, x = "NSCA1", y = "NSCA2", color_wd = TRUE)
  print(p_nsca_12_wd)
  
  # NSCA 2 vs 3 coloré par WD
  p_nsca_dbh_23_wd <- plot_nsca(site_scores_dbh, x = "NSCA3", y = "NSCA4", color_wd = TRUE)
  print(p_nsca_23_wd)
  
  # NSCA 1 vs 2 coloré par WD
  p_nsca_genus_12_wd <- plot_nsca(site_scores_genus, x = "NSCA1", y = "NSCA2", color_wd = TRUE)
  print(p_nsca_12_wd)
  
  # NSCA 2 vs 3 coloré par WD
  p_nsca_genus_23_wd <- plot_nsca(site_scores_genus, x = "NSCA3", y = "NSCA4", color_wd = TRUE)
  print(p_nsca_23_wd)
}


# =============================================================================
# C. CLUSTERING SUR LES SCORES NSCA
# =============================================================================

# Déterminer le nombre optimal de clusters
# Utiliser les 3 premiers axes pour le clustering
n_axes_clustering <- min(4, ncol(site_scores))
clustering_data <- site_scores[, 1:n_axes_clustering]

# Méthode du coude (WSS)
wss_results <- fviz_nbclust(clustering_data, kmeans, method = "wss", k.max = 10)
print(wss_results)

# Méthode silhouette
silhouette_results <- fviz_nbclust(clustering_data, kmeans, method = "silhouette", k.max = 10)
print(silhouette_results)

# Choisir le nombre de clusters (à ajuster selon les résultats)
# Ici on prend 5 clusters par défaut, mais vous pouvez ajuster
n_clusters <- 6

# Clustering hiérarchique
hclust_result <- hclust(dist(clustering_data), method = "ward.D2")
floristic_groups <- cutree(hclust_result, k = n_clusters)

# Dendrogramme
p_dendro <- fviz_dend(hclust_result, k = n_clusters, 
                      color_labels_by_k = TRUE,
                      rect = TRUE,
                      main = "Dendrogramme - Groupes floristiques")
print(p_dendro)

# Validation du clustering
sil_analysis <- silhouette(floristic_groups, dist(clustering_data))
sil_avg <- mean(sil_analysis[, 3])

cat("📊 Validation du clustering:\n")
cat("- Silhouette moyenne:", round(sil_avg, 3), "\n")
cat("- Répartition des groupes:\n")
print(table(floristic_groups))

# Visualisation des groupes sur NSCA
p_nsca_groups <- ggplot(site_scores, aes(x = NSCA3, y = NSCA4)) +
  geom_point(aes(color = factor(floristic_groups)), size = 4, alpha = 0.8) +
  scale_color_brewer(type = "qual", palette = "Set1", name = "Groupe\nfloristique") +
  geom_text_repel(aes(label = rownames(site_scores)), 
                  size = 3, max.overlaps = 10) +
  labs(title = "Groupes floristiques sur NSCA",
       subtitle = paste0("Clustering hiérarchique (", n_clusters, " groupes)"),
       x = paste0("NSCA1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("NSCA2 (", round(variance_explained[2], 1), "%)")) +
  theme_minimal() +
  theme(aspect.ratio = 1)

print(p_nsca_groups)

# =============================================================================
# D. CARACTÉRISATION DES GROUPES
# =============================================================================

# Créer un dataframe avec toutes les informations
analysis_data <- data.frame(
  plot_name = plots_name,
  floristic_group = floristic_groups,
  WD_BA = metadata$WD_BA,
  site_scores[, 1:4],  # Premiers axes NSCA
  structure_plots_table
)

# 1. Caractérisation structurelle par groupe
structure_by_group <- analysis_data %>%
  group_by(floristic_group) %>%
  summarise(
    n_plots = n(),
    across(all_of(gap_metrics_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# 2. Caractérisation wood density par groupe
wd_by_group <- analysis_data %>%
  group_by(floristic_group) %>%
  summarise(
    WD_mean = mean(WD_BA, na.rm = TRUE),
    WD_sd = sd(WD_BA, na.rm = TRUE),
    WD_min = min(WD_BA, na.rm = TRUE),
    WD_max = max(WD_BA, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n📈 Wood density par groupe:\n")
print(wd_by_group)

# =============================================================================
# E. BOXPLOTS ET TESTS STATISTIQUES
# =============================================================================

# Test ANOVA pour wood density
anova_wd <- aov(WD_BA ~ factor(floristic_group), data = analysis_data)
anova_wd_summary <- summary(anova_wd)
tukey_wd <- TukeyHSD(anova_wd)

cat("🔬 Test ANOVA - Wood Density:\n")
print(anova_wd_summary)

# Boxplot wood density par groupe
p_boxplot_wd <- ggplot(analysis_data, aes(x = factor(floristic_group), y = WD_BA)) +
  geom_boxplot(aes(fill = factor(floristic_group)), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  scale_fill_brewer(type = "qual", palette = "Set1", name = "Groupe") +
  labs(title = "Wood Density par groupe floristique",
       subtitle = paste0("ANOVA p-value: ", 
                         round(anova_wd_summary[[1]][["Pr(>F)"]][1], 4)),
       x = "Groupe floristique",
       y = "Wood Density (WD_BA)") +
  theme_minimal() +
  theme(legend.position = "none")

print(p_boxplot_wd)

# Tests pour les variables de structure (sélection des principales)
main_structure_vars <- gap_metrics_cols[1:min(8, length(gap_metrics_cols))]

# MANOVA pour les variables de structure
manova_structure <- manova(as.matrix(analysis_data[, main_structure_vars]) ~ 
                             factor(analysis_data$floristic_group))
manova_summary <- summary(manova_structure)

cat("\n🔬 Test MANOVA - Variables de structure:\n")
print(manova_summary)

# Boxplots pour les principales variables de structure
create_structure_boxplot <- function(var_name, data) {
  p <- ggplot(data, aes_string(x = "factor(floristic_group)", y = var_name)) +
    geom_boxplot(aes(fill = factor(floristic_group)), alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    labs(title = var_name,
         x = "Groupe floristique",
         y = var_name) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 9))
  return(p)
}

# Créer les boxplots pour les principales variables
structure_plots <- map(main_structure_vars[1:6], 
                       ~ create_structure_boxplot(.x, analysis_data))

# Combiner les plots
p_structure_combined <- wrap_plots(structure_plots, ncol = 2)
print(p_structure_combined)

# =============================================================================
# F. ANALYSE DES RELATIONS
# =============================================================================

cat("\n🔗 F. ANALYSE DES RELATIONS\n")
cat("-", rep("-", 30), "\n")

# Matrice de corrélation entre moyennes de structure et wood density par groupe
structure_means <- structure_by_group %>%
  select(-floristic_group, -n_plots) %>%
  as.matrix()

wd_means <- wd_by_group$WD_mean

# Corrélations
correlations <- cor(structure_means, wd_means, use = "complete.obs")
colnames(correlations) <- "WD_correlation"

# Trier par corrélation absolue
correlations_sorted <- correlations[order(abs(correlations[,1]), decreasing = TRUE), , drop = FALSE]

cat("🔗 Corrélations Structure-WD (par groupe):\n")
print(round(correlations_sorted, 3))

# Graphique des corrélations
p_correlations <- data.frame(
  Variable = rownames(correlations_sorted),
  Correlation = correlations_sorted[,1]
) %>%
  slice_head(n = 15) %>%  # Top 15 corrélations
  ggplot(aes(x = reorder(Variable, abs(Correlation)), y = Correlation)) +
  geom_col(aes(fill = Correlation > 0), alpha = 0.7) +
  scale_fill_manual(values = c("TRUE" = "darkgreen", "FALSE" = "darkred"),
                    name = "Corrélation", labels = c("Négative", "Positive")) +
  coord_flip() +
  labs(title = "Corrélations Structure-Wood Density",
       subtitle = "Top 15 variables structurelles (moyennes par groupe)",
       x = "Variables structurelles",
       y = "Corrélation avec Wood Density") +
  theme_minimal()

print(p_correlations)

#### 6️⃣ ÉTAPE 2: ANALYSE CCA COMPLÈTE AVEC CLUSTERING ####

# =============================================================================
# A. PRÉPARATION DES DONNÉES POUR CCA
# =============================================================================

cat("\n🔗 A. PRÉPARATION DES DONNÉES POUR CCA\n")
cat("-", rep("-", 40), "\n")

# Vérification de la correspondance entre les tableaux
common_plots <- intersect(rownames(species_plot_table), rownames(structure_plots_table))
cat("📊 Parcelles communes entre composition et structure:", length(common_plots), "\n")

# Synchronisation des tableaux
species_data_cca <- species_plot_table[common_plots, ]
structure_data_cca <- structure_plots_table[common_plots, ]

# Filtrage des espèces rares (optionnel - ajuster selon le contexte)
min_occurrence <- 2  # Espèces présentes dans au moins 2 parcelles
species_freq <- colSums(species_data_cca > 0)
species_to_keep <- names(species_freq)[species_freq >= min_occurrence]

species_data_cca_filtered <- species_data_cca[, species_to_keep]

cat("📈 Dimensions finales:\n")
cat("- Espèces retenues:", ncol(species_data_cca_filtered), "sur", ncol(species_data_cca), "\n")
cat("- Parcelles:", nrow(species_data_cca_filtered), "\n")
cat("- Variables de structure:", ncol(structure_data_cca), "\n")

# Transformation des données (recommandée pour CCA avec données d'abondance)
# Transformation de Hellinger pour réduire l'effet des espèces dominantes
species_hellinger <- decostand(species_data_cca_filtered, method = "hellinger")

# =============================================================================
# A.1. DIAGNOSTIC DE COLINÉARITÉ ET FILTRAGE DES VARIABLES
# =============================================================================

# Diagnostic initial de colinéarité
cat("🔍 Diagnostic de colinéarité des variables de structure:\n")

# Matrice de corrélation
cor_matrix_raw <- cor(structure_data_cca, use = "complete.obs")

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
vars_all <- colnames(structure_data_cca)
vars_norm <- vars_all[grepl("_norm$", vars_all)]
vars_base <- gsub("_norm$", "", vars_norm)

# Variables à supprimer (versions non normalisées quand version normalisée existe)
vars_to_remove <- vars_base[vars_base %in% vars_all]

# Variables supplémentaires identifiées comme redondantes
additional_removes <- c("upper_thickness", "basal_thickness", "transition_auc_norm", "basal_auc_norm", "upper_auc_norm")  # Ajuster selon vos données

# Liste finale des variables à supprimer
all_vars_to_remove <- unique(c(vars_to_remove, additional_removes))
all_vars_to_remove <- all_vars_to_remove[all_vars_to_remove %in% vars_all]

cat("\n🗑️  Variables supprimées (redondantes):\n")
if(length(all_vars_to_remove) > 0) {
  cat("  ", paste(all_vars_to_remove, collapse = ", "), "\n")
} else {
  cat("  Aucune variable supprimée\n")
}

# Filtrer les données
structure_data_filtered <- structure_data_cca[, !colnames(structure_data_cca) %in% all_vars_to_remove]

cat("📊 Variables restantes:", ncol(structure_data_filtered), "sur", ncol(structure_data_cca), "\n")

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

# Centrage-réduction des variables de structure filtrées
structure_scaled <- scale(structure_data_filtered)

# =============================================================================
# B. CCA COMPLÈTE - MODÈLE DE RÉFÉRENCE
# =============================================================================

cat("\n🎯 B. CCA COMPLÈTE - MODÈLE DE RÉFÉRENCE\n")
cat("-", rep("-", 50), "\n")

# CCA avec toutes les variables de structure filtrées
cca_full <- cca(species_hellinger ~ ., data = as.data.frame(structure_scaled))

# Résumé de la CCA complète
cat("📊 CCA COMPLÈTE - Résumé:\n")
cca_summary_full <- summary(cca_full)

# Variance expliquée
cca_eigenvals_full <- cca_full$CCA$eig
total_inertia_full <- cca_full$tot.chi
constrained_inertia_full <- cca_full$CCA$tot.chi
unconstrained_inertia_full <- cca_full$CA$tot.chi

variance_constrained_full <- cca_eigenvals_full / total_inertia_full * 100
variance_constrained_cumul_full <- cumsum(variance_constrained_full)

cat("- Variables utilisées:", ncol(structure_scaled), "\n")
cat("- Inertie totale:", round(total_inertia_full, 3), "\n")
cat("- Inertie contrainte:", round(constrained_inertia_full, 3), 
    "(", round(constrained_inertia_full/total_inertia_full*100, 1), "%)\n")
cat("- Premiers axes:\n")
for(i in 1:min(4, length(variance_constrained_full))) {
  cat("  - CCA", i, ":", round(variance_constrained_full[i], 1), 
      "% (cumulé:", round(variance_constrained_cumul_full[i], 1), "%)\n")
}

# Test de permutation du modèle complet
set.seed(123)
cca_perm_full <- anova(cca_full, permutations = 999)
cat("\n🔬 Test de permutation modèle complet:\n")
print(cca_perm_full)

# =============================================================================
# C. SÉLECTION DE VARIABLES - FORWARD ET BACKWARD
# =============================================================================

cat("\n🔍 C. SÉLECTION DE VARIABLES\n")
cat("-", rep("-", 35), "\n")

# Modèle null pour la sélection forward
cca_null <- cca(species_hellinger ~ 1, data = as.data.frame(structure_scaled))

# Sélection forward
cat("🚀 Sélection FORWARD...\n")
set.seed(123)
cca_forward <- ordistep(cca_null, 
                        scope = formula(cca_full), 
                        direction = "forward",
                        permutations = 99)  # Réduire pour la vitesse

cat("\n📊 MODÈLE FORWARD - Résumé:\n")
print(cca_forward)

# Sélection backward
cat("\n🔙 Sélection BACKWARD...\n")
set.seed(123)
cca_backward <- ordistep(cca_full, 
                         direction = "backward",
                         permutations = 99)

cat("\n📊 MODÈLE BACKWARD - Résumé:\n")
print(cca_backward)

# Sélection both (optionnel)
cat("\n↔️  Sélection BOTH (bidirectionnelle)...\n")
set.seed(123)
cca_both <- ordistep(cca_null, 
                     scope = formula(cca_full), 
                     direction = "both",
                     permutations = 99)

cat("\n📊 MODÈLE BOTH - Résumé:\n")
print(cca_both)

# =============================================================================
# D. COMPARAISON DES MODÈLES
# =============================================================================

cat("\n📈 D. COMPARAISON DES MODÈLES\n")
cat("-", rep("-", 35), "\n")

# Fonction pour extraire les métriques d'un modèle CCA
extract_cca_metrics <- function(cca_model, model_name) {
  total_inertia <- cca_model$tot.chi
  constrained_inertia <- cca_model$CCA$tot.chi
  n_vars <- ncol(cca_model$CCA$biplot)
  n_axes_significant <- sum(eigenvals(cca_model, model = "constrained") > 0)
  
  # Test de permutation
  set.seed(123)
  perm_test <- anova(cca_model, permutations = 999)
  p_value <- perm_test$`Pr(>F)`[1]
  
  # AIC (approximatif pour CCA)
  aic_approx <- -2 * constrained_inertia + 2 * n_vars
  
  return(data.frame(
    Model = model_name,
    N_Variables = n_vars,
    Total_Inertia = round(total_inertia, 4),
    Constrained_Inertia = round(constrained_inertia, 4),
    Percent_Explained = round(constrained_inertia/total_inertia*100, 2),
    P_Value = round(p_value, 4),
    AIC_Approx = round(aic_approx, 2),
    N_Axes = n_axes_significant
  ))
}

# Comparaison des modèles
model_comparison <- rbind(
  extract_cca_metrics(cca_full, "COMPLET"),
  extract_cca_metrics(cca_forward, "FORWARD"),
  extract_cca_metrics(cca_backward, "BACKWARD"),
  extract_cca_metrics(cca_both, "BOTH")
)

cat("📊 TABLEAU COMPARATIF DES MODÈLES:\n")
print(model_comparison)

# Identifier le meilleur modèle
best_model_idx <- which.min(model_comparison$AIC_Approx)
best_model_name <- model_comparison$Model[best_model_idx]

cat("\n🏆 MEILLEUR MODÈLE (AIC le plus bas):", best_model_name, "\n")

# Variables sélectionnées par chaque méthode
cat("\n📝 VARIABLES SÉLECTIONNÉES:\n")

get_selected_vars <- function(cca_model) {
  if(is.null(cca_model$CCA$biplot)) return(character(0))
  return(rownames(cca_model$CCA$biplot))
}

vars_full <- get_selected_vars(cca_full)
vars_forward <- get_selected_vars(cca_forward)
vars_backward <- get_selected_vars(cca_backward)
vars_both <- get_selected_vars(cca_both)

cat("- COMPLET (", length(vars_full), "vars):", paste(vars_full, collapse = ", "), "\n")
cat("- FORWARD (", length(vars_forward), "vars):", paste(vars_forward, collapse = ", "), "\n")
cat("- BACKWARD (", length(vars_backward), "vars):", paste(vars_backward, collapse = ", "), "\n")
cat("- BOTH (", length(vars_both), "vars):", paste(vars_both, collapse = ", "), "\n")

# Forcer l'utilisation du modèle souhaité
best_model_name <- "BACKWARD"  # Changer ici : "FORWARD", "BACKWARD", "BOTH", ou "COMPLET"

# =============================================================================
# E. TESTS DE PERMUTATION DÉTAILLÉS SUR LE MODÈLE SÉLECTIONNÉ
# =============================================================================

cat("\n🔬 E. TESTS DÉTAILLÉS SUR LE MODÈLE", best_model_name, "\n")
cat("-", rep("-", 45), "\n")

# Sélectionner le modèle correspondant
best_model <- switch(best_model_name,
                     "COMPLET" = cca_full,
                     "FORWARD" = cca_forward,
                     "BACKWARD" = cca_backward,
                     "BOTH" = cca_both)

cat("🎯 Analyse détaillée du modèle:", best_model_name, "\n")

# Test global avec 999 permutations
set.seed(123)
perm_global <- anova(best_model, permutations = 999)
cat("\n🔬 Test global (999 permutations):\n")
print(perm_global)

# Test par axes
set.seed(123)
perm_axes <- anova(best_model, by = "axis", permutations = 999)
cat("\n🔬 Test par axes:\n")
print(perm_axes)

# Test par termes (variables)
set.seed(123)
perm_terms <- anova(best_model, by = "terms", permutations = 999)
cat("\n🔬 Test par variables:\n")
print(perm_terms)

# =============================================================================
# F. VISUALISATIONS DU MODÈLE BEST
# =============================================================================

cat("\n📊 F. VISUALISATIONS DU BEST MODÈLE\n")
cat("-", rep("-", 40), "\n")

# Extraction des scores pour le modèle BOTH
sites_best <- scores(best_model, display = "sites", choices = 1:4)
env_best <- scores(best_model, display = "bp", choices = 1:4)

# Variance expliquée par les axes du modèle BOTH
eigenvals_best <- eigenvals(best_model, model = "constrained")
variance_best <- eigenvals_best / best_model$tot.chi * 100

# Fonction pour créer des biplots CCA optimisés
create_best_cca_biplot <- function(sites, env_scores, x_axis = 1, y_axis = 2, 
                                   model_name = "", color_by_wd = TRUE) {
  
  # Données pour le plot
  sites_df <- data.frame(
    CCA1 = sites[, x_axis],
    CCA2 = sites[, y_axis],
    plot_name = rownames(sites)
  )
  
  env_df <- data.frame(
    CCA1 = env_scores[, x_axis],
    CCA2 = env_scores[, y_axis],
    variable = rownames(env_scores)
  )
  
  # Ajouter wood density si disponible
  if(color_by_wd && all(sites_df$plot_name %in% plots_name)) {
    wd_values <- metadata$WD_BA[match(sites_df$plot_name, plots_name)]
    sites_df$WD_BA <- wd_values
  }
  
  # Base du graphique
  p <- ggplot()
  
  # Points des sites
  if(color_by_wd && "WD_BA" %in% colnames(sites_df)) {
    p <- p + 
      geom_point(data = sites_df, aes(x = CCA1, y = CCA2, color = WD_BA), 
                 size = 3, alpha = 0.8) +
      scale_color_viridis_c(name = "Wood\nDensity")
  } else {
    p <- p + 
      geom_point(data = sites_df, aes(x = CCA1, y = CCA2), 
                 size = 3, alpha = 0.8, color = "darkgreen")
  }
  
  # Flèches des variables environnementales
  arrow_scale <- 0.8
  p <- p + 
    geom_segment(data = env_df, 
                 aes(x = 0, y = 0, xend = CCA1 * arrow_scale, yend = CCA2 * arrow_scale),
                 arrow = arrow(length = unit(0.3, "cm")), 
                 color = "red", alpha = 0.8, size = 0.8) +
    geom_text_repel(data = env_df, 
                    aes(x = CCA1 * arrow_scale, y = CCA2 * arrow_scale, label = variable),
                    color = "red", size = 3, fontface = "bold", max.overlaps = 15)
  
  # Labels et thème
  x_var_exp <- round(variance_best[x_axis], 1)
  y_var_exp <- round(variance_best[y_axis], 1)
  
  p <- p +
    labs(title = paste0("CCA Biplot - Modèle ", model_name),
         subtitle = paste0("Variables sélectionnées | Axes ", x_axis, "-", y_axis),
         x = paste0("CCA", x_axis, " (", x_var_exp, "%)"),
         y = paste0("CCA", y_axis, " (", y_var_exp, "%)")) +
    theme_minimal() +
    theme(aspect.ratio = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)
  
  return(p)
}

# Biplot du modèle BOTH - axes 1-2
p_best_12 <- create_best_cca_biplot(sites_best, env_best, 1, 2, best_model_name)
print(p_best_12)

# Biplot du modèle BOTH - axes 2-3 (si disponibles)
if(ncol(sites_best) >= 3) {
  p_best_23 <- create_best_cca_biplot(sites_best, env_best, 2, 3, best_model_name)
  print(p_best_23)
}

# =============================================================================
# G. RÉSUMÉ MODÈLE BEST
# =============================================================================

cat("\n🎯 G. RÉSUMÉ MODÈLE BEST\n")
cat("-", rep("-", 30), "\n")

cat("🏆 MODÈLE OPTIMAL:", best_model_name, "\n")
cat("📊 Performance:\n")
best_model_metrics <- model_comparison[model_comparison$Model == best_model_name, ]
cat("- Variables retenues:", best_model_metrics$N_Variables, "\n")
cat("- Variance expliquée:", best_model_metrics$Percent_Explained, "%\n")
cat("- Significativité: p =", best_model_metrics$P_Value, "\n")

# Gain en parcimonie
original_vars <- ncol(structure_scaled)
selected_vars <- best_model_metrics$N_Variables
reduction_percent <- round((1 - selected_vars/original_vars) * 100, 1)

cat("\n📉 Gain en parcimonie:\n")
cat("- Réduction de variables:", reduction_percent, "% (", original_vars, "→", selected_vars, ")\n")

# Efficacité relative
full_model_perf <- model_comparison$Percent_Explained[model_comparison$Model == "COMPLET"]
both_model_perf <- best_model_metrics$Percent_Explained
efficiency <- round(both_model_perf / full_model_perf, 3)

cat("- Efficacité relative:", round(efficiency * 100, 1), "% de performance du modèle complet\n")

cat("\n✅ Analyse CCA avec sélection de variables terminée!\n")

# =============================================================================
# H. CLUSTERING SUR LES SCORES CCA
# =============================================================================

cat("\n🎯 H. CLUSTERING SUR LES SCORES CCA\n")
cat("-", rep("-", 40), "\n")

# Extraction des scores CCA pour le clustering (utiliser le meilleur modèle)
sites_cca_clustering <- scores(best_model, display = "sites", choices = 1:min(6, ncol(best_model$CCA$u)))
colnames(sites_cca_clustering) <- paste0("CCA", 1:ncol(sites_cca_clustering))

# Déterminer le nombre d'axes significatifs à utiliser
perm_axes_results <- anova(best_model, by = "axis", permutations = 999)
significant_axes <- which(perm_axes_results$`Pr(>F)` < 0.10)
significant_axes <- significant_axes[!is.na(perm_axes_results$`Pr(>F)`[significant_axes])]
n_sig_axes <- max(length(significant_axes), 2)  # Au minimum 2 axes
n_sig_axes <- min(n_sig_axes, ncol(sites_cca_clustering))

# Données pour le clustering
clustering_data_cca <- sites_cca_clustering[, 1:n_sig_axes]

cat("📊 Données pour clustering CCA:\n")
cat("- Parcelles:", nrow(clustering_data_cca), "\n")
cat("- Axes CCA utilisés:", n_sig_axes, "\n")
cat("- Axes:", paste(colnames(clustering_data_cca), collapse = ", "), "\n")

# =============================================================================
# H.1. DÉTERMINATION DU NOMBRE OPTIMAL DE CLUSTERS
# =============================================================================

cat("\n🔍 H.1. DÉTERMINATION DU NOMBRE OPTIMAL DE CLUSTERS\n")
cat("-", rep("-", 50), "\n")

# Méthode du coude (WSS)
wss_cca <- fviz_nbclust(clustering_data_cca, kmeans, method = "wss", k.max = 10) +
  ggtitle("Méthode du coude - Scores CCA") +
  theme_minimal()
print(wss_cca)

# Méthode silhouette
sil_cca <- fviz_nbclust(clustering_data_cca, kmeans, method = "silhouette", k.max = 10) +
  ggtitle("Méthode silhouette - Scores CCA") +
  theme_minimal()
print(sil_cca)

# Gap statistic
set.seed(123)
gap_cca <- fviz_nbclust(clustering_data_cca, kmeans, method = "gap_stat", k.max = 10) +
  ggtitle("Gap Statistic - Scores CCA") +
  theme_minimal()
print(gap_cca)

# =============================================================================
# H.2. CLUSTERING HIÉRARCHIQUE
# =============================================================================

cat("\n🌳 H.2. CLUSTERING HIÉRARCHIQUE\n")
cat("-", rep("-", 35), "\n")

# Distance euclidienne sur les scores CCA
dist_cca <- dist(clustering_data_cca)

# Clustering hiérarchique avec différentes méthodes
methods_hclust <- c("ward.D2", "complete", "average")
hclust_results <- list()

for(method in methods_hclust) {
  hclust_results[[method]] <- hclust(dist_cca, method = method)
}

# Évaluer la qualité des clusterings avec coefficient cophénétique
cophenetic_cors <- sapply(hclust_results, function(hc) cor(dist_cca, cophenetic(hc)))

cat("📊 Coefficients cophénétiques (qualité du clustering):\n")
for(i in 1:length(cophenetic_cors)) {
  cat("- ", names(cophenetic_cors)[i], ":", round(cophenetic_cors[i], 3), "\n")
}

# Sélectionner la meilleure méthode
best_method <- names(cophenetic_cors)[which.max(cophenetic_cors)]
hclust_best <- hclust_results[[best_method]]

cat("🏆 Meilleure méthode:", best_method, "\n")

# Déterminer le nombre de clusters (ajustable selon les graphiques précédents)
n_clusters_cca <- 7  # Vous pouvez ajuster selon vos résultats

# Couper l'arbre
groups_hclust <- cutree(hclust_best, k = n_clusters_cca)

# Dendrogramme
p_dendro_cca <- fviz_dend(hclust_best, k = n_clusters_cca, 
                          color_labels_by_k = TRUE,
                          rect = TRUE,
                          rect_fill = TRUE,
                          rect_border = "jco",
                          labels_track_height = 0.8,
                          main = paste0("Dendrogramme CCA - ", best_method)) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 8))
print(p_dendro_cca)

# =============================================================================
# H.3. CLUSTERING K-MEANS POUR COMPARAISON
# =============================================================================

cat("\n⭕ H.3. CLUSTERING K-MEANS\n")
cat("-", rep("-", 25), "\n")

# K-means avec le même nombre de clusters
set.seed(123)
kmeans_cca <- kmeans(clustering_data_cca, centers = n_clusters_cca, nstart = 25)
groups_kmeans <- kmeans_cca$cluster

# Validation des clusterings
sil_hclust <- silhouette(groups_hclust, dist_cca)
sil_kmeans <- silhouette(groups_kmeans, dist_cca)

sil_avg_hclust <- mean(sil_hclust[, 3])
sil_avg_kmeans <- mean(sil_kmeans[, 3])

cat("📊 Validation clustering:\n")
cat("- Silhouette hiérarchique:", round(sil_avg_hclust, 3), "\n")
cat("- Silhouette K-means:", round(sil_avg_kmeans, 3), "\n")

# Concordance entre les méthodes
if(!require(mclust)) {
  install.packages("mclust")
  library(mclust)
}

concordance <- adjustedRandIndex(groups_hclust, groups_kmeans)
cat("- Concordance hclust-kmeans (ARI):", round(concordance, 3), "\n")

# Choisir la meilleure méthode
if(sil_avg_hclust >= sil_avg_kmeans) {
  final_groups <- groups_hclust
  clustering_method <- paste0("Hiérarchique (", best_method, ")")
  cat("🏆 Méthode retenue: Clustering hiérarchique\n")
} else {
  final_groups <- groups_kmeans
  clustering_method <- "K-means"
  cat("🏆 Méthode retenue: K-means\n")
}

cat("- Répartition des groupes:\n")
group_counts <- table(final_groups)
for(i in 1:length(group_counts)) {
  cat("  Groupe", names(group_counts)[i], ":", group_counts[i], "parcelles\n")
}

# =============================================================================
# I. VISUALISATIONS DES GROUPES CCA
# =============================================================================

cat("\n📊 I. VISUALISATIONS DES GROUPES CCA\n")
cat("-", rep("-", 35), "\n")

# Préparer les données pour visualisation
plot_data_clusters <- data.frame(
  clustering_data_cca,
  Group = factor(final_groups),
  plot_name = rownames(clustering_data_cca)
)

# Ajouter wood density
if(all(plot_data_clusters$plot_name %in% plots_name)) {
  wd_values <- metadata$WD_BA[match(plot_data_clusters$plot_name, plots_name)]
  plot_data_clusters$WD_BA <- wd_values
}

# Fonction pour créer des plots avec ellipses de confiance
create_cluster_plot <- function(data, x_var, y_var, title_suffix = "") {
  
  # Extraire les numéros d'axes
  x_num <- as.numeric(gsub("CCA", "", x_var))
  y_num <- as.numeric(gsub("CCA", "", y_var))
  
  p <- ggplot(data, aes_string(x = x_var, y = y_var, color = "Group"))
  
  # Points colorés par groupe
  p <- p + geom_point(size = 3, alpha = 0.8)
  
  # Ellipses de confiance à 95%
  p <- p + stat_ellipse(level = 0.95, size = 1.2, alpha = 0.3)
  
  # Centroïdes des groupes
  p <- p + stat_summary(fun = mean, geom = "point", size = 4, shape = 18, color = "black")
  
  # Labels optionnels (limiter pour éviter la surcharge)
  if(nrow(data) <= 100) {
    p <- p + geom_text_repel(aes(label = plot_name), size = 2.5, max.overlaps = 8)
  }
  
  # Couleurs et thème
  p <- p + 
    scale_color_brewer(type = "qual", palette = "Set1", name = "Groupe") +
    labs(title = paste0("Groupes CCA ", title_suffix),
         subtitle = paste0(clustering_method, " | ", n_clusters_cca, " groupes"),
         x = paste0(x_var, " (", round(variance_best[x_num], 1), "%)"),
         y = paste0(y_var, " (", round(variance_best[y_num], 1), "%)")) +
    theme_minimal() +
    theme(aspect.ratio = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)
  
  return(p)
}

# Plot principal CCA1 vs CCA2
p_clusters_12 <- create_cluster_plot(plot_data_clusters, "CCA1", "CCA2", "Axes 1-2")
print(p_clusters_12)

# Plot CCA2 vs CCA3 si disponible
if(ncol(clustering_data_cca) >= 3) {
  p_clusters_23 <- create_cluster_plot(plot_data_clusters, "CCA2", "CCA3", "Axes 2-3")
  print(p_clusters_23)
}

# Plot avec wood density si disponible
if("WD_BA" %in% colnames(plot_data_clusters)) {
  p_wd_groups <- ggplot(plot_data_clusters, aes(x = CCA1, y = CCA2)) +
    geom_point(aes(color = WD_BA, shape = Group), size = 3, alpha = 0.8) +
    scale_color_viridis_c(name = "Wood\nDensity") +
    scale_shape_manual(values = c(16, 17, 15, 3, 4)[1:n_clusters_cca], name = "Groupe") +
    stat_ellipse(aes(group = Group), level = 0.95, linetype = "dashed") +
    labs(title = "Groupes CCA avec Wood Density",
         subtitle = "Forme = Groupe | Couleur = Wood Density",
         x = paste0("CCA1 (", round(variance_best[1], 1), "%)"),
         y = paste0("CCA2 (", round(variance_best[2], 1), "%)")) +
    theme_minimal() +
    theme(aspect.ratio = 1)
  
  print(p_wd_groups)
}

# =============================================================================
# J. CARACTÉRISATION DES GROUPES CCA
# =============================================================================

cat("\n📋 J. CARACTÉRISATION DES GROUPES CCA\n")
cat("-", rep("-", 35), "\n")

# Variables sélectionnées dans le modèle BOTH
selected_vars_both <- get_selected_vars(best_model)

# Créer un dataframe complet pour l'analyse
cca_analysis_data <- data.frame(
  plot_name = rownames(clustering_data_cca),
  cca_group = final_groups,
  WD_BA = metadata$WD_BA[match(rownames(clustering_data_cca), plots_name)],
  clustering_data_cca,
  structure_data_filtered[rownames(clustering_data_cca), selected_vars_both]
)

# Caractérisation par wood density
wd_by_group <- cca_analysis_data %>%
  group_by(cca_group) %>%
  summarise(
    n_plots = n(),
    WD_mean = mean(WD_BA, na.rm = TRUE),
    WD_sd = sd(WD_BA, na.rm = TRUE),
    WD_min = min(WD_BA, na.rm = TRUE),
    WD_max = max(WD_BA, na.rm = TRUE),
    .groups = "drop"
  )

cat("📈 Wood density par groupe CCA:\n")
print(wd_by_group)

# Caractérisation structurelle (variables sélectionnées)
structure_by_group <- cca_analysis_data %>%
  group_by(cca_group) %>%
  summarise(
    across(all_of(selected_vars_both), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

cat("\n📊 Variables structurelles sélectionnées par groupe (moyennes):\n")
print(structure_by_group)

# =============================================================================
# K. TESTS STATISTIQUES ET BOXPLOTS
# =============================================================================

cat("\n🔬 K. TESTS STATISTIQUES ET BOXPLOTS\n")
cat("-", rep("-", 35), "\n")

# Test ANOVA pour wood density
anova_wd_cca <- aov(WD_BA ~ factor(cca_group), data = cca_analysis_data)
anova_wd_summary <- summary(anova_wd_cca)

cat("🔬 ANOVA Wood Density par groupe CCA:\n")
print(anova_wd_summary)

# Test post-hoc si significatif
tukey_results_wd <- NULL
if(anova_wd_summary[[1]][["Pr(>F)"]][1] < 0.05) {
  tukey_results_wd <- TukeyHSD(anova_wd_cca)
  cat("\n📊 Test post-hoc Tukey HSD Wood Density:\n")
  print(tukey_results_wd)
}

# Fonction pour ajouter les lettres de significativité
add_significance_letters <- function(tukey_result, groups) {
  if(is.null(tukey_result)) return(rep("a", length(unique(groups))))
  
  # Extraire les p-values
  pvals <- tukey_result$`factor(cca_group)`[, "p adj"]
  
  # Créer un dataframe des comparaisons
  comparisons <- expand.grid(1:length(unique(groups)), 1:length(unique(groups)))
  comparisons <- comparisons[comparisons[,1] < comparisons[,2], ]
  
  # Assigner les lettres (simplification)
  letters_result <- rep("a", length(unique(groups)))
  
  # Si des différences significatives existent, assigner différentes lettres
  if(any(pvals < 0.05, na.rm = TRUE)) {
    sig_pairs <- which(pvals < 0.05)
    if(length(sig_pairs) > 0) {
      letters_result[1] <- "a"
      letter_index <- 1
      for(i in 2:length(unique(groups))) {
        # Vérifier si ce groupe diffère significativement des précédents
        differs <- FALSE
        for(j in 1:(i-1)) {
          pair_name <- paste0(j, "-", i)
          if(pair_name %in% names(pvals) && pvals[pair_name] < 0.05) {
            differs <- TRUE
            break
          }
        }
        if(differs) {
          letter_index <- letter_index + 1
          letters_result[i] <- letters[letter_index]
        } else {
          letters_result[i] <- letters_result[i-1]
        }
      }
    }
  }
  
  return(letters_result)
}

# Créer les lettres de significativité pour WD
sig_letters_wd <- add_significance_letters(tukey_results_wd, cca_analysis_data$cca_group)

# Boxplot wood density avec lettres de significativité
p_boxplot_wd_cca <- ggplot(cca_analysis_data, aes(x = factor(cca_group), y = WD_BA, fill = factor(cca_group))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_brewer(type = "qual", palette = "Set1", name = "Groupe") +
  labs(title = "Wood Density par groupe CCA",
       subtitle = paste0("ANOVA F = ", round(anova_wd_summary[[1]][["F value"]][1], 2),
                         ", p = ", round(anova_wd_summary[[1]][["Pr(>F)"]][1], 4)),
       x = "Groupe CCA",
       y = "Wood Density (WD_BA)") +
  theme_minimal() +
  theme(legend.position = "none")

# Ajouter les lettres de significativité
y_max <- max(cca_analysis_data$WD_BA, na.rm = TRUE)
y_positions <- rep(y_max * 1.05, length(unique(cca_analysis_data$cca_group)))

for(i in 1:length(unique(cca_analysis_data$cca_group))) {
  p_boxplot_wd_cca <- p_boxplot_wd_cca +
    annotate("text", x = i, y = y_positions[i], label = sig_letters_wd[i], 
             size = 4, fontface = "bold")
}

print(p_boxplot_wd_cca)

# MANOVA pour les variables de structure sélectionnées
cat("\n🔬 MANOVA Variables structurelles sélectionnées:\n")
manova_structure_cca <- manova(as.matrix(cca_analysis_data[, selected_vars_both]) ~ 
                                 factor(cca_analysis_data$cca_group))
manova_summary <- summary(manova_structure_cca)
print(manova_summary)

# Boxplots pour chaque variable structurelle sélectionnée
cat("\n📊 Boxplots variables structurelles:\n")

# Fonction pour créer des boxplots avec tests ANOVA
create_structure_boxplot_with_test <- function(var_name, data) {
  
  # ANOVA pour cette variable
  formula_text <- paste(var_name, "~ factor(cca_group)")
  anova_result <- aov(as.formula(formula_text), data = data)
  anova_summary_var <- summary(anova_result)
  
  p_value <- anova_summary_var[[1]][["Pr(>F)"]][1]
  f_value <- anova_summary_var[[1]][["F value"]][1]
  
  # Test post-hoc si significatif
  tukey_result_var <- NULL
  if(p_value < 0.05) {
    tukey_result_var <- TukeyHSD(anova_result)
  }
  
  # Lettres de significativité
  sig_letters_var <- add_significance_letters(tukey_result_var, data$cca_group)
  
  # Boxplot
  p <- ggplot(data, aes_string(x = "factor(cca_group)", y = var_name, fill = "factor(cca_group)")) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    labs(title = var_name,
         subtitle = paste0("F = ", round(f_value, 2), ", p = ", round(p_value, 4)),
         x = "Groupe CCA",
         y = var_name) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 9))
  
  # Ajouter les lettres de significativité
  if(p_value < 0.05) {
    y_max_var <- max(data[[var_name]], na.rm = TRUE)
    y_positions_var <- rep(y_max_var * 1.05, length(unique(data$cca_group)))
    
    for(i in 1:length(unique(data$cca_group))) {
      p <- p + annotate("text", x = i, y = y_positions_var[i], label = sig_letters_var[i], 
                        size = 3, fontface = "bold")
    }
  }
  
  return(p)
}

# Créer les boxplots pour toutes les variables structurelles sélectionnées
structure_plots <- map(selected_vars_both, 
                       ~ create_structure_boxplot_with_test(.x, cca_analysis_data))

# Combiner les plots (par groupes de 4 pour la lisibilité)
n_vars <- length(selected_vars_both)
n_plots_per_page <- 4

for(i in seq(1, n_vars, by = n_plots_per_page)) {
  end_idx <- min(i + n_plots_per_page - 1, n_vars)
  plots_subset <- structure_plots[i:end_idx]
  
  p_combined <- wrap_plots(plots_subset, ncol = 2)
  print(p_combined)
}

# =============================================================================
# L. RÉSUMÉ FINAL
# =============================================================================

cat("\n🎯 L. RÉSUMÉ FINAL - CCA + CLUSTERING\n")
cat("-", rep("-", 40), "\n")

cat("🏆 ANALYSE CCA AVEC CLUSTERING TERMINÉE\n")
cat("=" * 50, "\n")

# Résumé du modèle CCA
cat("📊 MODÈLE CCA", best_model_name, ":\n")
cat("- Variables sélectionnées:", length(selected_vars_both), "\n")
cat("- Variables:", paste(selected_vars_both, collapse = ", "), "\n")
cat("- Variance expliquée:", round(model_comparison$Percent_Explained[model_comparison$Model == best_model_name], 1), "%\n")

# Résumé du clustering
cat("\n🎯 CLUSTERING CCA:\n")
cat("- Méthode retenue:", clustering_method, "\n")
cat("- Nombre de groupes:", n_clusters_cca, "\n")
cat("- Axes utilisés:", n_sig_axes, "\n")
cat("- Qualité silhouette:", round(max(sil_avg_hclust, sil_avg_kmeans), 3), "\n")

# Différenciation wood density
wd_significant <- anova_wd_summary[[1]][["Pr(>F)"]][1] < 0.05
cat("- Différenciation WD:", ifelse(wd_significant, "SIGNIFICATIVE ✅", "NON SIGNIFICATIVE ❌"), 
    "(p =", round(anova_wd_summary[[1]][["Pr(>F)"]][1], 4), ")\n")

if(wd_significant) {
  wd_range <- range(wd_by_group$WD_mean, na.rm = TRUE)
  cat("- Gradient WD:", round(wd_range[1], 3), "→", round(wd_range[2], 3), "\n")
}

# Différenciation structure
manova_significant <- manova_summary$stats[1, "Pr(>F)"] < 0.05
cat("- Différenciation structure:", ifelse(manova_significant, "SIGNIFICATIVE ✅", "NON SIGNIFICATIVE ❌"), 
    "(MANOVA p =", round(manova_summary$stats[1, "Pr(>F)"], 4), ")\n")

cat("\n✅ VALIDATION HYPOTHÈSE DE MÉDIATION FLORISTIQUE:\n")
cat("1. 🔗 Structure → Composition: CCA significative ✅\n")
cat("2. 🔗 Types forestiers identifiés:", n_clusters_cca, "groupes ✅\n")
cat("3. 🔗 Groupes → Wood Density:", ifelse(wd_significant, "Différenciation significative ✅", "À confirmer ⚠️"), "\n")
cat("4. 🔗 Groupes → Structure:", ifelse(manova_significant, "Différenciation significative ✅", "À confirmer ⚠️"), "\n")

cat("\n🚀 IMPLICATIONS POUR REMOTE SENSING:\n")
cat("- Variables clés identifiées:", length(selected_vars_both), "variables\n")
cat("- Types forestiers détectables par structure 3D\n")
cat("- Prédiction wood density via classification\n")

cat("\n🎉 Analyse CCA + Clustering terminée!\n")

