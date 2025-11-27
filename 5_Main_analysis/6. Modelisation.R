#### 1️⃣ INITIALIZATION ####
# Clean environment and memory
rm(list = ls())
gc()

# Définition des packages nécessaires
pkgs <- c("tidyverse","FactoMineR", "factoextra", "corrplot", "patchwork", "ade4", "vegan", "broom", "here", "caret")

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
path_data <- file.path(path_output, "data/modelisation")
path_figures <- file.path(path_output, "figures/modelisation")
path_table <- file.path(path_output, "tables/modelisation")

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
gap_metrics_to_remove <- c("basal_thickness", "transition_thickness", "upper_thickness", 
                           "auc_0_25_norm", "auc_25_75_norm", "auc_75_100_norm",
                           "basal_auc_norm", "upper_auc_norm", "transition_auc_norm",
                           "basal_auc", "upper_auc", "transition_auc", 	
                           "height_min_d2", "min_d2", "height_max_d2", "max_d2")

gap_metrics_cols <- setdiff(names(plot_data_complete), c(metadata_cols, gap_metrics_to_remove))

# Metadata with wood density and succession
metadata <- plot_data_complete %>%
  select(all_of(metadata_cols))

plots_name <- metadata %>%
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

#### 5️⃣ MULTIVARIATE ANALYSIS (NSCA, CA, PCA) ####

# A. NSCA - Analyse de Correspondance Non-Symétrique

# Effectuer la NSCA sur les données floristiques
nsca_result <- dudi.nsc(species_plot_table, scannf = FALSE, nf = min(nrow(species_plot_table) - 1, ncol(species_plot_table) - 1))
nsca_dbh_result <- dudi.nsc(species_dbh_plot_table, scannf = FALSE, nf = min(nrow(species_dbh_plot_table) - 1, ncol(species_dbh_plot_table) - 1))
nsca_genus_result <- dudi.nsc(genus_plot_table, scannf = FALSE, nf = min(nrow(genus_plot_table) - 1, ncol(genus_plot_table) - 1))

# Scores des parcelles sur les axes NSCA
nsca_site_scores <- nsca_result$li
nsca_site_scores_dbh <- nsca_dbh_result$li
nsca_site_scores_genus <- nsca_genus_result$li

colnames(nsca_site_scores) <- paste0("NSCA", 1:ncol(nsca_site_scores))
colnames(nsca_site_scores_dbh) <- paste0("NSCA_DBH", 1:ncol(nsca_site_scores_dbh))
colnames(nsca_site_scores_genus) <- paste0("NSCA_GENUS", 1:ncol(nsca_site_scores_genus))

print(summary(nsca_result))
print(summary(nsca_dbh_result))
print(summary(nsca_genus_result))

# B. CA - Correspondance Analysis

# Effectuer la CA sur les données floristiques
ca_result <- dudi.coa(species_plot_table, scannf = F, nf = min(nrow(species_plot_table) - 1, ncol(species_plot_table) - 1))
ca_genus_result <- dudi.coa(genus_plot_table, scannf = F, nf = min(nrow(genus_plot_table) - 1, ncol(genus_plot_table) - 1))

# Scores des parcelles sur les axes CA
ca_site_scores <- ca_result$li
ca_site_scores_genus <- ca_genus_result$li

colnames(ca_site_scores) <- paste0("CA", 1:ncol(ca_site_scores))
colnames(ca_site_scores_genus) <- paste0("CA_GENUS", 1:ncol(ca_site_scores_genus))

print(summary(ca_result))
print(summary(ca_genus_result))

# C. PCA - Principal Component Analysis

# Réaliser la PCA sur les données de structure de canopée
pca_result <- dudi.pca(structure_plots_table, scannf = F, nf = min(nrow(structure_plots_table) - 1, ncol(structure_plots_table) - 1))

# Scores des parcelles sur les axes CA
pca_site_scores <- pca_result$li

colnames(pca_site_scores) <- paste0("PCA", 1:ncol(pca_site_scores))

print(summary(pca_result))

fviz_eig(pca_result)

fviz_pca_var(pca_result,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


# =============================================================================
# PARTIE 6 : SÉLECTION DE MÉTRIQUES ET MODÉLISATION WD ~ STRUCTURE
# =============================================================================

# Cette section poursuit l'analyse après la partie 5 (analyses multivariées)
# Objectifs :
# A. Sélection des métriques de structure les plus pertinentes
# B. Modèles prédictifs WD ~ Structure (objectif télédétection)
# C. Décomposition du rôle de la composition (partition de variance, path analysis)

#### 6️⃣ SÉLECTION DES MÉTRIQUES DE STRUCTURE ####

# Packages supplémentaires pour cette section
pkgs_model <- c("car", "vegan", "performance", "see", 
                "ggcorrplot", "lme4", "spaMM", "MuMIn")
to_install <- !pkgs_model %in% installed.packages()
if(any(to_install)) {renv::install(pkgs_model[to_install])}
inst <- lapply(pkgs_model, library, character.only = TRUE)

# =============================================================================
# 6A. Exploration des corrélations Structure - WD
# =============================================================================

# Identifier les colonnes de métriques de structure disponibles
str_cols_all <- names(structure_plots_table)
cat("Métriques de structure disponibles:", length(str_cols_all), "\n")
print(str_cols_all)

# Calculer les corrélations avec WD_BA
wd_data <- plot_data_complete %>% 
  select(plot_name, WD_BA)

structure_with_wd <- structure_plots_table %>%
  rownames_to_column("plot_name") %>%
  left_join(wd_data, by = "plot_name")

# Corrélations de Pearson avec WD
cor_with_wd <- structure_with_wd %>%
  select(-plot_name) %>%
  cor(use = "pairwise.complete.obs") %>%
  as.data.frame() %>%
  select(WD_BA) %>%
  rownames_to_column("metric") %>%
  filter(metric != "WD_BA") %>%
  arrange(desc(abs(WD_BA)))

cat("\n=== TOP 15 métriques corrélées avec WD_BA ===\n")
print(head(cor_with_wd, 15))

# Visualisation
p_cor_wd <- cor_with_wd %>%
  head(20) %>%
  mutate(metric = fct_reorder(metric, abs(WD_BA))) %>%
  ggplot(aes(x = metric, y = WD_BA, fill = WD_BA)) +
  geom_col() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  coord_flip() +
  labs(title = "Corrélation des métriques de structure avec WD_BA",
       x = NULL, y = "Corrélation de Pearson") +
  theme_minimal()

ggsave(file.path(path_figures, "correlation_structure_wd.jpg"), p_cor_wd,
       width = 10, height = 8, dpi = 300)

# =============================================================================
# 6B. Contributions PCA et sélection
# =============================================================================

# Récupérer les contributions des variables à PC1 et PC2
pca_contrib <- get_pca_var(pca_result)$contrib %>%
  as.data.frame() %>%
  rownames_to_column("metric")

# Corrélations des variables avec les axes
pca_cor <- get_pca_var(pca_result)$cor %>%
  as.data.frame() %>%
  rownames_to_column("metric") %>%
  rename(cor_PC1 = Dim.1, cor_PC2 = Dim.2)

# Combiner avec corrélations WD
selection_table <- cor_with_wd %>%
  rename(cor_WD = WD_BA) %>%
  left_join(pca_cor %>% select(metric, cor_PC1, cor_PC2), by = "metric") %>%
  left_join(pca_contrib %>% select(metric, Dim.1, Dim.2) %>% 
              rename(contrib_PC1 = Dim.1, contrib_PC2 = Dim.2), by = "metric") %>%
  mutate(abs_cor_WD = abs(cor_WD),
         total_contrib = contrib_PC1 + contrib_PC2)

cat("\n=== Table de sélection des métriques ===\n")
print(selection_table %>% arrange(desc(abs_cor_WD)) %>% head(15))

# =============================================================================
# 6C. Analyse de multicolinéarité et sélection finale
# =============================================================================

# Candidats potentiels basés sur interprétabilité et corrélation
# Groupes de métriques (éviter de prendre plusieurs du même groupe):
# - Groupe 1: Proportions d'ouverture (prop_at_Xm) 
# - Groupe 2: Percentiles de hauteur (hXX)
# - Groupe 3: Pentes (slope_*)
# - Groupe 4: Aires sous courbe (auc_*)

# Sélection manuelle de candidats représentatifs
candidates <- c("prop_at_20m", "prop_at_25m", "prop_at_30m",
                "h25", "h50", "h75", "h90",
                "auc_norm", "max_d1")

# Filtrer pour ne garder que les colonnes existantes
candidates <- candidates[candidates %in% str_cols_all]
cat("\nCandidats retenus:", paste(candidates, collapse = ", "), "\n")

# Matrice de corrélation entre candidats
cor_candidates <- structure_with_wd %>%
  select(all_of(candidates)) %>%
  cor(use = "pairwise.complete.obs")

# Visualisation
p_cor_matrix <- ggcorrplot(cor_candidates, 
                           type = "lower",
                           lab = TRUE,
                           lab_size = 3,
                           colors = c("blue", "white", "red"),
                           title = "Corrélations entre métriques candidates")

ggsave(file.path(path_figures, "correlation_matrix_candidates.jpg"), p_cor_matrix,
       width = 10, height = 8, dpi = 300)

# Test VIF pour une sélection de 3-4 métriques
# On commence par un modèle avec les meilleures candidates non redondantes

# Sélection initiale basée sur corrélation et groupes différents
selected_metrics <- c("prop_at_20m", "max_d1", "auc_norm")

# Vérifier que ces colonnes existent
selected_metrics <- selected_metrics[selected_metrics %in% str_cols_all]

if(length(selected_metrics) >= 2) {
  # Créer formule
  formula_vif <- as.formula(paste("WD_BA ~", paste(selected_metrics, collapse = " + ")))
  
  # Modèle pour VIF
  model_vif <- lm(formula_vif, data = structure_with_wd)
  
  cat("\n=== VIF pour les métriques sélectionnées ===\n")
  print(car::vif(model_vif))
  
  # Critère: VIF < 5 acceptable, < 10 tolérable
}

# =============================================================================
# 6D. Corrélation WD avec scores PCA
# =============================================================================

# Les scores PCA sont déjà calculés dans pca_site_scores
pca_wd_data <- data.frame(
  plot_name = rownames(pca_site_scores),
  PCA1 = pca_site_scores$PCA1,
  PCA2 = pca_site_scores$PCA2
) %>%
  left_join(wd_data, by = "plot_name")

cor_pca_wd <- cor.test(pca_wd_data$PCA1, pca_wd_data$WD_BA)
cat("\n=== Corrélation PCA1 - WD_BA ===\n")
cat("r =", round(cor_pca_wd$estimate, 3), ", p =", format.pval(cor_pca_wd$p.value, 3), "\n")

cor_pca2_wd <- cor.test(pca_wd_data$PCA2, pca_wd_data$WD_BA)
cat("r (PCA2) =", round(cor_pca2_wd$estimate, 3), ", p =", format.pval(cor_pca2_wd$p.value, 3), "\n")

#### 7️⃣ MODÈLES PRÉDICTIFS WD ~ STRUCTURE (SANS AUTOCORRÉLATION) ####

# =============================================================================
# 7A. Préparation des données pour modélisation
# =============================================================================

# Créer un dataframe complet pour la modélisation
model_data <- structure_with_wd %>%
  left_join(
    data.frame(
      plot_name = rownames(pca_site_scores),
      PCA1 = pca_site_scores$PCA1,
      PCA2 = pca_site_scores$PCA2
    ), by = "plot_name"
  ) %>%
  left_join(
    data.frame(
      plot_name = rownames(ca_site_scores),
      CA1 = ca_site_scores$CA1,
      CA2 = ca_site_scores$CA2
    ), by = "plot_name"
  ) %>%
  left_join(
    metadata %>% select(plot_name, prop_g_helio, prop_g_npld, prop_g_shade),
    by = "plot_name"
  )

# Vérifier les données complètes
cat("\n=== Dimensions du jeu de données ===\n")
cat("Observations:", nrow(model_data), "\n")
cat("Variables:", ncol(model_data), "\n")
cat("Données complètes:", sum(complete.cases(model_data)), "\n")

# =============================================================================
# 7B. Modèles simples : WD ~ métrique unique
# =============================================================================

# Tester chaque métrique individuellement
simple_models <- list()
simple_results <- data.frame()

for(metric in candidates) {
  if(metric %in% names(model_data)) {
    formula_i <- as.formula(paste("WD_BA ~", metric))
    model_i <- lm(formula_i, data = model_data)
    simple_models[[metric]] <- model_i
    
    # Extraire statistiques
    summary_i <- summary(model_i)
    simple_results <- rbind(simple_results, data.frame(
      metric = metric,
      R2 = summary_i$r.squared,
      R2_adj = summary_i$adj.r.squared,
      coef = coef(model_i)[2],
      pvalue = summary_i$coefficients[2, 4],
      AIC = AIC(model_i)
    ))
  }
}

simple_results <- simple_results %>% arrange(desc(R2_adj))
cat("\n=== Modèles simples WD ~ métrique unique ===\n")
print(simple_results)

# Sauvegarder le tableau
write_csv2(simple_results, file.path(path_table, "simple_models_results.csv"))

# =============================================================================
# 7C. Modèles multiples : combinaisons de métriques
# =============================================================================

# Modèle avec les 3 métriques sélectionnées
if(all(selected_metrics %in% names(model_data))) {
  formula_multi <- as.formula(paste("WD_BA ~", paste(selected_metrics, collapse = " + ")))
  model_multi <- lm(formula_multi, data = model_data)
  
  cat("\n=== Modèle multiple: WD ~ ", paste(selected_metrics, collapse = " + "), " ===\n")
  print(summary(model_multi))
  
  cat("\nVIF:\n")
  print(car::vif(model_multi))
}

# Modèle avec PCA1 (synthèse de la structure)
model_pca1 <- lm(WD_BA ~ PCA1, data = model_data)
cat("\n=== Modèle WD ~ PCA1 ===\n")
print(summary(model_pca1))

# Modèle avec PCA1 + PCA2
model_pca12 <- lm(WD_BA ~ PCA1 + PCA2, data = model_data)
cat("\n=== Modèle WD ~ PCA1 + PCA2 ===\n")
print(summary(model_pca12))

# =============================================================================
# 7D. Comparaison des modèles de structure
# =============================================================================

# Liste des modèles à comparer
structure_models <- list(
  "Meilleur_simple" = simple_models[[simple_results$metric[1]]],
  "Multi_3var" = model_multi,
  "PCA1" = model_pca1,
  "PCA1_PCA2" = model_pca12
)

# Tableau comparatif
comparison_structure <- data.frame(
  Modele = names(structure_models),
  R2_adj = sapply(structure_models, function(m) summary(m)$adj.r.squared),
  AIC = sapply(structure_models, AIC),
  BIC = sapply(structure_models, BIC),
  RMSE = sapply(structure_models, function(m) sqrt(mean(residuals(m)^2)))
) %>% arrange(AIC)

cat("\n=== Comparaison des modèles Structure -> WD ===\n")
print(comparison_structure)

write_csv2(comparison_structure, file.path(path_table, "comparison_structure_models.csv"))

#### 8️⃣ RÔLE DE LA COMPOSITION : PARTITION DE VARIANCE ####

# =============================================================================
# 8A. Modèles avec composition seule
# =============================================================================

# Modèle avec proportions de tempérament
model_temperament <- lm(WD_BA ~ prop_g_helio + prop_g_shade, data = model_data)
cat("\n=== Modèle WD ~ Tempérament (hélio + shade) ===\n")
print(summary(model_temperament))

# Modèle avec axes CA
model_ca <- lm(WD_BA ~ CA1 + CA2, data = model_data)
cat("\n=== Modèle WD ~ CA1 + CA2 ===\n")
print(summary(model_ca))

# =============================================================================
# 8B. Partition de variance (varpart)
# =============================================================================

# Préparer les matrices pour varpart
# X1 = Structure (métriques sélectionnées ou PCA)
# X2 = Composition (tempérament ou CA)

# Nettoyer les données (enlever NA)
varpart_data <- model_data %>%
  select(WD_BA, all_of(selected_metrics), prop_g_helio, prop_g_shade, CA1, CA2, PCA1, PCA2) %>%
  drop_na()

cat("\n=== Données pour varpart: n =", nrow(varpart_data), "===\n")

# Option 1: Structure (PCA1) vs Composition (tempérament)
if(nrow(varpart_data) >= 10) {
  vp1 <- varpart(varpart_data$WD_BA, 
                 ~ PCA1, 
                 ~ prop_g_helio + prop_g_shade, 
                 data = varpart_data)
  
  cat("\n=== Partition: Structure (PCA1) vs Tempérament ===\n")
  print(vp1)
  
  # Visualisation
  jpeg(file.path(path_figures, "varpart_structure_temperament.jpg"), 
       width = 800, height = 600, quality = 100)
  plot(vp1, digits = 2, Xnames = c("Structure\n(PCA1)", "Tempérament"),
       bg = c("lightblue", "lightgreen"))
  dev.off()
}

# Option 2: Structure (PCA1) vs Composition (CA)
if(nrow(varpart_data) >= 10) {
  vp2 <- varpart(varpart_data$WD_BA, 
                 ~ PCA1, 
                 ~ CA1 + CA2, 
                 data = varpart_data)
  
  cat("\n=== Partition: Structure (PCA1) vs CA ===\n")
  print(vp2)
  
  jpeg(file.path(path_figures, "varpart_structure_ca.jpg"), 
       width = 800, height = 600, quality = 100)
  plot(vp2, digits = 2, Xnames = c("Structure\n(PCA1)", "Composition\n(CA)"),
       bg = c("lightblue", "lightgreen"))
  dev.off()
}

# Option 3: Structure (métriques) vs Composition (tempérament)
if(nrow(varpart_data) >= 10 && all(selected_metrics %in% names(varpart_data))) {
  formula_struct <- as.formula(paste("~", paste(selected_metrics, collapse = " + ")))
  
  vp3 <- varpart(varpart_data$WD_BA, 
                 formula_struct, 
                 ~ prop_g_helio + prop_g_shade, 
                 data = varpart_data)
  
  cat("\n=== Partition: Structure (métriques) vs Tempérament ===\n")
  print(vp3)
  
  jpeg(file.path(path_figures, "varpart_metrics_temperament.jpg"), 
       width = 800, height = 600, quality = 100)
  plot(vp3, digits = 2, Xnames = c("Structure\n(métriques)", "Tempérament"),
       bg = c("lightblue", "lightgreen"))
  dev.off()
}

# =============================================================================
# 8C. Modèles combinés Structure + Composition
# =============================================================================

# Modèle complet : Structure + Tempérament
model_full_temp <- lm(WD_BA ~ PCA1 + prop_g_helio + prop_g_shade, data = model_data)
cat("\n=== Modèle complet: WD ~ PCA1 + Tempérament ===\n")
print(summary(model_full_temp))

# Modèle complet : Structure + CA
model_full_ca <- lm(WD_BA ~ PCA1 + CA1 + CA2, data = model_data)
cat("\n=== Modèle complet: WD ~ PCA1 + CA ===\n")
print(summary(model_full_ca))

# =============================================================================
# 8D. Comparaison : Structure seule vs Structure + Composition
# =============================================================================

# Tableau comparatif
comparison_all <- data.frame(
  Modele = c("Structure (PCA1)", "Tempérament", "CA", 
             "Structure + Tempérament", "Structure + CA"),
  Formule = c("WD ~ PCA1", "WD ~ helio + shade", "WD ~ CA1 + CA2",
              "WD ~ PCA1 + helio + shade", "WD ~ PCA1 + CA1 + CA2"),
  R2_adj = c(
    summary(model_pca1)$adj.r.squared,
    summary(model_temperament)$adj.r.squared,
    summary(model_ca)$adj.r.squared,
    summary(model_full_temp)$adj.r.squared,
    summary(model_full_ca)$adj.r.squared
  ),
  AIC = c(
    AIC(model_pca1),
    AIC(model_temperament),
    AIC(model_ca),
    AIC(model_full_temp),
    AIC(model_full_ca)
  )
) %>% arrange(AIC)

cat("\n=== Comparaison tous modèles ===\n")
print(comparison_all)

write_csv2(comparison_all, file.path(path_table, "comparison_all_models.csv"))

# Test: Est-ce que la composition améliore significativement le modèle?
cat("\n=== Test amélioration par ajout de la composition ===\n")
anova_temp <- anova(model_pca1, model_full_temp)
print(anova_temp)

anova_ca <- anova(model_pca1, model_full_ca)
print(anova_ca)

#### 9️⃣ SYNTHÈSE ET VISUALISATIONS ####

# =============================================================================
# 9A. Graphique principal : WD vs meilleur prédicteur structurel
# =============================================================================

best_metric <- simple_results$metric[1]

p_main <- ggplot(model_data, aes_string(x = best_metric, y = "WD_BA")) +
  geom_point(aes(color = prop_g_helio), size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_viridis_c(name = "Prop.\nhéliophiles", option = "plasma") +
  labs(
    title = paste("Relation WD - Structure forestière"),
    subtitle = paste0("R² = ", round(simple_results$R2[1], 3)),
    x = best_metric,
    y = "Densité du bois (g/cm³)"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(file.path(path_figures, "main_wd_structure_relation.jpg"), p_main,
       width = 10, height = 8, dpi = 300)

# =============================================================================
# 9B. Graphique WD vs PCA1 coloré par composition ----
# =============================================================================

p_pca_wd <- ggplot(model_data, aes(x = PCA1, y = WD_BA)) +
  geom_point(aes(color = prop_g_helio), size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_viridis_c(name = "Prop.\nhéliophiles", option = "plasma") +
  labs(
    title = "Densité du bois vs Structure de canopée (PCA1)",
    subtitle = paste0("R² = ", round(summary(model_pca1)$r.squared, 3)),
    x = "PCA1 (gradient structural)",
    y = "Densité du bois (g/cm³)"
  ) +
  theme_minimal()

ggsave(file.path(path_figures, "wd_vs_pca1.jpg"), p_pca_wd,
       width = 10, height = 8, dpi = 300)

# =============================================================================
# 9C. Résumé des résultats
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("                    SYNTHÈSE DES RÉSULTATS                      \n")
cat("================================================================\n")
cat("\n")
cat("1. MEILLEURE MÉTRIQUE SIMPLE:\n")
cat("   ", best_metric, "- R² =", round(simple_results$R2[1], 3), "\n")
cat("\n")
cat("2. MODÈLE STRUCTURE (PCA1):\n")
cat("   R² =", round(summary(model_pca1)$r.squared, 3), "\n")
cat("\n")
cat("3. PARTITION DE VARIANCE:\n")
cat("   (voir figures varpart_*.jpg)\n")
cat("\n")
cat("4. AMÉLIORATION PAR LA COMPOSITION:\n")
cat("   Structure seule:        R² =", round(summary(model_pca1)$adj.r.squared, 3), "\n")
cat("   Structure + Tempérament: R² =", round(summary(model_full_temp)$adj.r.squared, 3), "\n")
cat("   Structure + CA:          R² =", round(summary(model_full_ca)$adj.r.squared, 3), "\n")
cat("\n")
cat("================================================================\n")

# Sauvegarder l'environnement pour analyses ultérieures
save.image(file.path(path_data, "modelisation_results.RData"))

cat("\nAnalyse terminée. Résultats sauvegardés.\n")

