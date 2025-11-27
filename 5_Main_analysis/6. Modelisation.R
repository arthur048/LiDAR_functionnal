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
gaps_metrics <- read_csv2(file.path(path_input, "gaps_metrics.csv"))
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


#### 6️⃣ NOUVELLES MÉTRIQUES DE COURBE ----

# Packages supplémentaires
pkgs_new <- c("glmnet", "vegan", "lavaan", "spaMM", "sf", "nls2", "e1071")
to_install <- !pkgs_new %in% installed.packages()
if(any(to_install)) {install.packages(pkgs_new[to_install])}
inst <- lapply(pkgs_new, library, character.only = TRUE)

# 6.0 Création de model_data ----

# Créer model_data en joignant structure + WD + composition
model_data <- structure_plots_table %>%
  rownames_to_column("plot_name") %>%
  left_join(plot_metrics %>% select(plot_name, WD_BA, 
                                    prop_g_helio, prop_g_npld, prop_g_shade),
            by = "plot_name") %>%
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
  )

# Ajouter NSCA_DBH si disponible
if(exists("nsca_site_scores_dbh")) {
  model_data <- model_data %>%
    left_join(
      data.frame(
        plot_name = rownames(nsca_site_scores_dbh),
        NSCA_DBH1 = nsca_site_scores_dbh[,1],
        NSCA_DBH2 = nsca_site_scores_dbh[,2]
      ), by = "plot_name"
    )
}

cat("model_data créé: n =", nrow(model_data), "variables =", ncol(model_data), "\n")

# 6A. Fonction pour ajuster une sigmoïde et extraire les paramètres ----

fit_sigmoid <- function(height, proportion, plot_name = NA) {
  # Ajustement d'une fonction logistique: y = 1 / (1 + exp(-k*(x-x0)))
  # où k = pente au point d'inflexion, x0 = hauteur au point d'inflexion
  
  data_fit <- data.frame(x = height, y = proportion) %>%
    filter(!is.na(y) & !is.na(x))
  
  # Valeurs initiales estimées
  x0_init <- data_fit$x[which.min(abs(data_fit$y - 0.5))]
  k_init <- 0.1
  
  result <- tryCatch({
    # Ajustement NLS
    fit <- nls(y ~ 1 / (1 + exp(-k * (x - x0))),
               data = data_fit,
               start = list(k = k_init, x0 = x0_init),
               control = nls.control(maxiter = 100, warnOnly = TRUE))
    
    coefs <- coef(fit)
    residuals_fit <- residuals(fit)
    
    data.frame(
      sigmoid_k = coefs["k"],           # Pente (raideur de la transition)
      sigmoid_x0 = coefs["x0"],         # Point d'inflexion (hauteur 50%)
      sigmoid_rmse = sqrt(mean(residuals_fit^2)),  # Qualité ajustement
      sigmoid_residual_var = var(residuals_fit)    # Rugosité/hétérogénéité
    )
  }, error = function(e) {
    data.frame(
      sigmoid_k = NA_real_,
      sigmoid_x0 = NA_real_,
      sigmoid_rmse = NA_real_,
      sigmoid_residual_var = NA_real_
    )
  })
  
  return(result)
}

# 6B. Fonction pour calculer les nouvelles métriques ----

compute_new_metrics <- function(height, proportion) {
  # Largeur de transition (différence entre h75 et h25)
  h25 <- approx(proportion, height, xout = 0.25)$y
  h75 <- approx(proportion, height, xout = 0.75)$y
  transition_width <- h75 - h25
  
  # Ratio h25/h75 (proportionnalité des strates)
  ratio_h25_h75 <- ifelse(!is.na(h75) & h75 > 0, h25 / h75, NA)
  
  # Asymétrie (skewness) de la distribution de la dérivée
  # La dérivée représente la "densité" de fermeture par strate
  deriv <- diff(proportion) / diff(height)
  deriv_skewness <- ifelse(length(deriv) > 2, e1071::skewness(deriv, na.rm = TRUE), NA)
  
  # Kurtosis de la dérivée (forme de la distribution)
  deriv_kurtosis <- ifelse(length(deriv) > 2, e1071::kurtosis(deriv, na.rm = TRUE), NA)
  
  # Position relative du point d'inflexion par rapport à h50
  # (déjà calculé via max_d1 et height_max_d1, mais on peut normaliser)
  
  data.frame(
    transition_width = transition_width,
    ratio_h25_h75 = ratio_h25_h75,
    deriv_skewness = deriv_skewness,
    deriv_kurtosis = deriv_kurtosis
  )
}

# 6C. Appliquer les nouvelles métriques à toutes les parcelles ----

# Récupérer les données de courbes depuis gaps_metrics (buffer = 0)
curves_data <- gaps_metrics %>%
  filter(buffer == 0) %>%
  select(plot_name, height_aboveground, proportion) %>%
  arrange(plot_name, height_aboveground)

# Calculer les nouvelles métriques par parcelle
new_metrics_list <- list()

for(plot_i in unique(curves_data$plot_name)) {
  plot_curve <- curves_data %>% filter(plot_name == plot_i)
  
  # Métriques sigmoïdes
  sigmoid_metrics <- fit_sigmoid(plot_curve$height_aboveground, 
                                 plot_curve$proportion, 
                                 plot_i)
  
  # Autres nouvelles métriques
  other_metrics <- compute_new_metrics(plot_curve$height_aboveground,
                                       plot_curve$proportion)
  
  new_metrics_list[[plot_i]] <- cbind(
    data.frame(plot_name = plot_i),
    sigmoid_metrics,
    other_metrics
  )
}

new_metrics <- bind_rows(new_metrics_list)

cat("=== Nouvelles métriques calculées ===\n")
print(summary(new_metrics))

# Joindre aux données existantes
model_data <- model_data %>%
  left_join(new_metrics, by = "plot_name")

# Corrélations des nouvelles métriques avec WD
new_metric_cols <- c("sigmoid_k", "sigmoid_x0", "transition_width", 
                     "ratio_h25_h75", "deriv_skewness", "deriv_kurtosis")

cor_new_metrics <- model_data %>%
  select(WD_BA, all_of(new_metric_cols)) %>%
  cor(use = "pairwise.complete.obs")

cat("\n=== Corrélations nouvelles métriques - WD_BA ===\n")
print(round(cor_new_metrics[, "WD_BA"], 3))

#### 7️⃣ SÉLECTION DE VARIABLES PAR ELASTIC NET ----

# 7A. Préparation des données pour glmnet ----

# Sélectionner toutes les métriques de structure potentielles
all_structure_metrics <- c(
  # Métriques existantes
  "prop_at_5m", "prop_at_10m", "prop_at_15m", "prop_at_20m", "prop_at_25m",
  "prop_at_30m", "prop_at_35m", "prop_at_40m",
  "h10", "h25", "h50", "h75", "h90", "h100",
  "slope_basal", "slope_transition", "slope_upper",
  "auc", "auc_norm", "auc_0_25", "auc_25_75", "auc_75_100",
  "max_d1", "height_max_d1",
  # Nouvelles métriques
  "sigmoid_k", "sigmoid_x0", "transition_width", "ratio_h25_h75",
  "deriv_skewness", "deriv_kurtosis"
)

# Filtrer pour ne garder que les colonnes existantes
available_metrics <- all_structure_metrics[all_structure_metrics %in% names(model_data)]
cat("\nMétriques disponibles pour Elastic Net:", length(available_metrics), "\n")

# Créer la matrice X et le vecteur y
elastic_data <- model_data %>%
  select(WD_BA, all_of(available_metrics)) %>%
  drop_na()

X <- as.matrix(elastic_data %>% select(-WD_BA))
y <- elastic_data$WD_BA

cat("Dimensions: n =", nrow(X), ", p =", ncol(X), "\n")

# 7B. Cross-validation pour trouver lambda optimal ----

set.seed(42)

# Alpha = 0.5 pour Elastic Net (compromis LASSO/Ridge)
# On teste aussi alpha = 1 (LASSO pur) et alpha = 0 (Ridge)

cv_results <- list()
alphas <- c(0, 0.5, 1)

for(alpha_i in alphas) {
  cv_fit <- cv.glmnet(X, y, alpha = alpha_i, nfolds = 10)
  cv_results[[as.character(alpha_i)]] <- cv_fit
  cat("\nAlpha =", alpha_i, "- Lambda.min =", round(cv_fit$lambda.min, 5),
      "- MSE =", round(min(cv_fit$cvm), 5), "\n")
}

# Sélectionner le meilleur alpha
best_alpha <- alphas[which.min(sapply(cv_results, function(x) min(x$cvm)))]
cat("\n=== Meilleur alpha:", best_alpha, "===\n")

cv_best <- cv_results[[as.character(best_alpha)]]

# 7C. Extraire les coefficients et importance des variables ----

# Coefficients au lambda optimal
coef_best <- coef(cv_best, s = "lambda.min")
coef_df <- data.frame(
  variable = rownames(coef_best),
  coefficient = as.vector(coef_best)
) %>%
  filter(variable != "(Intercept)" & coefficient != 0) %>%
  mutate(abs_coef = abs(coefficient)) %>%
  arrange(desc(abs_coef))

cat("\n=== Variables sélectionnées par Elastic Net ===\n")
print(coef_df)

# Visualisation
p_elastic <- coef_df %>%
  mutate(variable = fct_reorder(variable, abs_coef)) %>%
  ggplot(aes(x = variable, y = coefficient, fill = coefficient > 0)) +
  geom_col() +
  scale_fill_manual(values = c("red", "blue"), 
                    labels = c("Négatif", "Positif"),
                    name = "Signe") +
  coord_flip() +
  labs(title = "Coefficients Elastic Net (lambda.min)",
       subtitle = paste("Alpha =", best_alpha),
       x = NULL, y = "Coefficient") +
  theme_minimal()

ggsave(file.path(path_figures, "elastic_net_coefficients.jpg"), p_elastic,
       width = 10, height = 8, dpi = 300)

# 7D. Modèle final avec variables sélectionnées ----

selected_vars_elastic <- coef_df$variable
cat("\nVariables sélectionnées:", paste(selected_vars_elastic, collapse = ", "), "\n")

# Modèle OLS avec ces variables pour interprétation
if(length(selected_vars_elastic) > 0) {
  formula_elastic <- as.formula(paste("WD_BA ~", paste(selected_vars_elastic, collapse = " + ")))
  model_elastic <- lm(formula_elastic, data = model_data)
  
  cat("\n=== Modèle avec variables Elastic Net ===\n")
  print(summary(model_elastic))
}

# Sauvegarder résultats
write_csv2(coef_df, file.path(path_table, "elastic_net_selection.csv"))

#### 8️⃣ PERMANOVA - DIFFÉRENCES STRUCTURELLES PAR CATÉGORIE WD ----

# 8A. Créer les catégories de WD ----

# Quartiles
model_data <- model_data %>%
  mutate(
    WD_quartile = cut(WD_BA, 
                      breaks = quantile(WD_BA, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
                      labels = c("Q1_low", "Q2", "Q3", "Q4_high"),
                      include.lowest = TRUE),
    # Catégorie simplifiée (3 groupes)
    WD_category = cut(WD_BA,
                      breaks = quantile(WD_BA, probs = c(0, 0.33, 0.66, 1), na.rm = TRUE),
                      labels = c("Low", "Medium", "High"),
                      include.lowest = TRUE)
  )

# 8B. Préparer la matrice de distance structurelle ----

# Variables structurelles pour PERMANOVA
permanova_vars <- available_metrics[available_metrics %in% names(model_data)]

# Enlever les variables avec trop de NA
permanova_data <- model_data %>%
  select(plot_name, WD_quartile, WD_category, all_of(permanova_vars)) %>%
  drop_na()

# Matrice de structure (standardisée)
struct_matrix <- permanova_data %>%
  select(all_of(permanova_vars)) %>%
  scale() %>%
  as.data.frame()

# Distance euclidienne
dist_struct <- dist(struct_matrix, method = "euclidean")

# 8C. PERMANOVA ----

set.seed(42)

# Test par quartiles
permanova_quartile <- adonis2(dist_struct ~ WD_quartile, 
                              data = permanova_data, 
                              permutations = 999)

cat("\n=== PERMANOVA: Structure ~ WD quartiles ===\n")
print(permanova_quartile)

# Test par catégorie (3 groupes)
permanova_category <- adonis2(dist_struct ~ WD_category, 
                              data = permanova_data, 
                              permutations = 999)

cat("\n=== PERMANOVA: Structure ~ WD catégorie ===\n")
print(permanova_category)

# 8D. Analyse de dispersion (betadisper) ----

# Vérifier l'homogénéité des dispersions
betadisp_quartile <- betadisper(dist_struct, permanova_data$WD_quartile)
permutest_disp <- permutest(betadisp_quartile, permutations = 999)

cat("\n=== Test d'homogénéité des dispersions ===\n")
print(permutest_disp)

# 8E. Pairwise comparisons (si significatif) ----

if(permanova_quartile$`Pr(>F)`[1] < 0.05) {
  # Comparaisons post-hoc simplifiées
  cat("\n=== Comparaisons par paires (à implémenter si besoin) ===\n")
  cat("PERMANOVA significatif - les groupes de WD diffèrent en structure\n")
}

# Sauvegarder résultats PERMANOVA
permanova_results <- data.frame(
  Test = c("WD_quartiles", "WD_category"),
  F_value = c(permanova_quartile$F[1], permanova_category$F[1]),
  R2 = c(permanova_quartile$R2[1], permanova_category$R2[1]),
  p_value = c(permanova_quartile$`Pr(>F)`[1], permanova_category$`Pr(>F)`[1])
)
write_csv2(permanova_results, file.path(path_table, "permanova_results.csv"))

#### 9️⃣ CHARGEMENT DES COORDONNÉES SPATIALES ----
# 9A. Charger les coordonnées depuis les gpkg ----

# Chemin vers les fichiers de plots
path_plots_gpkg <- "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/final/plots_unique"

# Lister tous les gpkg
gpkg_files <- list.files(path_plots_gpkg, pattern = "\\.gpkg$", full.names = TRUE)

# Fonction pour extraire le centroïde d'un gpkg
extract_centroid <- function(gpkg_path) {
  tryCatch({
    plot_sf <- st_read(gpkg_path, quiet = TRUE)
    centroid <- st_centroid(st_union(plot_sf))
    coords <- st_coordinates(centroid)
    
    # Extraire le nom du plot depuis le nom du fichier
    plot_name <- tools::file_path_sans_ext(basename(gpkg_path))
    
    data.frame(
      plot_name = plot_name,
      X_utm = coords[1],
      Y_utm = coords[2],
      epsg = st_crs(plot_sf)$epsg
    )
  }, error = function(e) {
    NULL
  })
}

# Extraire tous les centroïdes
cat("Extraction des coordonnées spatiales...\n")
coords_list <- lapply(gpkg_files, extract_centroid)
coords_df <- bind_rows(coords_list)

cat("Coordonnées extraites pour", nrow(coords_df), "plots\n")

# 9B. Joindre aux données de modélisation ----

# Filtrer pour ne garder que les plots de notre analyse
model_data <- model_data %>%
  left_join(coords_df, by = "plot_name")

# Vérifier les correspondances
n_with_coords <- sum(!is.na(model_data$X_utm))
cat("Plots avec coordonnées:", n_with_coords, "/", nrow(model_data), "\n")

# 9C. Transformer en coordonnées communes si nécessaire ----

# Si plusieurs EPSG, il faudra reprojeter
unique_epsg <- unique(model_data$epsg[!is.na(model_data$epsg)])
cat("EPSG uniques:", paste(unique_epsg, collapse = ", "), "\n")

# Pour spaMM, on peut utiliser les coordonnées UTM directement
# (même si EPSG différents, l'échelle est similaire en mètres)

#### 🔟 MODÈLES AVEC AUTOCORRÉLATION SPATIALE (spaMM) ----

# 10A. Modèle de base sans effet spatial ----

# Données complètes avec coordonnées
spatial_data <- model_data %>%
  filter(!is.na(X_utm) & !is.na(Y_utm) & !is.na(WD_BA))

cat("\nDonnées pour modèles spatiaux: n =", nrow(spatial_data), "\n")

# Meilleure variable simple (d'après analyses précédentes)
best_var <- coef_df$variable[1]  # Première variable Elastic Net
if(is.na(best_var)) best_var <- "prop_at_20m"  # Fallback

# Modèle OLS simple pour comparaison
formula_simple <- as.formula(paste("WD_BA ~", best_var))
model_ols <- lm(formula_simple, data = spatial_data)

cat("\n=== Modèle OLS (sans effet spatial) ===\n")
print(summary(model_ols))

# 10B. Modèle avec corrélation spatiale Matérn ----

# spaMM avec effet spatial gaussien
model_spatial <- fitme(
  formula = as.formula(paste("WD_BA ~", best_var, "+ Matern(1 | X_utm + Y_utm)")),
  data = spatial_data,
  method = "REML"
)

cat("\n=== Modèle spatial (Matérn) ===\n")
print(summary(model_spatial))

# 10C. Comparer AIC ----

aic_comparison <- data.frame(
  Modele = c("OLS", "Spatial_Matern"),
  AIC = c(AIC(model_ols), AIC(model_spatial)),
  logLik = c(logLik(model_ols), logLik(model_spatial))
)

cat("\n=== Comparaison AIC ===\n")
print(aic_comparison)

# 10D. Modèle spatial avec variables Elastic Net ----

if(length(selected_vars_elastic) > 0 && length(selected_vars_elastic) <= 5) {
  # Limiter à 5 variables max pour éviter surparamétrisation
  vars_to_use <- head(selected_vars_elastic, 5)
  formula_multi_spatial <- as.formula(
    paste("WD_BA ~", paste(vars_to_use, collapse = " + "), 
          "+ Matern(1 | X_utm + Y_utm)")
  )
  
  model_spatial_multi <- fitme(
    formula = formula_multi_spatial,
    data = spatial_data,
    method = "REML"
  )
  
  cat("\n=== Modèle spatial multivarié ===\n")
  print(summary(model_spatial_multi))
}

#### 1️⃣1️⃣ SEM AVEC LAVAAN - MODÈLE CAUSAL ----

# 11A. Préparer les données pour lavaan ----

# Variables pour le SEM
sem_data <- model_data %>%
  select(WD_BA, PCA1, CA1, CA2, prop_g_helio, prop_g_shade,
         NSCA_Genus_DBH1 = starts_with("NSCA_Genus_DBH1"),
         NSCA_Genus_DBH2 = starts_with("NSCA_Genus_DBH2")) %>%
  drop_na()

# Vérifier si NSCA_DBH existe dans les données
if(!"NSCA_Genus_DBH1" %in% names(sem_data)) {
  # Essayer d'autres noms possibles
  nsca_cols <- grep("NSCA.*DBH", names(model_data), value = TRUE)
  if(length(nsca_cols) >= 2) {
    sem_data <- model_data %>%
      select(WD_BA, PCA1, CA1, CA2, prop_g_helio, prop_g_shade,
             all_of(nsca_cols[1:2])) %>%
      drop_na()
    names(sem_data)[7:8] <- c("NSCA_DBH1", "NSCA_DBH2")
  }
}

cat("\nDonnées pour SEM: n =", nrow(sem_data), "\n")
cat("Variables:", paste(names(sem_data), collapse = ", "), "\n")

# 11B. Modèle 1: Composition (CA) comme cause commune ----

# Diagramme causal:
# CA1, CA2 --> PCA1 (structure)
# CA1, CA2 --> WD_BA
# PCA1 --> WD_BA (effet résiduel direct ?)

model_sem1 <- '
 # Composition prédit Structure
 PCA1 ~ CA1 + CA2
 
 # Composition prédit WD directement
 WD_BA ~ CA1 + CA2
 
 # Structure prédit WD (effet résiduel)
 WD_BA ~ PCA1
'

fit_sem1 <- sem(model_sem1, data = sem_data)

cat("\n=== SEM Modèle 1: CA comme cause commune ===\n")
summary(fit_sem1, standardized = TRUE, fit.measures = TRUE)

# 11C. Modèle 2: Tempérament comme cause commune ----

model_sem2 <- '
 # Tempérament prédit Structure
 PCA1 ~ prop_g_helio + prop_g_shade
 
 # Tempérament prédit WD directement  
 WD_BA ~ prop_g_helio + prop_g_shade
 
 # Structure prédit WD (effet résiduel)
 WD_BA ~ PCA1
'

fit_sem2 <- sem(model_sem2, data = sem_data)

cat("\n=== SEM Modèle 2: Tempérament comme cause commune ===\n")
summary(fit_sem2, standardized = TRUE, fit.measures = TRUE)

# 11D. Modèle 3: Sans effet direct Structure -> WD ----

model_sem3 <- '
 # Composition prédit Structure
 PCA1 ~ CA1 + CA2
 
 # Composition prédit WD directement (pas de lien Structure -> WD)
 WD_BA ~ CA1 + CA2
'

fit_sem3 <- sem(model_sem3, data = sem_data)

cat("\n=== SEM Modèle 3: Sans lien direct Structure->WD ===\n")
summary(fit_sem3, standardized = TRUE, fit.measures = TRUE)

# 11E. Comparaison des modèles SEM ----

sem_comparison <- data.frame(
  Modele = c("CA_cause_commune", "Temperament_cause_commune", "Sans_lien_direct"),
  ChiSq = c(fitmeasures(fit_sem1, "chisq"),
            fitmeasures(fit_sem2, "chisq"),
            fitmeasures(fit_sem3, "chisq")),
  df = c(fitmeasures(fit_sem1, "df"),
         fitmeasures(fit_sem2, "df"),
         fitmeasures(fit_sem3, "df")),
  CFI = c(fitmeasures(fit_sem1, "cfi"),
          fitmeasures(fit_sem2, "cfi"),
          fitmeasures(fit_sem3, "cfi")),
  RMSEA = c(fitmeasures(fit_sem1, "rmsea"),
            fitmeasures(fit_sem2, "rmsea"),
            fitmeasures(fit_sem3, "rmsea")),
  AIC = c(AIC(fit_sem1), AIC(fit_sem2), AIC(fit_sem3))
)

cat("\n=== Comparaison des modèles SEM ===\n")
print(sem_comparison)

write_csv2(sem_comparison, file.path(path_table, "sem_comparison.csv"))

# 11F. Extraire les coefficients standardisés du meilleur modèle ----

best_sem <- fit_sem1  # À ajuster selon résultats

sem_coefs <- parameterEstimates(best_sem, standardized = TRUE) %>%
  filter(op == "~") %>%
  select(lhs, rhs, est, std.all, pvalue)

cat("\n=== Coefficients standardisés (meilleur SEM) ===\n")
print(sem_coefs)

write_csv2(sem_coefs, file.path(path_table, "sem_coefficients.csv"))

#### 1️⃣2️⃣ VALIDATION CROISÉE LEAVE-ONE-SITE-OUT ----

# 12A. Identifier les sites ----

# Créer une variable "site" basée sur les préfixes des noms de plots
model_data <- model_data %>%
  mutate(
    site = case_when(
      grepl("^MAB", plot_name) ~ "Mabounie",
      grepl("^GIL", plot_name) ~ "Gilbertiodendron",
      grepl("^JEU|^Betamba", plot_name) ~ "Yangambi_secondary",
      grepl("^M2P", plot_name) ~ "Malebo",
      grepl("^PARAP", plot_name) ~ "Parapanga",
      grepl("^WWF|^Modele|^MIX|^Model", plot_name) ~ "WWF_sites",
      TRUE ~ "Other"
    )
  )

cat("\n=== Sites identifiés ===\n")
print(table(model_data$site))

# 12B. Leave-one-site-out CV ----

sites_unique <- unique(model_data$site)
n_sites <- length(sites_unique)

cv_results <- data.frame()

for(site_i in sites_unique) {
  # Split train/test
  train_data <- model_data %>% filter(site != site_i)
  test_data <- model_data %>% filter(site == site_i)
  
  if(nrow(test_data) == 0 || nrow(train_data) < 10) next
  
  # Modèle simple (meilleure variable)
  model_cv <- lm(formula_simple, data = train_data)
  
  # Prédiction
  pred <- predict(model_cv, newdata = test_data)
  
  # Métriques
  obs <- test_data$WD_BA
  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  mae <- mean(abs(obs - pred), na.rm = TRUE)
  r2 <- cor(obs, pred, use = "complete.obs")^2
  
  cv_results <- rbind(cv_results, data.frame(
    site_out = site_i,
    n_test = nrow(test_data),
    RMSE = rmse,
    MAE = mae,
    R2 = r2
  ))
}

cat("\n=== Résultats Leave-One-Site-Out CV ===\n")
print(cv_results)

# Résumé global
cat("\n=== Performance moyenne ===\n")
cat("RMSE moyen:", round(mean(cv_results$RMSE, na.rm = TRUE), 4), "\n")
cat("MAE moyen:", round(mean(cv_results$MAE, na.rm = TRUE), 4), "\n")
cat("R² moyen:", round(mean(cv_results$R2, na.rm = TRUE), 3), "\n")

write_csv2(cv_results, file.path(path_table, "cv_leave_site_out.csv"))

# 12C. Visualisation des prédictions CV ----

# Refaire les prédictions pour visualisation
all_cv_pred <- data.frame()

for(site_i in sites_unique) {
  train_data <- model_data %>% filter(site != site_i)
  test_data <- model_data %>% filter(site == site_i)
  
  if(nrow(test_data) == 0 || nrow(train_data) < 10) next
  
  model_cv <- lm(formula_simple, data = train_data)
  pred <- predict(model_cv, newdata = test_data)
  
  all_cv_pred <- rbind(all_cv_pred, data.frame(
    plot_name = test_data$plot_name,
    site = site_i,
    observed = test_data$WD_BA,
    predicted = pred
  ))
}

p_cv <- ggplot(all_cv_pred, aes(x = observed, y = predicted, color = site)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Validation croisée Leave-One-Site-Out",
    subtitle = paste("RMSE moyen =", round(mean(cv_results$RMSE, na.rm = TRUE), 3)),
    x = "WD observé (g/cm³)",
    y = "WD prédit (g/cm³)",
    color = "Site exclu"
  ) +
  theme_minimal() +
  coord_equal()

ggsave(file.path(path_figures, "cv_leave_site_out.jpg"), p_cv,
       width = 10, height = 8, dpi = 300)

#### 1️⃣3️⃣ SYNTHÈSE FINALE ----

cat("\n")
cat("================================================================\n")
cat("                  SYNTHÈSE COMPLÈTE DES ANALYSES                \n")
cat("================================================================\n")
cat("\n")

cat("1. NOUVELLES MÉTRIQUES:\n")
cat("   - Paramètres sigmoïdes (k, x0) calculés\n")
cat("   - Largeur de transition, asymétrie ajoutées\n")
cat("\n")

cat("2. SÉLECTION ELASTIC NET:\n")
cat("   Variables retenues:", paste(head(selected_vars_elastic, 5), collapse = ", "), "\n")
cat("\n")

cat("3. PERMANOVA:\n")
cat("   R² (quartiles WD):", round(permanova_quartile$R2[1], 3), "\n")
cat("   p-value:", permanova_quartile$`Pr(>F)`[1], "\n")
cat("\n")

cat("4. MODÈLE SPATIAL:\n")
cat("   AIC OLS:", round(aic_comparison$AIC[1], 1), "\n")
cat("   AIC Spatial:", round(aic_comparison$AIC[2], 1), "\n")
cat("\n")

cat("5. SEM:\n")
cat("   Meilleur modèle:", sem_comparison$Modele[which.min(sem_comparison$AIC)], "\n")
cat("\n")

cat("6. VALIDATION CROISÉE:\n")
cat("   RMSE moyen:", round(mean(cv_results$RMSE, na.rm = TRUE), 4), "\n")
cat("   R² moyen:", round(mean(cv_results$R2, na.rm = TRUE), 3), "\n")
cat("\n")
cat("================================================================\n")

# Sauvegarder l'environnement
save.image(file.path(path_data, "modelisation_complete.RData"))

cat("\nAnalyse complète terminée. Tous les résultats sauvegardés.\n")

