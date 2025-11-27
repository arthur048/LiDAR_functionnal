# =================================================================
# ANALYSE CONSOLIDÉE : WD ~ STRUCTURE FORESTIÈRE
# Script propre avec export structuré des résultats
# =================================================================

#### 1️⃣ INITIALISATION ----

rm(list = ls())
gc()

pkgs <- c("tidyverse", "here", "FactoMineR", "factoextra", "vegan",
          "lavaan", "spaMM", "sf", "e1071", "ade4")
to_install <- !pkgs %in% installed.packages()
if(any(to_install)) install.packages(pkgs[to_install])
lapply(pkgs, library, character.only = TRUE)

#### 2️⃣ CHEMINS ----

# Chemins
project_dir <- here()
path_main <- file.path(project_dir, "5_Main_analysis")
path_input <- file.path(path_main, "input")
path_output <- file.path(path_main, "output")
path_figures <- file.path(path_output, "figures", "modelisation")
path_tables <- file.path(path_output, "tables", "modelisation")
path_data <- file.path(path_output, "data", "modelisation")

dir.create(path_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(path_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(path_data, recursive = TRUE, showWarnings = FALSE)

# Fichier pour exporter les résultats clés
results_file <- file.path(path_tables, "results_summary.txt")
cat("=== RÉSULTATS D'ANALYSE - WD ~ STRUCTURE ===\n", file = results_file)
cat("Date:", as.character(Sys.time()), "\n\n", file = results_file, append = TRUE)

# Fonction helper pour exporter
export_result <- function(title, content, file = results_file) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n", file = file, append = TRUE)
  cat(title, "\n", file = file, append = TRUE)
  cat(paste(rep("=", 60), collapse = ""), "\n\n", file = file, append = TRUE)
  capture.output(print(content), file = file, append = TRUE)
  cat("\n", file = file, append = TRUE)
}

#### 3️⃣ CHARGEMENT DES DONNÉES ----

gaps_metrics <- read_csv2(file.path(path_input, "gaps_metrics.csv"))
plot_metrics <- read_csv2(file.path(path_input, "plot_metrics.csv"))
inventory_temperament <- read_csv2(file.path(path_input, "inventory_with_temperament.csv"))

cat("Données chargées:\n")
cat("  - gaps_metrics: n =", nrow(gaps_metrics), "\n")
cat("  - plot_metrics: n =", nrow(plot_metrics), "\n")
cat("  - inventory_temperament: n =", nrow(inventory_temperament), "\n")

#### 4️⃣ PRÉPARATION DES DONNÉES ----

# 4A. Filtrer les parcelles avec données de tempérament ----
plots_with_temperament <- inventory_temperament %>%
  pull(plot_ref) %>%
  unique()

plot_data_complete <- plot_metrics %>%
  filter(plot_name %in% plots_with_temperament)

# 4B. Définir les colonnes ----
metadata_cols <- c("plot_name", "WD_BA", "prop_g_helio", "prop_g_npld", "prop_g_shade",
                   "prop_ind_helio", "prop_ind_npld", "prop_ind_shade")

# Métriques de trouées à exclure
gap_metrics_to_remove <- c("basal_thickness", "transition_thickness", "upper_thickness",
                           "auc_0_25_norm", "auc_25_75_norm", "auc_75_100_norm",
                           "basal_auc_norm", "upper_auc_norm", "transition_auc_norm",
                           "basal_auc", "upper_auc", "transition_auc",
                           "height_min_d2", "min_d2", "height_max_d2", "max_d2")

gap_metrics_cols <- setdiff(names(plot_data_complete), c(metadata_cols, gap_metrics_to_remove))

# Metadata avec densité de bois et succession
metadata <- plot_data_complete %>%
  select(all_of(metadata_cols))

plots_name <- metadata %>%
  pull(plot_name)

# 4C. Nettoyer les données taxonomiques ----
inventory_clean <- inventory_temperament %>%
  filter(!is.na(species_full) &
           species_full != "" &
           species_full != "NA NA" &
           species_full != "sp.2 NA" &
           species_full != "indet indet " &
           species_full != "bois mort" &
           species_full != "Unknown NA"
  )

cat("\nInventaire nettoyé: n =", nrow(inventory_clean), "\n")

# 4D. Construire tables de contingence ----

# Table espèces × parcelles (abondance)
species_plot_table <- inventory_clean %>%
  count(plot_ref, species_full) %>%
  pivot_wider(names_from = species_full,
              values_from = n,
              values_fill = 0) %>%
  column_to_rownames("plot_ref")

# Table espèces × parcelles (surface terrière)
species_dbh_plot_table <- inventory_clean %>%
  group_by(plot_ref, species_full) %>%
  summarise(basal_area = sum(G), .groups = "drop") %>%
  pivot_wider(names_from = species_full,
              values_from = basal_area,
              values_fill = 0) %>%
  column_to_rownames("plot_ref")

# Table genre × parcelles
genus_plot_table <- inventory_clean %>%
  count(plot_ref, genus) %>%
  pivot_wider(names_from = genus,
              values_from = n,
              values_fill = 0) %>%
  column_to_rownames("plot_ref")

# Table espèce × densité du bois
species_wd_table <- inventory_clean %>%
  group_by(species_full) %>%
  reframe(wd = min(meanWD, na.rm = TRUE)) %>%
  distinct(species_full, .keep_all = TRUE) %>%
  column_to_rownames("species_full")

# 4E. Table structure × parcelles ----
structure_plots_table <- plot_data_complete %>%
  mutate(h100 = ifelse(is.na(h100), 1, h100)) %>%
  select(plot_name, all_of(gap_metrics_cols)) %>%
  column_to_rownames("plot_name")

cat("\nTables de contingence créées:\n")
cat("  - species_plot_table:", nrow(species_plot_table), "×", ncol(species_plot_table), "\n")
cat("  - species_dbh_plot_table:", nrow(species_dbh_plot_table), "×", ncol(species_dbh_plot_table), "\n")
cat("  - genus_plot_table:", nrow(genus_plot_table), "×", ncol(genus_plot_table), "\n")
cat("  - structure_plots_table:", nrow(structure_plots_table), "×", ncol(structure_plots_table), "\n")

#### 5️⃣ ANALYSES MULTIVARIÉES (NSCA, CA, PCA) ----

# 5A. NSCA - Analyse de Correspondance Non-Symétrique ----

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

cat("\n=== NSCA Summary ===\n")
cat("Variance expliquée axes 1-2:\n")
cat("  NSCA (abondance):", round(sum(nsca_result$eig[1:2])/sum(nsca_result$eig)*100, 1), "%\n")
cat("  NSCA DBH:", round(sum(nsca_dbh_result$eig[1:2])/sum(nsca_dbh_result$eig)*100, 1), "%\n")
cat("  NSCA Genus:", round(sum(nsca_genus_result$eig[1:2])/sum(nsca_genus_result$eig)*100, 1), "%\n")

# 5B. CA - Analyse de Correspondance ----

ca_result <- dudi.coa(species_plot_table, scannf = FALSE, nf = min(nrow(species_plot_table) - 1, ncol(species_plot_table) - 1))
ca_genus_result <- dudi.coa(genus_plot_table, scannf = FALSE, nf = min(nrow(genus_plot_table) - 1, ncol(genus_plot_table) - 1))

# Scores des parcelles sur les axes CA
ca_site_scores <- ca_result$li
ca_site_scores_genus <- ca_genus_result$li

colnames(ca_site_scores) <- paste0("CA", 1:ncol(ca_site_scores))
colnames(ca_site_scores_genus) <- paste0("CA_GENUS", 1:ncol(ca_site_scores_genus))

cat("\n=== CA Summary ===\n")
cat("Variance expliquée axes 1-2:\n")
cat("  CA (espèces):", round(sum(ca_result$eig[1:2])/sum(ca_result$eig)*100, 1), "%\n")
cat("  CA (genre):", round(sum(ca_genus_result$eig[1:2])/sum(ca_genus_result$eig)*100, 1), "%\n")

# 5C. PCA - Analyse en Composantes Principales ----

pca_result <- dudi.pca(structure_plots_table, scannf = FALSE, nf = min(nrow(structure_plots_table) - 1, ncol(structure_plots_table) - 1))

# Scores des parcelles sur les axes PCA
pca_site_scores <- pca_result$li

colnames(pca_site_scores) <- paste0("PCA", 1:ncol(pca_site_scores))

cat("\n=== PCA Summary ===\n")
cat("Variance expliquée axes 1-2:\n")
cat("  PCA1:", round(pca_result$eig[1]/sum(pca_result$eig)*100, 1), "%\n")
cat("  PCA2:", round(pca_result$eig[2]/sum(pca_result$eig)*100, 1), "%\n")

# 5D. Visualisations ----
p_pca_eig <- fviz_eig(pca_result, addlabels = TRUE)
ggsave(file.path(path_figures, "pca_screeplot.jpg"), p_pca_eig, width = 8, height = 6, dpi = 300)

p_pca_var <- fviz_pca_var(pca_result,
                          col.var = "contrib",
                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                          repel = TRUE)
ggsave(file.path(path_figures, "pca_variables.jpg"), p_pca_var, width = 10, height = 8, dpi = 300)

#### 6️⃣ CRÉATION DE model_data ----

# Joindre structure + WD + composition + scores multivariés
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
  ) %>%
  left_join(
    data.frame(
      plot_name = rownames(nsca_site_scores_dbh),
      NSCA_DBH1 = nsca_site_scores_dbh[,1],
      NSCA_DBH2 = nsca_site_scores_dbh[,2]
    ), by = "plot_name"
  )

cat("\nmodel_data créé: n =", nrow(model_data), "variables =", ncol(model_data), "\n")

#### 8️⃣ CORRÉLATIONS AVEC WD ----

# 8A. Corrélations de Pearson et Spearman ----
structure_cols <- c("prop_at_5m", "prop_at_10m", "prop_at_15m", "prop_at_20m",
                    "prop_at_25m", "prop_at_30m", "prop_at_35m", "prop_at_40m",
                    "h10", "h25", "h50", "h75", "h90", "h100",
                    "slope_basal", "slope_transition", "slope_upper",
                    "auc", "auc_norm", "auc_0_25", "auc_25_75", "auc_75_100",
                    "max_d1", "height_max_d1")

# Filtrer les colonnes existantes
structure_cols <- structure_cols[structure_cols %in% names(model_data)]

cor_results <- data.frame()

for(var in structure_cols) {
  if(var %in% names(model_data)) {
    x <- model_data[[var]]
    y <- model_data$WD_BA
    valid <- complete.cases(x, y)

    if(sum(valid) >= 10) {
      pearson <- cor.test(x[valid], y[valid], method = "pearson")
      spearman <- cor.test(x[valid], y[valid], method = "spearman")

      cor_results <- rbind(cor_results, data.frame(
        variable = var,
        r_pearson = pearson$estimate,
        p_pearson = pearson$p.value,
        rho_spearman = spearman$estimate,
        p_spearman = spearman$p.value,
        n = sum(valid)
      ))
    }
  }
}

cor_results <- cor_results %>%
  mutate(abs_r = abs(r_pearson),
         abs_rho = abs(rho_spearman)) %>%
  arrange(desc(abs_r))

cat("\n=== TOP 10 corrélations avec WD_BA ===\n")
print(head(cor_results %>% select(variable, r_pearson, rho_spearman, p_pearson), 10))

export_result("CORRÉLATIONS STRUCTURE - WD",
              cor_results %>% select(variable, r_pearson, rho_spearman, p_pearson, n))

# 8B. Visualisation ----
p_cor <- cor_results %>%
  head(15) %>%
  mutate(variable = fct_reorder(variable, abs_r)) %>%
  pivot_longer(cols = c(r_pearson, rho_spearman),
               names_to = "method", values_to = "correlation") %>%
  ggplot(aes(x = variable, y = correlation, fill = method)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("steelblue", "coral"),
                    labels = c("Pearson", "Spearman")) +
  coord_flip() +
  labs(title = "Corrélations Structure - WD",
       x = NULL, y = "Coefficient de corrélation") +
  theme_minimal()

ggsave(file.path(path_figures, "correlations_wd.jpg"), p_cor,
       width = 10, height = 8, dpi = 300)

write_csv2(cor_results, file.path(path_tables, "correlations_structure_wd.csv"))

#### 9️⃣ SÉLECTION DE VARIABLES (MANUELLE) ----

# 9A. Top variables par corrélation ----
top_vars <- cor_results %>%
  head(10) %>%
  pull(variable)

cat("\n=== TOP 10 variables à considérer ===\n")
print(cor_results %>%
        head(10) %>%
        select(variable, r_pearson, p_pearson, abs_r))

# 9B. Matrice de corrélation entre top variables ----
cor_matrix <- model_data %>%
  select(all_of(top_vars)) %>%
  cor(use = "pairwise.complete.obs")

cat("\n=== Corrélations entre top 10 variables ===\n")
cat("(pour identifier la colinéarité)\n")
print(round(cor_matrix, 2))

# Visualisation de la matrice de corrélation
p_cor_matrix <- cor_matrix %>%
  as.data.frame() %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "correlation") %>%
  ggplot(aes(x = var1, y = var2, fill = correlation)) +
  geom_tile() +
  geom_text(aes(label = round(correlation, 2)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Matrice de corrélation - Top 10 variables",
       x = NULL, y = NULL)

ggsave(file.path(path_figures, "correlation_matrix_top10.jpg"), p_cor_matrix,
       width = 10, height = 9, dpi = 300)

# 9C. SÉLECTION MANUELLE DES VARIABLES ----
# À compléter par l'utilisateur après analyse des corrélations
# Décommenter et modifier la ligne suivante avec vos variables choisies:

selected_vars <- c("prop_at_20m",
                   "h50",
                   "max_d1")  # EXEMPLE - À MODIFIER

# Par défaut, si aucune sélection manuelle, on prend la meilleure variable
if(!exists("selected_vars") || length(selected_vars) == 0) {
  selected_vars <- c(cor_results$variable[1])  # Meilleure variable uniquement
  cat("\n⚠️  Aucune sélection manuelle - utilisation de la meilleure variable uniquement\n")
  cat("Pour sélectionner plusieurs variables, décommentez et modifiez la ligne:\n")
  cat('selected_vars <- c("var1", "var2", "var3")\n')
}

cat("\n=== Variables sélectionnées pour les modèles ===\n")
cat(paste(selected_vars, collapse = ", "), "\n")

# 9D. Test VIF sur les variables sélectionnées ----
if(length(selected_vars) >= 2) {
  formula_vif <- as.formula(paste("WD_BA ~", paste(selected_vars, collapse = " + ")))
  model_vif <- lm(formula_vif, data = model_data)
  vif_values <- car::vif(model_vif)

  cat("\n=== VIF - Variables sélectionnées ===\n")
  print(vif_values)

  if(any(vif_values > 10)) {
    cat("\n⚠️  ATTENTION: VIF > 10 détecté (multicolinéarité forte)\n")
    cat("Variables problématiques:", names(vif_values)[vif_values > 10], "\n")
    cat("Conseil: Retirez une des variables fortement corrélées\n")
  } else if(any(vif_values > 5)) {
    cat("\n⚠️  VIF > 5 détecté (multicolinéarité modérée)\n")
  } else {
    cat("\n✓ Pas de problème de multicolinéarité (VIF < 5)\n")
  }

  export_result("VIF - VARIABLES SÉLECTIONNÉES", vif_values)
} else {
  cat("\n(Une seule variable sélectionnée - pas de test VIF nécessaire)\n")
}

#### 🔟 MODÈLES PRÉDICTIFS ----

# 10A. Modèles simples (une variable à la fois) ----
# Modèles pour chaque variable sélectionnée
simple_models <- list()
simple_results <- data.frame()

cat("\n=== Modèles simples sur variables sélectionnées ===\n")

for(var in selected_vars) {
  formula_i <- as.formula(paste("WD_BA ~", var))
  model_i <- lm(formula_i, data = model_data)
  simple_models[[var]] <- model_i

  s <- summary(model_i)
  simple_results <- rbind(simple_results, data.frame(
    variable = var,
    R2 = s$r.squared,
    R2_adj = s$adj.r.squared,
    coef = coef(model_i)[2],
    pvalue = s$coefficients[2,4],
    AIC = AIC(model_i),
    RMSE = sqrt(mean(residuals(model_i)^2))
  ))
}

simple_results <- simple_results %>% arrange(desc(R2_adj))

cat("\n=== Résultats des modèles simples ===\n")
print(simple_results)
export_result("MODÈLES SIMPLES WD ~ variable", simple_results)

# 10B. Meilleur modèle simple parmi les variables sélectionnées ----
best_var <- simple_results$variable[1]
best_model <- simple_models[[best_var]]

cat("\n=== Meilleur modèle simple ===\n")
cat("Variable:", best_var, "\n")
print(summary(best_model))
export_result("MEILLEUR MODÈLE SIMPLE (détail)", summary(best_model))

# 10C. Modèle multiple (variables sélectionnées) ----
if(length(selected_vars) >= 2) {
  formula_multi <- as.formula(paste("WD_BA ~", paste(selected_vars, collapse = " + ")))
  model_multi <- lm(formula_multi, data = model_data)

  cat("\n=== Modèle multiple ===\n")
  print(summary(model_multi))
  export_result("MODÈLE MULTIPLE", summary(model_multi))
}

# 10D. Modèle PCA1 ----
model_pca1 <- lm(WD_BA ~ PCA1, data = model_data)

cat("\n=== Modèle WD ~ PCA1 ===\n")
print(summary(model_pca1))
export_result("MODÈLE WD ~ PCA1", summary(model_pca1))

# 10E. Tableau comparatif ----
comparison_models <- data.frame(
  Modele = c(paste0("Simple (", best_var, ")"),
             "Multiple",
             "PCA1"),
  R2_adj = c(summary(best_model)$adj.r.squared,
             if(exists("model_multi")) summary(model_multi)$adj.r.squared else NA,
             summary(model_pca1)$adj.r.squared),
  AIC = c(AIC(best_model),
          if(exists("model_multi")) AIC(model_multi) else NA,
          AIC(model_pca1)),
  RMSE = c(sqrt(mean(residuals(best_model)^2)),
           if(exists("model_multi")) sqrt(mean(residuals(model_multi)^2)) else NA,
           sqrt(mean(residuals(model_pca1)^2)))
) %>%
  filter(!is.na(R2_adj)) %>%
  arrange(AIC)

cat("\n=== Comparaison des modèles ===\n")
print(comparison_models)
export_result("COMPARAISON MODÈLES PRÉDICTIFS", comparison_models)

write_csv2(simple_results, file.path(path_tables, "simple_models_comparison.csv"))
write_csv2(comparison_models, file.path(path_tables, "all_models_comparison.csv"))

#### 1️⃣1️⃣ PARTITION DE VARIANCE ----

# Préparer données complètes
varpart_data <- model_data %>%
  select(WD_BA, PCA1, prop_g_helio, prop_g_shade,
         starts_with("CA"), starts_with("NSCA")) %>%
  drop_na()

cat("\nDonnées pour varpart: n =", nrow(varpart_data), "\n")

# 11A. Structure (PCA1) vs Tempérament ----
if(all(c("prop_g_helio", "prop_g_shade") %in% names(varpart_data))) {
  vp_temp <- varpart(varpart_data$WD_BA,
                     ~ PCA1,
                     ~ prop_g_helio + prop_g_shade,
                     data = varpart_data)

  cat("\n=== Partition: Structure vs Tempérament ===\n")
  print(vp_temp)
  export_result("VARPART: Structure (PCA1) vs Tempérament", vp_temp)

  # Figure
  jpeg(file.path(path_figures, "varpart_temperament.jpg"),
       width = 800, height = 600, quality = 100)
  plot(vp_temp, digits = 2,
       Xnames = c("Structure", "Tempérament"),
       bg = c("lightblue", "lightgreen"))
  dev.off()
}

# 11B. Structure (PCA1) vs CA ----
if(all(c("CA1", "CA2") %in% names(varpart_data))) {
  vp_ca <- varpart(varpart_data$WD_BA,
                   ~ PCA1,
                   ~ CA1 + CA2,
                   data = varpart_data)

  cat("\n=== Partition: Structure vs CA ===\n")
  print(vp_ca)
  export_result("VARPART: Structure (PCA1) vs CA", vp_ca)

  # Figure
  jpeg(file.path(path_figures, "varpart_ca.jpg"),
       width = 800, height = 600, quality = 100)
  plot(vp_ca, digits = 2,
       Xnames = c("Structure", "CA"),
       bg = c("lightblue", "lightyellow"))
  dev.off()
}

# 11C. Structure vs NSCA_DBH (si disponible) ----
nsca_cols <- grep("NSCA.*DBH", names(varpart_data), value = TRUE)
if(length(nsca_cols) >= 2) {
  formula_nsca <- as.formula(paste("~", paste(nsca_cols[1:2], collapse = " + ")))

  vp_nsca <- varpart(varpart_data$WD_BA,
                     ~ PCA1,
                     formula_nsca,
                     data = varpart_data)

  cat("\n=== Partition: Structure vs NSCA_DBH ===\n")
  print(vp_nsca)
  export_result("VARPART: Structure (PCA1) vs NSCA_DBH", vp_nsca)
}

#### 1️⃣2️⃣ SEM - MODÈLES CAUSAUX ----

# 12A. Préparer données ----
sem_data <- model_data %>%
  select(WD_BA, PCA1, prop_g_helio, prop_g_shade) %>%
  drop_na() %>%
  mutate(across(everything(), scale))  # Standardiser

cat("\nDonnées pour SEM: n =", nrow(sem_data), "\n")

# 12B. Modèle 1: Avec lien Structure -> WD ----
model_sem_full <- '
  # Composition prédit Structure
  PCA1 ~ a1*prop_g_helio + a2*prop_g_shade

  # Composition prédit WD
  WD_BA ~ b1*prop_g_helio + b2*prop_g_shade

  # Structure prédit WD (lien résiduel)
  WD_BA ~ c*PCA1
'

fit_full <- sem(model_sem_full, data = sem_data)

# 12C. Modèle 2: Sans lien Structure -> WD ----
model_sem_constrained <- '
  # Composition prédit Structure
  PCA1 ~ a1*prop_g_helio + a2*prop_g_shade

  # Composition prédit WD (pas de lien via Structure)
  WD_BA ~ b1*prop_g_helio + b2*prop_g_shade

  # Contrainte: pas de lien direct Structure -> WD
  WD_BA ~ 0*PCA1
'

fit_constrained <- sem(model_sem_constrained, data = sem_data)

# 12D. Comparaison des modèles ----
comparison_sem <- anova(fit_constrained, fit_full)

cat("\n=== Comparaison SEM: avec vs sans lien Structure->WD ===\n")
print(comparison_sem)
export_result("SEM: Comparaison modèles (test du lien Structure->WD)", comparison_sem)

# 12E. Coefficients du modèle complet ----
sem_params <- parameterEstimates(fit_full, standardized = TRUE) %>%
  filter(op == "~") %>%
  select(lhs, op, rhs, est, std.all, pvalue)

cat("\n=== Coefficients SEM (modèle complet) ===\n")
print(sem_params)
export_result("SEM: Coefficients standardisés", sem_params)

write_csv2(sem_params, file.path(path_tables, "sem_coefficients.csv"))

#### 1️⃣3️⃣ DIFFÉRENCES STRUCTURELLES PAR WD ----

# 13A. PERMANOVA ----
permanova_data <- model_data %>%
  select(plot_name, WD_BA, all_of(structure_cols)) %>%
  drop_na()

# Créer quartiles
permanova_data <- permanova_data %>%
  mutate(WD_quartile = cut(WD_BA,
                           breaks = quantile(WD_BA, c(0, 0.25, 0.5, 0.75, 1)),
                           labels = c("Q1", "Q2", "Q3", "Q4"),
                           include.lowest = TRUE))

# Matrice de distance
dist_matrix <- permanova_data %>%
  select(all_of(structure_cols)) %>%
  scale() %>%
  dist()

set.seed(42)
permanova_result <- adonis2(dist_matrix ~ WD_quartile,
                            data = permanova_data,
                            permutations = 999)

cat("\n=== PERMANOVA ===\n")
print(permanova_result)
export_result("PERMANOVA: Structure ~ WD quartiles", permanova_result)

# 13B. Test de Mantel (corrélation continue) ----
dist_wd <- dist(permanova_data$WD_BA)
mantel_result <- mantel(dist_matrix, dist_wd, method = "pearson", permutations = 999)

cat("\n=== Test de Mantel ===\n")
print(mantel_result)
export_result("MANTEL: Distance structurelle vs Distance WD", mantel_result)

# 13C. Statistiques descriptives par quartile ----
key_metrics <- head(cor_results$variable, 8)

desc_by_quartile <- permanova_data %>%
  group_by(WD_quartile) %>%
  summarise(
    n = n(),
    WD_mean = mean(WD_BA),
    across(all_of(key_metrics[key_metrics %in% names(.)]),
           list(mean = ~mean(., na.rm = TRUE),
                sd = ~sd(., na.rm = TRUE))),
    .groups = "drop"
  )

cat("\n=== Statistiques par quartile WD ===\n")
print(desc_by_quartile)
export_result("DESCRIPTION: Métriques par quartile WD", desc_by_quartile)

write_csv2(desc_by_quartile, file.path(path_tables, "stats_by_wd_quartile.csv"))

# 13D. Tests de différence par métrique ----
metric_tests <- data.frame()

for(metric in key_metrics) {
  if(metric %in% names(permanova_data)) {
    kruskal <- kruskal.test(permanova_data[[metric]] ~ permanova_data$WD_quartile)

    metric_tests <- rbind(metric_tests, data.frame(
      metric = metric,
      chi_squared = kruskal$statistic,
      df = kruskal$parameter,
      p_value = kruskal$p.value
    ))
  }
}

metric_tests <- metric_tests %>% arrange(p_value)

cat("\n=== Tests Kruskal-Wallis par métrique ===\n")
print(metric_tests)
export_result("KRUSKAL-WALLIS: Différences par quartile WD", metric_tests)

write_csv2(metric_tests, file.path(path_tables, "kruskal_wallis_by_metric.csv"))

#### 1️⃣4️⃣ AUTOCORRÉLATION SPATIALE ----

# 14A. Charger les coordonnées ----
path_plots_gpkg <- file.path(dirname(project_dir), "0_Inventories_plot_preparation", "final", "plots_unique")

if(dir.exists(path_plots_gpkg)) {
  gpkg_files <- list.files(path_plots_gpkg, pattern = "\\.gpkg$", full.names = TRUE)

  coords_list <- lapply(gpkg_files, function(f) {
    tryCatch({
      sf_obj <- st_read(f, quiet = TRUE)
      centroid <- st_centroid(st_union(sf_obj))
      coords <- st_coordinates(centroid)
      data.frame(
        plot_name = tools::file_path_sans_ext(basename(f)),
        X_utm = coords[1],
        Y_utm = coords[2]
      )
    }, error = function(e) NULL)
  })

  coords_df <- bind_rows(coords_list)
  model_data <- model_data %>% left_join(coords_df, by = "plot_name")

  cat("\nCoordonnées chargées pour", sum(!is.na(model_data$X_utm)), "plots\n")
}

# 14B. Modèle spatial (si coordonnées disponibles) ----
if("X_utm" %in% names(model_data) && sum(!is.na(model_data$X_utm)) > 50) {

  spatial_data <- model_data %>%
    filter(!is.na(X_utm) & !is.na(Y_utm) & !is.na(WD_BA) & !is.na(!!sym(best_var)))

  # Modèle OLS
  formula_best <- as.formula(paste("WD_BA ~", best_var))
  model_ols <- lm(formula_best, data = spatial_data)

  # Modèle spatial Matérn
  model_spatial <- fitme(
    as.formula(paste("WD_BA ~", best_var, "+ Matern(1 | X_utm + Y_utm)")),
    data = spatial_data,
    method = "REML"
  )

  cat("\n=== Modèle OLS vs Spatial ===\n")
  cat("AIC OLS:", AIC(model_ols), "\n")
  cat("AIC Spatial (marginal):", AIC(model_spatial, which = "marginal"), "\n")

  spatial_comparison <- data.frame(
    Modele = c("OLS", "Spatial_Matern"),
    AIC = c(AIC(model_ols), AIC(model_spatial, which = "marginal"))
  )
  export_result("COMPARAISON: OLS vs Spatial", spatial_comparison)

  write_csv2(spatial_comparison, file.path(path_tables, "spatial_models_comparison.csv"))
}

#### 1️⃣5️⃣ LEAVE-ONE-OUT CROSS-VALIDATION ----

# LOOCV simple
loocv_data <- model_data %>%
  select(plot_name, WD_BA, all_of(best_var)) %>%
  drop_na()

predictions <- numeric(nrow(loocv_data))

for(i in 1:nrow(loocv_data)) {
  train <- loocv_data[-i, ]
  test <- loocv_data[i, ]

  model_cv <- lm(as.formula(paste("WD_BA ~", best_var)), data = train)
  predictions[i] <- predict(model_cv, newdata = test)
}

loocv_results <- data.frame(
  observed = loocv_data$WD_BA,
  predicted = predictions,
  residual = loocv_data$WD_BA - predictions
)

loocv_metrics <- data.frame(
  RMSE = sqrt(mean(loocv_results$residual^2)),
  MAE = mean(abs(loocv_results$residual)),
  R2 = cor(loocv_results$observed, loocv_results$predicted)^2
)

cat("\n=== LOOCV Résultats ===\n")
print(loocv_metrics)
export_result("LOOCV: Métriques de validation", loocv_metrics)

# Visualisation
p_loocv <- ggplot(loocv_results, aes(x = observed, y = predicted)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Leave-One-Out Cross-Validation",
       subtitle = paste("RMSE =", round(loocv_metrics$RMSE, 4),
                        "| R² =", round(loocv_metrics$R2, 3)),
       x = "WD observé", y = "WD prédit") +
  theme_minimal() +
  coord_equal()

ggsave(file.path(path_figures, "loocv_validation.jpg"), p_loocv,
       width = 8, height = 8, dpi = 300)

write_csv2(loocv_results, file.path(path_tables, "loocv_predictions.csv"))
write_csv2(loocv_metrics, file.path(path_tables, "loocv_metrics.csv"))

#### 1️⃣6️⃣ SYNTHÈSE FINALE ----

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("                         SYNTHÈSE DES RÉSULTATS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("1. MEILLEUR PRÉDICTEUR SIMPLE:\n")
cat("   Variable:", best_var, "\n")
cat("   R² ajusté:", round(summary(best_model)$adj.r.squared, 3), "\n")
cat("   RMSE:", round(sqrt(mean(residuals(best_model)^2)), 4), "\n\n")

cat("2. PARTITION DE VARIANCE (Structure vs Tempérament):\n")
if(exists("vp_temp")) {
  cat("   Structure unique:", round(vp_temp$part$indfract$Adj.R.squared[1], 3), "\n")
  cat("   Tempérament unique:", round(vp_temp$part$indfract$Adj.R.squared[2], 3), "\n")
  cat("   Partagée:", round(vp_temp$part$indfract$Adj.R.squared[3], 3), "\n\n")
}

cat("3. SEM - Lien résiduel Structure->WD:\n")
if(exists("comparison_sem")) {
  pval_sem <- comparison_sem$`Pr(>Chisq)`[2]
  cat("   Test χ² p-value:", ifelse(is.na(pval_sem), "NA", format.pval(pval_sem)), "\n")
  cat("   Interprétation:", ifelse(pval_sem < 0.05, "Lien significatif", "Lien non significatif"), "\n\n")
}

cat("4. PERMANOVA:\n")
cat("   R²:", round(permanova_result$R2[1], 3), "\n")
cat("   p-value:", permanova_result$`Pr(>F)`[1], "\n\n")

cat("5. VALIDATION (LOOCV):\n")
cat("   RMSE:", round(loocv_metrics$RMSE, 4), "\n")
cat("   R²:", round(loocv_metrics$R2, 3), "\n\n")

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Résultats exportés dans:", results_file, "\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Export synthèse
export_result("SYNTHÈSE COMPLÈTE", capture.output({
  cat("1. Meilleur prédicteur:", best_var, "\n")
  cat("2. PERMANOVA R²:", round(permanova_result$R2[1], 3), "\n")
  cat("3. Mantel r:", round(mantel_result$statistic, 3), "\n")
  if(exists("vp_temp")) {
    cat("4. Varpart Structure unique:", round(vp_temp$part$indfract$Adj.R.squared[1], 3), "\n")
  }
  cat("5. LOOCV RMSE:", round(loocv_metrics$RMSE, 4), "\n")
}))

# Sauvegarder
save.image(file.path(path_data, "analysis_complete.RData"))
write_csv2(model_data, file.path(path_data, "model_data_complete.csv"))

cat("\nAnalyse terminée.\n")
cat("  - Données:", file.path(path_data, "analysis_complete.RData"), "\n")
cat("  - Résultats:", results_file, "\n")
cat("  - Figures:", path_figures, "\n")
cat("  - Tableaux:", path_tables, "\n\n")
