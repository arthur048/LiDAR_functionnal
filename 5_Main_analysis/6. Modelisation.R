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
if(any(to_install)) renv::install(pkgs[to_install])
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

p_cor

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

p_cor_matrix

ggsave(file.path(path_figures, "correlation_matrix_top10.jpg"), p_cor_matrix,
       width = 10, height = 9, dpi = 300)

# 9C. SÉLECTION MANUELLE DES VARIABLES ----
# À compléter par l'utilisateur après analyse des corrélations
# Décommenter et modifier la ligne suivante avec vos variables choisies:

selected_vars <- c("prop_at_20m",
                   "auc_norm",
                   "max_d1",
                   "h90")  # EXEMPLE - À MODIFIER

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
  select(WD_BA, PCA1, prop_g_helio, prop_g_npld, prop_g_shade,
         starts_with("CA"), starts_with("NSCA")) %>%
  drop_na()

cat("\nDonnées pour varpart: n =", nrow(varpart_data), "\n")

# 11A. Structure (PCA1) vs Tempérament (3 catégories) ----
# INTERPRÉTATION: Compare la variance de WD expliquée par la structure verticale
# vs celle expliquée par les proportions de tempéraments des espèces
if(all(c("prop_g_helio", "prop_g_npld", "prop_g_shade") %in% names(varpart_data))) {
  vp_temp <- varpart(varpart_data$WD_BA,
                     ~ PCA1,
                     ~ prop_g_helio + prop_g_npld + prop_g_shade,
                     data = varpart_data)

  cat("\n=== Partition: Structure vs Tempérament (3 catégories) ===\n")
  cat("INTERPRÉTATION:\n")
  cat("  [a] = Variance WD expliquée UNIQUEMENT par Structure (effet pur)\n")
  cat("  [b] = Variance WD expliquée UNIQUEMENT par Tempérament (effet pur)\n")
  cat("  [c] = Variance PARTAGÉE (covariation Structure-Tempérament)\n")
  cat("  [d] = Variance RÉSIDUELLE non expliquée\n\n")
  print(vp_temp)
  export_result("VARPART: Structure vs Tempérament (3 catégories)", vp_temp)

  # Figure
  jpeg(file.path(path_figures, "varpart_temperament.jpg"),
       width = 800, height = 600, quality = 100)
  plot(vp_temp, digits = 2,
       Xnames = c("Structure", "Tempérament"),
       bg = c("lightblue", "lightgreen"))
  dev.off()
}

# 11B. Structure (PCA1) vs Composition floristique (CA) ----
# INTERPRÉTATION: CA capture les différences de composition en espèces
if(all(c("CA1", "CA2") %in% names(varpart_data))) {
  vp_ca <- varpart(varpart_data$WD_BA,
                   ~ PCA1,
                   ~ CA1 + CA2,
                   data = varpart_data)

  cat("\n=== Partition: Structure vs Composition floristique (CA) ===\n")
  cat("INTERPRÉTATION: CA = axes de variation de la composition en espèces\n\n")
  print(vp_ca)
  export_result("VARPART: Structure vs CA", vp_ca)

  # Figure
  jpeg(file.path(path_figures, "varpart_ca.jpg"),
       width = 800, height = 600, quality = 100)
  plot(vp_ca, digits = 2,
       Xnames = c("Structure", "CA"),
       bg = c("lightblue", "lightyellow"))
  dev.off()
}

# 11C. Structure vs NSCA_DBH (composition pondérée par surface terrière) ----
nsca_cols <- grep("NSCA.*DBH", names(varpart_data), value = TRUE)
if(length(nsca_cols) >= 2) {
  formula_nsca <- as.formula(paste("~", paste(nsca_cols[1:2], collapse = " + ")))

  vp_nsca <- varpart(varpart_data$WD_BA,
                     ~ PCA1,
                     formula_nsca,
                     data = varpart_data)

  cat("\n=== Partition: Structure vs NSCA_DBH ===\n")
  cat("INTERPRÉTATION: NSCA_DBH = composition pondérée par dominance (surface terrière)\n\n")
  print(vp_nsca)
  export_result("VARPART: Structure vs NSCA_DBH", vp_nsca)

  # Figure
  jpeg(file.path(path_figures, "varpart_nsca_dbh.jpg"),
       width = 800, height = 600, quality = 100)
  plot(vp_nsca, digits = 2,
       Xnames = c("Structure", "NSCA_DBH"),
       bg = c("lightblue", "orange"))
  dev.off()
}

# 11D. Structure vs TOUTES les dimensions compositionnelles regroupées ----
# INTERPRÉTATION: Test le plus conservateur - effet structure au-delà de TOUTE la composition
comp_cols <- c()
if("CA1" %in% names(varpart_data)) comp_cols <- c(comp_cols, "CA1", "CA2")
if("NSCA_DBH1" %in% names(varpart_data)) comp_cols <- c(comp_cols, "NSCA_DBH1", "NSCA_DBH2")
if("prop_g_helio" %in% names(varpart_data)) comp_cols <- c(comp_cols, "prop_g_helio", "prop_g_npld", "prop_g_shade")

if(length(comp_cols) >= 2) {
  formula_all_comp <- as.formula(paste("~", paste(comp_cols, collapse = " + ")))

  vp_all_comp <- varpart(varpart_data$WD_BA,
                         ~ PCA1,
                         formula_all_comp,
                         data = varpart_data)

  cat("\n=== Partition: Structure vs TOUTE la Composition (CA+NSCA+Temp) ===\n")
  cat("INTERPRÉTATION: Test conservateur - effet structure RÉSIDUEL après contrôle de:\n")
  cat("  - Composition en espèces (CA)\n")
  cat("  - Composition pondérée (NSCA_DBH)\n")
  cat("  - Tempéraments fonctionnels\n")
  cat("\nVariables composition:", paste(comp_cols, collapse = ", "), "\n\n")
  print(vp_all_comp)
  export_result("VARPART: Structure vs Composition complète (CA+NSCA+Temp)", vp_all_comp)

  # Figure
  jpeg(file.path(path_figures, "varpart_all_composition.jpg"),
       width = 800, height = 600, quality = 100)
  plot(vp_all_comp, digits = 2,
       Xnames = c("Structure", "Composition\n(CA+NSCA+Temp)"),
       bg = c("lightblue", "coral"))
  dev.off()
}

# 11E. Test de significativité des fractions varpart ----
# OBJECTIF: Tester si la fraction unique de la structure (2.7%) est significative
# via permutations sur le modèle RDA partiel

cat("\n=== Test RDA - Significativité des fractions ===\n")

# RDA partielle : Structure | Composition
rda_struct_partial <- rda(varpart_data$WD_BA ~ PCA1 + 
                            Condition(prop_g_helio + prop_g_shade + CA1 + CA2 + NSCA_DBH1 + NSCA_DBH2),
                          data = varpart_data)

set.seed(42)
rda_test <- anova(rda_struct_partial, permutations = 999)

cat("\nTest permutationnel - Structure | Composition :\n")
print(rda_test)
export_result("RDA: Test fraction unique Structure", rda_test)

# Interpréter
if(rda_test$`Pr(>F)`[1] < 0.05) {
  cat("\n✓ La fraction unique de la structure (", 
      round(vp_all_comp$part$indfract$Adj.R.squared[1] * 100, 1), 
      "%) est SIGNIFICATIVE\n")
} else {
  cat("\n✗ La fraction unique de la structure n'est PAS significative\n")
}
#### 1️⃣2️⃣ SEM - MODÈLES CAUSAUX ----
# INTERPRÉTATION: Les SEM testent les relations causales entre composition, structure et WD
# Question clé: Y a-t-il un effet DIRECT de la structure sur WD, ou tout passe par la composition?

# 12A. Préparer données avec TOUTES les variables compositionnelles ----
sem_data_full <- model_data %>%
  select(WD_BA, PCA1, prop_g_helio, prop_g_npld, prop_g_shade,
         starts_with("CA"), starts_with("NSCA_DBH")) %>%
  drop_na()

cat("\nDonnées pour SEM: n =", nrow(sem_data_full), "\n")
cat("Variables:", paste(names(sem_data_full), collapse = ", "), "\n\n")

# ==============================================================================
# SÉRIE 1: SEM avec Tempérament (3 catégories)
# ==============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SÉRIE 1: SEM avec Tempérament (héliophile, tolérant)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Données standardisées pour SEM tempérament
sem_data_temp <- sem_data_full %>%
  select(WD_BA, PCA1, prop_g_helio, prop_g_shade) %>%
  drop_na() %>%
  mutate(across(everything(), scale))

# 12A.1. Modèle AVEC lien direct Structure -> WD ----
model_sem_temp_full <- '
  # Tempérament prédit Structure
  PCA1 ~ a1*prop_g_helio + a2*prop_g_shade

  # Tempérament ET Structure prédisent WD (tout sur une ligne!)
  WD_BA ~ b1*prop_g_helio + b2*prop_g_shade + c*PCA1
'

fit_temp_full <- sem(model_sem_temp_full, data = sem_data_temp)

# 12A.2. Modèle SANS lien direct Structure -> WD (contraint à 0) ----
model_sem_temp_constrained <- '
  # Tempérament prédit Structure
  PCA1 ~ a1*prop_g_helio + a2*prop_g_shade

  # Tempérament prédit WD (pas de lien via Structure)
  WD_BA ~ b1*prop_g_helio + b2*prop_g_shade

  # Contrainte: PAS de lien direct Structure -> WD
  WD_BA ~ 0*PCA1
'

fit_temp_constrained <- sem(model_sem_temp_constrained, data = sem_data_temp)

# 12A.3. Test du lien Structure->WD ----
comparison_sem_temp <- anova(fit_temp_constrained, fit_temp_full)

cat("\n=== SEM Tempérament: Test du lien Structure->WD ===\n")
cat("H0: Le lien Structure->WD est nul (c = 0)\n")
cat("H1: Il existe un lien résiduel Structure->WD (c ≠ 0)\n\n")
print(comparison_sem_temp)
export_result("SEM Tempérament: Test lien Structure->WD", comparison_sem_temp)

# Coefficients
sem_temp_params <- parameterEstimates(fit_temp_full, standardized = TRUE) %>%
  filter(op == "~") %>%
  select(lhs, op, rhs, est, std.all, pvalue)

cat("\n=== Coefficients SEM Tempérament (standardisés) ===\n")
print(sem_temp_params)
export_result("SEM Tempérament: Coefficients", sem_temp_params)

# ==============================================================================
# SÉRIE 2: SEM avec CA (Composition floristique)
# ==============================================================================

if(all(c("CA1", "CA2") %in% names(sem_data_full))) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("SÉRIE 2: SEM avec CA (axes de composition en espèces)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  sem_data_ca <- sem_data_full %>%
    select(WD_BA, PCA1, CA1, CA2) %>%
    mutate(across(everything(), scale))

  # 12B.1. Modèle AVEC lien Structure -> WD ----
  model_sem_ca_full <- '
    PCA1 ~ a1*CA1 + a2*CA2
    WD_BA ~ b1*CA1 + b2*CA2 + c*PCA1
  '

  fit_ca_full <- sem(model_sem_ca_full, data = sem_data_ca)

  # 12B.2. Modèle SANS lien Structure -> WD ----
  model_sem_ca_constrained <- '
    PCA1 ~ a1*CA1 + a2*CA2
    WD_BA ~ b1*CA1 + b2*CA2 + 0*PCA1
  '

  fit_ca_constrained <- sem(model_sem_ca_constrained, data = sem_data_ca)

  # Test
  comparison_sem_ca <- anova(fit_ca_constrained, fit_ca_full)

  cat("\n=== SEM CA: Test du lien Structure->WD ===\n")
  print(comparison_sem_ca)
  export_result("SEM CA: Test lien Structure->WD", comparison_sem_ca)

  sem_ca_params <- parameterEstimates(fit_ca_full, standardized = TRUE) %>%
    filter(op == "~") %>%
    select(lhs, op, rhs, est, std.all, pvalue)

  cat("\n=== Coefficients SEM CA ===\n")
  print(sem_ca_params)
  export_result("SEM CA: Coefficients", sem_ca_params)
}

# ==============================================================================
# SÉRIE 3: SEM avec NSCA_DBH (Composition pondérée)
# ==============================================================================

nsca_dbh_cols <- grep("NSCA_DBH[12]", names(sem_data_full), value = TRUE)
if(length(nsca_dbh_cols) >= 2) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("SÉRIE 3: SEM avec NSCA_DBH (composition pondérée par dominance)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  sem_data_nsca <- sem_data_full %>%
    select(WD_BA, PCA1, all_of(nsca_dbh_cols[1:2])) %>%
    mutate(across(everything(), scale))

  names(sem_data_nsca)[3:4] <- c("NSCA1", "NSCA2")

  # 12C.1. Modèle AVEC lien ----
  model_sem_nsca_full <- '
    PCA1 ~ a1*NSCA1 + a2*NSCA2
    WD_BA ~ b1*NSCA1 + b2*NSCA2 + c*PCA1
  '

  fit_nsca_full <- sem(model_sem_nsca_full, data = sem_data_nsca)

  # 12C.2. Modèle SANS lien ----
  model_sem_nsca_constrained <- '
    PCA1 ~ a1*NSCA1 + a2*NSCA2
    WD_BA ~ b1*NSCA1 + b2*NSCA2 + 0*PCA1
  '

  fit_nsca_constrained <- sem(model_sem_nsca_constrained, data = sem_data_nsca)

  # Test
  comparison_sem_nsca <- anova(fit_nsca_constrained, fit_nsca_full)

  cat("\n=== SEM NSCA_DBH: Test du lien Structure->WD ===\n")
  print(comparison_sem_nsca)
  export_result("SEM NSCA_DBH: Test lien Structure->WD", comparison_sem_nsca)

  sem_nsca_params <- parameterEstimates(fit_nsca_full, standardized = TRUE) %>%
    filter(op == "~") %>%
    select(lhs, op, rhs, est, std.all, pvalue)

  cat("\n=== Coefficients SEM NSCA_DBH ===\n")
  print(sem_nsca_params)
  export_result("SEM NSCA_DBH: Coefficients", sem_nsca_params)
}

# ==============================================================================
# SÉRIE 4: SEM COMPLET avec TOUTES les dimensions compositionnelles
# ==============================================================================

# Vérifier quelles variables sont disponibles
comp_vars_sem <- c()
if("CA1" %in% names(sem_data_full)) comp_vars_sem <- c(comp_vars_sem, "CA1", "CA2")
if("NSCA_DBH1" %in% names(sem_data_full)) comp_vars_sem <- c(comp_vars_sem, "NSCA_DBH1", "NSCA_DBH2")
if("prop_g_helio" %in% names(sem_data_full)) comp_vars_sem <- c(comp_vars_sem, "prop_g_helio", "prop_g_shade")

if(length(comp_vars_sem) >= 3) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("SÉRIE 4: SEM COMPLET (CA + NSCA + Tempérament)\n")
  cat("Test CONSERVATEUR: effet structure après contrôle de TOUTE la composition\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  sem_data_complete <- sem_data_full %>%
    select(WD_BA, PCA1, all_of(comp_vars_sem)) %>%
    mutate(across(everything(), scale))

  # Construire formules dynamiquement
  formula_pca <- paste("PCA1 ~", paste(comp_vars_sem, collapse = " + "))
  formula_wd_full <- paste("WD_BA ~", paste(comp_vars_sem, collapse = " + "), "+ c*PCA1")
  formula_wd_constrained <- paste("WD_BA ~", paste(comp_vars_sem, collapse = " + "), "+ 0*PCA1")

  model_sem_complete_full <- paste(formula_pca, formula_wd_full, sep = "\n  ")
  model_sem_complete_constrained <- paste(formula_pca, formula_wd_constrained, sep = "\n  ")

  fit_complete_full <- sem(model_sem_complete_full, data = sem_data_complete)
  fit_complete_constrained <- sem(model_sem_complete_constrained, data = sem_data_complete)

  # Test
  comparison_sem_complete <- anova(fit_complete_constrained, fit_complete_full)

  cat("\n=== SEM COMPLET: Test du lien résiduel Structure->WD ===\n")
  cat("Variables composition contrôlées:", paste(comp_vars_sem, collapse = ", "), "\n\n")
  print(comparison_sem_complete)
  export_result("SEM COMPLET: Test lien résiduel Structure->WD", comparison_sem_complete)

  sem_complete_params <- parameterEstimates(fit_complete_full, standardized = TRUE) %>%
    filter(op == "~") %>%
    select(lhs, op, rhs, est, std.all, pvalue)

  cat("\n=== Coefficients SEM COMPLET ===\n")
  print(sem_complete_params)
  export_result("SEM COMPLET: Coefficients", sem_complete_params)
}

# ==============================================================================
# SYNTHÈSE DES RÉSULTATS SEM
# ==============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SYNTHÈSE SEM: Lien résiduel Structure->WD selon les modèles\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

sem_summary <- data.frame(
  Modele = character(),
  p_value_lien = numeric(),
  Significatif = character(),
  stringsAsFactors = FALSE
)

if(exists("comparison_sem_temp")) {
  sem_summary <- rbind(sem_summary, data.frame(
    Modele = "Tempérament (3 cat.)",
    p_value_lien = comparison_sem_temp$`Pr(>Chisq)`[2],
    Significatif = ifelse(comparison_sem_temp$`Pr(>Chisq)`[2] < 0.05, "OUI", "NON")
  ))
}

if(exists("comparison_sem_ca")) {
  sem_summary <- rbind(sem_summary, data.frame(
    Modele = "CA (composition)",
    p_value_lien = comparison_sem_ca$`Pr(>Chisq)`[2],
    Significatif = ifelse(comparison_sem_ca$`Pr(>Chisq)`[2] < 0.05, "OUI", "NON")
  ))
}

if(exists("comparison_sem_nsca")) {
  sem_summary <- rbind(sem_summary, data.frame(
    Modele = "NSCA_DBH (dominance)",
    p_value_lien = comparison_sem_nsca$`Pr(>Chisq)`[2],
    Significatif = ifelse(comparison_sem_nsca$`Pr(>Chisq)`[2] < 0.05, "OUI", "NON")
  ))
}

if(exists("comparison_sem_complete")) {
  sem_summary <- rbind(sem_summary, data.frame(
    Modele = "COMPLET (toutes dim.)",
    p_value_lien = comparison_sem_complete$`Pr(>Chisq)`[2],
    Significatif = ifelse(comparison_sem_complete$`Pr(>Chisq)`[2] < 0.05, "OUI", "NON")
  ))
}

cat("Résultats:\n")
print(sem_summary)
cat("\nINTERPRÉTATION:\n")
cat("  - Si Significatif = OUI: Lien résiduel structure->WD existe\n")
cat("  - Si Significatif = NON: Tout l'effet passe par la composition\n\n")

export_result("SEM: SYNTHÈSE - Lien résiduel Structure->WD", sem_summary)
write_csv2(sem_summary, file.path(path_tables, "sem_summary_all_models.csv"))

#### 1️⃣3️⃣ DIFFÉRENCES STRUCTURELLES PAR WD ----

# 13A. PERMANOVA ----
# OBJECTIF: Tester si les parcelles groupées par quartiles de WD présentent
# des structures verticales significativement différentes (approche multivariée)
#
# INTERPRÉTATION:
# - R² élevé = Les quartiles de WD expliquent une forte variance dans la structure
# - p < 0.05 = Différence significative de structure entre quartiles de WD
# - Complète le test de Mantel en testant les différences de groupe plutôt que la corrélation continue

permanova_data <- model_data %>%
  select(plot_name, WD_BA, all_of(structure_cols)) %>%
  drop_na()

# Créer quartiles de WD
permanova_data <- permanova_data %>%
  mutate(WD_quartile = cut(WD_BA,
                           breaks = quantile(WD_BA, c(0, 0.25, 0.5, 0.75, 1)),
                           labels = c("Q1", "Q2", "Q3", "Q4"),
                           include.lowest = TRUE))

# Matrice de distance euclidienne sur variables structurelles standardisées
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
# OBJECTIF: Tester si les distances entre parcelles dans l'espace structural
# sont corrélées avec les distances dans l'espace WD (approche continue)
#
# INTERPRÉTATION:
# - r positif = Les parcelles similaires en WD ont aussi des structures similaires
# - r négatif = Relation inverse (rare dans ce contexte)
# - p < 0.05 = Corrélation significative entre distances structurelles et distances de WD
# - Complète PERMANOVA: Mantel teste la corrélation continue, PERMANOVA teste les différences de groupes

dist_wd <- dist(permanova_data$WD_BA)
mantel_result <- mantel(dist_matrix, dist_wd, method = "pearson", permutations = 999)

cat("\n=== Test de Mantel ===\n")
print(mantel_result)
export_result("MANTEL: Distance structurelle vs Distance WD", mantel_result)

# 13C. Statistiques descriptives par quartile ----
# Examiner les valeurs moyennes des métriques clés dans chaque quartile de WD

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
# OBJECTIF: Identifier quelles métriques structurelles individuelles diffèrent
# significativement entre quartiles de WD (tests univariés)
#
# INTERPRÉTATION:
# - p < 0.05 = Cette métrique varie significativement selon les quartiles de WD
# - Métriques avec p faible = Candidats importants pour expliquer la relation WD-structure
# - Complète PERMANOVA: Tests univariés vs test multivarié global

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
# OBJECTIF: Tester si la prise en compte de l'autocorrélation spatiale améliore
# le modèle de prédiction de WD par la structure verticale
#
# INTERPRÉTATION:
# - AIC spatial < AIC OLS = L'autocorrélation spatiale améliore le modèle
# - Différence AIC > 2 = Amélioration substantielle
# - Si modèle spatial meilleur = Nécessité de contrôler l'effet spatial dans les analyses

# 14A. Charger les coordonnées ----
path_plots_gpkg <- file.path(dirname(project_dir), "LiDAR_functionnal","0_Inventories_plot_preparation", "final", "plots_unique")

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
# 14B. Modèle spatial (si coordonnées disponibles) ----
if("X_utm" %in% names(model_data) && sum(!is.na(model_data$X_utm)) > 50) {
  
  spatial_data <- model_data %>%
    filter(!is.na(X_utm) & !is.na(Y_utm) & !is.na(WD_BA) & !is.na(!!sym(best_var)))
  
  # Modèle OLS classique
  formula_best <- as.formula(paste("WD_BA ~", best_var))
  model_ols <- lm(formula_best, data = spatial_data)
  
  # Modèle spatial avec corrélation Matérn
  model_spatial <- fitme(
    as.formula(paste("WD_BA ~", best_var, "+ Matern(1 | X_utm + Y_utm)")),
    data = spatial_data,
    method = "REML"
  )
  
  # Extraire les AIC correctement
  aic_ols <- AIC(model_ols)
  aic_spatial_all <- AIC(model_spatial)
  aic_spatial <- aic_spatial_all[1]  # Premier élément = marginal AIC (mAIC)
  
  cat("\n=== Modèle OLS vs Spatial ===\n")
  cat("AIC OLS:", round(aic_ols, 2), "\n")
  cat("AIC Spatial (mAIC):", round(aic_spatial, 2), "\n")
  cat("Différence (OLS - Spatial):", round(aic_ols - aic_spatial, 2), "\n")
  
  if(aic_spatial < aic_ols) {
    cat("→ Le modèle spatial est meilleur (ΔAIC =", round(aic_ols - aic_spatial, 1), ")\n")
  } else {
    cat("→ Le modèle OLS est suffisant\n")
  }
  
  spatial_comparison <- data.frame(
    Modele = c("OLS", "Spatial_Matern"),
    AIC = c(aic_ols, aic_spatial),
    Delta_AIC = c(0, aic_spatial - aic_ols)
  )
  
  cat("\nTableau de comparaison:\n")
  print(spatial_comparison)
  
  export_result("COMPARAISON: OLS vs Spatial", spatial_comparison)
  export_result("MODÈLE SPATIAL (détail)", summary(model_spatial))
  
  write_csv2(spatial_comparison, file.path(path_tables, "spatial_models_comparison.csv"))
}

#### 1️⃣5️⃣  CROSS-VALIDATION ----

# 15B. VALIDATION CROISÉE LOOCV ----

# OBJECTIF: Évaluer la capacité prédictive réelle du modèle en testant sur des données
# non utilisées pour l'entraînement (validation croisée rigoureuse)
#
# INTERPRÉTATION:
# - RMSE faible = Bonnes prédictions en validation
# - R² proche de celui du modèle complet = Modèle stable et généralisable
# - Pente proche de 1 et intercept proche de 0 = Prédictions non biaisées
# - Pente < 1 = Sous-estimation des valeurs élevées, sur-estimation des faibles

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

# Calculer les métriques de qualité incluant pente et intercept
lm_obs_pred <- lm(predicted ~ observed, data = loocv_results)

loocv_metrics <- data.frame(
  RMSE = sqrt(mean(loocv_results$residual^2)),
  MAE = mean(abs(loocv_results$residual)),
  R2 = cor(loocv_results$observed, loocv_results$predicted)^2,
  Slope = coef(lm_obs_pred)[2],
  Intercept = coef(lm_obs_pred)[1]
)

cat("\n=== LOOCV Résultats ===\n")
print(loocv_metrics)
export_result("LOOCV: Métriques de validation", loocv_metrics)

# Visualisation publication-ready avec toutes les métriques centrées
metrics_text <- paste0(
  "RMSE = ", round(loocv_metrics$RMSE, 4), "\n",
  "R² = ", round(loocv_metrics$R2, 3), "\n",
  "MAE = ", round(loocv_metrics$MAE, 4), "\n",
  "Slope = ", round(loocv_metrics$Slope, 3), "\n",
  "Intercept = ", round(loocv_metrics$Intercept, 3)
)

# Calculer les limites des axes pour coordonnées égales
range_vals <- range(c(loocv_results$observed, loocv_results$predicted))
axis_limits <- c(floor(range_vals[1] * 20) / 20, ceiling(range_vals[2] * 20) / 20)

p_loocv <- ggplot(loocv_results, aes(x = observed, y = predicted)) +
  geom_point(alpha = 0.6, size = 3, color = "#2C3E50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "#3498DB", linewidth = 1) +
  annotate("text", x = mean(axis_limits), y = axis_limits[1] + 0.85 * diff(axis_limits),
           label = metrics_text, hjust = 0.5, vjust = 0.5,
           size = 4.5, fontface = "bold", family = "sans") +
  labs(title = "Leave-One-Out Cross-Validation",
       subtitle = paste0("Prédiction de WD par ", best_var),
       x = "WD observé (g/cm³)",
       y = "WD prédit (g/cm³)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  coord_equal(xlim = axis_limits, ylim = axis_limits)

ggsave(file.path(path_figures, "loocv_validation.jpg"), p_loocv,
       width = 10, height = 10, dpi = 600)

write_csv2(loocv_results, file.path(path_tables, "loocv_predictions.csv"))
write_csv2(loocv_metrics, file.path(path_tables, "loocv_metrics.csv"))

# 15B. VALIDATION CROISÉE SPATIALE (Leave-Site-Out) ----
# OBJECTIF: Évaluer la capacité de généralisation à de NOUVEAUX sites
# Plus conservateur que LOOCV car teste la transférabilité géographique

cat("\n=== Leave-Site-Out Cross-Validation ===\n")

# Extraire le site depuis le nom du plot (ex: "MAB01h" -> "MAB")
# Extraire le site avec la correspondance fournie
loocv_spatial_data <- model_data %>%
  select(plot_name, WD_BA, all_of(best_var), X_utm, Y_utm) %>%
  drop_na() %>%
  mutate(site = case_when(
    # Malebo : chiffres simples + M2P*
    plot_name %in% c("123", "148", "162", "166", "180", "184", "189", "199", "41", "53") ~ "Malebo",
    str_detect(plot_name, "^M2P") ~ "Malebo",
    
    # Yangambi : GIL*, JEU*, MIX*
    str_detect(plot_name, "^GIL|^JEU|^MIX") ~ "Yangambi",
    
    # Monkoto : Lokofa*, Betamba*
    str_detect(plot_name, "^Lokofa|^Betamba") ~ "Monkoto",
    
    # Mabounie : MAB*
    str_detect(plot_name, "^MAB") ~ "Mabounie",
    
    # Lac Mai Ndombe : FRM 289, 373, 408
    plot_name %in% c("FRM 289", "FRM 373", "FRM 408") ~ "Lac_Mai_Ndombe",
    
    # Mbala : FRM 546, 627
    plot_name %in% c("FRM 546", "FRM 627") ~ "Mbala",
    
    # Kiri : FRM 162, 227
    plot_name %in% c("FRM 162", "FRM 227") ~ "Kiri",
    
    # Bongimba : FRM 15, 39, 80, 102
    plot_name %in% c("FRM 15", "FRM 39", "FRM 80", "FRM 102") ~ "Bongimba",
    
    # Kabambare
    plot_name %in% c("PARAP_Maniema_178", "Modele_4503_Principal") ~ "Kabambare",
    
    # Autres = sites individuels (nom du plot = site)
    TRUE ~ plot_name
  ))

table(loocv_spatial_data$site)

# Vérifier les sites
sites <- unique(loocv_spatial_data$site)
cat("Sites identifiés:", length(sites), "\n")
cat(paste(sites, collapse = ", "), "\n\n")

# Leave-Site-Out CV
predictions_spatial <- numeric(nrow(loocv_spatial_data))

for(s in sites) {
  train_idx <- loocv_spatial_data$site != s
  test_idx <- loocv_spatial_data$site == s
  
  if(sum(train_idx) >= 10 && sum(test_idx) >= 1) {
    model_cv <- lm(as.formula(paste("WD_BA ~", best_var)), 
                   data = loocv_spatial_data[train_idx, ])
    predictions_spatial[test_idx] <- predict(model_cv, 
                                             newdata = loocv_spatial_data[test_idx, ])
  }
}

loocv_spatial_results <- data.frame(
  observed = loocv_spatial_data$WD_BA,
  predicted = predictions_spatial,
  site = loocv_spatial_data$site,
  residual = loocv_spatial_data$WD_BA - predictions_spatial
)

# Métriques
lm_spatial <- lm(predicted ~ observed, data = loocv_spatial_results)

loocv_spatial_metrics <- data.frame(
  Method = "Leave-Site-Out",
  RMSE = sqrt(mean(loocv_spatial_results$residual^2, na.rm = TRUE)),
  MAE = mean(abs(loocv_spatial_results$residual), na.rm = TRUE),
  R2 = cor(loocv_spatial_results$observed, loocv_spatial_results$predicted, 
           use = "complete.obs")^2,
  Slope = coef(lm_spatial)[2],
  Intercept = coef(lm_spatial)[1]
)

# Comparaison avec LOOCV simple
loocv_comparison <- rbind(
  data.frame(Method = "LOOCV simple", loocv_metrics),
  loocv_spatial_metrics
)

cat("\n=== Comparaison LOOCV simple vs Leave-Site-Out ===\n")
print(loocv_comparison)
export_result("LOOCV: Comparaison simple vs spatial", loocv_comparison)

# Figure
p_loocv_spatial <- ggplot(loocv_spatial_results, aes(x = observed, y = predicted)) +
  geom_point(aes(color = site), alpha = 0.7, size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  labs(title = "Leave-Site-Out Cross-Validation",
       subtitle = paste0("RMSE = ", round(loocv_spatial_metrics$RMSE, 4),
                         ", R² = ", round(loocv_spatial_metrics$R2, 3)),
       x = "WD observé (g/cm³)",
       y = "WD prédit (g/cm³)",
       color = "Site") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right") +
  coord_equal()

ggsave(file.path(path_figures, "loocv_spatial_validation.jpg"), p_loocv_spatial,
       width = 12, height = 10, dpi = 600)

write_csv2(loocv_comparison, file.path(path_tables, "loocv_comparison.csv"))


#### 1️⃣6️⃣ANALYSE COMPLÉMENTAIRES ----
# 16A. IDENTIFICATION DES OUTLIERS ----
# OBJECTIF: Identifier les plots où la prédiction est mauvaise
# Potentiellement Gilbertiodendron ou forêts secondaires extrêmes

cat("\n=== Identification des outliers ===\n")

# Calculer résidus standardisés
outlier_data <- model_data %>%
  select(plot_name, WD_BA, all_of(best_var)) %>%
  drop_na() %>%
  mutate(
    predicted = predict(best_model, newdata = .),
    residual = WD_BA - predicted,
    residual_std = scale(residual)[,1]
  ) %>%
  arrange(desc(abs(residual_std)))

# Outliers = résidus standardisés > 2
outliers <- outlier_data %>%
  filter(abs(residual_std) > 2)

cat("\nPlots avec résidus standardisés > 2:\n")
print(outliers %>% select(plot_name, WD_BA, predicted, residual, residual_std))
export_result("OUTLIERS: Plots mal prédits", outliers)

# Catégoriser les outliers
outliers_summary <- outliers %>%
  mutate(
    type = case_when(
      residual > 0 ~ "WD sous-estimé (observé > prédit)",
      residual < 0 ~ "WD sur-estimé (observé < prédit)"
    )
  )

cat("\nRésumé des outliers:\n")
print(table(outliers_summary$type))

# Sauvegarder
write_csv2(outlier_data, file.path(path_tables, "residuals_all_plots.csv"))
write_csv2(outliers, file.path(path_tables, "outliers_identified.csv"))

# 16B. MODÈLE SPATIAL AVEC COMPOSITION ----
# OBJECTIF: Tester si l'autocorrélation spatiale disparaît quand on ajoute la composition
# Si oui → la structure spatiale vient de la composition

if("X_utm" %in% names(model_data) && sum(!is.na(model_data$X_utm)) > 50) {
  
  cat("\n=== Modèle spatial avec composition ===\n")
  
  spatial_data_full <- model_data %>%
    filter(!is.na(X_utm) & !is.na(Y_utm) & !is.na(WD_BA) & 
             !is.na(!!sym(best_var)) & !is.na(CA1) & !is.na(CA2))
  
  # Modèle 1: Structure seule + Matérn
  model_spatial_struct <- fitme(
    as.formula(paste("WD_BA ~", best_var, "+ Matern(1 | X_utm + Y_utm)")),
    data = spatial_data_full, method = "REML"
  )
  
  # Modèle 2: Structure + Composition + Matérn
  model_spatial_full <- fitme(
    as.formula(paste("WD_BA ~", best_var, "+ CA1 + CA2 + Matern(1 | X_utm + Y_utm)")),
    data = spatial_data_full, method = "REML"
  )
  
  # Modèle 3: Structure + Composition SANS Matérn (pour comparer)
  model_ols_full <- lm(
    as.formula(paste("WD_BA ~", best_var, "+ CA1 + CA2")),
    data = spatial_data_full
  )
  
  # Extraire lambda (variance spatiale)
  lambda_struct <- VarCorr(model_spatial_struct)$`X_utm + Y.`$`(Intercept)`
  lambda_full <- VarCorr(model_spatial_full)$`X_utm + Y.`$`(Intercept)`
  
  spatial_evolution <- data.frame(
    Modele = c("Structure + Matérn", "Structure + Compo + Matérn", "Structure + Compo (OLS)"),
    AIC = c(AIC(model_spatial_struct)[1], AIC(model_spatial_full)[1], AIC(model_ols_full)),
    Lambda_spatial = c(lambda_struct, lambda_full, NA)
  )
  
  cat("\nÉvolution de la variance spatiale:\n")
  print(spatial_evolution)
  export_result("SPATIAL: Évolution avec composition", spatial_evolution)
  
  reduction_pct <- (1 - lambda_full / lambda_struct) * 100
  cat("\nRéduction de la variance spatiale avec composition:", round(reduction_pct, 1), "%\n")
  
  if(reduction_pct > 50) {
    cat("→ La structure spatiale est largement expliquée par la composition\n")
  } else {
    cat("→ La structure spatiale persiste même après contrôle de la composition\n")
  }
  
  write_csv2(spatial_evolution, file.path(path_tables, "spatial_models_evolution.csv"))
}

#### SYNTHÈSE FINALE ----

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
