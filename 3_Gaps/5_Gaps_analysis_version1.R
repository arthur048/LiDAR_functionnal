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
pkgs <- c("tidyverse", "rio", "ggpubr", "viridis", "corrplot", "gridExtra", "mgcv", "patchwork", "here")

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
#'
#' @param data Données d'un plot filtré à buffer 0
#' @param threshold Seuil de proportion (entre 0 et 1)
#' @return Hauteur à laquelle la proportion atteint le seuil
gap_height_at_threshold <- function(data, threshold) {
  # Vérifier que les données ont au moins 2 lignes
  if (nrow(data) < 2) return(NA)
  
  # Ordonner les données par hauteur
  data <- data %>% arrange(height_aboveground)
  
  # Rechercher où la proportion passe au-dessus du seuil (et non en dessous)
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
#'
#' @param data Données d'un plot filtré à buffer 0
#' @return Pente de la régression linéaire (taux de décroissance)
calculate_decay_rate <- function(data) {
  # Vérifier que les données ont au moins 2 lignes
  if (nrow(data) < 3) return(NA)
  
  # Ajuster un modèle linéaire
  model <- lm(proportion ~ height_aboveground, data = data)
  
  # Récupérer la pente (taux de décroissance)
  decay_rate <- coef(model)[2]
  return(decay_rate)
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
buffer0_data <- gaps_metrics %>%
  filter(buffer == 0)

# Joindre avec les infos des plots
if ("WD_mean" %in% colnames(plots_info)) {
  combined_data <- buffer0_data %>%
    left_join(plots_info %>% select(plot_name, WD_mean), by = "plot_name") %>%
    filter(!is.na(WD_mean))  # Garder uniquement les plots avec WD_mean
} else {
  stop("La colonne WD_mean n'est pas disponible dans le fichier plots_info.")
}

#### 5️⃣ ANALYSE ####
# Calculer diverses métriques pour chaque plot
plot_metrics <- tibble(
  plot_name = character(),
  WD_mean = numeric(),
  prop_at_10m = numeric(),
  prop_at_20m = numeric(),
  prop_at_30m = numeric(),
  height_at_50pct = numeric(),
  height_at_25pct = numeric(),
  height_at_10pct = numeric(),
  decay_rate = numeric(),
  auc = numeric()
)

for (plot_i in plot_names) {
  # Filtrer les données pour ce plot
  plot_data <- combined_data %>%
    filter(plot_name == plot_i)
  
  if (nrow(plot_data) < 3) next  # Sauter les plots avec trop peu de données
  
  # Récupérer WD_mean
  wd_mean_val <- plot_data$WD_mean[1]
  
  # Proportion à hauteurs spécifiques
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
  
  # Hauteurs à seuils spécifiques
  height_50pct <- gap_height_at_threshold(plot_data, 0.5)
  height_25pct <- gap_height_at_threshold(plot_data, 0.25)
  height_10pct <- gap_height_at_threshold(plot_data, 0.1)
  
  # Taux de décroissance
  decay <- calculate_decay_rate(plot_data)
  
  # Aire sous la courbe
  auc_val <- calculate_auc(plot_data)
  
  # Ajouter au dataframe de métriques
  plot_metrics <- plot_metrics %>%
    bind_rows(tibble(
      plot_name = plot_i,
      WD_mean = wd_mean_val,
      prop_at_10m = prop_10m,
      prop_at_20m = prop_20m,
      prop_at_30m = prop_30m,
      height_at_50pct = height_50pct,
      height_at_25pct = height_25pct,
      height_at_10pct = height_10pct,
      decay_rate = decay,
      auc = auc_val
    ))
}

# Exporter les métriques calculées
rio::export(
  plot_metrics,
  paste0(path_analysis, "plot_gap_metrics_buffer0.csv"),
  sep = ";",
  dec = ",",
  append = FALSE
)

#### 6️⃣ TESTS DE CORRÉLATION ####
# Calculer les corrélations entre métriques et WD_mean
correlation_vars <- c("prop_at_10m", "prop_at_20m", "prop_at_30m", 
                      "height_at_50pct", "height_at_25pct", "height_at_10pct",
                      "decay_rate", "auc")

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
  prop_at_10m = "Proportion de trouées à 10m",
  prop_at_20m = "Proportion de trouées à 20m",
  prop_at_30m = "Proportion de trouées à 30m",
  height_at_50pct = "Hauteur à 50% de trouées",
  height_at_25pct = "Hauteur à 25% de trouées",
  height_at_10pct = "Hauteur à 10% de trouées",
  decay_rate = "Taux de décroissance",
  auc = "Aire sous la courbe"
)

# Ajouter les noms lisibles
correlations$metric_label <- metric_labels[correlations$metric]

# Exporter les corrélations
rio::export(
  correlations,
  paste0(path_analysis, "gap_wd_correlations_buffer0.csv"),
  sep = ";",
  dec = ",",
  append = FALSE
)

#### 7️⃣ VISUALISATION ####
# 1. Graphique de corrélation
corr_plot <- ggplot(correlations %>% filter(method == "Pearson"), 
                    aes(x = metric_label, y = correlation, fill = p_value < 0.05)) +
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
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold")
  )

# Sauvegarder le graphique de corrélation
ggsave(
  paste0(path_analysis, "gap_wd_correlations.jpg"),
  corr_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# 2. Sélectionner les 2-3 métriques avec les corrélations les plus fortes
top_correlations <- correlations %>%
  filter(method == "Pearson") %>%
  arrange(desc(abs(correlation))) %>%
  slice(1:3)

top_metrics <- top_correlations$metric

scatter_plots <- list()
for (i in seq_along(top_metrics)) {
  metric <- top_metrics[i]
  
  # Créer un graphique de dispersion
  p <- ggplot(plot_metrics, aes_string(x = "WD_mean", y = metric)) +
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
  mod <- lm(reformulate("WD_mean", metric), data = plot_metrics)
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2, 
                   list(a = format(coef(mod)[1], digits = 3),
                        b = format(coef(mod)[2], digits = 3),
                        r2 = format(summary(mod)$r.squared, digits = 3)))
  
  p <- p + annotate("text", x = max(plot_metrics$WD_mean, na.rm = TRUE) * 0.8, 
                    y = max(plot_metrics[[metric]], na.rm = TRUE) * 0.9, 
                    label = as.character(as.expression(eq)), parse = TRUE, size = 3)
  
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
combined_scatter <- gridExtra::grid.arrange(grobs = scatter_plots, ncol = 2)

# Sauvegarder
ggsave(
  paste0(path_analysis, "combined_scatter_plots.jpg"),
  combined_scatter,
  width = 14,
  height = 10,
  dpi = 300
)

# 3. Examiner comment la courbe complète varie avec WD
# Créer des catégories de WD pour faciliter la visualisation
wd_quantiles <- quantile(plot_metrics$WD_mean, probs = seq(0, 1, 0.25), na.rm = TRUE)
combined_data <- combined_data %>%
  mutate(wd_category = cut(WD_mean, 
                           breaks = wd_quantiles, 
                           include.lowest = TRUE,
                           labels = c("Q1 (WD bas)", "Q2", "Q3", "Q4 (WD élevé)")))

# Visualiser les courbes moyennes par catégorie de WD
curve_plot <- ggplot(combined_data, aes(x = height_aboveground, y = proportion, color = wd_category, group = interaction(plot_name, wd_category))) +
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
  paste0(path_analysis, "gap_curves_by_wd_category.jpg"),
  curve_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# 4. Heatmap des corrélations entre toutes les métriques
cor_matrix <- cor(plot_metrics %>% select(-plot_name), use = "pairwise.complete.obs")

# Visualiser la matrice de corrélation
corrplot::corrplot(cor_matrix, method = "color", type = "upper", 
                   order = "hclust", tl.col = "black", tl.srt = 45,
                   addCoef.col = "black", number.cex = 0.7,
                   col = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200),
                   title = "Matrice de corrélation entre métriques de trouées et WD")

# Sauvegarder
jpeg(paste0(path_analysis, "correlation_matrix.jpg"), width = 10, height = 8, units = "in", res = 300)
corrplot::corrplot(cor_matrix, method = "color", type = "upper", 
                   order = "hclust", tl.col = "black", tl.srt = 45,
                   addCoef.col = "black", number.cex = 0.7,
                   col = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200),
                   title = "Matrice de corrélation entre métriques de trouées et WD")
dev.off()

# 5. Visualisation des métriques par boxplot selon les catégories de WD
# Préparer les données au format long pour les boxplots
plot_metrics_long <- plot_metrics %>%
  mutate(wd_category = cut(WD_mean, 
                           breaks = wd_quantiles, 
                           include.lowest = TRUE,
                           labels = c("Q1 (WD bas)", "Q2", "Q3", "Q4 (WD élevé)"))) %>%
  pivot_longer(cols = c(prop_at_10m, prop_at_20m, prop_at_30m, 
                        height_at_50pct, height_at_25pct, height_at_10pct,
                        decay_rate, auc),
               names_to = "metric", 
               values_to = "value") %>%
  mutate(metric_label = metric_labels[metric])

# Créer des boxplots pour les métriques principales
box_plots <- list()

for (i in seq_along(top_metrics)) {
  metric <- top_metrics[i]
  
  box_p <- ggplot(plot_metrics_long %>% filter(metric == !!metric), 
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
  paste0(path_analysis, "boxplots_by_wd_category.jpg"),
  combined_boxes,
  width = 12,
  height = 8,
  dpi = 300
)

# 6. Visualisation 3D : GAM pour représenter la relation entre hauteur, WD et proportion
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
  paste0(path_analysis, "contour_plot_wd_height_proportion.jpg"),
  contour_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# Résumé des résultats
cat("\n===== RÉSUMÉ DES ANALYSES =====\n\n")
cat("Métriques calculées et exportées dans:", paste0(path_analysis, "plot_gap_metrics_buffer0.csv"), "\n")
cat("Corrélations calculées et exportées dans:", paste0(path_analysis, "gap_wd_correlations_buffer0.csv"), "\n\n")

cat("Top 3 métriques les plus corrélées avec WD_mean:\n")
print(top_correlations %>% select(metric_label, correlation, p_value))

cat("\nGraphiques générés dans le dossier:", path_analysis, "\n")
