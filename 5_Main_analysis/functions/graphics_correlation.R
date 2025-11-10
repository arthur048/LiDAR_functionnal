# ===================================================================
# FONCTIONS GRAPHIQUES POUR L'ANALYSE DE CORRÉLATION
# Objectif : Centraliser les fonctions de visualisation pour les analyses
#            de corrélation, avec une approche modulaire et extensible
# ===================================================================

#' @title Créer un graphique en barres des corrélations triées
#' @description Génère un graphique en barres des corrélations entre métriques et une variable cible, 
#'             triées par valeur mais sélectionnées par valeur absolue
#' @param corr_data Tableau de données contenant les corrélations
#' @param target_var Variable cible (ex: WD_BA, WD_BA)
#' @param corr_method Méthode de corrélation ("Pearson" ou "Spearman")
#' @param top_n Nombre de métriques à inclure (les plus corrélées en valeur absolue)
#' @param use_full_names Utiliser les noms complets des métriques (TRUE) ou les abréviations (FALSE)
#' @return Objet ggplot
plot_correlation_bars <- function(corr_data, target_var, corr_method, top_n = 20, use_full_names = TRUE) {
  # Filtrer les données et sélectionner top_n métriques en valeur absolue
  filtered_data <- corr_data %>%
    filter(target_variable == target_var, method == corr_method, !is.na(correlation), !is.na(p_value)) %>%
    mutate(abs_corr = abs(correlation)) %>%
    arrange(desc(abs_corr)) %>%
    slice(1:min(top_n, n())) %>%
    # Réordonner par valeur réelle pour l'affichage
    arrange(correlation)
  
  # Choisir la colonne pour les étiquettes
  label_col <- if(use_full_names) "metric_label" else "metric"
  
  # Créer catégories de significativité
  filtered_data <- filtered_data %>%
    mutate(
      significance = case_when(
        p_value >= 0.05 ~ "p ≥ 0.05",
        p_value < 0.05 & p_value >= 0.01 ~ "p < 0.05",
        p_value < 0.01 & p_value >= 0.001 ~ "p < 0.01",
        p_value < 0.001 ~ "p < 0.001"
      ),
      significance = factor(significance, levels = c("p ≥ 0.05", "p < 0.05", "p < 0.01", "p < 0.001"))
    )
  
  # Calculer la plage pour les limites d'axe
  max_abs_corr <- max(abs(filtered_data$correlation), na.rm = TRUE)
  axis_limit <- ceiling(max_abs_corr * 10) / 10
  
  # Palette de couleurs plus adaptée à la publication (bleu vers rouge)
  significance_colors <- c(
    "p ≥ 0.05" = "#cccccc",       # Gris pour non significatif
    "p < 0.05" = "#6baed6",       # Bleu clair 
    "p < 0.01" = "#3182bd",       # Bleu moyen
    "p < 0.001" = "#08519c"       # Bleu foncé
  )
  
  # Créer le graphique
  p <- ggplot(
    filtered_data,
    aes(x = reorder(!!sym(label_col), correlation), y = correlation, fill = significance)
  ) +
    geom_bar(stat = "identity", color = "darkblue", linewidth = 0.2) +
    coord_flip() +
    scale_fill_manual(
      values = significance_colors,
      name = "Significativité"
    ) +
    scale_y_continuous(limits = c(-axis_limit, axis_limit)) +
    labs(
      title = paste0("Top ", nrow(filtered_data), " des corrélations (", corr_method, ") avec ", target_var),
      x = "",
      y = paste("Coefficient de corrélation de", corr_method)
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  
  return(p)
}


#' @title Créer un scatterplot entre une métrique et une variable cible
#' @description Génère un graphique de dispersion entre une métrique et une variable cible,
#'              incluant régression linéaire et statistiques
#' @param data Jeu de données
#' @param metric_name Nom de la métrique
#' @param target_var Variable cible (ex: WD_BA, WD_BA)
#' @param metric_labels Dictionnaire des étiquettes des métriques
#' @return Liste contenant l'objet ggplot et les résultats du modèle linéaire
plot_correlation_scatter <- function(data, metric_name, target_var, metric_labels) {
  # Utiliser l'étiquette complète si disponible
  metric_label <- NULL
  if (metric_name %in% names(metric_labels)) {
    metric_label <- metric_labels[[metric_name]]
  } else {
    metric_label <- metric_name
  }
  
  # Données valides uniquement
  valid_data <- data %>% 
    filter(!is.na(.data[[metric_name]]), !is.na(.data[[target_var]]))
  
  # Calculer le modèle linéaire
  formula_str <- paste(target_var, "~", metric_name)
  lm_model <- lm(as.formula(formula_str), data = valid_data)
  lm_summary <- summary(lm_model)
  
  # Extraire les coefficients et statistiques
  intercept <- round(coef(lm_model)[1], 4)
  slope <- round(coef(lm_model)[2], 4)
  r_squared <- round(lm_summary$r.squared, 4)
  adj_r_squared <- round(lm_summary$adj.r.squared, 4)
  f_stat <- round(lm_summary$fstatistic[1], 2)
  p_value <- round(lm_summary$coefficients[2, 4], 4)
  correlation <- round(cor(valid_data[[metric_name]], valid_data[[target_var]], use = "complete.obs"), 4)
  
  # Trouver les min/max pour le scaling des axes
  x_range <- range(valid_data[[metric_name]], na.rm = TRUE)
  y_range <- range(valid_data[[target_var]], na.rm = TRUE)
  
  # Ajouter une marge de 5% aux limites
  x_margin <- diff(x_range) * 0.05
  y_margin <- diff(y_range) * 0.05
  x_limits <- c(x_range[1] - x_margin, x_range[2] + x_margin)
  y_limits <- c(y_range[1] - y_margin, y_range[2] + y_margin)
  
  # Équation de la régression
  eq_text <- paste0("y = ", slope, "x + ", intercept)
  stat_text <- paste0("R² = ", r_squared, ", p = ", p_value)
  
  # Créer le graphique
  p <- ggplot(valid_data, aes(x = .data[[metric_name]], y = .data[[target_var]])) +
    geom_point(alpha = 0.8, size = 3, color = "darkblue") +
    geom_smooth(method = "lm", se = TRUE, color = "firebrick", alpha = 0.2) +
    labs(
      title = paste("Relation entre", metric_label, "et", target_var),
      subtitle = paste(eq_text, "\n", stat_text),
      x = metric_label,
      y = target_var
    ) +
    xlim(x_limits) +
    ylim(y_limits) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, color = "darkblue", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  # Créer une liste avec le graphique et les statistiques
  result <- list(
    plot = p,
    stats = data.frame(
      metric = metric_name,
      metric_label = metric_label,
      target_var = target_var,
      intercept = intercept,
      slope = slope,
      r_squared = r_squared,
      adj_r_squared = adj_r_squared,
      f_statistic = f_stat,
      p_value = p_value,
      correlation = correlation,
      n_observations = nrow(valid_data)
    )
  )
  
  return(result)
}

#' @title Générer des visualisations batch pour les métriques les plus corrélées
#' @description Génère des graphiques de dispersion pour les métriques les plus corrélées
#' @param data Jeu de données
#' @param corr_data Tableau de données contenant les corrélations
#' @param target_var Variable cible (ex: WD_BA, WD_BA)
#' @param metric_labels Dictionnaire des étiquettes des métriques
#' @param base_dir Répertoire de base pour les sorties
#' @param path_target_var Répertoire pour la variable cible
#' @param corr_method Méthode de corrélation ("Pearson" ou "Spearman")
#' @param top_n Nombre de métriques à inclure (les plus corrélées en valeur absolue)
#' @param dpi Résolution en DPI
#' @return Dataframe avec les statistiques des modèles linéaires
generate_correlation_scatterplots <- function(data, corr_data, target_var, metric_labels, 
                                              base_dir, path_target_var, corr_method = "Pearson", 
                                              top_n = 10, dpi = 600) {
  # Filtrer les métriques les plus corrélées en valeur absolue
  filtered_corr <- corr_data %>%
    filter(target_variable == target_var, method == corr_method, !is.na(correlation)) %>%
    mutate(abs_corr = abs(correlation)) %>%
    arrange(desc(abs_corr)) %>%
    slice(1:min(top_n, n()))
  
  top_metrics <- filtered_corr$metric
  
  # Créer un objet combinant tous les graphiques
  scatter_plots <- list()
  all_stats <- data.frame()
  
  # Générer chaque graphique
  for (i in seq_along(top_metrics)) {
    metric <- top_metrics[i]
    result <- plot_correlation_scatter(data, metric, target_var, metric_labels)
    scatter_plots[[i]] <- result$plot
    all_stats <- bind_rows(all_stats, result$stats)
    
    # Sauvegarder individuellement dans le dossier de la variable cible
    ggsave(
      file.path(path_target_var, paste0("scatter_", metric, ".jpg")),
      result$plot,
      width = 8,
      height = 6,
      dpi = dpi
    )
  }
  
  # Sauvegarder les statistiques dans le dossier de la variable cible
  write_csv2(
    all_stats,
    file.path(path_target_var, paste0("lm_stats_top", top_n, "_", corr_method, ".csv"))
  )
  
  # Créer un graphique combiné
  n_rows <- ceiling(length(scatter_plots) / 2)
  combined_plot <- wrap_plots(scatter_plots, ncol = 2, nrow = n_rows)
  
  # Sauvegarder le graphique combiné dans le dossier de la variable cible
  ggsave(
    file.path(path_target_var, paste0("combined_scatter_top", top_n, "_", corr_method, ".jpg")),
    combined_plot,
    width = 16,
    height = n_rows * 5,
    dpi = dpi,
    limitsize = FALSE
  )
  
  return(all_stats)
}

#' @title Comparer les corrélations entre deux variables cibles
#' @description Crée un graphique de comparaison des corrélations entre deux variables cibles
#' @param corr_data Tableau de données contenant les corrélations
#' @param var1 Première variable cible
#' @param var2 Deuxième variable cible
#' @param metric_labels Dictionnaire des étiquettes des métriques
#' @param corr_method Méthode de corrélation ("Pearson" ou "Spearman")
#' @param top_n Nombre de métriques à inclure (les plus corrélées en valeur absolue avec var1)
#' @param use_full_names Utiliser les noms complets des métriques (TRUE) ou les abréviations (FALSE)
#' @return Liste contenant le graphique et le dataframe de comparaison
compare_correlations <- function(corr_data, var1, var2, metric_labels, corr_method = "Pearson", 
                                 top_n = NULL, use_full_names = TRUE) {
  # Vérifier si les variables existent dans les données
  if (!var1 %in% unique(corr_data$target_variable)) {
    warning(paste("Variable", var1, "non trouvée dans les données de corrélation"))
    return(NULL)
  }
  
  if (!var2 %in% unique(corr_data$target_variable)) {
    warning(paste("Variable", var2, "non trouvée dans les données de corrélation"))
    return(NULL)
  }
  
  # Filtrer les métriques pour les deux variables
  var1_corr <- corr_data %>%
    filter(target_variable == var1, method == corr_method) %>%
    dplyr::select(metric, correlation, p_value, metric_label)
  
  var2_corr <- corr_data %>%
    filter(target_variable == var2, method == corr_method) %>%
    dplyr::select(metric, correlation, p_value, metric_label)
  
  # Vérifier si les dataframes contiennent des données
  if (nrow(var1_corr) == 0) {
    warning(paste("Aucune donnée de corrélation trouvée pour", var1))
    return(NULL)
  }
  
  if (nrow(var2_corr) == 0) {
    warning(paste("Aucune donnée de corrélation trouvée pour", var2))
    return(NULL)
  }
  
  # Préparer les données pour la comparaison
  comparison_data <- full_join(
    var1_corr, 
    var2_corr, 
    by = c("metric", "metric_label"), 
    suffix = c("_var1", "_var2")
  )
  
  # Vérifier si le join a produit des données
  if (nrow(comparison_data) == 0) {
    warning(paste("Aucune métrique commune entre", var1, "et", var2))
    return(NULL)
  }
  
  # Gérer les valeurs manquantes potentielles
  comparison_data <- comparison_data %>%
    filter(!is.na(correlation_var1) & !is.na(correlation_var2)) %>%
    mutate(
      abs_corr_var1 = abs(correlation_var1),
      abs_corr_var2 = abs(correlation_var2),
      diff = correlation_var1 - correlation_var2,
      abs_diff = abs(diff),
      # Direction (pour la visualisation)
      direction = case_when(
        correlation_var1 >= 0 & correlation_var2 >= 0 ~ "Les deux positives",
        correlation_var1 < 0 & correlation_var2 < 0 ~ "Les deux négatives",
        correlation_var1 >= 0 & correlation_var2 < 0 ~ "Var1+ / Var2-",
        correlation_var1 < 0 & correlation_var2 >= 0 ~ "Var1- / Var2+",
        TRUE ~ "NA"
      )
    ) %>%
    arrange(desc(abs_corr_var1))
  
  # Vérifier à nouveau si nous avons des données après filtrage
  if (nrow(comparison_data) == 0) {
    warning("Aucune donnée valide après filtrage des valeurs manquantes")
    return(NULL)
  }
  
  # Limiter aux top_n métriques si spécifié
  if (!is.null(top_n)) {
    comparison_data <- comparison_data %>% 
      slice(1:min(top_n, n()))
  }
  
  # Déterminer quelle colonne utiliser pour les labels
  label_col <- if(use_full_names) "metric_label" else "metric"
  
  # Déterminer les limites des axes pour une symétrie
  min_corr <- min(c(comparison_data$abs_corr_var1, comparison_data$abs_corr_var2), na.rm = TRUE)
  max_corr <- max(c(comparison_data$abs_corr_var1, comparison_data$abs_corr_var2), na.rm = TRUE)
  
  axis_limit_min <- floor(min_corr * 10) / 10  # Arrondir à la décimale inférieure
  axis_limit_max <- ceiling(max_corr * 10) / 10  # Arrondir à la décimale supérieure
  
  # Graphique de comparaison des valeurs absolues
  p <- ggplot(
    comparison_data,
    aes(x = abs_corr_var1, y = abs_corr_var2, color = direction, shape = direction)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_point(size = 3, alpha = 0.8) +
    geom_text_repel(
      aes(label = !!sym(label_col)),
      size = 3,
      max.overlaps = 30,
      box.padding = 0.5
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_shape_manual(values = c(16, 17, 15, 18)) +
    scale_x_continuous(limits = c(axis_limit_min, axis_limit_max)) +
    scale_y_continuous(limits = c(axis_limit_min, axis_limit_max)) +
    labs(
      title = paste("Comparaison des corrélations", corr_method, "absolues entre", var1, "et", var2),
      subtitle = paste("Top", nrow(comparison_data), "métriques les plus corrélées avec", var1),
      x = paste("Corrélation absolue avec", var1),
      y = paste("Corrélation absolue avec", var2),
      color = "Direction",
      shape = "Direction"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  return(list(plot = p, data = comparison_data))
}

#' @title Créer un résumé des métriques les plus corrélées
#' @description Génère un résumé des métriques les plus fortement corrélées avec une variable cible
#' @param corr_data Tableau de données contenant les corrélations
#' @param target_var Variable cible
#' @param corr_method Méthode de corrélation ("Pearson" ou "Spearman")
#' @param top_n Nombre de métriques à inclure (en valeur absolue)
#' @return Dataframe avec le résumé des corrélations
summarize_top_correlations <- function(corr_data, target_var, corr_method = "Pearson", top_n = 10) {
  # Filtrer les corrélations pour la variable et méthode cible
  filtered_data <- corr_data %>%
    filter(target_variable == target_var, method == corr_method, !is.na(correlation), !is.na(p_value))
  
  # Sélectionner les top n en valeur absolue
  top_abs <- filtered_data %>%
    mutate(abs_corr = abs(correlation)) %>%
    arrange(desc(abs_corr)) %>%
    slice(1:min(top_n, n())) %>%
    mutate(
      significance = case_when(
        p_value >= 0.05 ~ "",
        p_value < 0.05 & p_value >= 0.01 ~ "*",
        p_value < 0.01 & p_value >= 0.001 ~ "**",
        p_value < 0.001 ~ "***"
      ),
      corr_formatted = paste0(round(correlation, 3), significance),
      direction = ifelse(correlation > 0, "positive", "negative"),
      rank = row_number()
    )
  
  return(top_abs)
}

#' @title Créer un corrplot pour les corrélations avec une variable cible
#' @description Génère un graphique corrplot pour visualiser les corrélations entre une variable cible 
#'              et les principales métriques, en utilisant la significativité des coefficients
#' @param data Jeu de données
#' @param target_var Variable cible
#' @param exclude_vars Variables à exclure de l'analyse
#' @param metric_labels Dictionnaire des étiquettes des métriques
#' @param path_target_var Répertoire pour la variable cible
#' @param method Méthode de corrélation ("Pearson" ou "Spearman")
#' @param top_n Nombre de métriques à inclure (NULL = toutes les métriques)
#' @param dpi Résolution en DPI
#' @return Invisiblement, la matrice de corrélation
create_target_corrplot <- function(data, target_var, exclude_vars, metric_labels, path_target_var, 
                                   method = "Pearson", top_n = NULL, dpi = 600) {
  # Vérifions si la variable cible existe dans les données
  if (!target_var %in% names(data)) {
    warning(paste("La variable cible", target_var, "n'existe pas dans les données."))
    return(NULL)
  }
  
  # Définir les variables pour l'analyse
  all_vars <- setdiff(names(data), exclude_vars)
  
  # Vérifier qu'il reste des variables après l'exclusion
  if (length(all_vars) == 0) {
    warning("Aucune variable disponible après exclusion.")
    return(NULL)
  }
  
  # S'assurer que la variable cible n'est pas exclue
  all_vars <- unique(c(target_var, all_vars))
  
  # Créer un dataframe avec uniquement les colonnes numériques
  numeric_cols <- sapply(data[, all_vars, drop = FALSE], is.numeric)
  if (sum(numeric_cols) == 0) {
    warning("Aucune variable numérique disponible pour l'analyse.")
    return(NULL)
  }
  
  # Filtrer pour ne garder que les colonnes numériques
  selected_vars <- all_vars[all_vars %in% names(data)[numeric_cols]]
  
  # Vérifier que nous avons suffisamment de variables pour créer un corrplot
  if (length(selected_vars) < 2) {
    warning("Pas assez de variables numériques pour créer un corrplot (minimum 2 requis).")
    return(NULL)
  }
  
  # Si top_n est spécifié, sélectionner les variables les plus corrélées au target_var
  if (!is.null(top_n) && top_n < length(selected_vars)) {
    # S'assurer que la variable cible est numérique
    if (!is.numeric(data[[target_var]])) {
      warning(paste("La variable cible", target_var, "n'est pas numérique."))
      return(NULL)
    }
    
    # Calculer les corrélations avec la variable cible
    corr_values <- numeric(length(selected_vars))
    names(corr_values) <- selected_vars
    
    for (i in seq_along(selected_vars)) {
      var <- selected_vars[i]
      if (var != target_var) {
        valid_data <- !is.na(data[[var]]) & !is.na(data[[target_var]])
        if (sum(valid_data) >= 3) {
          tryCatch({
            corr_values[i] <- cor(data[[var]][valid_data], data[[target_var]][valid_data], 
                                  method = tolower(method), use = "complete.obs")
          }, error = function(e) {
            corr_values[i] <- 0
            warning(paste("Erreur lors du calcul de la corrélation pour", var, ":", e$message))
          })
        } else {
          corr_values[i] <- 0
        }
      } else {
        corr_values[i] <- 1  # Corrélation parfaite avec soi-même
      }
    }
    
    # Trier par corrélation décroissante (en valeur absolue)
    sorted_indices <- order(abs(corr_values), decreasing = TRUE)
    sorted_vars <- selected_vars[sorted_indices]
    
    # S'assurer que la variable cible est toujours incluse et en première position
    sorted_vars <- unique(c(target_var, setdiff(sorted_vars, target_var)))
    selected_vars <- sorted_vars[1:min(top_n + 1, length(sorted_vars))]  # +1 pour inclure la variable cible
  } else {
    # Pour le cas où nous utilisons toutes les variables, placer la cible en premier
    selected_vars <- c(target_var, setdiff(selected_vars, target_var))
  }
  
  # Vérifier si certaines variables ont un écart-type nul et les retirer
  valid_vars <- sapply(selected_vars, function(var) {
    !all(is.na(data[[var]])) && sd(data[[var]], na.rm = TRUE) > 0
  })
  
  if (!all(valid_vars)) {
    invalid_vars <- selected_vars[!valid_vars]
    warning(paste("Variables avec écart-type nul ou toutes NA, exclues de l'analyse:", 
                  paste(invalid_vars, collapse = ", ")))
    selected_vars <- selected_vars[valid_vars]
  }
  
  # S'assurer que nous avons toujours assez de variables
  if (length(selected_vars) < 2) {
    warning("Pas assez de variables valides pour créer un corrplot après filtrage (minimum 2 requis).")
    return(NULL)
  }
  
  # Sélectionner les données pour les variables choisies
  selected_data <- data[, selected_vars, drop = FALSE]
  
  # Vérifier s'il y a des données manquantes et en informer l'utilisateur
  na_counts <- colSums(is.na(selected_data))
  if (any(na_counts > 0)) {
    message("Attention : certaines variables contiennent des valeurs manquantes:")
    for (var in names(na_counts)[na_counts > 0]) {
      message(paste(" -", var, ":", na_counts[var], "valeurs manquantes"))
    }
  }
  
  # Calculer la matrice de corrélation
  tryCatch({
    cors <- cor(selected_data, use = "pairwise.complete.obs", method = tolower(method))
  }, error = function(e) {
    warning(paste("Erreur lors du calcul de la matrice de corrélation:", e$message))
    return(NULL)
  })
  
  # Calculer la matrice des p-values
  p_matrix <- matrix(NA, nrow = ncol(selected_data), ncol = ncol(selected_data))
  colnames(p_matrix) <- colnames(selected_data)
  rownames(p_matrix) <- colnames(selected_data)
  
  for (i in 1:ncol(selected_data)) {
    for (j in 1:ncol(selected_data)) {
      if (i != j) {
        valid_data <- !is.na(selected_data[, i]) & !is.na(selected_data[, j])
        if (sum(valid_data) >= 3) {
          tryCatch({
            if (tolower(method) == "pearson") {
              test_result <- cor.test(selected_data[valid_data, i] %>% pull(), 
                                      selected_data[valid_data, j] %>% pull(), 
                                      method = "pearson")
            } else {
              test_result <- suppressWarnings(
                cor.test(selected_data[valid_data, i] %>% pull(), 
                         selected_data[valid_data, j] %>% pull(), 
                         method = "spearman", 
                         exact = FALSE)
              )
            }
            p_matrix[i, j] <- test_result$p.value
          }, error = function(e) {
            p_matrix[i, j] <- NA
            warning(paste("Erreur lors du calcul de la p-value pour", 
                          colnames(selected_data)[i], "et", 
                          colnames(selected_data)[j], ":", e$message))
          })
        }
      } else {
        p_matrix[i, j] <- 0  # Diagonale: p-value = 0
      }
    }
  }
  
  # Réorganiser pour avoir la corrélation avec la variable cible en ordre décroissant
  # Mais conserver la variable cible en première position
  if (length(selected_vars) > 2) {
    target_index <- which(colnames(cors) == target_var)
    other_indices <- setdiff(1:ncol(cors), target_index)
    
    # Trier les autres variables par leur corrélation avec la variable cible
    other_cors <- cors[target_index, other_indices]
    sorted_indices <- order(other_cors, decreasing = TRUE)
    
    # Nouvel ordre avec variable cible en premier, puis autres variables triées
    new_order <- c(target_index, other_indices[sorted_indices])
    cors <- cors[new_order, new_order]
    p_matrix <- p_matrix[new_order, new_order]
  }
  
  # Vérifier que toutes les lignes/colonnes des matrices ont des noms
  if (is.null(colnames(cors)) || is.null(rownames(cors))) {
    warning("La matrice de corrélation n'a pas de noms de lignes/colonnes.")
    colnames(cors) <- rownames(cors) <- colnames(selected_data)
  }
  
  if (is.null(colnames(p_matrix)) || is.null(rownames(p_matrix))) {
    warning("La matrice de p-value n'a pas de noms de lignes/colonnes.")
    colnames(p_matrix) <- rownames(p_matrix) <- colnames(selected_data)
  }
  
  # S'assurer que p_matrix et cors ont les mêmes dimensions et noms
  if (!identical(dim(p_matrix), dim(cors)) || 
      !identical(rownames(p_matrix), rownames(cors)) || 
      !identical(colnames(p_matrix), colnames(cors))) {
    warning("Les dimensions ou les noms des matrices de corrélation et de p-values ne correspondent pas.")
    return(NULL)
  }
  
  # Vérifier si metric_labels est fourni
  if (is.null(metric_labels)) {
    # Utiliser les noms des variables comme étiquettes si aucun dictionnaire n'est fourni
    metric_labels <- setNames(as.list(names(data)), names(data))
  }
  
  # Créer un vecteur de labels complets
  labels <- character(length(colnames(cors)))
  names(labels) <- colnames(cors)
  
  for (i in seq_along(colnames(cors))) {
    var_name <- colnames(cors)[i]
    if (!is.null(metric_labels) && var_name %in% names(metric_labels)) {
      labels[i] <- metric_labels[[var_name]]
    } else {
      labels[i] <- var_name
    }
  }
  
  rownames(cors) <- labels
  colnames(cors) <- labels
  rownames(p_matrix) <- labels
  colnames(p_matrix) <- labels
  
  # Fichier de sortie
  suffix <- if (is.null(top_n)) "all" else paste0("top", top_n)
  output_file <- file.path(path_target_var, paste0("corrplot_", method, "_", suffix, ".jpg"))
  
  # Calculer la taille d'étiquette appropriée en fonction du nombre de variables
  label_size <- max(0.8, min(1.5, 10 / length(colnames(cors))))
  
  # Déterminer la taille de l'image en fonction du nombre de variables
  img_size <- max(10, min(20, length(colnames(cors)) * 0.6))
  
  # Ouvrir le périphérique graphique avec gestion d'erreur
  tryCatch({
    jpeg(output_file, width = img_size, height = img_size, units = "in", res = dpi)
    
    # Créer le corrplot
    corrplot(
      cors, 
      p.mat = p_matrix, 
      method = 'circle', 
      type = 'lower', 
      insig = 'blank',
      addCoef.col = if(is.null(top_n) || length(colnames(cors)) > 25) NULL else 'black', # Ne pas afficher les valeurs pour les grands corrplots
      number.cex = max(0.5, min(0.9, 12 / length(colnames(cors)))),  # Ajuster la taille des nombres
      order = 'original',  # Utiliser l'ordre que nous avons défini
      diag = FALSE,
      tl.col = "black",
      tl.srt = 45,
      tl.cex = label_size,
      cl.cex = 1.0,  # Augmenter la taille de la légende
      cl.ratio = 0.2,  # Augmenter la taille de la barre de couleur
      col = colorRampPalette(c("#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", 
                               "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4"))(100),
      mar = c(1, 1, 3, 1),  # Augmenter les marges
      main = paste("Corrélations", method, "avec", target_var)
    )
    
    # Fermer le périphérique graphique
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()  # Fermer le périphérique en cas d'erreur
    warning(paste("Erreur lors de la création du corrplot:", e$message))
  })
  
  # Retourner invisiblement la matrice de corrélation
  invisible(cors)
}

#' @title Créer un scatterplot amélioré entre une métrique et une variable cible
#' @description Génère un graphique de dispersion entre une métrique et une variable cible,
#'              avec formes et couleurs différentes selon le pays et le site
#' @param data Jeu de données
#' @param metric_name Nom de la métrique
#' @param target_var Variable cible (ex: WD_BA, WD_BA)
#' @param metric_labels Dictionnaire des étiquettes des métriques
#' @param plots_info Dataframe contenant les colonnes plot_name, site, pays
#' @return Liste contenant l'objet ggplot et les résultats du modèle linéaire
plot_correlation_scatter_enhanced <- function(data, metric_name, target_var, metric_labels, plots_info) {
  # Utiliser l'étiquette complète si disponible
  metric_label <- NULL
  if (metric_name %in% names(metric_labels)) {
    metric_label <- metric_labels[[metric_name]]
  } else {
    metric_label <- metric_name
  }
  
  # Vérifier que data contient la colonne plot_name pour joindre avec plots_info
  if (!"plot_name" %in% names(data)) {
    stop("La colonne 'plot_name' est requise dans le jeu de données pour fusionner avec les informations de sites/pays.")
  }
  
  # Données valides uniquement (sans les NA pour la métrique et la variable cible)
  valid_data <- data %>% 
    filter(!is.na(.data[[metric_name]]), !is.na(.data[[target_var]]))
  
  # Joindre avec les informations de sites et pays
  valid_data <- valid_data %>%
    left_join(plots_info %>% select(plot_name, site, pays), by = "plot_name")
  
  # Vérifier que la jointure a fonctionné
  if (!"site" %in% names(valid_data) || !"pays" %in% names(valid_data)) {
    warning("Impossible de joindre les informations de sites/pays avec les données.")
    valid_data$site <- "Inconnu"
    valid_data$pays <- "Inconnu"
  }
  
  # Calculer le modèle linéaire
  formula_str <- paste(target_var, "~", metric_name)
  lm_model <- lm(as.formula(formula_str), data = valid_data)
  lm_summary <- summary(lm_model)
  
  # Extraire les coefficients et statistiques
  intercept <- round(coef(lm_model)[1], 4)
  slope <- round(coef(lm_model)[2], 4)
  r_squared <- round(lm_summary$r.squared, 4)
  adj_r_squared <- round(lm_summary$adj.r.squared, 4)
  f_stat <- round(lm_summary$fstatistic[1], 2)
  p_value <- round(lm_summary$coefficients[2, 4], 4)
  correlation <- round(cor(valid_data[[metric_name]], valid_data[[target_var]], use = "complete.obs"), 4)
  
  # Trouver les min/max pour le scaling des axes
  x_range <- range(valid_data[[metric_name]], na.rm = TRUE)
  y_range <- range(valid_data[[target_var]], na.rm = TRUE)
  
  # Ajouter une marge de 5% aux limites
  x_margin <- diff(x_range) * 0.05
  y_margin <- diff(y_range) * 0.05
  x_limits <- c(x_range[1] - x_margin, x_range[2] + x_margin)
  y_limits <- c(y_range[1] - y_margin, y_range[2] + y_margin)
  
  # Équation de la régression
  eq_text <- paste0("y = ", slope, "x + ", intercept)
  stat_text <- paste0("R² = ", r_squared, ", p = ", p_value)
  
  # Définir le nombre de formes différentes disponibles dans ggplot2
  shapes_list <- c(16, 17, 15, 18, 8, 3, 4, 5, 6, 0, 1, 2)
  
  # Créer une palette de couleurs avec assez de couleurs distinctes
  # Utiliser une palette personnalisée pour les sites
  site_colors <- scales::hue_pal()(length(unique(valid_data$site)))
  names(site_colors) <- unique(valid_data$site)
  
  # Créer le graphique
  p <- ggplot(valid_data, aes(x = .data[[metric_name]], y = .data[[target_var]])) +
    geom_smooth(method = "lm", se = TRUE, color = "grey60", alpha = 0.2) +
    geom_point(aes(shape = pays, color = site), size = 3, alpha = 0.9) +
    scale_shape_manual(values = shapes_list[1:length(unique(valid_data$pays))]) +
    scale_color_viridis_d(option = color_palette_discrete) +
    labs(
      title = paste("Relation entre", metric_label, "et", target_var),
      subtitle = paste(eq_text, "\n", stat_text),
      x = metric_label,
      y = target_var,
      shape = "Pays",
      color = "Site"
    ) +
    xlim(x_limits) +
    ylim(y_limits) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, color = "darkblue", hjust = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 8)
    )
  
  # Créer une liste avec le graphique et les statistiques
  result <- list(
    plot = p,
    stats = data.frame(
      metric = metric_name,
      metric_label = metric_label,
      target_var = target_var,
      intercept = intercept,
      slope = slope,
      r_squared = r_squared,
      adj_r_squared = adj_r_squared,
      f_statistic = f_stat,
      p_value = p_value,
      correlation = correlation,
      n_observations = nrow(valid_data)
    )
  )
  
  return(result)
}

#' @title Générer des visualisations batch améliorées pour les métriques les plus corrélées
#' @description Génère des graphiques de dispersion pour les métriques les plus corrélées avec formes et couleurs selon pays/sites
#' @param data Jeu de données
#' @param corr_data Tableau de données contenant les corrélations
#' @param target_var Variable cible (ex: WD_BA, WD_BA)
#' @param metric_labels Dictionnaire des étiquettes des métriques
#' @param plots_info Dataframe contenant les colonnes plot_name, site, pays
#' @param path_target_var Répertoire pour la variable cible
#' @param corr_method Méthode de corrélation ("Pearson" ou "Spearman")
#' @param top_n Nombre de métriques à inclure (les plus corrélées en valeur absolue)
#' @param dpi Résolution en DPI
#' @return Dataframe avec les statistiques des modèles linéaires
generate_enhanced_scatterplots <- function(data, corr_data, target_var, metric_labels, plots_info,
                                           path_target_var, corr_method = "Pearson", 
                                           top_n = 10, dpi = 600) {
  # Filtrer les métriques les plus corrélées en valeur absolue
  filtered_corr <- corr_data %>%
    filter(target_variable == target_var, method == corr_method, !is.na(correlation)) %>%
    mutate(abs_corr = abs(correlation)) %>%
    arrange(desc(abs_corr)) %>%
    slice(1:min(top_n, n()))
  
  top_metrics <- filtered_corr$metric
  
  # Créer un objet combinant tous les graphiques
  scatter_plots <- list()
  all_stats <- data.frame()
  
  # Préparer une légende commune
  # Nous allons créer un graphique spécial juste pour récupérer la légende
  
  # Générer chaque graphique
  for (i in seq_along(top_metrics)) {
    metric <- top_metrics[i]
    result <- plot_correlation_scatter_enhanced(data, metric, target_var, metric_labels, plots_info)
    scatter_plots[[i]] <- result$plot
    all_stats <- bind_rows(all_stats, result$stats)
    
    # Sauvegarder individuellement dans le dossier de la variable cible
    ggsave(
      file.path(path_target_var, paste0("scatter_", metric, ".jpg")),
      result$plot,
      width = 9,  # Un peu plus large pour accommoder la légende
      height = 7,
      dpi = dpi
    )
  }
  
  # Créer un graphique combiné avec tous les graphiques
  # mais en extrayant la légende pour la mettre une seule fois
  # Récupérer une légende commune du premier graphique
  if (length(scatter_plots) > 0) {
    # Fonction pour extraire la légende d'un graphique
    get_legend <- function(p) {
      tmp <- ggplot_gtable(ggplot_build(p))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      if (length(leg) > 0) {
        legend <- tmp$grobs[[leg]]
        return(legend)
      }
      return(NULL)
    }
    
    # Extraire la légende du premier graphique
    legend <- get_legend(scatter_plots[[1]])
    
    # Supprimer les légendes de tous les graphiques
    scatter_plots_no_legend <- lapply(scatter_plots, function(p) p + theme(legend.position = "none"))
    
    # Créer une mise en page pour les graphiques sans légende
    n_rows <- ceiling(length(scatter_plots_no_legend) / 2)
    combined_plot <- wrap_plots(scatter_plots_no_legend, ncol = 2, nrow = n_rows)
    
    # Sauvegarder le graphique combiné
    ggsave(
      file.path(path_target_var, paste0("combined_scatter_top", top_n, "_", corr_method, ".jpg")),
      combined_plot,
      width = 16,
      height = n_rows * 5,
      dpi = dpi,
      limitsize = FALSE
    )
    
    # Créer un deuxième graphique combiné, mais cette fois avec la légende
    if (!is.null(legend)) {
      # Combiner les graphiques sans légende et la légende commune
      final_plot <- wrap_plots(
        combined_plot, 
        ggdraw() + draw_grob(legend),
        ncol = 2,
        widths = c(0.8, 0.2)
      )
      
      # Sauvegarder le graphique combiné avec légende commune
      ggsave(
        file.path(path_target_var, paste0("combined_scatter_with_legend_top", top_n, "_", corr_method, ".jpg")),
        final_plot,
        width = 18,  # Un peu plus large pour accommoder la légende
        height = n_rows * 5,
        dpi = dpi,
        limitsize = FALSE
      )
    }
  }
  
  # Sauvegarder les statistiques dans le dossier de la variable cible
  write_csv2(
    all_stats,
    file.path(path_target_var, paste0("lm_stats_top", top_n, "_", corr_method, ".csv"))
  )
  
  return(all_stats)
}
