# =================================================================
# FONCTIONS DE VISUALISATION DES PROFILS DE TROUÉES
# Objectif : Fournir des fonctions pour créer des visualisations avancées
#            des profils de trouées en relation avec la densité du bois
# =================================================================

#### FONCTIONS D'ANALYSE DES PROFILS PAR CATÉGORIE de WD ####

#' Crée un graphique des courbes de proportion de trouées par catégorie de WD
#'
#' @param data Données combinées avec colonnes height_aboveground, proportion, wd_category
#' @param wd_palette Palette de couleurs pour les catégories de WD
#' @param title Titre du graphique
#' @param smooth Booléen indiquant si on veut lisser les courbes moyennes
#' @return Un objet ggplot
plot_gap_curves_by_wd <- function(data, wd_palette, 
                                  title = "Courbes de proportion de trouées par catégorie de densité de bois",
                                  smooth = TRUE) {
  
  # Vérifier que les colonnes nécessaires sont présentes
  required_cols <- c("height_aboveground", "proportion", "wd_category", "plot_name")
  if (!all(required_cols %in% colnames(data))) {
    stop("Les données doivent contenir les colonnes: ", paste(required_cols, collapse = ", "))
  }
  
  # Créer le graphique
  p <- ggplot(data, 
              aes(x = height_aboveground, y = proportion, 
                  color = wd_category, group = interaction(plot_name, wd_category))) +
    geom_line(alpha = 0.3) +
    scale_color_manual(name = "Catégorie WD", values = wd_palette) +
    labs(
      title = title,
      x = "Hauteur au-dessus du sol (m)",
      y = "Proportion de trouées"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold")
    )
  
  # Ajouter les courbes moyennes
  if (smooth) {
    p <- p + stat_summary(aes(group = wd_category), fun = mean, geom = "line", linewidth = 1.5)
  } else {
    # Calculer les moyennes manuellement pour plus de contrôle
    avg_data <- data %>%
      group_by(wd_category, height_aboveground) %>%
      summarise(prop_mean = mean(proportion, na.rm = TRUE), .groups = "drop")
    
    p <- p + geom_line(data = avg_data, 
                       aes(x = height_aboveground, y = prop_mean, group = wd_category),
                       linewidth = 1.5)
  }
  
  return(p)
}

#' Crée une heatmap des profils de trouées ordonnés par densité du bois
#'
#' @param data Dataframe contenant les données de profils
#' @param plot_names Noms des plots à inclure
#' @param color_palette Palette de couleurs pour la heatmap (défaut: "turbo")
#' @param height_max Hauteur maximale à représenter
#' @param height_step Pas d'incrémentation pour la grille de hauteurs
#' @param title Titre du graphique
#' @param n_wd_ticks Nombre de valeurs pivots de WD à afficher (défaut: 5)
#' @return Un objet ggplot
create_profiles_heatmap <- function(data, plot_names, color_palette = "turbo",
                                    height_max = 60, height_step = 1,
                                    title = "Heatmap des profils de trouées ordonnés par densité du bois croissante",
                                    n_wd_ticks = 5) {
  
  # Créer une grille régulière de hauteurs
  height_grid_heatmap <- seq(0, height_max, by = height_step)
  
  # Créer un tableau pour stocker les profils interpolés
  plot_profiles <- tibble(
    height = rep(height_grid_heatmap, times = length(plot_names)),
    plot_name = rep(plot_names, each = length(height_grid_heatmap)),
    proportion = NA_real_
  )
  
  # Pour chaque plot, interpoler le profil
  for (plot_i in plot_names) {
    plot_data <- data %>%
      filter(plot_name == plot_i) %>%
      arrange(height_aboveground)
    
    # Vérifier suffisamment de points
    if (nrow(plot_data) >= 3) {
      # Interpoler avec une spline
      spline_fit <- smooth.spline(plot_data$height_aboveground, plot_data$proportion)
      interp_values <- predict(spline_fit, height_grid_heatmap)$y
      
      # Limiter entre 0 et 1
      interp_values <- pmin(pmax(interp_values, 0), 1)
      
      # Ajouter au dataframe
      plot_profiles$proportion[plot_profiles$plot_name == plot_i] <- interp_values
    }
  }
  
  # Joindre avec les valeurs de WD
  plot_profiles <- plot_profiles %>%
    left_join(data %>% 
                select(plot_name, WD_BA) %>% 
                distinct(), 
              by = "plot_name") %>%
    filter(!is.na(proportion), !is.na(WD_BA)) %>%
    arrange(WD_BA)
  
  # Liste ordonnée des plots par WD
  ordered_plots_df <- plot_profiles %>%
    group_by(plot_name) %>%
    summarise(WD_BA = first(WD_BA)) %>%
    arrange(WD_BA)
  
  ordered_plots <- ordered_plots_df %>% pull(plot_name)
  
  # Identifier les plots pivots pour afficher les valeurs WD
  n_plots <- length(ordered_plots)
  if (n_wd_ticks > n_plots) n_wd_ticks <- n_plots
  
  pivot_indices <- round(seq(1, n_plots, length.out = n_wd_ticks))
  pivot_plots <- ordered_plots[pivot_indices]
  pivot_wds <- ordered_plots_df$WD_BA[pivot_indices]
  
  # Créer un vecteur d'étiquettes qui affiche la valeur WD pour les plots pivots
  # et des chaînes vides pour les autres plots
  x_labels <- rep("", n_plots)
  names(x_labels) <- ordered_plots
  x_labels[pivot_plots] <- sprintf("WD: %.2f", pivot_wds)
  
  # Créer la heatmap
  p <- ggplot(plot_profiles, aes(x = reorder(plot_name, WD_BA), y = height, fill = proportion)) +
    geom_tile() +
    scale_fill_viridis_c(option = color_palette, name = "Proportion\nde trouées") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(labels = x_labels) +  # Utiliser nos étiquettes personnalisées
    labs(
      title = title,
      x = "Plots (WD croissant →)",
      y = "Hauteur au-dessus du sol (m)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # Rotation des étiquettes
      panel.grid = element_blank(),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  return(p)
}

#' Crée un graphique des profils types avec zones d'incertitude optionnelles
#'
#' @param data Données combinées avec colonnes height_aboveground, proportion
#' @param wd_palette Palette de couleurs pour les catégories de WD
#' @param category_column Nom de la colonne contenant les catégories à utiliser
#' @param show_deviation Booléen pour afficher ou non les zones d'incertitude
#' @param title Titre du graphique
#' @param height_max Hauteur maximale pour la grille
#' @param height_step Pas de la grille de hauteur
#' @return Un objet ggplot
create_typical_profiles <- function(data, wd_palette, 
                                    category_column = "wd_category_quartile",
                                    show_deviation = TRUE,
                                    title = "Profils types de proportion de trouées par catégorie de densité du bois",
                                    height_max = 60, height_step = 0.5) {
  
  # Vérifier que la colonne de catégorie existe
  if (!(category_column %in% colnames(data))) {
    stop("La colonne de catégorie spécifiée n'existe pas dans les données.")
  }
  
  # Créer une grille de hauteurs commune
  height_grid <- seq(0, height_max, by = height_step)
  
  # Dataframe pour stocker les profils types
  prototypes <- data.frame(height = height_grid)
  
  categories <- data %>%
    group_by(!!sym(category_column)) %>%
    summarise(mean_wd = mean(WD_BA, na.rm=TRUE)) %>%
    arrange(mean_wd) %>%
    pull(!!sym(category_column))
  
  # Ne garder que les catégories qui existent réellement dans les données
  categories <- categories[categories %in% unique(data[[category_column]])]
  
  for (cat in categories) {
    # Filtrer les plots de cette catégorie
    cat_plots <- data %>%
      filter(!!sym(category_column) == cat) %>%
      pull(plot_name) %>%
      unique()
    
    # Préparer un dataframe pour stocker les valeurs interpolées
    all_interp <- data.frame(height = height_grid)
    
    # Pour chaque plot, interpoler sur la grille commune
    for (plot_i in cat_plots) {
      plot_data <- data %>%
        filter(plot_name == plot_i) %>%
        arrange(height_aboveground)
      
      # Vérifier suffisamment de points
      if (nrow(plot_data) >= 3) {
        # Interpoler avec une spline
        spline_fit <- smooth.spline(plot_data$height_aboveground, plot_data$proportion)
        interp_values <- predict(spline_fit, height_grid)$y
        
        # Limiter entre 0 et 1
        interp_values <- pmin(pmax(interp_values, 0), 1)
        
        # Ajouter au dataframe
        all_interp[[plot_i]] <- interp_values
      }
    }
    
    # Calculer la moyenne et écart-type pour cette catégorie
    if (ncol(all_interp) > 1) {
      all_interp$mean <- rowMeans(all_interp[, -1], na.rm = TRUE)
      all_interp$sd <- apply(all_interp[, -1], 1, sd, na.rm = TRUE)
      
      # Ajouter au dataframe des prototypes
      prototypes[[paste0("mean_", cat)]] <- all_interp$mean
      prototypes[[paste0("sd_", cat)]] <- all_interp$sd
    }
  }
  
  # Vérifier qu'on a au moins une colonne de moyenne
  mean_cols <- grep("^mean_", names(prototypes), value = TRUE)
  if (length(mean_cols) == 0) {
    stop("Aucune catégorie n'a suffisamment de données pour créer un profil type")
  }
  
  # Convertir en format long pour ggplot
  prototypes_long <- prototypes %>%
    pivot_longer(
      cols = all_of(mean_cols),
      names_to = "category",
      values_to = "mean_value"
    ) %>%
    mutate(
      category = gsub("mean_", "", category),
      sd_col = paste0("sd_", category)
    ) %>%
    rowwise() %>%
    mutate(
      sd_value = prototypes[[sd_col]][which(prototypes$height == height)]
    ) %>%
    ungroup()
  
  # S'assurer que la catégorie est un facteur avec le bon ordre
  prototypes_long$category <- factor(prototypes_long$category, levels = categories)
  
  # Créer la visualisation de base
  p <- ggplot(prototypes_long, aes(x = height, y = mean_value, color = category)) +
    geom_line(size = 1.5) +
    scale_color_manual(values = wd_palette) +
    labs(
      title = title,
      x = "Hauteur au-dessus du sol (m)",
      y = "Proportion de trouées",
      color = "Catégorie WD"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold")
    )
  
  # Ajouter les zones d'incertitude si demandé
  if (show_deviation) {
    p <- p + 
      geom_ribbon(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value, fill = category), alpha = 0.2) +
      scale_fill_manual(values = wd_palette) +
      labs(fill = "Catégorie WD")
  }
  
  return(p)
}

#' Crée un graphique des signatures différentielles avec zones d'incertitude
#'
#' @param data Données combinées avec colonnes height_aboveground, proportion, wd_category
#' @param wd_palette Palette de couleurs pour les catégories de WD
#' @param show_deviation Booléen pour afficher ou non les zones d'incertitude
#' @param title Titre du graphique
#' @param confidence_interval Type d'intervalle de confiance (sd, se)
#' @param smooth_factor Facteur de lissage pour le calcul des dérivées
#' @return Un objet ggplot
create_differential_signatures <- function(data, wd_palette,
                                           title = "Signatures différentielles des profils par catégorie de densité du bois",
                                           confidence_interval = "sd", smooth_factor = 0.4,
                                           show_deviation = TRUE) {
  
  # Vérifier les colonnes nécessaires
  required_cols <- c("height_aboveground", "proportion", "wd_category", "plot_name")
  if (!all(required_cols %in% colnames(data))) {
    stop("Les données doivent contenir les colonnes: ", paste(required_cols, collapse = ", "))
  }
  
  # S'assurer que wd_category est un facteur et, si possible, ordonné selon la palette fournie
  if (!is.factor(data$wd_category)) {
    if (!is.null(wd_palette) && !is.null(names(wd_palette))) {
      data$wd_category <- factor(data$wd_category, levels = names(wd_palette))
    } else if ("WD_BA" %in% colnames(data)) {
      categories <- data %>%
        dplyr::group_by(wd_category) %>%
        dplyr::summarise(mean_wd = mean(WD_BA, na.rm = TRUE)) %>%
        dplyr::arrange(mean_wd) %>%
        dplyr::pull(wd_category)
      data$wd_category <- factor(data$wd_category, levels = categories)
    } else {
      data$wd_category <- factor(data$wd_category)
    }
  }
  
  # Créer une grille de hauteurs commune
  height_grid <- seq(min(data$height_aboveground), max(data$height_aboveground), length.out = 200)
  
  # Fonction pour calculer la dérivée première
  calculate_first_derivative <- function(heights, props, smooth_span = smooth_factor) {
    # S'assurer que les données sont triées par hauteur
    ordered <- order(heights)
    heights <- heights[ordered]
    props <- props[ordered]
    
    # Ajuster une spline pour le lissage
    spline_fit <- smooth.spline(heights, props, spar = smooth_span)
    
    # Prédire sur la grille régulière
    pred_heights <- height_grid
    pred_props <- predict(spline_fit, pred_heights)$y
    
    # Calculer la dérivée première
    d1 <- diff(pred_props) / diff(pred_heights)
    d1_heights <- (pred_heights[-1] + pred_heights[-length(pred_heights)]) / 2
    
    return(list(heights = d1_heights, derivatives = d1))
  }
  
  # Dataframe pour stocker les signatures différentielles
  derivatives_data <- data.frame()
  
  # Pour chaque catégorie de WD
  for (cat in levels(data$wd_category)) {
    # Filtrer les plots de cette catégorie
    cat_data <- data %>% 
      dplyr::filter(wd_category == cat)
    
    # Pour chaque plot, calculer la dérivée
    plot_derivatives <- list()
    for (plot_i in unique(cat_data$plot_name)) {
      plot_data <- cat_data %>% 
        dplyr::filter(plot_name == plot_i) %>%
        dplyr::arrange(height_aboveground)
      
      # Vérifier qu'il y a suffisamment de points
      if (nrow(plot_data) >= 5) {
        # Calculer la dérivée
        deriv_result <- calculate_first_derivative(
          plot_data$height_aboveground, 
          plot_data$proportion, 
          smooth_span = smooth_factor
        )
        
        # Ajouter au dataframe
        plot_derivatives[[plot_i]] <- data.frame(
          height = deriv_result$heights,
          derivative = deriv_result$derivatives,
          plot_name = plot_i,
          wd_category = cat
        )
      }
    }
    
    # Combiner tous les plots de cette catégorie
    if (length(plot_derivatives) > 0) {
      cat_derivatives <- dplyr::bind_rows(plot_derivatives)
      derivatives_data <- dplyr::bind_rows(derivatives_data, cat_derivatives)
    }
  }
  
  # Calculer les moyennes et intervalles de confiance par catégorie et hauteur
  if (nrow(derivatives_data) > 0) {
    summary_derivatives <- derivatives_data %>%
      dplyr::group_by(wd_category, height) %>%
      dplyr::summarise(
        mean_derivative = mean(derivative, na.rm = TRUE),
        sd_derivative = sd(derivative, na.rm = TRUE),
        se_derivative = sd_derivative / sqrt(dplyr::n()),
        n = dplyr::n(),
        .groups = "drop"
      )
    
    # Choisir l'intervalle de confiance
    if (confidence_interval == "sd") {
      summary_derivatives$ci <- summary_derivatives$sd_derivative
    } else if (confidence_interval == "se") {
      summary_derivatives$ci <- summary_derivatives$se_derivative
    } else {
      stop("confidence_interval doit être 'sd' ou 'se'")
    }
    
    # Créer la visualisation
    p <- ggplot2::ggplot(summary_derivatives, ggplot2::aes(x = height, y = mean_derivative, color = wd_category, fill = wd_category))
    
    # Ajouter la zone d'incertitude seulement si show_deviation est TRUE
    if (show_deviation) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = mean_derivative - ci, ymax = mean_derivative + ci), alpha = 0.2)
    }
    
    p <- p + 
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::scale_color_manual(values = wd_palette) +
      ggplot2::scale_fill_manual(values = wd_palette) +
      ggplot2::labs(
        title = title,
        x = "Hauteur au-dessus du sol (m)",
        y = "Dérivée première (dP/dh)",
        color = "Catégorie WD",
        fill = "Catégorie WD"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(size = 14, face = "bold")
      )
    
    return(p)
  } else {
    stop("Pas assez de données pour calculer les signatures différentielles")
  }
}

#### FONCTIONS D'ANALYSE AVANCÉE DES PROFILS ####

#' Décompose une courbe de proportion de trouées en segments fonctionnels
#'
#' @param height Vecteur des hauteurs
#' @param prop Vecteur des proportions de trouées
#' @param plot_name Nom du plot pour le titre
#' @param wd_value Valeur WD pour le titre
#' @param threshold_low Seuil inférieur pour la zone basale (défaut: 0.25)
#' @param threshold_high Seuil supérieur pour la zone de transition (défaut: 0.75)
#' @return Un ggplot avec la courbe décomposée en segments
plot_decompose_curve_segments <- function(height, prop, plot_name, wd_value,
                                          threshold_low = 0.25, threshold_high = 0.75) {
  # S'assurer que les données sont triées par hauteur
  data <- tibble(height = height, prop = prop) %>%
    arrange(height)
  
  # Ajuster une spline pour l'interpolation
  spline_fit <- smooth.spline(data$height, data$prop)
  
  # Prédire sur une grille régulière pour faciliter l'analyse
  height_grid <- seq(min(data$height), max(data$height), length.out = 200)
  prop_grid <- predict(spline_fit, height_grid)$y
  
  # Trouver les hauteurs correspondant aux seuils
  h_low <- approx(prop_grid, height_grid, xout = threshold_low)$y
  h_high <- approx(prop_grid, height_grid, xout = threshold_high)$y
  
  # Diviser la courbe en segments
  segment_basal <- tibble(
    height = height_grid[height_grid <= h_low],
    prop = prop_grid[height_grid <= h_low],
    segment = "basal"
  )
  
  segment_transition <- tibble(
    height = height_grid[height_grid > h_low & height_grid <= h_high],
    prop = prop_grid[height_grid > h_low & height_grid <= h_high],
    segment = "transition"
  )
  
  segment_upper <- tibble(
    height = height_grid[height_grid > h_high],
    prop = prop_grid[height_grid > h_high],
    segment = "upper"
  )
  
  # Combiner les segments
  segments <- bind_rows(segment_basal, segment_transition, segment_upper)
  
  # Utiliser une palette viridis pour plus de clarté
  segment_colors <- viridis(3, option = "viridis")
  names(segment_colors) <- c("basal", "transition", "upper")
  
  # Créer le graphique
  p <- ggplot() +
    # Données originales
    geom_line(data = data.frame(height = height, prop = prop),
              aes(x = height, y = prop), color = "gray50", size = 0.7, alpha = 0.7) +
    # Segments colorés
    geom_line(data = segments, aes(x = height, y = prop, color = segment), size = 1.2) +
    scale_color_manual(values = segment_colors,
                       labels = c("Basale", "Transition", "Supérieure")) +
    labs(
      title = "Décomposition fonctionnelle",
      subtitle = paste0(plot_name, " (WD = ", round(wd_value, 3), ")"),
      x = "Hauteur au-dessus du sol (m)",
      y = "Proportion de trouées",
      color = "Zone"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

#' Calcule et affiche les dérivées d'une courbe de proportion de trouées
#'
#' @param height Vecteur des hauteurs
#' @param prop Vecteur des proportions de trouées
#' @param plot_name Nom du plot pour le titre
#' @param wd_value Valeur de WD pour le titre
#' @param smooth_span Paramètre de lissage pour la spline (défaut: 0.4)
#' @return Liste avec deux graphiques (courbe+dérivée et dérivées)
plot_derivatives <- function(height, prop, plot_name, wd_value, smooth_span = 0.4) {
  # S'assurer que les données sont triées par hauteur
  data <- tibble(height = height, prop = prop) %>%
    arrange(height)
  
  # Ajuster une spline pour le lissage
  spline_fit <- smooth.spline(data$height, data$prop, spar = smooth_span)
  
  # Prédire sur une grille régulière
  height_grid <- seq(min(data$height), max(data$height), length.out = 200)
  prop_smooth <- predict(spline_fit, height_grid)$y
  
  # Calculer la dérivée première
  d1 <- diff(prop_smooth) / diff(height_grid)
  height_d1 <- (height_grid[-1] + height_grid[-length(height_grid)]) / 2
  
  # Calculer la dérivée seconde
  d2 <- diff(d1) / diff(height_d1)
  height_d2 <- (height_d1[-1] + height_d1[-length(height_d1)]) / 2
  
  # Trouver les points caractéristiques
  max_d1_idx <- which.max(d1)
  max_d2_idx <- which.max(d2)
  min_d2_idx <- which.min(d2)
  
  height_max_d1 <- height_d1[max_d1_idx]
  height_max_d2 <- height_d2[max_d2_idx]
  height_min_d2 <- height_d2[min_d2_idx]
  
  # Créer un dataframe pour les résultats
  derivatives <- tibble(
    height = height_grid,
    prop = prop_smooth
  )
  
  # Ajouter les dérivées avec interpolation pour aligner sur la grille originale
  derivatives$d1 <- c(NA, approx(height_d1, d1, xout = height_grid[-1])$y)
  derivatives$d2 <- c(NA, NA, approx(height_d2, d2, xout = height_grid[-c(1,2)])$y)
  
  # Échelles ajustées pour mieux voir les détails
  scale_d1 <- 5   # Échelle plus petite pour la dérivée première
  scale_d2 <- 50  # Échelle plus petite pour la dérivée seconde
  
  # Graphique 1: Courbe originale avec dérivée première
  p1 <- ggplot() +
    geom_line(data = derivatives, aes(x = height, y = prop), color = "black", size = 0.8) +
    geom_line(data = derivatives, aes(x = height, y = d1 * scale_d1), color = "red", size = 0.8) +
    geom_vline(xintercept = height_max_d1, linetype = "dashed", color = "red") +
    scale_y_continuous(
      name = "Proportion de trouées",
      sec.axis = sec_axis(~./scale_d1, name = "Dérivée première (dP/dh)")
    ) +
    labs(
      title = "Analyse différentielle : dérivée première",
      x = "Hauteur (m)"
    ) +
    theme_minimal() +
    theme(
      axis.title.y.right = element_text(color = "red"),
      axis.text.y.right = element_text(color = "red")
    ) +
    annotate("text", x = height_max_d1 + 2, y = 0.8, label = "Point d'inflexion", color = "red")
  
  # Graphique 2: Dérivées première et seconde
  p2 <- ggplot() +
    geom_line(data = derivatives, aes(x = height, y = d1 * scale_d1), color = "red", size = 0.8) +
    geom_line(data = derivatives, aes(x = height, y = d2 * scale_d2), color = "blue", size = 0.8) +
    geom_vline(xintercept = height_max_d2, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = height_min_d2, linetype = "dotted", color = "blue") +
    scale_y_continuous(
      name = "Dérivée première (dP/dh)",
      sec.axis = sec_axis(~./scale_d2*scale_d1, name = "Dérivée seconde (d²P/dh²)")
    ) +
    labs(
      title = "Analyse différentielle : dérivée seconde",
      x = "Hauteur (m)"
    ) +
    theme_minimal() +
    theme(
      axis.title.y.left = element_text(color = "red"),
      axis.text.y.left = element_text(color = "red"),
      axis.title.y.right = element_text(color = "blue"),
      axis.text.y.right = element_text(color = "blue")
    )
  
  return(list(p1 = p1, p2 = p2))
}

#' Réalise une analyse multi-échelle d'une courbe de proportion de trouées
#'
#' @param height Vecteur des hauteurs
#' @param prop Vecteur des proportions de trouées
#' @param plot_name Nom du plot pour le titre
#' @param wd_value Valeur de WD pour le titre
#' @param scales Échelles à analyser en mètres
#' @return Liste avec deux graphiques (courbes lissées et résidus)
plot_multiscale_analysis <- function(height, prop, plot_name, wd_value, scales = c(5, 10, 15, 20)) {
  # S'assurer que les données sont triées par hauteur
  data <- tibble(height = height, prop = prop) %>%
    arrange(height)
  
  # Ajuster une spline pour avoir des points équidistants
  spline_fit <- smooth.spline(data$height, data$prop)
  
  # Prédire sur une grille régulière
  height_grid <- seq(min(data$height), max(data$height), length.out = 512)
  prop_grid <- predict(spline_fit, height_grid)$y
  
  # Utiliser une fonction simple de moyennage pour différentes échelles
  scale_results <- list()
  
  for (scale in scales) {
    # Convertir l'échelle en nombre de points
    window_size <- max(3, round(scale / (diff(range(height_grid)) / length(height_grid))))
    if (window_size %% 2 == 0) window_size <- window_size + 1  # Assurer un nombre impair
    
    # Appliquer un filtre de moyennage
    prop_smooth <- zoo::rollmean(prop_grid, window_size, fill = "extend")
    
    # Calculer les différences entre la courbe originale et la courbe lissée
    residuals <- prop_grid - prop_smooth
    
    # Stocker les résultats
    scale_results[[as.character(scale)]] <- tibble(
      height = height_grid,
      prop_original = prop_grid,
      prop_smooth = prop_smooth,
      residuals = residuals,
      scale = paste0(scale, "m")
    )
  }
  
  # Combiner tous les résultats
  multi_df <- bind_rows(scale_results)
  
  # Définir une palette de couleurs turbo
  color_scale <- scales::viridis_pal(option = "viridis")(length(scales))
  names(color_scale) <- paste0(scales, "m")
  
  # Courbes originales et lissées
  p1 <- ggplot() +
    # D'abord les courbes lissées
    geom_line(data = multi_df, aes(x = height, y = prop_smooth, color = scale), size = 0.8) +
    # Puis la courbe originale en noir pointillé par-dessus pour la rendre plus visible
    geom_line(data = multi_df %>% filter(scale == paste0(scales[1], "m")), 
              aes(x = height, y = prop_original), color = "red", size = 1, linetype = "dotted") +
    scale_color_manual(values = color_scale) +
    labs(
      title = "Analyse multi-échelle : lissage",
      subtitle = "Filtrage de la courbe originale (pointillés) à différentes échelles spatiales",
      x = "Hauteur (m)",
      y = "Proportion",
      color = "Échelle"
    ) +
    theme_minimal()
  
  # Résidus aux différentes échelles
  p2 <- ggplot(multi_df, aes(x = height, y = residuals, color = scale)) +
    geom_line(size = 0.8) +
    scale_color_manual(values = color_scale) +
    labs(
      title = "Analyse multi-échelle : résidus",
      subtitle = "Variations fines détectées à chaque échelle d'analyse",
      x = "Hauteur (m)",
      y = "Résidus",
      color = "Échelle"
    ) +
    theme_minimal()
  
  return(list(p1 = p1, p2 = p2))
}
