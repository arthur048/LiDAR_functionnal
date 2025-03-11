# =================================================================
# SCRIPT DE VISUALISATION DES MÉTRIQUES DE TROUÉES
# Objectif : Créer des visualisations de l'évolution de la proportion
#            de trouées en fonction de la hauteur pour différents plots
# =================================================================

#### 1️⃣ INITIALISATION ####
# Nettoyage de l'environnement et de la mémoire
rm(list = ls())
gc()

# Définition des packages nécessaires
pkgs <- c("terra", "tidyverse", "sf", "rio", "foreach", "doParallel", "ggpubr", "viridis")

# Installation et chargement automatique des packages manquants
to_install <- !pkgs %in% installed.packages()
if(any(to_install)) {install.packages(pkgs[to_install])}
inst <- lapply(pkgs, library, character.only = TRUE)

# Définition du répertoire de travail
path0 <- "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/3_Gaps/"
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

#' Crée un graphique de proportion de trouées vs hauteur pour un buffer spécifique
#'
#' @param data Données filtrées pour un plot et un buffer
#' @param buffer Valeur du buffer
#' @param title Titre à afficher sur le graphique
#' @return Un objet ggplot
create_gap_plot <- function(data, buffer, title) {
  ggplot(data, aes(x = height_aboveground, y = proportion)) +
    geom_line(size = 0.7, color = "darkblue") +
    geom_point(size = 2, color = "darkblue", alpha = 0.7) +
    theme_minimal() +
    labs(
      title = paste0("Buffer ", buffer, "m"),
      x = "Height Above Ground (m)",
      y = "Gap Proportion"
    ) +
    theme(
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 9),
      axis.title = element_text(color = "black", size = 10, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90"),
      panel.border = element_rect(color = "black", fill = NA),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(0, 60, by = 5))
}

#### 3️⃣ CHARGEMENT DES DONNÉES ####
# Chemins des fichiers
path_output <- paste0(path0, "output/")
path_plots_info <- "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/3_Gaps/plots_infos.csv"

# Chemin vers les données d'élévation
path_elevation = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/2_ElevationData/"

# Chemin vers les CHM
path_chm_folder <- file.path(path_elevation, "CHM_final")
path_chm <- list.files(path = path_chm_folder, 
                       pattern = "\\.tif$", 
                       full.names = TRUE)

# Extraction des noms des plots
plot_name_chm <- basename(gsub("\\.tif$", "", path_chm))

plots_info <- read_csv2(path_plots_info) %>% 
  as_tibble() %>%
  filter(plot_ref %in% plot_name_chm) %>%
  rename_with(~ ifelse(tolower(.) == "plot_ref", "plot_name", .)) 

plot_name = plots_info$plot_name

# Lecture des métriques de trouées
gaps_metrics <- read_csv2(paste0(path_output, "results_gaps_metrics.csv")) %>%
  filter(plot_name %in% plots_info$plot_name)

#### 4️⃣ PRÉTRAITEMENT ####
# Joindre les informations des plots si la colonne WD_mean existe

combined_data <- gaps_metrics %>%
  left_join(plots_info %>% dplyr::select(plot_name, WD_mean), by = "plot_name")

# Calculer les statistiques pour WD_mean
has_wd_data <- any(!is.na(combined_data$WD_mean))

if (has_wd_data) {
  wd_stats <- plots_info %>%
    filter(!is.na(WD_mean)) %>%
    summarise(
      min_wd = min(WD_mean, na.rm = TRUE),
      max_wd = max(WD_mean, na.rm = TRUE)
    )
  
  # Créer des catégories basées sur WD_mean
  min_wd <- floor(wd_stats$min_wd * 20) / 20  # Arrondi à 0.05 en-dessous
  max_wd <- ceiling(wd_stats$max_wd * 20) / 20  # Arrondi à 0.05 au-dessus
  wd_breaks <- seq(min_wd, max_wd, by = 0.05)
  
  # Créer un vecteur pour les labels numériques
  wd_labels <- 1:length(wd_breaks[-1])
  
  # Ajouter la catégorie aux données
  combined_data <- combined_data %>%
    mutate(wd_category = cut(WD_mean, 
                             breaks = wd_breaks, 
                             include.lowest = TRUE, 
                             labels = wd_labels))
} else {
  combined_data$wd_category <- NA
}
  
#### 5️⃣ ANALYSE ####
# Création des dossiers de sortie
path_graphics <- paste0(path0, "graphics/")
path_buffer_0 <- paste0(path_graphics, "buffer_0/")
path_buffer_20 <- paste0(path_graphics, "buffer_20/")
path_buffer_50 <- paste0(path_graphics, "buffer_50/")
path_buffer_100 <- paste0(path_graphics, "buffer_100/")
path_buffer_200 <- paste0(path_graphics, "buffer_200/")
path_combined <- paste0(path_graphics, "combined/")

lapply(list(path_graphics, path_buffer_0, path_buffer_20, path_buffer_50, 
            path_buffer_100, path_buffer_200, path_combined), create_dir)

# Liste des plots uniques
plot_names <- unique(combined_data$plot_name)
buffers <- c(0, 20, 50, 100, 200)

#### 6️⃣ VISUALISATION ####

#### 6.1 GRAPHIQUES INDIVIDUELS PAR PLOT ####
# Boucle sur chaque plot
for (plot_name_i in plot_names) {
  # Filtrer les données pour le plot actuel
  plot_data <- combined_data %>%
    filter(plot_name == plot_name_i)
  
  # Extraire les infos pour le titre
  if (has_wd_data) {
    plot_info <- plots_info %>% 
      filter(plot_name == plot_name_i) %>% 
      slice(1)
    
    # Obtenir WD_mean en gérant de manière sécurisée
    wd_mean_val <- if ("WD_mean" %in% names(plot_info) && length(plot_info$WD_mean) > 0 && !is.na(plot_info$WD_mean[1])) {
      plot_info$WD_mean[1]
    } else {
      NA
    }
    
    # Obtenir la catégorie WD pour le nom de fichier
    wd_cat_suffix <- ""
    if (!is.na(wd_mean_val)) {
      wd_category_val <- cut(wd_mean_val, 
                             breaks = wd_breaks,
                             include.lowest = TRUE,
                             labels = wd_labels)
      
      wd_cat_suffix <- paste0("_WDcat_", wd_category_val)
    }
    
    # Créer le titre avec les informations disponibles
    if (!is.na(wd_mean_val)) {
      title_info <- paste0("Plot: ", plot_name_i, ", WD Mean: ", round(wd_mean_val, 3))
    } else {
      title_info <- paste0("Plot: ", plot_name_i)
    }
  } else {
    title_info <- paste0("Plot: ", plot_name_i)
    wd_cat_suffix <- ""
  }
  
  # Créer un graphique pour chaque buffer
  buffer_plots <- list()
  
  for (buf in buffers) {
    buffer_data <- plot_data %>%
      filter(buffer == buf)
    
    if (nrow(buffer_data) > 0) {
      p <- create_gap_plot(buffer_data, buf, title_info)
      
      # Sauvegarder le graphique dans le dossier approprié
      path_buffer <- paste0(path_graphics, "buffer_", buf, "/")
      ggsave(
        paste0(path_buffer, plot_name_i, "_buffer_", buf, wd_cat_suffix, ".jpg"), 
        p, 
        width = 8, 
        height = 6,
        dpi = 300,
        device = "jpeg"
      )
      
      buffer_plots[[as.character(buf)]] <- p
    }
  }
  
  #### 6.2 GRAPHIQUES COMBINÉS PAR PLOT ####
  # Créer un graphique combiné pour tous les buffers disponibles
  if (length(buffer_plots) > 0) {
    # Déterminer le nombre de colonnes pour l'arrangement
    n_cols <- min(length(buffer_plots), 5)
    
    # Créer le graphique combiné
    combined_plot <- ggpubr::ggarrange(
      plotlist = buffer_plots[as.character(sort(as.numeric(names(buffer_plots))))],
      ncol = n_cols, 
      common.legend = TRUE
    )
    
    # Ajouter le titre global
    combined_plot <- ggpubr::annotate_figure(
      combined_plot,
      top = ggpubr::text_grob(title_info, size = 14, face = "bold")
    )
    
    # Sauvegarder le graphique combiné
    ggsave(
      paste0(path_combined, plot_name_i, "_combined", wd_cat_suffix, ".jpg"), 
      combined_plot, 
      width = min(n_cols * 5, 25), 
      height = 6,
      dpi = 300,
      device = "jpeg"
    )
  }
  
  # Afficher la progression
  cat("Processed plot:", plot_name_i, "\n")
}

#### 6.3 GRAPHIQUE COMPARATIF PAR DENSITÉ DE BOIS ####
# Créer un dossier pour les graphiques comparatifs
path_comparative <- paste0(path_graphics, "comparative/")
create_dir(path_comparative)

# Pour chaque buffer, créer un graphique comparatif avec gradient de couleur par WD
for (buf in buffers) {
  # Filtrer les données pour le buffer et uniquement les plots avec valeur WD
  buffer_data <- combined_data %>%
    filter(buffer == buf, !is.na(WD_mean))
  
  if (nrow(buffer_data) > 0) {
    # Définir les limites min/max pour l'échelle de couleur
    min_wd <- min(buffer_data$WD_mean, na.rm = TRUE)
    max_wd <- max(buffer_data$WD_mean, na.rm = TRUE)
    
    # Créer un graphique avec toutes les lignes colorées par WD
    p_comparative <- ggplot(buffer_data, aes(x = height_aboveground, y = proportion, group = plot_name, color = WD_mean)) +
      geom_line(size = 1, alpha = 0.8) +
      theme_minimal() +
      labs(
        title = paste0("Distribution des trouées par densité de bois - Buffer ", buf, "m"),
        x = "Hauteur au-dessus du sol (m)",
        y = "Proportion de trouées",
        color = "Densité du bois\n(WD Mean)"
      ) +
      theme(
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        legend.key.height = unit(1.5, "cm")
      ) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      scale_x_continuous(breaks = seq(0, 60, by = 5)) +
      scale_color_viridis_c(
        option = "plasma",  # Palette similaire à celle utilisée pour la topographie (bleu->violet->rouge)
        direction = -1,     # Inverse les couleurs pour avoir bleu (bas) -> rouge (haut)
        limits = c(min_wd, max_wd),
        breaks = seq(min_wd, max_wd, length.out = 5),
        labels = function(x) sprintf("%.3f", x),
        guide = guide_colorbar(
          barwidth = 1,
          barheight = 10,
          ticks.linewidth = 1.5,
          frame.colour = "black",
          ticks.colour = "black"
        )
      )
    
    # Sauvegarder le graphique comparatif
    ggsave(
      paste0(path_comparative, "comparative_WD_buffer_", buf, ".jpg"), 
      p_comparative, 
      width = 12, 
      height = 8,
      dpi = 300,
      device = "jpeg"
    )
    
    # Version sans légende pour insertion dans les publications
    p_no_legend <- p_comparative + 
      theme(legend.position = "none") +
      labs(title = paste0("Buffer ", buf, "m"))
    
    ggsave(
      paste0(path_comparative, "comparative_WD_buffer_", buf, "_nolegend.jpg"), 
      p_no_legend, 
      width = 10, 
      height = 8,
      dpi = 300,
      device = "jpeg"
    )
    
    # Créer une version avec légende seule pour pouvoir l'utiliser séparément
    legend_only <- cowplot::get_legend(
      p_comparative + 
        guides(color = guide_colorbar(title = "Densité du bois (WD Mean)", 
                                      barwidth = 20, 
                                      barheight = 2)) +
        theme(legend.position = "bottom")
    )
    
    legend_plot <- cowplot::ggdraw() + 
      cowplot::draw_grob(legend_only)
    
    ggsave(
      paste0(path_comparative, "legend_WD.jpg"), 
      legend_plot, 
      width = 8, 
      height = 2,
      dpi = 300,
      device = "jpeg"
    )
  }
}

print("Graphics generation completed!")