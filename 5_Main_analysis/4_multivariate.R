# =================================================================
# ANALYSE MULTIVARIÉE - MÉTRIQUES DE TROUÉES ET COMPOSITION TAXONOMIQUE
# Objectif: Explorer les relations entre structure de canopée, traits fonctionnels
#           et composition spécifique par tempérament
# =================================================================

#### 1️⃣ INITIALISATION ####
# Nettoyage de l'environnement et de la mémoire
rm(list = ls())
gc()

# Définition des packages nécessaires
pkgs <- c("tidyverse", "FactoMineR", "factoextra", "ggrepel", "viridis", 
          "corrplot", "patchwork", "ade4", "vegan", "broom", "here")

# Installation et chargement automatique des packages manquants
to_install <- !pkgs %in% installed.packages()
if(any(to_install)) {install.packages(pkgs[to_install])}
inst <- lapply(pkgs, library, character.only = TRUE)

# Configuration de l'environnement
options(scipen = 999)        # Éviter la notation scientifique
set.seed(42)                 # Fixer la graine pour reproductibilité

# Configuration du thème général pour tous les graphiques
theme_set(theme_minimal() + 
            theme(plot.title = element_text(face = "bold"),
                  axis.title = element_text(face = "bold")))

#### 2️⃣ FONCTIONS ####

#### 2.1 Fonctions de visualisation ####

#' Crée un graphique d'ordination avec coloration par variable
#'
#' @param data Tibble contenant les données d'ordination et les métadonnées
#' @param x_col Nom de la colonne pour l'axe X
#' @param y_col Nom de la colonne pour l'axe Y
#' @param color_var Nom de la variable pour la coloration
#' @param title Titre du graphique
#' @param subtitle Sous-titre du graphique
#' @param eigenvalues_df Dataframe contenant les valeurs propres pour les pourcentages
#' @param point_size Taille des points
#' @param alpha Transparence des points
#' @param max_overlaps Nombre maximal de chevauchements pour les étiquettes
#' @param label_size Taille des étiquettes
#' @param id_col Nom de la colonne contenant les identifiants des points
#' @param palette Option de palette de couleur (turbo par défaut)
#' @return Objet ggplot
create_ordination_plot <- function(data, x_col, y_col, color_var, title, subtitle, 
                                   eigenvalues_df, point_size = 3, alpha = 0.8, 
                                   max_overlaps = 30, label_size = 3, id_col = "plot_ref",
                                   palette = "turbo") {
  x_var_name <- sym(x_col)
  y_var_name <- sym(y_col)
  color_var_name <- sym(color_var)
  id_var_name <- sym(id_col)
  
  # Extraire le numéro d'axe à partir du nom de colonne (ex: "NSCA1" -> 1)
  x_axis_num <- as.numeric(gsub("[^0-9]", "", x_col))
  y_axis_num <- as.numeric(gsub("[^0-9]", "", y_col))
  
  # Récupérer les pourcentages d'explication
  x_percent <- round(eigenvalues_df$Percentage[eigenvalues_df$Axis == x_axis_num], 1)
  y_percent <- round(eigenvalues_df$Percentage[eigenvalues_df$Axis == y_axis_num], 1)
  
  # Créer le graphique
  ggplot(data, aes(x = !!x_var_name, y = !!y_var_name, color = !!color_var_name)) +
    geom_point(size = point_size, alpha = alpha) +
    scale_color_viridis_c(option = palette) +
    geom_text_repel(aes(label = !!id_var_name), size = label_size, max.overlaps = max_overlaps, 
                    segment.alpha = 0.6, min.segment.length = 0.2) +
    labs(title = title,
         subtitle = subtitle,
         x = paste0(x_col, " (", x_percent, "%)"),
         y = paste0(y_col, " (", y_percent, "%)")) +
    theme_minimal()
}

#' Crée un graphique de projection des espèces contributives sans les sites
#'
#' @param species_coord Tibble avec les coordonnées des espèces
#' @param species_contrib Tibble avec les contributions des espèces
#' @param x_axis Nom de la colonne pour l'axe X
#' @param y_axis Nom de la colonne pour l'axe Y
#' @param n_species Nombre d'espèces contributives à afficher
#' @param eigenvalues_df Dataframe contenant les valeurs propres pour les pourcentages
#' @param title Titre du graphique optionnel
#' @param species_color Couleur pour les espèces
#' @param species_col Nom de la colonne contenant les identifiants d'espèces
#' @param max_overlaps Nombre maximal de chevauchements pour les étiquettes
#' @param label_size Taille des étiquettes
#' @return Objet ggplot
create_species_projection <- function(species_coord, species_contrib, x_axis, y_axis, 
                                      n_species = 15, eigenvalues_df, title = NULL,
                                      species_color = "darkred", species_col = "species",
                                      max_overlaps = 30, label_size = 3) {
  # Obtenir les numéros d'axes à partir des chaînes d'entrée
  x_axis_num <- as.numeric(gsub("[^0-9]", "", x_axis))
  y_axis_num <- as.numeric(gsub("[^0-9]", "", y_axis))
  
  # Extraire le préfixe des axes (tout ce qui précède le dernier chiffre)
  x_prefix <- gsub("[0-9]+$", "", x_axis)
  
  # Identifier les espèces les plus contributives pour chacun des deux axes
  top_species_x <- get_top_species(species_contrib, x_axis_num, contrib_prefix = x_prefix, taxon_col = species_col)
  top_species_y <- get_top_species(species_contrib, y_axis_num, contrib_prefix = x_prefix, taxon_col = species_col)
  
  # Combiner les deux listes d'espèces et éliminer les doublons
  combined_top_species <- unique(c(top_species_x, top_species_y))
  
  # Pourcentages d'explication des axes
  x_percent <- round(eigenvalues_df$Percentage[eigenvalues_df$Axis == x_axis_num], 1)
  y_percent <- round(eigenvalues_df$Percentage[eigenvalues_df$Axis == y_axis_num], 1)
  
  # Définir un titre par défaut si non fourni
  if (is.null(title)) {
    title <- paste0("Espèces contributives - Axes ", x_axis_num, "-", y_axis_num)
  }
  
  # Créer le graphique
  projection_plot <- ggplot() +
    # Axes
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", size = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 0.3) +
    # Ajouter les espèces contributives
    geom_point(data = species_coord %>% filter(!!sym(species_col) %in% combined_top_species), 
               aes(x = !!sym(x_axis), y = !!sym(y_axis)), shape = 17, size = 3, color = species_color) +
    geom_text_repel(data = species_coord %>% filter(!!sym(species_col) %in% combined_top_species),
                    aes(x = !!sym(x_axis), y = !!sym(y_axis), label = !!sym(species_col)), 
                    color = species_color, size = label_size, max.overlaps = max_overlaps,
                    segment.alpha = 0.7, min.segment.length = 0.2) +
    labs(title = title,
         x = paste0(x_axis, " (", x_percent, "%)"),
         y = paste0(y_axis, " (", y_percent, "%)")) +
    theme_minimal()
  
  return(projection_plot)
}

#' Crée un vrai biplot avec sites et projection des espèces contributives
#'
#' @param plot_coord Tibble avec les coordonnées des parcelles
#' @param species_coord Tibble avec les coordonnées des espèces
#' @param species_contrib Tibble avec les contributions des espèces
#' @param x_axis Nom de la colonne pour l'axe X
#' @param y_axis Nom de la colonne pour l'axe Y
#' @param color_var Nom de la variable pour la coloration des sites
#' @param n_species Nombre d'espèces contributives à afficher
#' @param eigenvalues_df Dataframe contenant les valeurs propres pour les pourcentages
#' @param title Titre du graphique optionnel
#' @param subtitle Sous-titre du graphique optionnel
#' @param id_col Nom de la colonne contenant les identifiants des sites
#' @param species_col Nom de la colonne contenant les identifiants d'espèces
#' @param species_color Couleur pour les lignes de projection des espèces
#' @param palette Option de palette de couleur
#' @param show_site_labels Afficher les étiquettes des sites
#' @return Objet ggplot
create_biplot <- function(plot_coord, species_coord, species_contrib, x_axis, y_axis, 
                          color_var, n_species = 15, eigenvalues_df, 
                          title = NULL, subtitle = NULL, id_col = "plot_ref", 
                          species_col = "species", species_color = "darkred", 
                          palette = "turbo", show_site_labels = FALSE) {
  
  # Obtenir les numéros d'axes
  x_axis_num <- as.numeric(gsub("[^0-9]", "", x_axis))
  y_axis_num <- as.numeric(gsub("[^0-9]", "", y_axis))
  
  # Extraire le préfixe des axes (tout ce qui précède le dernier chiffre)
  x_prefix <- gsub("[0-9]+$", "", x_axis)
  
  # Identifier les espèces les plus contributives pour chacun des deux axes
  top_species_x <- get_top_species(species_contrib, x_axis_num, n = n_species, 
                                   contrib_prefix = x_prefix, taxon_col = species_col)
  top_species_y <- get_top_species(species_contrib, y_axis_num, n = n_species, 
                                   contrib_prefix = x_prefix, taxon_col = species_col)
  
  # Combiner les deux listes d'espèces et éliminer les doublons
  combined_top_species <- unique(c(top_species_x, top_species_y))
  
  # Pourcentages d'explication des axes
  x_percent <- round(eigenvalues_df$Percentage[eigenvalues_df$Axis == x_axis_num], 1)
  y_percent <- round(eigenvalues_df$Percentage[eigenvalues_df$Axis == y_axis_num], 1)
  
  # Définir titre et sous-titre par défaut si non fournis
  if (is.null(title)) {
    title <- paste0("Biplot - Axes ", x_axis_num, "-", y_axis_num)
  }
  
  if (is.null(subtitle)) {
    subtitle <- "Sites et espèces contributives"
  }
  
  # Créer le biplot
  biplot <- ggplot() +
    # Axes
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", size = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 0.3) +
    # Sites
    geom_point(data = plot_coord, 
               aes(x = !!sym(x_axis), y = !!sym(y_axis), color = !!sym(color_var)), 
               size = 3, alpha = 0.7) +
    scale_color_viridis_c(option = palette) +
    # Lignes de projection des espèces
    geom_segment(data = species_coord %>% filter(!!sym(species_col) %in% combined_top_species),
                 aes(x = 0, y = 0, xend = !!sym(x_axis), yend = !!sym(y_axis)),
                 arrow = arrow(length = unit(0.2, "cm")), color = species_color, alpha = 0.5) +
    # Étiquettes des espèces
    geom_text_repel(data = species_coord %>% filter(!!sym(species_col) %in% combined_top_species),
                    aes(x = !!sym(x_axis), y = !!sym(y_axis), label = !!sym(species_col)), 
                    color = species_color, size = 3, max.overlaps = 30,
                    segment.alpha = 0) +
    # Étiquettes des sites (conditionnelles)
    {if(show_site_labels) 
      geom_text_repel(data = plot_coord, aes(x = !!sym(x_axis), y = !!sym(y_axis), label = !!sym(id_col)), 
                      size = 3, max.overlaps = 30, segment.alpha = 0.3)
    } +
    labs(title = title,
         subtitle = subtitle,
         x = paste0(x_axis, " (", x_percent, "%)"),
         y = paste0(y_axis, " (", y_percent, "%)")) +
    theme_minimal()
  
  return(biplot)
}

#' Identifie les espèces/genres les plus contributifs à un axe donné
#'
#' @param contrib_data Tibble avec les contributions
#' @param axis Numéro de l'axe
#' @param n Nombre d'espèces à retourner
#' @param contrib_prefix Préfixe pour les colonnes de contribution (défaut = "NSCA")
#' @param taxon_col Nom de la colonne contenant les identifiants de taxons
#' @return Vecteur des n taxons les plus contributifs
get_top_species <- function(contrib_data, axis, n = 15, 
                            contrib_prefix = NULL, taxon_col = "species") {
  
  # Si aucun préfixe n'est spécifié, extraire le préfixe du nom des colonnes
  if (is.null(contrib_prefix)) {
    # Prendre la première colonne qui n'est pas taxon_col et extraire le préfixe
    col_names <- names(contrib_data)
    col_names <- col_names[col_names != taxon_col]
    if (length(col_names) > 0) {
      # Extraire le préfixe (tout ce qui précède le premier chiffre)
      prefix_match <- regexpr("[0-9]", col_names[1])[1]
      if (prefix_match > 1) {
        contrib_prefix <- substr(col_names[1], 1, prefix_match - 1)
      } else {
        contrib_prefix <- "NSCA" # Fallback au préfixe par défaut
      }
    } else {
      contrib_prefix <- "NSCA" # Fallback au préfixe par défaut
    }
  }
  
  # Convertir l'axe en nom de colonne
  axis_col <- sym(paste0(contrib_prefix, axis))
  
  # Trier les taxons par contribution décroissante et sélectionner les n premiers
  top_taxa <- contrib_data %>%
    arrange(desc(!!axis_col)) %>%
    slice_head(n = n) %>%
    pull(!!sym(taxon_col))
  
  return(top_taxa)
}

#' Crée un scree plot pour visualiser les valeurs propres
#'
#' @param eigenvalues Vecteur des valeurs propres
#' @param title Titre du graphique
#' @param bar_color Couleur des barres
#' @param min_percent Pourcentage minimal pour inclure un axe
#' @param prefix Préfixe pour le nom des axes dans l'étiquette
#' @return Tibble avec les valeurs propres formatées et un objet ggplot
create_scree_plot <- function(eigenvalues, title, bar_color = "#4C72B0", 
                              min_percent = 2, prefix = "") {
  
  # Formater les valeurs propres
  eigenvalues_df <- tibble(
    Axis = 1:length(eigenvalues),
    Eigenvalue = eigenvalues,
    Percentage = eigenvalues / sum(eigenvalues) * 100
  ) %>%
    # Ne garder que les axes expliquant au moins min_percent
    filter(Percentage >= min_percent)
  
  # Créer le graphique
  scree_plot <- ggplot(eigenvalues_df, aes(x = factor(Axis), y = Percentage)) +
    geom_bar(stat = "identity", fill = bar_color, color = "white") +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)), vjust = -0.5, size = 3.5) +
    labs(title = title, 
         x = "Dimensions", 
         y = "% de variance") +
    theme_minimal()
  
  # Retourner à la fois les données et le graphique
  return(list(
    data = eigenvalues_df,
    plot = scree_plot
  ))
}

#' Crée un dashboard d'ordinations avec différentes colorations
#'
#' @param ordinations Liste d'objets ggplot d'ordination
#' @param species_plot Objet ggplot de projection des espèces
#' @param title Titre du dashboard
#' @param layout Mise en page (ex: c(2,2) pour 2x2)
#' @param title_size Taille du titre
#' @return Objet patchwork combinant les graphiques
create_dashboard <- function(ordinations, species_plot, title, 
                             layout = NULL, title_size = 16) {
  
  # Créer un arrangement par défaut avec tous les graphiques + projection des espèces
  if (is.null(layout)) {
    n_plots <- length(ordinations)
    if (n_plots == 3) {
      dashboard <- (ordinations[[1]] + ordinations[[2]]) / 
        (ordinations[[3]] + species_plot)
    } else if (n_plots == 2) {
      dashboard <- (ordinations[[1]] + ordinations[[2]]) / 
        (species_plot + plot_spacer())
    } else {
      # Convertir la liste en grille avec patchwork
      all_plots <- c(ordinations, list(species_plot))
      dashboard <- wrap_plots(all_plots, ncol = 2)
    }
  } else {
    # Utiliser la mise en page personnalisée
    dashboard <- wrap_plots(c(ordinations, list(species_plot)), 
                            ncol = layout[1], nrow = layout[2])
  }
  
  # Ajouter l'annotation
  result <- dashboard + 
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(size = title_size, face = "bold"))
    )
  
  return(result)
}

#' Crée un dashboard comparatif de plusieurs ensembles d'axes
#'
#' @param ordination_list Liste des plots d'ordination
#' @param species_list Liste des plots de projection des espèces
#' @param title Titre du dashboard
#' @param subtitle Sous-titre du dashboard
#' @return Objet patchwork combinant les graphiques
create_comparative_dashboard <- function(ordination_list, species_list, 
                                         title, subtitle = NULL) {
  
  # Vérifier que les listes ont la même longueur
  if (length(ordination_list) != length(species_list)) {
    stop("Les listes d'ordination et de projection doivent avoir la même longueur")
  }
  
  # Assembler les paires d'ordination et de projection
  pairs <- lapply(1:length(ordination_list), function(i) {
    ordination_list[[i]] | species_list[[i]]
  })
  
  # Combiner toutes les paires verticalement
  dashboard <- Reduce(`/`, pairs)
  
  # Ajouter l'annotation
  result <- dashboard +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14)
      )
    ) &
    theme(plot.margin = margin(5, 5, 5, 5))
  
  return(result)
}

#### 2.2 Fonctions d'analyse ####


#' Exécute une analyse multivariée et prépare les résultats
#'
#' @param data_table Table de contingence pour l'analyse
#' @param n_axes Nombre d'axes à retenir
#' @param scannf Demander à l'utilisateur le nombre d'axes à conserver
#' @param prefix Préfixe pour les noms des colonnes dans les résultats (NSCA, CA, etc.)
#' @param bar_color Couleur pour le scree plot
#' @param min_percent Pourcentage minimal pour garder un axe dans le scree plot
#' @param id_col Nom de la colonne d'identifiants pour les observations
#' @param taxon_col Nom de la colonne d'identifiants pour les taxons
#' @param metadata Données supplémentaires à joindre aux coordonnées des observations
#' @param title_scree Titre du scree plot
#' @return Liste avec tous les résultats de l'analyse
run_multivariate_analysis <- function(data_table, n_axes = 5, scannf = FALSE, 
                                      prefix = "", bar_color = "#4C72B0", 
                                      min_percent = 2, id_col = "plot_ref", 
                                      taxon_col = "species", metadata = NULL,
                                      title_scree = NULL, analysis_type = "") {
  
  # Déterminer la fonction d'analyse à utiliser
  analysis_func <- switch(
    tolower(analysis_type),
    "nsca" = dudi.nsc,
    "ca" = dudi.coa,
    "pca" = dudi.pca,
    stop("Type d'analyse non reconnu. Utiliser 'nsca', 'ca' ou 'pca'.")
  )
  
  # Exécuter l'analyse
  result <- analysis_func(data_table, scannf = scannf, nf = n_axes)
  
  # Titre par défaut pour le scree plot
  if (is.null(title_scree)) {
    title_scree <- paste("Valeurs propres", toupper(analysis_type))
  }
  
  # Extraire et formater les résultats
  # 1. Scree plot et données des valeurs propres
  scree_results <- create_scree_plot(
    eigenvalues = result$eig,
    title = title_scree,
    bar_color = bar_color,
    min_percent = min_percent,
    prefix = prefix
  )
  
  # 2. Coordonnées des observations (sites/parcelles)
  plot_coord <- as.data.frame(result$li) %>%
    setNames(paste0(prefix, 1:n_axes)) %>%
    as_tibble()
  
  # Pour l'ACP, les noms de lignes ne correspondent pas aux identifiants des parcelles
  # Il faut ajouter les identifiants manuellement
  if (analysis_type == "pca") {
    # Récupérer les identifiants de parcelles à partir des métadonnées
    if (!is.null(metadata)) {
      plot_ids <- metadata[[id_col]]
      if (length(plot_ids) == nrow(plot_coord)) {
        plot_coord <- plot_coord %>%
          mutate(!!sym(id_col) := plot_ids)
      } else {
        warning("Le nombre d'identifiants de parcelles ne correspond pas au nombre de lignes dans les résultats de l'ACP.")
      }
    }
  } else {
    # Pour les autres analyses, utiliser les noms de lignes
    plot_coord <- plot_coord %>%
      rownames_to_column(id_col)
  }
  
  # Joindre les métadonnées si fournies
  if (!is.null(metadata)) {
    plot_coord <- plot_coord %>%
      left_join(metadata, by = id_col)
  }
  
  # 3. Coordonnées des variables (espèces/genres)
  taxon_coord <- as.data.frame(result$co) %>%
    setNames(paste0(prefix, 1:n_axes)) %>%
    rownames_to_column(taxon_col) %>%
    as_tibble()
  
  # 4. Contributions des variables
  taxon_contrib <- inertia.dudi(result, col.inertia = TRUE)$col.abs %>%
    as.data.frame() %>%
    setNames(paste0(prefix, 1:min(n_axes, ncol(.)))) %>%
    rownames_to_column(taxon_col) %>%
    as_tibble()
  
  # Retourner tous les résultats dans une liste
  return(list(
    result = result,
    eigenvalues = scree_results$data,
    scree_plot = scree_results$plot,
    plot_coord = plot_coord,
    taxon_coord = taxon_coord,
    taxon_contrib = taxon_contrib,
    total_inertia = sum(result$eig)
  ))
}

#' Crée un dataframe de métadonnées pour les parcelles avec les variables appropriées
#'
#' @param metadata Dataframe de métadonnées source
#' @param ponderation Type de pondération ('G' pour surface terrière, 'ind' pour individus)
#' @param plot_col Nom de la colonne contenant les identifiants de parcelle
#' @param target_col_name Nom de la colonne pour les identifiants dans le résultat
#' @return Tibble avec les métadonnées formatées
create_plot_metadata <- function(metadata, ponderation = "G", 
                                 plot_col = "plot_name", 
                                 target_col_name = "plot_ref") {
  
  # Valider le type de pondération
  if (!ponderation %in% c("G", "ind")) {
    stop("Le type de pondération doit être 'G' (surface terrière) ou 'ind' (individus)")
  }
  
  # Choix des colonnes en fonction du type de pondération
  wd_col <- if(ponderation == "G") "WD_BA" else "WD_mean"
  prop_prefix <- if(ponderation == "G") "prop_g_" else "prop_ind_"
  
  # Créer le dataframe avec les colonnes appropriées
  plot_metadata <- metadata %>%
    rename(!!target_col_name := !!sym(plot_col)) %>%
    select(
      !!target_col_name, 
      WD = !!sym(wd_col), 
      prop_helio = !!sym(paste0(prop_prefix, "helio")), 
      prop_npld = !!sym(paste0(prop_prefix, "npld")), 
      prop_shade = !!sym(paste0(prop_prefix, "shade"))
    )
  
  return(plot_metadata)
}

#### 2.3 Fonctions utilitaires ####

#' Crée un chemin relatif au répertoire du projet
#'
#' @param ... Composants du chemin à joindre
#' @return Chemin complet
project_path <- function(...) {
  here(...)
}

#' S'assure que les répertoires existent, les crée si nécessaire
#'
#' @param dirs Vecteur de chemins de répertoires à vérifier/créer
#' @param verbose Afficher un message lors de la création
#' @return Invisible NULL
ensure_directories <- function(dirs, verbose = TRUE) {
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      if (verbose) {
        cat("Répertoire créé:", dir, "\n")
      }
    }
  }
  return(invisible(NULL))
}

#' Exporte les résultats d'analyse en fichiers CSV
#'
#' @param results Liste des résultats d'analyse
#' @param base_path Chemin de base pour l'export
#' @param analysis_name Nom de l'analyse pour les noms de fichiers
#' @param export_contrib Exporter les contributions
#' @param export_coords Exporter les coordonnées
#' @return Invisible NULL
export_analysis_results <- function(results, base_path, analysis_name,
                                    export_contrib = TRUE, export_coords = TRUE) {

  # Exporter les contributions si demandé
  if (export_contrib) {
    contrib_path <- file.path(base_path, paste0(analysis_name, "_contributions.csv"))
    write_csv(results$taxon_contrib, contrib_path)
  }

  # Exporter les coordonnées si demandé
  if (export_coords) {
    # Coordonnées des sites
    sites_path <- file.path(base_path, paste0(analysis_name, "_site_coordinates.csv"))
    write_csv(results$plot_coord, sites_path)

    # Coordonnées des taxons
    taxon_path <- file.path(base_path, paste0(analysis_name, "_taxon_coordinates.csv"))
    write_csv(results$taxon_coord, taxon_path)
  }

  return(invisible(NULL))
}

#### 3️⃣ CHARGEMENT DES DONNÉES ####
# Définition du répertoire de base (à adapter selon votre environnement)
project_dir <- here()

# Définition des chemins essentiels
path_main <- project_path("5_Main_analysis")
path_input <- file.path(path_main, "input")
path_output <- file.path(path_main, "output")

# Chemins pour les résultats
path_figures <- file.path(path_output, "figures", "multivariate")
path_data <- file.path(path_output, "data", "multivariate")

# Sous-répertoires pour les figures
path_figures_taxa <- file.path(path_figures, "taxonomic")
path_figures_pca <- file.path(path_figures, "pca")
path_figures_combined <- file.path(path_figures, "combined")

# Sous-répertoires par type d'analyse
path_figures_species <- file.path(path_figures_taxa, "species")
path_figures_genus <- file.path(path_figures_taxa, "genus")

# S'assurer que tous les répertoires existent
dirs_to_check <- c(
  path_figures_taxa, path_figures_pca, path_figures_combined,
  path_figures_species, path_figures_genus, path_data
)
ensure_directories(dirs_to_check)

# Chargement des données
plot_metrics <- read_csv2(file.path(path_input, "plot_metrics.csv"))
inventory_temperament <- read_csv2(file.path(path_input, "inventory_with_temperament.csv"))

#### 4️⃣ PRÉPARATION DES DONNÉES ####

# Filtrer les parcelles qui ont des données de tempérament
plots_with_temperament <- inventory_temperament %>%
  pull(plot_ref) %>%
  unique()

# Filtrer plot_metrics pour ne garder que ces parcelles
plot_data_complete <- plot_metrics %>%
  filter(plot_name %in% plots_with_temperament)

# Définir les colonnes pour l'ACP et les métadonnées
metadata_cols = c("plot_name", "WD_BA", "prop_g_helio", "prop_g_npld", "prop_g_shade", 
                  "prop_ind_helio", "prop_ind_npld", "prop_ind_shade")

# Métriques issues des sigmoïdes et des analyses structurelles de la canopée
gap_metrics_cols <- names(plot_data_complete)[!names(plot_data_complete) %in% metadata_cols]

# Séparer les données d'ACP et les métadonnées
data_for_pca <- plot_data_complete %>%
  select(all_of(gap_metrics_cols))

# Métadonnées pour interprétation avec toutes les variables d'intérêt
metadata <- plot_data_complete %>%
  select(all_of(metadata_cols))

# Créer les métadonnées formatées pour les deux types de pondération
metadata_G <- create_plot_metadata(metadata, ponderation = "G")
metadata_ind <- create_plot_metadata(metadata, ponderation = "ind")

# 1. Nettoyage des données taxonomiques
# Retirer les espèces avec des noms NA ou problématiques
inventory_clean <- inventory_temperament %>%
  filter(!is.na(species_full) & species_full != "" & species_full != "NA NA")

# 2. Construction des tables de contingence
# 2.1 Table espèces × parcelles
species_plot_table <- inventory_clean %>%
  count(plot_ref, species_full) %>%
  # Convertir au format large (contingence)
  pivot_wider(names_from = species_full, 
              values_from = n,
              values_fill = 0) %>%
  # Conserver plot_ref comme identifiant de ligne
  column_to_rownames("plot_ref")

# 2.2 Table genres × parcelles
genus_plot_table <- inventory_clean %>%
  count(plot_ref, genus) %>%
  pivot_wider(names_from = genus, 
              values_from = n,
              values_fill = 0) %>%
  column_to_rownames("plot_ref")

# 2.3 Table espèces × parcelles pondérée par DBH
species_plot_dbh_table <- inventory_clean %>%
  # Somme des sections par espèce et parcelle (DBH en cm)
  group_by(plot_ref, species_full) %>%
  summarise(basal_area = sum(G), .groups = "drop") %>% # Surface terrière en m²
  # Convertir au format large (contingence)
  pivot_wider(names_from = species_full, 
              values_from = basal_area,
              values_fill = 0) %>%
  # Conserver plot_ref comme identifiant de ligne
  column_to_rownames("plot_ref")

# 2.4 Table genres × parcelles pondérée par DBH
genus_plot_dbh_table <- inventory_clean %>%
  # Somme des sections par genre et parcelle
  group_by(plot_ref, genus) %>%
  summarise(basal_area = sum(G), .groups = "drop") %>%
  # Convertir au format large
  pivot_wider(names_from = genus, 
              values_from = basal_area,
              values_fill = 0) %>%
  column_to_rownames("plot_ref")

# 4. Calcul des indices de diversité pour chaque parcelle
diversity_indices <- species_plot_table %>%
  vegan::diversity(index = "simpson") %>%
  enframe(name = "plot_ref", value = "simpson_diversity") %>%
  mutate(
    # Richesse spécifique (nombre d'espèces)
    species_richness = apply(species_plot_table, 1, function(x) sum(x > 0)),
    # Indice de Shannon
    shannon_diversity = vegan::diversity(species_plot_table, index = "shannon"),
    # Équitabilité de Pielou (Shannon/log(richesse))
    pielou_evenness = shannon_diversity / log(species_richness),
    # Indice de Hill (nombre effectif d'espèces basé sur Shannon)
    hill_diversity = exp(shannon_diversity)
  )

diversity_indices_dbh <- species_plot_dbh_table %>%
  vegan::diversity(index = "simpson") %>%
  enframe(name = "plot_ref", value = "simpson_diversity") %>%
  mutate(
    # Richesse spécifique (nombre d'espèces)
    species_richness = apply(species_plot_table, 1, function(x) sum(x > 0)),
    # Indice de Shannon
    shannon_diversity = vegan::diversity(species_plot_table, index = "shannon"),
    # Équitabilité de Pielou (Shannon/log(richesse))
    pielou_evenness = shannon_diversity / log(species_richness),
    # Indice de Hill (nombre effectif d'espèces basé sur Shannon)
    hill_diversity = exp(shannon_diversity)
  )

#### 5️⃣ ANALYSES TAXONOMIQUES ####

#### 5.1 NSCA pour espèces (non pondérée) ####

# Exécution de l'analyse NSCA pour les espèces (non pondérée)
nsca_species <- run_multivariate_analysis(
  data_table = species_plot_table,
  n_axes = 5,
  prefix = "NSCA",
  bar_color = "#4C72B0",
  metadata = metadata_G,
  title_scree = "Valeurs propres NSCA-Simpson (espèces)",
  analysis_type = "nsca"
)

# Afficher l'inertie totale
cat("Diversité de Simpson (NSCA espèces) :", nsca_species$total_inertia, "\n")

# Créer les ordinations NSCA pour les axes 1-2 avec différentes colorations
nsca_ord_1_2_wd <- create_ordination_plot(
  nsca_species$plot_coord, "NSCA1", "NSCA2", "WD",
  "NSCA-Simpson - Axes 1-2", "Coloration par densité du bois",
  nsca_species$eigenvalues
)

nsca_ord_1_2_helio <- create_ordination_plot(
  nsca_species$plot_coord, "NSCA1", "NSCA2", "prop_helio",
  "NSCA-Simpson - Axes 1-2", "Coloration par proportion d'héliophiles",
  nsca_species$eigenvalues
)

nsca_ord_1_2_shade <- create_ordination_plot(
  nsca_species$plot_coord, "NSCA1", "NSCA2", "prop_shade",
  "NSCA-Simpson - Axes 1-2", "Coloration par proportion d'espèces d'ombre",
  nsca_species$eigenvalues
)

# Créer les projections des espèces contributives pour les axes 1-2
nsca_species_proj_1_2 <- create_species_projection(
  nsca_species$taxon_coord, nsca_species$taxon_contrib,
  "NSCA1", "NSCA2", n_species = 15, 
  nsca_species$eigenvalues, 
  title = "Espèces contributives NSCA-Simpson - Axes 1-2"
)

# Créer le biplot pour les axes 1-2
nsca_biplot_1_2 <- create_biplot(
  nsca_species$plot_coord, nsca_species$taxon_coord, nsca_species$taxon_contrib,
  "NSCA1", "NSCA2", "WD", n_species = 15,
  nsca_species$eigenvalues,
  title = "Biplot NSCA-Simpson - Axes 1-2",
  subtitle = "Sites colorés par WD et projections d'espèces"
)

# Créer le dashboard pour les axes 1-2
nsca_dashboard_1_2 <- create_dashboard(
  list(nsca_ord_1_2_wd, nsca_ord_1_2_helio, nsca_ord_1_2_shade),
  nsca_species_proj_1_2,
  "NSCA-Simpson: Axes 1-2 avec différentes colorations et espèces contributives"
)

# Répéter pour les autres paires d'axes (3-4)
# Créer les ordinations NSCA pour les axes 3-4 avec différentes colorations
nsca_ord_3_4_wd <- create_ordination_plot(
  nsca_species$plot_coord, "NSCA3", "NSCA4", "WD",
  "NSCA-Simpson - Axes 3-4", "Coloration par densité du bois",
  nsca_species$eigenvalues
)

nsca_ord_3_4_helio <- create_ordination_plot(
  nsca_species$plot_coord, "NSCA3", "NSCA4", "prop_helio",
  "NSCA-Simpson - Axes 3-4", "Coloration par proportion d'héliophiles",
  nsca_species$eigenvalues
)

nsca_ord_3_4_shade <- create_ordination_plot(
  nsca_species$plot_coord, "NSCA3", "NSCA4", "prop_shade",
  "NSCA-Simpson - Axes 3-4", "Coloration par proportion d'espèces d'ombre",
  nsca_species$eigenvalues
)

# Créer les projections des espèces contributives pour les axes 3-4
nsca_species_proj_3_4 <- create_species_projection(
  nsca_species$taxon_coord, nsca_species$taxon_contrib,
  "NSCA3", "NSCA4", n_species = 15, 
  nsca_species$eigenvalues, 
  title = "Espèces contributives NSCA-Simpson - Axes 3-4"
)

# Créer le biplot pour les axes 3-4
nsca_biplot_3_4 <- create_biplot(
  nsca_species$plot_coord, nsca_species$taxon_coord, nsca_species$taxon_contrib,
  "NSCA3", "NSCA4", "WD", n_species = 15,
  nsca_species$eigenvalues,
  title = "Biplot NSCA-Simpson - Axes 3-4",
  subtitle = "Sites colorés par WD et projections d'espèces"
)

# Créer le dashboard pour les axes 3-4
nsca_dashboard_3_4 <- create_dashboard(
  list(nsca_ord_3_4_wd, nsca_ord_3_4_helio, nsca_ord_3_4_shade),
  nsca_species_proj_3_4,
  "NSCA-Simpson: Axes 3-4 avec différentes colorations et espèces contributives"
)

# Sauvegarde des dashboards
ggsave(file.path(path_figures_species, "nsca_dashboard_1_2.pdf"), 
       nsca_dashboard_1_2, width = 12, height = 12)
ggsave(file.path(path_figures_species, "nsca_dashboard_3_4.pdf"), 
       nsca_dashboard_3_4, width = 12, height = 12)
ggsave(file.path(path_figures_species, "nsca_biplot_1_2.pdf"), 
       nsca_biplot_1_2, width = 10, height = 8)
ggsave(file.path(path_figures_species, "nsca_biplot_3_4.pdf"), 
       nsca_biplot_3_4, width = 10, height = 8)

#### 5.2 NSCA pour espèces (pondérée par DBH) ####

# Exécution de l'analyse NSCA pour les espèces (pondérée par DBH)
nsca_species_dbh <- run_multivariate_analysis(
  data_table = species_plot_dbh_table,
  n_axes = 5,
  prefix = "NSCA_DBH",
  bar_color = "#B04C72",
  metadata = metadata_G,
  title_scree = "Valeurs propres NSCA-Simpson pondérée par DBH (espèces)",
  analysis_type = "nsca"
)

# Afficher l'inertie totale
cat("Diversité de Simpson pondérée par DBH (NSCA espèces) :", nsca_species_dbh$total_inertia, "\n")

# Créer les ordinations NSCA DBH pour les axes 1-2 avec différentes colorations
nsca_dbh_ord_1_2_wd <- create_ordination_plot(
  nsca_species_dbh$plot_coord, "NSCA_DBH1", "NSCA_DBH2", "WD",
  "NSCA pondérée par DBH - Axes 1-2", "Coloration par densité du bois",
  nsca_species_dbh$eigenvalues
)

nsca_dbh_ord_1_2_helio <- create_ordination_plot(
  nsca_species_dbh$plot_coord, "NSCA_DBH1", "NSCA_DBH2", "prop_helio",
  "NSCA pondérée par DBH - Axes 1-2", "Coloration par proportion d'héliophiles",
  nsca_species_dbh$eigenvalues
)

nsca_dbh_ord_1_2_shade <- create_ordination_plot(
  nsca_species_dbh$plot_coord, "NSCA_DBH1", "NSCA_DBH2", "prop_shade",
  "NSCA pondérée par DBH - Axes 1-2", "Coloration par proportion d'espèces d'ombre",
  nsca_species_dbh$eigenvalues
)

# Créer les projections des espèces contributives pour les axes 1-2
nsca_dbh_species_proj_1_2 <- create_species_projection(
  nsca_species_dbh$taxon_coord, nsca_species_dbh$taxon_contrib,
  "NSCA_DBH1", "NSCA_DBH2", n_species = 15, 
  nsca_species_dbh$eigenvalues, 
  title = "Espèces contributives NSCA pondérée par DBH - Axes 1-2",
  species_color = "darkred"
)

# Créer le biplot pour les axes 1-2
nsca_dbh_biplot_1_2 <- create_biplot(
  nsca_species_dbh$plot_coord, nsca_species_dbh$taxon_coord, nsca_species_dbh$taxon_contrib,
  "NSCA_DBH1", "NSCA_DBH2", "WD", n_species = 15,
  nsca_species_dbh$eigenvalues,
  title = "Biplot NSCA pondérée par DBH - Axes 1-2",
  subtitle = "Sites colorés par WD et projections d'espèces",
  species_color = "darkred"
)

# Créer le dashboard pour les axes 1-2
nsca_dbh_dashboard_1_2 <- create_dashboard(
  list(nsca_dbh_ord_1_2_wd, nsca_dbh_ord_1_2_helio, nsca_dbh_ord_1_2_shade),
  nsca_dbh_species_proj_1_2,
  "NSCA pondérée par DBH: Axes 1-2 avec différentes colorations et espèces contributives"
)

# Répéter pour les axes 3-4
# Créer les ordinations NSCA DBH pour les axes 3-4 avec différentes colorations
nsca_dbh_ord_3_4_wd <- create_ordination_plot(
  nsca_species_dbh$plot_coord, "NSCA_DBH3", "NSCA_DBH4", "WD",
  "NSCA pondérée par DBH - Axes 3-4", "Coloration par densité du bois",
  nsca_species_dbh$eigenvalues
)

nsca_dbh_ord_3_4_helio <- create_ordination_plot(
  nsca_species_dbh$plot_coord, "NSCA_DBH3", "NSCA_DBH4", "prop_helio",
  "NSCA pondérée par DBH - Axes 3-4", "Coloration par proportion d'héliophiles",
  nsca_species_dbh$eigenvalues
)

nsca_dbh_ord_3_4_shade <- create_ordination_plot(
  nsca_species_dbh$plot_coord, "NSCA_DBH3", "NSCA_DBH4", "prop_shade",
  "NSCA pondérée par DBH - Axes 3-4", "Coloration par proportion d'espèces d'ombre",
  nsca_species_dbh$eigenvalues
)

# Créer les projections des espèces contributives pour les axes 3-4
nsca_dbh_species_proj_3_4 <- create_species_projection(
  nsca_species_dbh$taxon_coord, nsca_species_dbh$taxon_contrib,
  "NSCA_DBH3", "NSCA_DBH4", n_species = 15, 
  nsca_species_dbh$eigenvalues, 
  title = "Espèces contributives NSCA pondérée par DBH - Axes 3-4",
  species_color = "darkred"
)

# Créer le biplot pour les axes 3-4
nsca_dbh_biplot_3_4 <- create_biplot(
  nsca_species_dbh$plot_coord, nsca_species_dbh$taxon_coord, nsca_species_dbh$taxon_contrib,
  "NSCA_DBH3", "NSCA_DBH4", "WD", n_species = 15,
  nsca_species_dbh$eigenvalues,
  title = "Biplot NSCA pondérée par DBH - Axes 3-4",
  subtitle = "Sites colorés par WD et projections d'espèces",
  species_color = "darkred"
)

# Créer le dashboard pour les axes 3-4
nsca_dbh_dashboard_3_4 <- create_dashboard(
  list(nsca_dbh_ord_3_4_wd, nsca_dbh_ord_3_4_helio, nsca_dbh_ord_3_4_shade),
  nsca_dbh_species_proj_3_4,
  "NSCA pondérée par DBH: Axes 3-4 avec différentes colorations et espèces contributives"
)

# Sauvegarde des dashboards
ggsave(file.path(path_figures_species, "nsca_dbh_dashboard_1_2.pdf"), 
       nsca_dbh_dashboard_1_2, width = 12, height = 12)
ggsave(file.path(path_figures_species, "nsca_dbh_dashboard_3_4.pdf"), 
       nsca_dbh_dashboard_3_4, width = 12, height = 12)
ggsave(file.path(path_figures_species, "nsca_dbh_biplot_1_2.pdf"), 
       nsca_dbh_biplot_1_2, width = 10, height = 8)
ggsave(file.path(path_figures_species, "nsca_dbh_biplot_3_4.pdf"), 
       nsca_dbh_biplot_3_4, width = 10, height = 8)

#### 5.3 NSCA pour genres (non pondérée) ####

# Exécution de l'analyse NSCA pour les genres (non pondérée)
nsca_genus <- run_multivariate_analysis(
  data_table = genus_plot_table,
  n_axes = 5,
  prefix = "NSCA_Genus",
  bar_color = "#4C72B0",
  metadata = metadata_G,
  title_scree = "Valeurs propres NSCA-Simpson (genres)",
  analysis_type = "nsca",
  taxon_col = "genus"
)

# Afficher l'inertie totale
cat("Diversité de Simpson (NSCA genres) :", nsca_genus$total_inertia, "\n")

# Créer les ordinations NSCA pour les axes 1-2 avec différentes colorations
nsca_genus_ord_1_2_wd <- create_ordination_plot(
  nsca_genus$plot_coord, "NSCA_Genus1", "NSCA_Genus2", "WD",
  "NSCA-Simpson Genres - Axes 1-2", "Coloration par densité du bois",
  nsca_genus$eigenvalues
)

nsca_genus_ord_1_2_helio <- create_ordination_plot(
  nsca_genus$plot_coord, "NSCA_Genus1", "NSCA_Genus2", "prop_helio",
  "NSCA-Simpson Genres - Axes 1-2", "Coloration par proportion d'héliophiles",
  nsca_genus$eigenvalues
)

nsca_genus_ord_1_2_shade <- create_ordination_plot(
  nsca_genus$plot_coord, "NSCA_Genus1", "NSCA_Genus2", "prop_shade",
  "NSCA-Simpson Genres - Axes 1-2", "Coloration par proportion d'espèces d'ombre",
  nsca_genus$eigenvalues
)

# Créer les projections des genres contributifs pour les axes 1-2
nsca_genus_proj_1_2 <- create_species_projection(
  nsca_genus$taxon_coord, nsca_genus$taxon_contrib,
  "NSCA_Genus1", "NSCA_Genus2", n_species = 15, 
  nsca_genus$eigenvalues, 
  title = "Genres contributifs NSCA-Simpson - Axes 1-2",
  species_color = "darkblue",
  species_col = "genus"
)

# Créer le biplot pour les axes 1-2
nsca_genus_biplot_1_2 <- create_biplot(
  nsca_genus$plot_coord, nsca_genus$taxon_coord, nsca_genus$taxon_contrib,
  "NSCA_Genus1", "NSCA_Genus2", "WD", n_species = 15,
  nsca_genus$eigenvalues,
  title = "Biplot NSCA-Simpson Genres - Axes 1-2",
  subtitle = "Sites colorés par WD et projections de genres",
  species_color = "darkblue",
  species_col = "genus"
)

# Créer le dashboard pour les axes 1-2
nsca_genus_dashboard_1_2 <- create_dashboard(
  list(nsca_genus_ord_1_2_wd, nsca_genus_ord_1_2_helio, nsca_genus_ord_1_2_shade),
  nsca_genus_proj_1_2,
  "NSCA-Simpson Genres: Axes 1-2 avec différentes colorations et genres contributifs"
)

# Répéter pour les axes 3-4
# Créer les ordinations NSCA pour les axes 3-4 avec différentes colorations
nsca_genus_ord_3_4_wd <- create_ordination_plot(
  nsca_genus$plot_coord, "NSCA_Genus3", "NSCA_Genus4", "WD",
  "NSCA-Simpson Genres - Axes 3-4", "Coloration par densité du bois",
  nsca_genus$eigenvalues
)

nsca_genus_ord_3_4_helio <- create_ordination_plot(
  nsca_genus$plot_coord, "NSCA_Genus3", "NSCA_Genus4", "prop_helio",
  "NSCA-Simpson Genres - Axes 3-4", "Coloration par proportion d'héliophiles",
  nsca_genus$eigenvalues
)

nsca_genus_ord_3_4_shade <- create_ordination_plot(
  nsca_genus$plot_coord, "NSCA_Genus3", "NSCA_Genus4", "prop_shade",
  "NSCA-Simpson Genres - Axes 3-4", "Coloration par proportion d'espèces d'ombre",
  nsca_genus$eigenvalues
)

# Créer les projections des genres contributifs pour les axes 3-4
nsca_genus_proj_3_4 <- create_species_projection(
  nsca_genus$taxon_coord, nsca_genus$taxon_contrib,
  "NSCA_Genus3", "NSCA_Genus4", n_species = 15, 
  nsca_genus$eigenvalues, 
  title = "Genres contributifs NSCA-Simpson - Axes 3-4",
  species_color = "darkblue",
  species_col = "genus"
)

# Créer le biplot pour les axes 3-4
nsca_genus_biplot_3_4 <- create_biplot(
  nsca_genus$plot_coord, nsca_genus$taxon_coord, nsca_genus$taxon_contrib,
  "NSCA_Genus3", "NSCA_Genus4", "WD", n_species = 15,
  nsca_genus$eigenvalues,
  title = "Biplot NSCA-Simpson Genres - Axes 3-4",
  subtitle = "Sites colorés par WD et projections de genres",
  species_color = "darkblue",
  species_col = "genus"
)

# Créer le dashboard pour les axes 3-4
nsca_genus_dashboard_3_4 <- create_dashboard(
  list(nsca_genus_ord_3_4_wd, nsca_genus_ord_3_4_helio, nsca_genus_ord_3_4_shade),
  nsca_genus_proj_3_4,
  "NSCA-Simpson Genres: Axes 3-4 avec différentes colorations et genres contributifs"
)

# Sauvegarde des dashboards
ggsave(file.path(path_figures_genus, "nsca_genus_dashboard_1_2.pdf"), 
       nsca_genus_dashboard_1_2, width = 12, height = 12)
ggsave(file.path(path_figures_genus, "nsca_genus_dashboard_3_4.pdf"), 
       nsca_genus_dashboard_3_4, width = 12, height = 12)
ggsave(file.path(path_figures_genus, "nsca_genus_biplot_1_2.pdf"), 
       nsca_genus_biplot_1_2, width = 10, height = 8)
ggsave(file.path(path_figures_genus, "nsca_genus_biplot_3_4.pdf"), 
       nsca_genus_biplot_3_4, width = 10, height = 8)

#### 5.4 NSCA pour genres (pondérée par DBH) ####

# Exécution de l'analyse NSCA pour les genres (pondérée par DBH)
nsca_genus_dbh <- run_multivariate_analysis(
  data_table = genus_plot_dbh_table,
  n_axes = 5,
  prefix = "NSCA_Genus_DBH",
  bar_color = "#B04C72",
  metadata = metadata_G,
  title_scree = "Valeurs propres NSCA-Simpson pondérée par DBH (genres)",
  analysis_type = "nsca",
  taxon_col = "genus"
)

# Afficher l'inertie totale
cat("Diversité de Simpson pondérée par DBH (NSCA genres) :", nsca_genus_dbh$total_inertia, "\n")

# Créer les ordinations NSCA DBH pour les axes 1-2 avec différentes colorations
nsca_genus_dbh_ord_1_2_wd <- create_ordination_plot(
  nsca_genus_dbh$plot_coord, "NSCA_Genus_DBH1", "NSCA_Genus_DBH2", "WD",
  "NSCA Genres pondérée par DBH - Axes 1-2", "Coloration par densité du bois",
  nsca_genus_dbh$eigenvalues
)

nsca_genus_dbh_ord_1_2_helio <- create_ordination_plot(
  nsca_genus_dbh$plot_coord, "NSCA_Genus_DBH1", "NSCA_Genus_DBH2", "prop_helio",
  "NSCA Genres pondérée par DBH - Axes 1-2", "Coloration par proportion d'héliophiles",
  nsca_genus_dbh$eigenvalues
)

nsca_genus_dbh_ord_1_2_shade <- create_ordination_plot(
  nsca_genus_dbh$plot_coord, "NSCA_Genus_DBH1", "NSCA_Genus_DBH2", "prop_shade",
  "NSCA Genres pondérée par DBH - Axes 1-2", "Coloration par proportion d'espèces d'ombre",
  nsca_genus_dbh$eigenvalues
)

# Créer les projections des genres contributifs pour les axes 1-2
nsca_genus_dbh_proj_1_2 <- create_species_projection(
  nsca_genus_dbh$taxon_coord, nsca_genus_dbh$taxon_contrib,
  "NSCA_Genus_DBH1", "NSCA_Genus_DBH2", n_species = 15, 
  nsca_genus_dbh$eigenvalues, 
  title = "Genres contributifs NSCA pondérée par DBH - Axes 1-2",
  species_color = "darkmagenta",
  species_col = "genus"
)

# Créer le biplot pour les axes 1-2
nsca_genus_dbh_biplot_1_2 <- create_biplot(
  nsca_genus_dbh$plot_coord, nsca_genus_dbh$taxon_coord, nsca_genus_dbh$taxon_contrib,
  "NSCA_Genus_DBH1", "NSCA_Genus_DBH2", "WD", n_species = 15,
  nsca_genus_dbh$eigenvalues,
  title = "Biplot NSCA Genres pondérée par DBH - Axes 1-2",
  subtitle = "Sites colorés par WD et projections de genres",
  species_color = "darkmagenta",
  species_col = "genus"
)

# Créer le dashboard pour les axes 1-2
nsca_genus_dbh_dashboard_1_2 <- create_dashboard(
  list(nsca_genus_dbh_ord_1_2_wd, nsca_genus_dbh_ord_1_2_helio, nsca_genus_dbh_ord_1_2_shade),
  nsca_genus_dbh_proj_1_2,
  "NSCA Genres pondérée par DBH: Axes 1-2 avec différentes colorations et genres contributifs"
)

# Répéter pour les axes 3-4
# Créer les ordinations NSCA DBH pour les axes 3-4 avec différentes colorations
nsca_genus_dbh_ord_3_4_wd <- create_ordination_plot(
  nsca_genus_dbh$plot_coord, "NSCA_Genus_DBH3", "NSCA_Genus_DBH4", "WD",
  "NSCA Genres pondérée par DBH - Axes 3-4", "Coloration par densité du bois",
  nsca_genus_dbh$eigenvalues
)

nsca_genus_dbh_ord_3_4_helio <- create_ordination_plot(
  nsca_genus_dbh$plot_coord, "NSCA_Genus_DBH3", "NSCA_Genus_DBH4", "prop_helio",
  "NSCA Genres pondérée par DBH - Axes 3-4", "Coloration par proportion d'héliophiles",
  nsca_genus_dbh$eigenvalues
)

nsca_genus_dbh_ord_3_4_shade <- create_ordination_plot(
  nsca_genus_dbh$plot_coord, "NSCA_Genus_DBH3", "NSCA_Genus_DBH4", "prop_shade",
  "NSCA Genres pondérée par DBH - Axes 3-4", "Coloration par proportion d'espèces d'ombre",
  nsca_genus_dbh$eigenvalues
)

# Créer les projections des genres contributifs pour les axes 3-4
nsca_genus_dbh_proj_3_4 <- create_species_projection(
  nsca_genus_dbh$taxon_coord, nsca_genus_dbh$taxon_contrib,
  "NSCA_Genus_DBH3", "NSCA_Genus_DBH4", n_species = 15, 
  nsca_genus_dbh$eigenvalues, 
  title = "Genres contributifs NSCA pondérée par DBH - Axes 3-4",
  species_color = "darkmagenta",
  species_col = "genus"
)

# Créer le biplot pour les axes 3-4
nsca_genus_dbh_biplot_3_4 <- create_biplot(
  nsca_genus_dbh$plot_coord, nsca_genus_dbh$taxon_coord, nsca_genus_dbh$taxon_contrib,
  "NSCA_Genus_DBH3", "NSCA_Genus_DBH4", "WD", n_species = 15,
  nsca_genus_dbh$eigenvalues,
  title = "Biplot NSCA Genres pondérée par DBH - Axes 3-4",
  subtitle = "Sites colorés par WD et projections de genres",
  species_color = "darkmagenta",
  species_col = "genus"
)

# Créer le dashboard pour les axes 3-4
nsca_genus_dbh_dashboard_3_4 <- create_dashboard(
  list(nsca_genus_dbh_ord_3_4_wd, nsca_genus_dbh_ord_3_4_helio, nsca_genus_dbh_ord_3_4_shade),
  nsca_genus_dbh_proj_3_4,
  "NSCA Genres pondérée par DBH: Axes 3-4 avec différentes colorations et genres contributifs"
)

# Sauvegarde des dashboards
ggsave(file.path(path_figures_genus, "nsca_genus_dbh_dashboard_1_2.pdf"), 
       nsca_genus_dbh_dashboard_1_2, width = 12, height = 12)
ggsave(file.path(path_figures_genus, "nsca_genus_dbh_dashboard_3_4.pdf"), 
       nsca_genus_dbh_dashboard_3_4, width = 12, height = 12)
ggsave(file.path(path_figures_genus, "nsca_genus_dbh_biplot_1_2.pdf"), 
       nsca_genus_dbh_biplot_1_2, width = 10, height = 8)
ggsave(file.path(path_figures_genus, "nsca_genus_dbh_biplot_3_4.pdf"), 
       nsca_genus_dbh_biplot_3_4, width = 10, height = 8)

#### 5.5 CA pour espèces (richesse) ####

# Exécution de l'analyse CA pour les espèces
ca_species <- run_multivariate_analysis(
  data_table = species_plot_table,
  n_axes = 5,
  prefix = "CA",
  bar_color = "#72B04C",
  metadata = metadata_G,
  title_scree = "Valeurs propres CA-Richness (espèces)",
  analysis_type = "ca"
)

# Afficher l'inertie totale
cat("Richesse spécifique - 1 (CA) :", ca_species$total_inertia, "\n")

# Créer les ordinations CA pour les axes 1-2 avec différentes colorations
ca_ord_1_2_wd <- create_ordination_plot(
  ca_species$plot_coord, "CA1", "CA2", "WD",
  "CA-Richness - Axes 1-2", "Coloration par densité du bois",
  ca_species$eigenvalues
)

ca_ord_1_2_helio <- create_ordination_plot(
  ca_species$plot_coord, "CA1", "CA2", "prop_helio",
  "CA-Richness - Axes 1-2", "Coloration par proportion d'héliophiles",
  ca_species$eigenvalues
)

ca_ord_1_2_shade <- create_ordination_plot(
  ca_species$plot_coord, "CA1", "CA2", "prop_shade",
  "CA-Richness - Axes 1-2", "Coloration par proportion d'espèces d'ombre",
  ca_species$eigenvalues
)

# Créer les projections des espèces contributives pour les axes 1-2
ca_species_proj_1_2 <- create_species_projection(
  ca_species$taxon_coord, ca_species$taxon_contrib,
  "CA1", "CA2", n_species = 15, 
  ca_species$eigenvalues, 
  title = "Espèces contributives CA-Richness - Axes 1-2",
  species_color = "darkgreen"
)

# Créer le biplot pour les axes 1-2
ca_biplot_1_2 <- create_biplot(
  ca_species$plot_coord, ca_species$taxon_coord, ca_species$taxon_contrib,
  "CA1", "CA2", "WD", n_species = 15,
  ca_species$eigenvalues,
  title = "Biplot CA-Richness - Axes 1-2",
  subtitle = "Sites colorés par WD et projections d'espèces",
  species_color = "darkgreen"
)

# Créer le dashboard pour les axes 1-2
ca_dashboard_1_2 <- create_dashboard(
  list(ca_ord_1_2_wd, ca_ord_1_2_helio, ca_ord_1_2_shade),
  ca_species_proj_1_2,
  "CA-Richness: Axes 1-2 avec différentes colorations et espèces contributives"
)

# Répéter pour les axes 3-4
# Créer les ordinations CA pour les axes 3-4 avec différentes colorations
ca_ord_3_4_wd <- create_ordination_plot(
  ca_species$plot_coord, "CA3", "CA4", "WD",
  "CA-Richness - Axes 3-4", "Coloration par densité du bois",
  ca_species$eigenvalues
)

ca_ord_3_4_helio <- create_ordination_plot(
  ca_species$plot_coord, "CA3", "CA4", "prop_helio",
  "CA-Richness - Axes 3-4", "Coloration par proportion d'héliophiles",
  ca_species$eigenvalues
)

ca_ord_3_4_shade <- create_ordination_plot(
  ca_species$plot_coord, "CA3", "CA4", "prop_shade",
  "CA-Richness - Axes 3-4", "Coloration par proportion d'espèces d'ombre",
  ca_species$eigenvalues
)

# Créer les projections des espèces contributives pour les axes 3-4
ca_species_proj_3_4 <- create_species_projection(
  ca_species$taxon_coord, ca_species$taxon_contrib,
  "CA3", "CA4", n_species = 15, 
  ca_species$eigenvalues, 
  title = "Espèces contributives CA-Richness - Axes 3-4",
  species_color = "darkgreen"
)

# Créer le biplot pour les axes 3-4
ca_biplot_3_4 <- create_biplot(
  ca_species$plot_coord, ca_species$taxon_coord, ca_species$taxon_contrib,
  "CA3", "CA4", "WD", n_species = 15,
  ca_species$eigenvalues,
  title = "Biplot CA-Richness - Axes 3-4",
  subtitle = "Sites colorés par WD et projections d'espèces",
  species_color = "darkgreen"
)

# Créer le dashboard pour les axes 3-4
ca_dashboard_3_4 <- create_dashboard(
  list(ca_ord_3_4_wd, ca_ord_3_4_helio, ca_ord_3_4_shade),
  ca_species_proj_3_4,
  "CA-Richness: Axes 3-4 avec différentes colorations et espèces contributives"
)

# Sauvegarde des dashboards
ggsave(file.path(path_figures_species, "ca_dashboard_1_2.pdf"), 
       ca_dashboard_1_2, width = 12, height = 12)
ggsave(file.path(path_figures_species, "ca_dashboard_3_4.pdf"), 
       ca_dashboard_3_4, width = 12, height = 12)
ggsave(file.path(path_figures_species, "ca_biplot_1_2.pdf"), 
       ca_biplot_1_2, width = 10, height = 8)
ggsave(file.path(path_figures_species, "ca_biplot_3_4.pdf"), 
       ca_biplot_3_4, width = 10, height = 8)

#### 5.6 Dashboards comparatifs ####

# 1. Comparaison des valeurs propres
scree_comparison <- (nsca_species$scree_plot + ca_species$scree_plot) / 
  (nsca_species_dbh$scree_plot + nsca_genus$scree_plot) /
  (nsca_genus_dbh$scree_plot + plot_spacer()) +
  plot_annotation(
    title = "Comparaison des valeurs propres entre les différentes analyses",
    subtitle = "Différentes approches pour l'analyse de la composition taxonomique",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

# 2. Comparaison des biplots axes 1-2
biplot_comparison_1_2 <- (nsca_biplot_1_2 + ca_biplot_1_2) / 
  (nsca_dbh_biplot_1_2 + nsca_genus_biplot_1_2) +
  plot_annotation(
    title = "Comparaison des biplots entre les différentes analyses - Axes 1-2",
    subtitle = "Sites colorés par WD et projections de taxons",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

# 3. Comparaison des biplots axes 3-4
biplot_comparison_3_4 <- (nsca_biplot_3_4 + ca_biplot_3_4) / 
  (nsca_dbh_biplot_3_4 + nsca_genus_biplot_3_4) +
  plot_annotation(
    title = "Comparaison des biplots entre les différentes analyses - Axes 3-4",
    subtitle = "Sites colorés par WD et projections de taxons",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

# 4. Dashboard comparatif pour toutes les combinaisons d'axes
nsca_dashboard_comparatif <- create_comparative_dashboard(
  list(nsca_ord_1_2_wd, nsca_ord_3_4_wd, ca_ord_1_2_wd, ca_ord_3_4_wd),
  list(nsca_species_proj_1_2, nsca_species_proj_3_4, ca_species_proj_1_2, ca_species_proj_3_4),
  "Comparaison des axes et taxons contributifs",
  "Analyses taxonomiques avec coloration par densité du bois (WD)"
)

# Sauvegarde des dashboards comparatifs
ggsave(file.path(path_figures_taxa, "scree_comparison.pdf"), 
       scree_comparison, width = 12, height = 12)
ggsave(file.path(path_figures_taxa, "biplot_comparison_1_2.pdf"), 
       biplot_comparison_1_2, width = 12, height = 12)
ggsave(file.path(path_figures_taxa, "biplot_comparison_3_4.pdf"), 
       biplot_comparison_3_4, width = 12, height = 12)
ggsave(file.path(path_figures_taxa, "nsca_dashboard_comparatif.pdf"), 
       nsca_dashboard_comparatif, width = 12, height = 24)

#### 6️⃣ INTÉGRATION TAXONOMIE-TRAITS ####

# Préparer un dataframe combinant tous les axes des analyses taxonomiques
all_axes_data <- nsca_species$plot_coord %>%
  select(plot_ref, NSCA1, NSCA2, NSCA3, NSCA4) %>%
  left_join(
    nsca_species_dbh$plot_coord %>%
      select(plot_ref, NSCA_DBH1, NSCA_DBH2, NSCA_DBH3, NSCA_DBH4),
    by = "plot_ref"
  ) %>%
  left_join(
    ca_species$plot_coord %>%
      select(plot_ref, CA1, CA2, CA3, CA4),
    by = "plot_ref"
  ) %>%
  left_join(
    nsca_genus$plot_coord %>%
      select(plot_ref, NSCA_Genus1, NSCA_Genus2, NSCA_Genus3, NSCA_Genus4),
    by = "plot_ref"
  ) %>%
  left_join(
    nsca_genus_dbh$plot_coord %>%
      select(plot_ref, NSCA_Genus_DBH1, NSCA_Genus_DBH2, NSCA_Genus_DBH3, NSCA_Genus_DBH4),
    by = "plot_ref"
  )

# Joindre les indices de diversité
all_axes_data <- all_axes_data %>%
  left_join(diversity_indices, by = "plot_ref")

# Joindre toutes les données de traits fonctionnels
taxa_with_traits <- all_axes_data %>%
  left_join(metadata, by = c("plot_ref" = "plot_name"))

# Définir les variables d'axes taxonomiques pour la corrélation
taxa_vars <- c(
  # NSCA espèces
  "NSCA1", "NSCA2", "NSCA3", "NSCA4",
  # NSCA espèces pondérée par DBH
  "NSCA_DBH1", "NSCA_DBH2", "NSCA_DBH3", "NSCA_DBH4",
  # CA espèces
  "CA1", "CA2", "CA3", "CA4",
  # NSCA genres
  "NSCA_Genus1", "NSCA_Genus2", "NSCA_Genus3", "NSCA_Genus4",
  # NSCA genres pondérée par DBH
  "NSCA_Genus_DBH1", "NSCA_Genus_DBH2", "NSCA_Genus_DBH3", "NSCA_Genus_DBH4"
)

# Définir les variables de traits fonctionnels pour la corrélation
trait_vars <- c(
  "WD_mean", "WD_BA", 
  "prop_ind_helio", "prop_ind_npld", "prop_ind_shade",
  "prop_g_helio", "prop_g_npld", "prop_g_shade",
  "simpson_diversity", "shannon_diversity", "species_richness"
)

# Calculer les corrélations entre axes taxonomiques et traits fonctionnels
taxa_trait_cor <- cor(
  taxa_with_traits %>% select(all_of(taxa_vars)),
  taxa_with_traits %>% select(all_of(trait_vars)),
  use = "pairwise.complete.obs"
)

# Visualiser les corrélations
taxa_trait_cor_plot <- corrplot(
  taxa_trait_cor, 
  method = "color", 
  addCoef.col = "black", 
  is.corr = FALSE,
  col = viridis(100, option = "turbo"),
  title = "Corrélations entre axes taxonomiques et traits fonctionnels"
)

# Sauvegarder la visualisation
ggsave(file.path(path_figures_combined, "taxa_trait_correlations.pdf"), 
       ggplotify::as.ggplot(function() {
         corrplot(
           taxa_trait_cor, 
           method = "color", 
           addCoef.col = "black", 
           is.corr = FALSE,
           col = viridis(100, option = "turbo"),
           title = "Corrélations entre axes taxonomiques et traits fonctionnels"
         )
       }), 
       width = 16, height = 14)

#### 7️⃣ ANALYSE EN COMPOSANTES PRINCIPALES ####

#### 7.1 Sélection optimisée des variables ####

# 1. Examiner les corrélations entre les métriques et les variables d'intérêt
# Créer un dataframe combinant les métriques et les variables d'intérêt
correlation_data <- plot_data_complete %>%
  select(all_of(gap_metrics_cols), WD_mean, WD_BA, prop_g_helio, prop_g_shade, prop_ind_helio, prop_ind_shade)

# Calculer la matrice de corrélation
cor_matrix <- cor(correlation_data, use = "pairwise.complete.obs")

# Visualiser les corrélations avec les variables d'intérêt
target_vars <- c("WD_mean", "WD_BA", "prop_g_helio", "prop_g_shade", "prop_ind_helio", "prop_ind_shade")
cor_with_targets <- cor_matrix[!rownames(cor_matrix) %in% target_vars, target_vars, drop = FALSE]

# Tri des variables par force de corrélation absolue avec les cibles
cor_strength <- apply(abs(cor_with_targets), 1, mean)
cor_ordered <- sort(cor_strength, decreasing = TRUE)

# Visualiser les 15 variables les plus corrélées avec nos variables d'intérêt
top_metrics <- names(cor_ordered)[1:15]

# Heatmap des corrélations pour les variables les plus importantes
cor_subset <- cor_matrix[c(top_metrics, target_vars), c(top_metrics, target_vars)]
cor_heatmap <- corrplot(cor_subset, method = "color", type = "upper", order = "hclust", 
                        tl.col = "black", tl.srt = 80, addCoef.col = "black", 
                        col = viridis(100, option = "turbo"))

# 2. Clustering hiérarchique des variables pour identifier les groupes redondants
# Distance basée sur la corrélation (1 - |cor|)
dist_cor <- as.dist(1 - abs(cor_matrix[gap_metrics_cols, gap_metrics_cols]))

# Clustering hiérarchique
hc <- hclust(dist_cor, method = "ward.D2")

# Visualiser le dendrogramme
plot(hc, hang = -1, cex = 0.6, main = "Clustering des variables")

# Découper en k clusters (à ajuster selon le dendrogramme)
k <- 8  # À ajuster après visualisation
clusters <- cutree(hc, k = k)

# Créer un dataframe avec les variables et leurs clusters
var_clusters <- data.frame(
  variable = names(clusters),
  cluster = as.factor(clusters)
)

# 3. Sélection des représentants de chaque cluster en fonction de leur corrélation avec les variables cibles
# Pour chaque cluster, trouver la variable avec la plus forte corrélation moyenne avec les cibles
best_representatives <- sapply(1:k, function(i) {
  cluster_vars <- var_clusters$variable[var_clusters$cluster == i]
  if(length(cluster_vars) == 1) return(cluster_vars)
  
  # Calculer la corrélation moyenne avec les variables cibles
  mean_cor <- sapply(cluster_vars, function(var) {
    mean(abs(cor_matrix[var, target_vars]))
  })
  
  # Retourner la variable avec la plus forte corrélation moyenne
  return(cluster_vars[which.max(mean_cor)])
})

# Afficher les variables sélectionnées
selected_metrics <- unname(best_representatives)
cat("Variables sélectionnées (les plus représentatives de chaque cluster) :\n")
print(selected_metrics)

#### 7.2 ACP avec G (WD_BA, prop_g_*) ####

# Réaliser une ACP avec uniquement les variables sélectionnées
selected_data <- plot_data_complete %>%
  select(all_of(selected_metrics)) %>%
  drop_na()

# Réaliser l'ACP
pca_G <- run_multivariate_analysis(
  data_table = selected_data,
  n_axes = 5,
  prefix = "PC",
  bar_color = "#4C72B0",
  metadata = metadata_G,
  title_scree = "Valeurs propres ACP (pondération G)",
  analysis_type = "pca"
)

# Calculer les contributions des variables aux axes
pca_G_var_contrib <- as.data.frame(pca_G$result$co) %>%
  setNames(paste0("PC", 1:ncol(.))) %>%
  rownames_to_column("variable") %>%
  as_tibble()

# Créer le biplot standard sans coloration
pca_G_biplot_standard <- fviz_pca_biplot(pca_G$result,
                                         label = "var",
                                         col.var = "black",
                                         col.ind = "grey70",
                                         repel = TRUE) +
  theme_minimal() +
  ggtitle("ACP - Projection des parcelles et variables (pondération G)")

# Créer des visualisations colorées pour les axes 1-2
pca_G_1_2_wd <- create_ordination_plot(
  pca_G$plot_coord, "PC1", "PC2", "WD",
  "ACP - Axes 1-2", "Coloration par densité du bois (WD_BA)",
  pca_G$eigenvalues
)

pca_G_1_2_helio <- create_ordination_plot(
  pca_G$plot_coord, "PC1", "PC2", "prop_helio",
  "ACP - Axes 1-2", "Coloration par proportion d'héliophiles (G)",
  pca_G$eigenvalues
)

pca_G_1_2_shade <- create_ordination_plot(
  pca_G$plot_coord, "PC1", "PC2", "prop_shade",
  "ACP - Axes 1-2", "Coloration par proportion d'espèces d'ombre (G)",
  pca_G$eigenvalues
)

# Créer des visualisations colorées pour les axes 3-4
pca_G_3_4_wd <- create_ordination_plot(
  pca_G$plot_coord, "PC3", "PC4", "WD",
  "ACP - Axes 3-4", "Coloration par densité du bois (WD_BA)",
  pca_G$eigenvalues
)

pca_G_3_4_helio <- create_ordination_plot(
  pca_G$plot_coord, "PC3", "PC4", "prop_helio",
  "ACP - Axes 3-4", "Coloration par proportion d'héliophiles (G)",
  pca_G$eigenvalues
)

pca_G_3_4_shade <- create_ordination_plot(
  pca_G$plot_coord, "PC3", "PC4", "prop_shade",
  "ACP - Axes 3-4", "Coloration par proportion d'espèces d'ombre (G)",
  pca_G$eigenvalues
)

# Créer le dashboard pour les axes 1-2
pca_G_dashboard_1_2 <- (pca_G_1_2_wd + pca_G_1_2_helio) /
  (pca_G_1_2_shade + pca_G_biplot_standard) +
  plot_annotation(
    title = "ACP - Axes 1-2 (pondération G)",
    subtitle = "Métriques de trouées et traits fonctionnels",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14)
    )
  )

# Créer le dashboard pour les axes 3-4
pca_G_dashboard_3_4 <- (pca_G_3_4_wd + pca_G_3_4_helio) /
  (pca_G_3_4_shade + pca_G$scree_plot) +
  plot_annotation(
    title = "ACP - Axes 3-4 (pondération G)",
    subtitle = "Métriques de trouées et traits fonctionnels",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14)
    )
  )

# Sauvegarde des figures
ggsave(file.path(path_figures_pca, "pca_G_dashboard_1_2.pdf"), 
       pca_G_dashboard_1_2, width = 16, height = 12)
ggsave(file.path(path_figures_pca, "pca_G_dashboard_3_4.pdf"), 
       pca_G_dashboard_3_4, width = 16, height = 12)
ggsave(file.path(path_figures_pca, "pca_G_biplot_standard.pdf"), 
       pca_G_biplot_standard, width = 10, height = 8)

#### 7.3 ACP avec ind (WD_mean, prop_ind_*) ####

# Réaliser une ACP avec les mêmes variables mais avec métadonnées ind
pca_ind <- run_multivariate_analysis(
  data_table = selected_data,
  n_axes = 5,
  prefix = "PC",
  bar_color = "#72B04C",
  metadata = metadata_ind,
  title_scree = "Valeurs propres ACP (pondération ind)",
  analysis_type = "pca"
)

# Créer le biplot standard sans coloration
pca_ind_biplot_standard <- fviz_pca_biplot(pca_ind$result,
                                           label = "var",
                                           col.var = "black",
                                           col.ind = "grey70",
                                           repel = TRUE) +
  theme_minimal() +
  ggtitle("ACP - Projection des parcelles et variables (pondération ind)")

# Créer des visualisations colorées pour les axes 1-2
pca_ind_1_2_wd <- create_ordination_plot(
  pca_ind$plot_coord, "PC1", "PC2", "WD",
  "ACP - Axes 1-2", "Coloration par densité du bois (WD_mean)",
  pca_ind$eigenvalues
)

pca_ind_1_2_helio <- create_ordination_plot(
  pca_ind$plot_coord, "PC1", "PC2", "prop_helio",
  "ACP - Axes 1-2", "Coloration par proportion d'héliophiles (ind)",
  pca_ind$eigenvalues
)

pca_ind_1_2_shade <- create_ordination_plot(
  pca_ind$plot_coord, "PC1", "PC2", "prop_shade",
  "ACP - Axes 1-2", "Coloration par proportion d'espèces d'ombre (ind)",
  pca_ind$eigenvalues
)

# Créer des visualisations colorées pour les axes 3-4
pca_ind_3_4_wd <- create_ordination_plot(
  pca_ind$plot_coord, "PC3", "PC4", "WD",
  "ACP - Axes 3-4", "Coloration par densité du bois (WD_mean)",
  pca_ind$eigenvalues
)

pca_ind_3_4_helio <- create_ordination_plot(
  pca_ind$plot_coord, "PC3", "PC4", "prop_helio",
  "ACP - Axes 3-4", "Coloration par proportion d'héliophiles (ind)",
  pca_ind$eigenvalues
)

pca_ind_3_4_shade <- create_ordination_plot(
  pca_ind$plot_coord, "PC3", "PC4", "prop_shade",
  "ACP - Axes 3-4", "Coloration par proportion d'espèces d'ombre (ind)",
  pca_ind$eigenvalues
)

# Créer le dashboard pour les axes 1-2
pca_ind_dashboard_1_2 <- (pca_ind_1_2_wd + pca_ind_1_2_helio) /
  (pca_ind_1_2_shade + pca_ind_biplot_standard) +
  plot_annotation(
    title = "ACP - Axes 1-2 (pondération ind)",
    subtitle = "Métriques de trouées et traits fonctionnels",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14)
    )
  )

# Créer le dashboard pour les axes 3-4
pca_ind_dashboard_3_4 <- (pca_ind_3_4_wd + pca_ind_3_4_helio) /
  (pca_ind_3_4_shade + pca_ind$scree_plot) +
  plot_annotation(
    title = "ACP - Axes 3-4 (pondération ind)",
    subtitle = "Métriques de trouées et traits fonctionnels",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14)
    )
  )

# Sauvegarde des figures
ggsave(file.path(path_figures_pca, "pca_ind_dashboard_1_2.pdf"), 
       pca_ind_dashboard_1_2, width = 16, height = 12)
ggsave(file.path(path_figures_pca, "pca_ind_dashboard_3_4.pdf"), 
       pca_ind_dashboard_3_4, width = 16, height = 12)
ggsave(file.path(path_figures_pca, "pca_ind_biplot_standard.pdf"), 
       pca_ind_biplot_standard, width = 10, height = 8)

#### 7.4 Corrélations entre composantes et traits ####

# Combiner les scores des composantes principales de l'ACP G avec les traits fonctionnels
pc_with_traits <- pca_G$plot_coord %>%
  select(plot_ref, starts_with("PC")) %>%
  left_join(metadata, by = c("plot_ref" = "plot_name"))

# Calculer les corrélations entre composantes principales et traits pour G
pc_trait_cor <- cor(
  pc_G_with_traits %>% select(PC1, PC2, PC3, PC4, PC5),
  pc_G_with_traits %>% select(WD_mean, WD_BA, prop_g_helio, prop_g_shade, prop_ind_helio, prop_ind_shade),
  use = "pairwise.complete.obs"
)

# Visualiser les corrélations pour G
pc_cor_plot <- corrplot(
  pc_trait_cor, 
  method = "color", 
  addCoef.col = "black", 
  col = viridis(100, option = "turbo"),
  title = "Corrélations entre composantes principales et traits - Pondération G"
)


# Sauvegarde des figures
ggsave(file.path(path_figures_pca, "pc_correlations.pdf"), 
       ggplotify::as.ggplot(function() {
         corrplot(
           pc_trait_cor, 
           method = "color", 
           addCoef.col = "black", 
           col = viridis(100, option = "turbo"),
           title = "Corrélations entre composantes principales et traits - Pondération G"
         )
       }), 
       width = 12, height = 6)

# Afficher les corrélations pour les 3 premières composantes
cat("Corrélations entre les composantes principales et les traits (pondération G) :\n")
print(round(pc_G_trait_cor[1:3, ], 2))

#### 8️⃣ EXPORT DES RÉSULTATS ####

# Créer un dataframe complet avec tous les résultats
all_results <- all_axes_data %>%
  # Ajouter les coordonnées de l'ACP
  left_join(
    pca_G$plot_coord %>%
      select(plot_ref, all_of(paste0("PC", 1:5))),
    by = "plot_ref"
  ) %>%
  # Ajouter les métadonnées
  left_join(metadata, by = c("plot_ref" = "plot_name"))

# Exporter les résultats globaux
write_csv(all_results, file.path(path_data, "multivariate_combined_results.csv"))

# Exporter les contributions des variables aux composantes principales
write_csv(pca_G_var_contrib, file.path(path_data, "pca_variable_contributions.csv"))

# Exporter les matrices de corrélation
# Pour l'ACP G
write.csv(
  as.data.frame(pc_G_trait_cor),
  file.path(path_data, "pc_trait_correlations.csv")
)

# Pour les analyses taxonomiques
write.csv(
  as.data.frame(taxa_trait_cor),
  file.path(path_data, "taxa_trait_correlations.csv")
)

# Exporter les contributions des analyses taxonomiques
export_analysis_results(
  nsca_species,
  path_data,
  "nsca_species",
  export_contrib = TRUE,
  export_coords = TRUE
)

export_analysis_results(
  nsca_species_dbh,
  path_data,
  "nsca_species_dbh",
  export_contrib = TRUE,
  export_coords = TRUE
)

export_analysis_results(
  ca_species,
  path_data,
  "ca_species",
  export_contrib = TRUE,
  export_coords = TRUE
)

export_analysis_results(
  nsca_genus,
  path_data,
  "nsca_genus",
  export_contrib = TRUE,
  export_coords = TRUE
)

export_analysis_results(
  nsca_genus_dbh,
  path_data,
  "nsca_genus_dbh",
  export_contrib = TRUE,
  export_coords = TRUE
)

# Message de fin
cat("\n===== ANALYSE MULTIVARIÉE TERMINÉE =====\n")
cat("Tous les résultats ont été exportés dans:", path_data, "\n")
cat("Toutes les figures ont été exportées dans:", path_figures, "\n")
