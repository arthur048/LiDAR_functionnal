rm(list=ls())
gc()

# Initialisation ----------------------------------------------------------
pkgs <- c("tidyverse", "rio", "BIOMASS", "foreach", "doParallel", "sf", "terra", "lidR", "gridExtra", "grid")

dir.create(file.path(tempdir(), "BIOMASS"), showWarnings = FALSE, recursive = TRUE)

to_install <- !pkgs %in% installed.packages()
if(any(to_install)) install.packages(pkgs[to_install])
inst <- lapply(pkgs, library, character.only = TRUE)

# Définition du chemin de base
path0 <- "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/4_Plots_metrics/"

# Lecture des données -------------------------------------------------------------------------
plot_name = read.csv2("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/final_plot_name.csv") %>% pull(plot_name)

plots_data <- rio::import("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/final/plots_inventories.csv", 
                          dec = ",", sep = ";") %>%
  dplyr::as_tibble() %>%
  dplyr::filter(dbh > 0 & plot_ref %in% plot_name) %>%  # Garde uniquement les arbres avec dbh valide
  dplyr::mutate(sp = str_replace_all(sp, " ", "_")) %>%
  separate(sp, into = c("genus", "species"), sep = "_") %>%  # Séparation du nom d'espèce
  dplyr::select(plot_ref, site, genus, species, dbh, h)

plots_afrisar_wd <- rio::import("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/final/plots_AfriSAR_WD.csv", 
                                dec = ",", sep = ";") %>%
  as_tibble() %>%
  dplyr::filter(plot_ref %in% plot_name)

# Lecture du fichier plots_info pour les métadonnées
plots_info <- rio::import(
  "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/final/plots_info.csv",
  dec = ",", sep = ";") %>%
  as_tibble() %>%
  filter(plot_ref %in% plot_name)

# =================================================================
# SCRIPT DE CALCUL DES MÉTRIQUES À L'ÉCHELLE DE L'ARBRE ET DE LA PARCELLE
# 
# Ce script calcule différentes métriques forestières à partir de données d'inventaire :
# 1. Correction taxonomique et complétion des noms d'espèces
# 2. Attribution des wood density (WD) avec priorité aux données locales (CoForTraits)
# 3. Estimation des hauteurs manquantes avec modèle local ou régional
# 4. Calcul des métriques (surface terrière, biomasse) à l'échelle de l'arbre
#
# Données requises :
# - plots_data : données d'inventaire (colonnes : plot_ref, site, sp, dbh, h)
# - cofortraits.csv : base de données locale de wood density
# 
# Packages requis :
# - tidyverse : manipulation des données
# - BIOMASS : fonctions de correction taxonomique et calcul de biomasse
# - rio : import/export de données
# =================================================================

# 1. Préparation des données CoForTraits ----
# CoForTraits est une base de données locale de wood density.
# On l'utilise en priorité car elle contient des mesures locales plus précises
# que la base de données globale (Global Wood Density Database)
CoFor_raw <- rio::import(
  "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/cofortraits.csv",
  sep = ';', dec = ',', encoding = "UTF-8"
) %>%
  as_tibble() %>%
  tidyr::separate(Species, into = c("genus", "species"), sep = " ") %>%
  filter(!is.na(id_taxon_name))

# Correction taxonomique et calcul des WD moyens par espèce
# Cette étape assure l'uniformité des noms d'espèces et calcule une WD moyenne
# quand plusieurs mesures sont disponibles pour une même espèce
meanWD_by_species_CoFor <- CoFor_raw %>%
  mutate(
    # Correction des noms selon les standards taxonomiques actuels
    taxo_info = correctTaxo(genus, species, useCache = TRUE, verbose = FALSE),
    genus = taxo_info$genusCorrected,
    species = taxo_info$speciesCorrected,
    species_full = paste(genus, species)
  ) %>%
  # Calcul des WD moyens
  filter(trait_concept_name == "densité du bois") %>%
  mutate(numeric_value = as.numeric(str_extract(sample_value, "^[0-9]+\\.[0-9]+"))) %>%
  group_by(species_full) %>%
  summarise(
    WD_CoFor = mean(numeric_value, na.rm = TRUE),
    n_samples_CoFor = n()
  ) %>%
  filter(!is.na(WD_CoFor))

# 2. Traitement des données d'inventaire ----
# Traitement principal des données
field_inventories <- plots_data %>%
  # 2.1 Correction taxonomique
  mutate(
    taxo_info = correctTaxo(genus, species, useCache = TRUE, verbose = FALSE),
    genus = taxo_info$genusCorrected,
    species = taxo_info$speciesCorrected,
    species_full = paste(genus, species),
    correctTaxo_BIOMASS = taxo_info$nameModified
  ) %>%
  # 2.2 Ajout des niveaux taxonomiques supérieurs (famille et ordre)
  # Nécessaire pour la recherche de WD quand l'espèce n'est pas connue
  bind_cols(
    getTaxonomy(.$genus, findOrder = TRUE) %>%
      as_tibble() %>%
      dplyr::select(fam = family, order = order)
  ) %>%
  # 2.3 Estimation des wood density (avec priorité à CoForTraits)
  left_join(meanWD_by_species_CoFor, by = "species_full") %>%
  mutate(
    # Récupération des données Global Wood Density
    GWD_info = getWoodDensity(family = fam, genus = genus, species = species, stand = plot_ref),
    # Application de la priorité
    source_WD = if_else(!is.na(WD_CoFor), "CoForTrait", "GlobalWD"),
    meanWD = coalesce(WD_CoFor, GWD_info$meanWD),
    levelWD = case_when(
      !is.na(WD_CoFor) ~ "species",
      TRUE ~ GWD_info$levelWD
    )
  ) %>%
  # 2.4 Estimation des hauteurs manquantes
  # Utilise d'abord le modèle local, sinon le modèle régional pour l'Afrique Centrale
  mutate(
    h = if_else(
      is.na(h),
      retrieveH(D = dbh, region = "CAfrica")$H,
      h
    )
  ) %>%
  # 2.5 Calcul des métriques par parcelle
  group_by(plot_ref) %>%
  mutate(
    G = pi * (0.01 * dbh/2)^2,  # Surface terrière (m²)
    G_rel = G/sum(G),  # Surface terrière relative
    WD_pond = meanWD * G_rel,  # WD pondérée par la surface terrière
    AGB_Chave = computeAGB(D = dbh, WD = meanWD, H = h)  # Biomasse selon Chave et al. 2014
  ) %>%
  ungroup() %>%
  # 2.6 Sélection et organisation des colonnes finales
  dplyr::select(
    # Informations de base
    plot_ref, site, species_full, order, fam, genus, species,
    # Mesures et estimations
    dbh, h, meanWD, G, G_rel, WD_pond, levelWD, source_WD,
    # Résultats finaux et métadonnées
    AGB_Chave, correctTaxo_BIOMASS, n_samples_CoFor
  )

# 3. Export des résultats ----
rio::export(
  meanWD_by_species_CoFor,
  file = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/meanWD_by_species_CoFor.csv",
  sep = ";", dec = ",")

rio::export(
  field_inventories,
  file = paste0(path0, "FieldDataInventories_computation.csv"),
  sep = ";", dec = ",")

# =================================================================
# ANALYSE DE L'IMPACT DU NIVEAU D'IDENTIFICATION DE LA WOOD DENSITY
# 
# Cette partie analyse et visualise l'impact du niveau taxonomique 
# d'identification de la wood density (espèce, genre, famille, parcelle)
# sur l'estimation de la biomasse (AGB), surface terrière (G) et nombre d'arbres
# =================================================================

# 1. Calcul des métriques par parcelle et niveau taxonomique ----
summary_by_plot <- field_inventories %>%
  group_by(plot_ref, site) %>%
  summarise(
    # Biomasse (AGB) par niveau taxonomique et total
    AGB_Chave_species = sum(AGB_Chave[levelWD == 'species'], na.rm = TRUE),
    AGB_Chave_genus = sum(AGB_Chave[levelWD == 'genus'], na.rm = TRUE),
    AGB_Chave_family = sum(AGB_Chave[levelWD == 'family'], na.rm = TRUE),
    AGB_Chave_plot = sum(AGB_Chave[levelWD == 'plot'], na.rm = TRUE),
    AGB_Chave_Total = sum(AGB_Chave, na.rm = TRUE),
    
    # Surface terrière (G) par niveau taxonomique et total
    G_species = sum(G[levelWD == 'species'], na.rm = TRUE),
    G_genus = sum(G[levelWD == 'genus'], na.rm = TRUE),
    G_family = sum(G[levelWD == 'family'], na.rm = TRUE),
    G_plot = sum(G[levelWD == 'plot'], na.rm = TRUE),
    G_Total = sum(G, na.rm = TRUE),
    
    # Nombre d'arbres par niveau taxonomique et total
    Count_species = sum(levelWD == 'species'),
    Count_genus = sum(levelWD == 'genus'),
    Count_family = sum(levelWD == 'family'),
    Count_plot = sum(levelWD == 'plot'),
    Total_Count = n(),
    
    # Calcul des ratios
    AGB_Ratio_species = AGB_Chave_species / AGB_Chave_Total,
    AGB_Ratio_genus = AGB_Chave_genus / AGB_Chave_Total,
    AGB_Ratio_family = AGB_Chave_family / AGB_Chave_Total,
    AGB_Ratio_plot = AGB_Chave_plot / AGB_Chave_Total,
    
    G_Ratio_species = G_species / G_Total,
    G_Ratio_genus = G_genus / G_Total,
    G_Ratio_family = G_family / G_Total,
    G_Ratio_plot = G_plot / G_Total,
    
    Count_Ratio_species = Count_species / Total_Count,
    Count_Ratio_genus = Count_genus / Total_Count,
    Count_Ratio_family = Count_family / Total_Count,
    Count_Ratio_plot = Count_plot / Total_Count,
    
    .groups = 'drop'
  )

# 2. Création des visualisations ----
# Fonction générique pour créer les boxplots
create_ratio_plot <- function(data, ratio_type, ylabel, show_legend = FALSE) {
  data %>%
    pivot_longer(
      cols = starts_with(paste0(ratio_type, "_Ratio_")), 
      names_to = "Metric", 
      values_to = "Value"
    ) %>%
    mutate(
      Metric = factor(
        Metric,
        levels = paste0(ratio_type, "_Ratio_", c("species", "genus", "family", "plot"))
      ),
      Metric = fct_recode(
        Metric,
        Species = paste0(ratio_type, "_Ratio_species"),
        Genus = paste0(ratio_type, "_Ratio_genus"),
        Family = paste0(ratio_type, "_Ratio_family"),
        Plot = paste0(ratio_type, "_Ratio_plot")
      )
    ) %>%
    ggplot(aes(x = Metric, y = Value, fill = Metric)) +
    geom_boxplot() +
    theme_minimal() +
    labs(
      x = "Niveau taxonomique",
      y = ylabel,
      title = paste("Distribution des ratios de", ylabel)
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.position = if(show_legend) "right" else "none"
    ) +
    if(show_legend) guides(fill = guide_legend(title = "Niveau")) else guides(fill = "none")
}

# Création des trois graphiques
p1 <- create_ratio_plot(summary_by_plot, "AGB", "Biomasse (%)", show_legend = FALSE)
p2 <- create_ratio_plot(summary_by_plot, "G", "Surface terrière (%)", show_legend = FALSE)
p3 <- create_ratio_plot(summary_by_plot, "Count", "Nombre d'arbres (%)", show_legend = TRUE)


# Extraction de la légende du dernier plot
shared_legend <- get_legend(p3)
p3 <- create_ratio_plot(summary_by_plot, "Count", "Nombre d'arbres (%)", show_legend = FALSE)

# 3. Assemblage et sauvegarde de la figure ----
plot_combined <- grid.arrange(
  arrangeGrob(
    p1, p2, p3, 
    nrow = 1,
    left = textGrob("Proportion", rot = 90, gp = gpar(fontsize = 12))
  ),
  shared_legend,
  widths = c(8, 1),
  top = textGrob(
    "Impact du niveau d'identification sur la détermination de Wood Density",
    gp = gpar(fontsize = 14, fontface = "bold")
  )
)

# Sauvegarde de la figure
ggsave(
  paste0(path0, "ratio_niveaux_WD_impacts.svg"),
  plot = plot_combined,
  width = 15,
  height = 6,
  device = "svg"
)
