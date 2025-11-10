rm(list=ls())
gc()

# Définition des packages nécessaires
pkgs = c("rio","foreach", "doParallel","lidR", "terra", "sf", "tidyverse", "ggplot2", "BIOMASS", "here")

to_install = !pkgs %in% installed.packages() ; if(any(to_install)) {install.packages(pkgs[to_install])} ; inst = lapply(pkgs, library, character.only = TRUE) # load them

path0 = here("4_Plots_metrics")

path_output = file.path(path0, "output_temperament")
if(!dir.exists(path_output)){dir.create(path_output)}

setwd(path_output)

# =========================================================================
# CHARGEMENT DES DONNÉES BRUTES
# =========================================================================
# Chargement de la base de données CoforTraits
CoFor_raw = rio::import(here("cofortraits.csv"),
                        sep = ';', dec = ',', encoding = "UTF-8") %>%
  as_tibble() %>%
  tidyr::separate(Species, into = c("genus", "species"), sep = " ") %>%
  filter(!is.na(id_taxon_name))

plot_names <- read.csv2(here("final_plot_name.csv")) %>%
  pull(plot_name)

# Chargement des données d'inventaire terrain
field_inventories = rio::import(
  paste0(path0, "FieldDataInventories_computation.csv"),
  sep = ";",
  dec = ","
) %>% 
  as_tibble() %>%
  filter(plot_ref %in% plot_names)

# =========================================================================
# DÉFINITION DES FONCTIONS DE TEMPÉRAMENT
# =========================================================================

#' Détermine le tempérament consensus basé sur les occurrences des catégories
#'
#' Cette fonction applique les règles de décision pour déterminer le tempérament
#' dominant en se basant sur les occurrences de chaque catégorie.
#'
#' @param a Occurrences de "catégorie inconnue"
#' @param b Occurrences de "tolérante à l'ombrage"
#' @param c Occurrences de "héliophile non pionnière"
#' @param d Occurrences de "pionnière" (ou "héliophile" selon le niveau)
#' @param e Occurrences de "other" (catégories non standard)
#' @param subdivision Nombre de catégories de tempérament (2 ou 3)
#' @return Tempérament consensus
determine_temperament <- function(a, b, c, d, e, subdivision) {
  if (subdivision == 2) {
    case_when(
      a > b & a > c & a > d ~ "catégorie inconnue",
      b >= a & b >= c & b > d ~ "tolérante à l'ombrage",
      (c > b & c >= d & c >= a) | (d >= c & d >= a & d > b) ~ "héliophile",
      b == d & c > 0 ~ "héliophile",
      b == d & (a > 0 | e > 0) ~ "no-consensus",
      b == 0 & d == 0 & c == 0 & (a > 0 | e > 0) ~ "no-consensus",
      TRUE ~ "need check"
    )
  } else if (subdivision == 3) {
    case_when(
      a > b & a > c & a > d ~ "catégorie inconnue",
      b >= a & b >= c & b > d ~ "tolérante à l'ombrage",
      c > b & c >= d & c >= a ~ "héliophile non pionnière",
      d > c & d >= a & d > b ~ "héliophile",
      b == d & c > 0 ~ "héliophile non pionnière",
      c == d & d > 0 ~ "héliophile non pionnière",
      b == d & (a > 0 | e > 0) ~ "no-consensus",
      b == 0 & d == 0 & c == 0 & (a > 0 | e > 0) ~ "no-consensus",
      TRUE ~ "need check"
    )
  } else {
    "Subdivision doit être égal à 2 ou 3"
  }
}




#' Détermine les tempéraments à tous les niveaux taxonomiques
#'
#' Cette fonction combine l'extraction des tempéraments au niveau espèce,
#' puis l'agrégation au niveau genre et famille en une seule opération.
#' Elle suit exactement la même logique que les fonctions CoFor.UniqueTemperament et
#' CoFor.UniqueTemperament_taxonomic_consensus originales.
#'
#' @param CoFor Base de données CoforTraits
#' @param subdivision Nombre de catégories de tempérament (2 ou 3)
#' @return Liste de 3 dataframes: species, genus, family avec leurs tempéraments
get_all_temperaments <- function(CoFor, subdivision = 2) {
  # Vérification du paramètre
  if (!subdivision %in% c(2, 3)) {
    stop("Le paramètre subdivision doit être 2 ou 3")
  }
  
  # 1. Extraction des informations taxonomiques
  taxo <- CoFor %>%
    distinct(species_full, Genus, Family)
  
  # 2. Traitement des données de tempérament au niveau espèce
  # (reproduit exactement le code de CoFor.UniqueTemperament)
  species_db <- CoFor %>%
    # Filtrer pour ne garder que les données de tempérament
    filter(trait_concept_name == 'tempérament (guilde de régénération)') %>% 
    dplyr::select(species_full, attributed_category_with_thesaurus, SampleValueMRM) %>%
    # Utiliser attributed_category_with_thesaurus si SampleValueMRM est vide
    mutate(SampleValueMRM = case_when(SampleValueMRM == "" ~ attributed_category_with_thesaurus, TRUE ~ SampleValueMRM)) %>%
    dplyr::select(species_full, Temperament = SampleValueMRM) %>%
    # Standardiser les catégories de tempérament
    mutate(Temperament = case_when(
      Temperament %in% c("catégorie inconnue", "tolérante à l'ombrage", "héliophile non pionnière", "pionnière") ~ Temperament,
      TRUE ~ "other"
    )) %>%
    # Compter les occurrences par espèce et tempérament
    count(species_full, Temperament, name = "count") %>%
    # Transformer en format large
    pivot_wider(names_from = Temperament, values_from = count, values_fill = list(count = 0)) %>%
    # Joindre avec les informations taxonomiques
    left_join(taxo, by = "species_full") %>%
    mutate(n_subdiv = subdivision)
  
  # Gestion des colonnes manquantes (comme dans le code original)
  if (!'other' %in% colnames(species_db)) {
    species_db <- species_db %>%
      mutate(a = `catégorie inconnue`,
             b = `tolérante à l'ombrage`,
             c = `héliophile non pionnière`,
             d = `pionnière`,
             e = 0)
  } else {
    species_db <- species_db %>%
      mutate(a = `catégorie inconnue`,
             b = `tolérante à l'ombrage`,
             c = `héliophile non pionnière`,
             d = `pionnière`,
             e = `other`)
  }
  
  # Détermination du tempérament au niveau espèce
  species_db <- species_db %>%
    mutate(max_temp = determine_temperament(a, b, c, d, e, subdivision)) %>%
    dplyr::select(-a, -b, -c, -d, -e) # Suppression des colonnes temporaires
  
  # 3. Agrégation au niveau du genre
  # (reproduit exactement le code de CoFor.UniqueTemperament_taxonomic_consensus)
  genus_db <- species_db %>%
    group_by(Genus) %>%
    summarize(
      `catégorie inconnue` = sum(`catégorie inconnue`),
      `tolérante à l'ombrage` = sum(`tolérante à l'ombrage`),
      `héliophile non pionnière` = sum(`héliophile non pionnière`),
      pionnière = sum(pionnière),
      other = sum(other),
      .groups = 'drop'
    ) %>%
    mutate(
      a = `catégorie inconnue`,
      b = `tolérante à l'ombrage`,
      c = `héliophile non pionnière`,
      d = pionnière,
      e = other
    )
  
  # Détermination du tempérament au niveau genre
  genus_db <- genus_db %>%
    mutate(max_temp = determine_temperament(a, b, c, d, e, subdivision)) %>%
    dplyr::select(-a, -b, -c, -d, -e) # Suppression des colonnes temporaires
  
  # 4. Agrégation au niveau de la famille
  # (reproduit exactement le code de CoFor.UniqueTemperament_taxonomic_consensus)
  family_db <- species_db %>%
    group_by(Family) %>%
    summarize(
      `catégorie inconnue` = sum(max_temp == "catégorie inconnue"),
      `tolérante à l'ombrage` = sum(max_temp == "tolérante à l'ombrage"),
      `héliophile non pionnière` = sum(max_temp == "héliophile non pionnière"),
      pionnière = sum(max_temp == "héliophile" | max_temp == "pionnière"),
      other = sum(max_temp == "no-consensus" | max_temp == "need check"),
      .groups = 'drop'
    ) %>%
    mutate(
      a = `catégorie inconnue`,
      b = `tolérante à l'ombrage`,
      c = `héliophile non pionnière`,
      d = pionnière,
      e = other
    )
  
  # Détermination du tempérament au niveau famille
  family_db <- family_db %>%
    mutate(max_temp = determine_temperament(a, b, c, d, e, subdivision)) %>%
    dplyr::select(-a, -b, -c, -d, -e) # Suppression des colonnes temporaires
  
  # Retourner tous les résultats
  return(list(
    species = species_db,
    genus = genus_db,
    family = family_db
  ))
}

#' Attribue le tempérament de manière hiérarchique (espèce > genre > famille)
#'
#' @param data Données d'inventaire
#' @param sp_temp Tempéraments au niveau espèce
#' @param genus_temp Tempéraments au niveau genre
#' @param family_temp Tempéraments au niveau famille
#' @param suffix Suffixe pour les noms de colonnes (optionnel)
#' @return Données avec tempéraments et niveau taxonomique attribués
assign_hierarchical_temperament <- function(data, sp_temp, genus_temp, family_temp, suffix = "") {
  # Noms des colonnes de sortie
  temp_col <- paste0("Temperament", suffix)
  level_col <- paste0("level_Temperament", suffix)
  
  # Noms des colonnes temporaires
  sp_col <- paste0("Temperament.sp", suffix)
  genus_col <- paste0("Temperament.genus", suffix)
  family_col <- paste0("Temperament.family", suffix)
  
  # 1. Jointure avec les différents niveaux taxonomiques
  result <- data %>%
    left_join(
      sp_temp %>% dplyr::select(species_full, !!sym(sp_col) := max_temp),
      by = "species_full"
    ) %>%
    left_join(
      genus_temp %>% dplyr::select(Genus, !!sym(genus_col) := max_temp),
      by = c("genus" = "Genus")
    ) %>%
    left_join(
      family_temp %>% dplyr::select(Family, !!sym(family_col) := max_temp),
      by = c("fam" = "Family")
    )
  
  # 2. Attribution hiérarchique (espèce > genre > famille)
  result <- result %>%
    mutate(
      # Tempérament - prendre le premier disponible dans l'ordre hiérarchique
      !!sym(temp_col) := case_when(
        !is.na(!!sym(sp_col)) ~ !!sym(sp_col),
        !is.na(!!sym(genus_col)) ~ !!sym(genus_col),
        !is.na(!!sym(family_col)) ~ !!sym(family_col),
        TRUE ~ NA_character_
      ),
      
      # Niveau taxonomique utilisé pour l'attribution
      !!sym(level_col) := case_when(
        !is.na(!!sym(sp_col)) ~ "species",
        !is.na(!!sym(genus_col)) ~ "genus", 
        !is.na(!!sym(family_col)) ~ "family",
        TRUE ~ NA_character_
      )
    ) %>%
    # Supprimer les colonnes temporaires
    dplyr::select(-!!sym(sp_col), -!!sym(genus_col), -!!sym(family_col))
  
  return(result)
}


# =========================================================================
# CORRECTION TAXONOMIQUE
# =========================================================================

# Correction taxonomique des noms d'espèces
taxo = correctTaxo(genus = CoFor_raw$genus, species = CoFor_raw$species, useCache = F, verbose = F) %>% as_tibble()

# Application des corrections taxonomiques à la base CoforTraits
CoFor = CoFor_raw %>%
  mutate(genus = taxo$genusCorrected) %>% 
  mutate(species = taxo$speciesCorrected) %>%
  mutate(species_full = paste(genus, species))

# Identification des espèces présentes/absentes dans CoforTraits
field_inventories = field_inventories %>% 
  mutate(sp_in_CoFor = species_full %in% CoFor$species_full) 

# Liste des espèces absentes de CoforTraits
liste_sp_out_CoFor = field_inventories %>%
  filter(!sp_in_CoFor) %>%
  dplyr::select(species_full) %>% 
  unique()

# Liste des espèces présentes dans CoforTraits
liste_sp_in_CoFor = field_inventories %>%
  filter(sp_in_CoFor) %>%
  dplyr::select(species_full) %>% 
  unique()

# =========================================================================
# CALCUL DES TEMPÉRAMENTS 
# =========================================================================

# Calcul des tempéraments pour tous les niveaux taxonomiques en une fois

temperaments_3cat <- get_all_temperaments(CoFor, subdivision = 3)

CoFor.Temp <- temperaments_3cat$species
CoFor.Temp_genus_consensus <- temperaments_3cat$genus
CoFor.Temp_family_consensus <- temperaments_3cat$family

# =========================================================================
# FILTRAGE ET FUSION DES DONNÉES
# =========================================================================

# Filtrage des tempéraments pour les espèces présentes dans l'inventaire

# (3 subdivisions)
field_temps_3cat <- list(
  species = CoFor.Temp %>% filter(species_full %in% field_inventories$species_full),
  genus = CoFor.Temp_genus_consensus %>% filter(Genus %in% field_inventories$genus),
  family = CoFor.Temp_family_consensus %>% filter(Family %in% field_inventories$fam)
)

# Attribution hiérarchique des tempéraments aux données d'inventaire

field_inventories <- assign_hierarchical_temperament(
  field_inventories,
  field_temps_3cat$species,
  field_temps_3cat$genus,
  field_temps_3cat$family,
)

# 3. Remplacement des NA par "Inconnu"
field_inventories <- field_inventories %>%
  mutate(
    Temperament = replace_na(Temperament, "Inconnu"),
    level_Temperament = replace_na(level_Temperament, "Inconnu"),
    Temperament.3subdiv = replace_na(Temperament, "Inconnu"),
    level_Temperament.3subdiv = replace_na(level_Temperament, "Inconnu")
  )

# =========================================================================
# EXPORT DES RÉSULTATS
# =========================================================================

rio::export(
  field_inventories,
  paste0(path0, "field_data_inventories_with_temperament.csv"),
  sep = ";",
  dec = ","
)


