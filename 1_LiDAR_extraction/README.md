# 1. Extraction LiDAR

## Vue d'ensemble

Ce module extrait les nuages de points LiDAR pour chaque parcelle d'inventaire à partir des fichiers bruts de la campagne NASA 2015 en RDC. Il découpe les grands fichiers LAS selon les limites des parcelles et les prépare pour l'étape suivante de traitement.

## Objectif

Créer des fichiers LAS individuels pour chaque parcelle, découpés précisément selon leurs limites spatiales, à partir des données LiDAR brutes organisées en transects et lignes de ferry.

## Données d'entrée

### Données LiDAR brutes
**Localisation**: `E:/Arthur/Doctorat_DataAnnexe/LiDAR_RDC_2015/00-raw-data/`

**Structure**:
- `Plots/{UTM}/{transect_name}/Las/*.las` - Données des transects
- `Ferry Data/{UTM}/{ferry_folder}/*.las` - Données des lignes de ferry

### Données spatiales
- `plots_limites_NASA2015.gpkg` - Limites des parcelles (du module 0)
- `plots_info_NASA2015.csv` - Métadonnées des parcelles (du module 0)
- `DRC01_lidar_transects.shp` - Polygones des transects
- `DRC02_ferry_lines.shp` - Polygones des lignes de ferry

## Scripts principaux

### `LiDAR_extraction.R`

**Fonction**: Extrait et découpe les fichiers LAS pour chaque parcelle

**Workflow**:
1. **Lecture des métadonnées** des parcelles
2. **Identification de la source LiDAR**:
   - Transects (plots le long de lignes droites)
   - Ferry lines (plots le long de trajectoires de bateau)
3. **Localisation des fichiers LAS** correspondants
4. **Découpage spatial** selon les limites de parcelle
5. **Export** des fichiers LAS découpés

**Fonctions utilitaires**:
```r
extract_utm_plot(name)           # Extrait UTM et nom du plot
format_ferry_folder(plot_name)   # Formate le nom du dossier ferry
get_lidar_path(source, utm, folder_name)  # Construit le chemin LiDAR
get_epsg_from_utm(utm_zone)      # Convertit zone UTM en code EPSG
```

## Configuration

### Parallélisation
Le script utilise la moitié des cœurs disponibles:
```r
nb.coeurs <- parallel::detectCores()
set_lidr_threads(nb.coeurs / 2)
```

### Zones UTM
```r
utm_epsg <- c(
  "UTM33S" = 32733,
  "UTM34N" = 32634,
  "UTM34S" = 32734,
  "UTM35N" = 32635,
  "UTM35S" = 32735,
  "UTM36N" = 32636
)
```

## Fichiers de sortie

### `output/{plot_ref}/`
Dossier créé pour chaque parcelle contenant:
- `{plot_ref}.las` - Fichier LAS découpé pour la parcelle
- Métadonnées préservées (CRS, attributs des points)

**Caractéristiques des fichiers LAS**:
- Points classifiés (sol, végétation, etc.)
- Coordonnées en UTM (selon la zone)
- Attributs: X, Y, Z, Intensity, ReturnNumber, Classification

## Traitement des données

### Découpage spatial
Le script effectue un découpage précis basé sur:
1. **Intersection** entre la parcelle et le transect/ferry line
2. **Clippage** du nuage de points selon les limites exactes de la parcelle
3. **Vérification** que les fichiers de sortie contiennent des points

### Gestion des erreurs
- Vérifie l'existence des dossiers sources
- Teste la présence de fichiers LAS
- Gère les CRS manquants ou incorrects
- Enregistre les erreurs pour analyse

## Dépendances R

```r
library(tidyverse)  # Manipulation de données
library(sf)         # Données spatiales
library(lidR)       # Traitement LiDAR
```

## Utilisation

```r
# Depuis la racine du projet
source("1_LiDAR_extraction/LiDAR_extraction.R")
```

**⚠️ Attention**:
- Nécessite ~50-100 Go d'espace disque temporaire
- Temps de traitement: 1-3h selon le nombre de parcelles
- Les chemins vers les données brutes doivent être ajustés selon votre configuration

## Validation des données

Après l'extraction, vérifier:
```r
# Liste des fichiers générés
output_files <- list.files("1_LiDAR_extraction/output",
                           pattern = "\\.las$",
                           recursive = TRUE)

# Nombre de parcelles traitées
length(unique(dirname(output_files)))

# Vérifier qu'un fichier contient des points
library(lidR)
test_las <- readLAS("1_LiDAR_extraction/output/{plot_ref}/{plot_ref}.las")
summary(test_las)
```

## Problèmes courants

### Fichier LAS vide
**Cause**: La parcelle ne se superpose pas avec les données LiDAR brutes
**Solution**: Vérifier l'alignement spatial des shapefiles et des transects

### Erreur de CRS
**Cause**: Projection non définie ou incorrecte
**Solution**: S'assurer que tous les shapefiles sont en UTM approprié

### Mémoire insuffisante
**Cause**: Fichiers LAS trop volumineux
**Solution**: Réduire `set_lidr_threads()` ou traiter par lots
