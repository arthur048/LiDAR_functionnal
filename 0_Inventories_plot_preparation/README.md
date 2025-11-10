# 0. Préparation des Inventaires de Parcelles

## Vue d'ensemble

Ce module prépare les données d'inventaires forestiers et les limites spatiales des parcelles pour l'ensemble du workflow d'analyse. Il combine plusieurs sources de données provenant de campagnes AfriSAR et de la mission NASA 2015 en Afrique centrale.

## Objectif

Consolider les données d'inventaires de terrain provenant de multiples sources (Rabi, Mondah, Mabounie, Lope, Yangambi, Malebo) et créer des fichiers géospatiaux standardisés pour toutes les parcelles à analyser.

## Données d'entrée

### Shapefiles des campagnes
- `shapefile/Rabi_AfriSAR_1ha.shp`
- `shapefile/Mondah_AfriSAR_1ha.shp`
- `shapefile/Mabounie_AfriSAR_1ha.shp`
- `shapefile/Lope_AfriSAR_1ha.shp`
- `shapefile/plots_under_lidar_without_ituri_final.shp`

### Données d'inventaires de terrain
- `field_data/NASA_JPL/db_Africa_NASA.xlsx` - Inventaires NASA
- `field_data/RDC_others/MaleboFieldPlots.csv` - Parcelles Malebo
- `field_data/RDC_others/Yangambi_Nestor_FieldPlots.csv` - Parcelles Yangambi (Nestor)
- `field_data/RDC_others/YangambiFieldPlots.xlsx` - Parcelles Yangambi
- `field_data/PlotsMabounie.RData` - Données Mabounie

### Zones LiDAR
- `shapefile/lidar/DRC01_lidar_transects.shp` - Transects LiDAR
- `shapefile/lidar/DRC02_ferry_lines.shp` - Lignes de ferry

## Scripts principaux

### `plots_preparation.R`

**Fonction**: Script principal qui consolide toutes les données

**Workflow**:
1. **Lecture des shapefiles** pour chaque site
2. **Import des inventaires de terrain** avec standardisation des colonnes
3. **Filtrage des données**:
   - DBH ≥ 10 cm
   - Hauteur ≤ 75 m (si disponible)
4. **Attribution des zones UTM** selon la localisation géographique
5. **Harmonisation des références de parcelles** (codes de parcelle)
6. **Création des fichiers de sortie** standardisés

**Zones UTM supportées**:
- UTM 32N/S (EPSG: 32632/32732)
- UTM 33S (EPSG: 32733)
- UTM 34N/S (EPSG: 32634/32734)
- UTM 35N/S (EPSG: 32635/32735)
- UTM 36N (EPSG: 32636)

## Fichiers de sortie

### `final/plots_inventories.csv`
Inventaire consolidé de tous les arbres (DBH ≥ 10 cm):
- **Colonnes**: `plot_ref`, `site`, `sp` (nom d'espèce), `dbh` (cm), `h` (hauteur, m)
- **Format**: CSV avec `;` comme délimiteur, `,` comme séparateur décimal
- **Usage**: Entrée principale pour le calcul des métriques de parcelles

### `final/plots_info.csv` / `plots_info_NASA2015.csv`
Métadonnées géographiques des parcelles:
- **Colonnes**: `plot_ref`, `source` (transects/ferry), `utm`, `transect_name`, `folder_name`
- **Usage**: Localisation des données LiDAR brutes et informations de projection

### `final/plots_limites.gpkg` / `plots_limites_NASA2015.gpkg`
Limites spatiales des parcelles:
- **Format**: GeoPackage (GPKG)
- **Géométrie**: Polygones des parcelles de 1 ha
- **Usage**: Masquage et découpage des données raster

### `final/plots_unique/`
Dossier contenant un fichier GPKG par parcelle:
- **Nomenclature**: `{plot_ref}.gpkg`
- **Usage**: Traitement individuel des parcelles

### `final/plots_AfriSAR_WD.csv`
Densités du bois pour les parcelles AfriSAR:
- **Source**: Données préexistantes de densité du bois
- **Usage**: Comparaison et validation

## Notes importantes

### Harmonisation des noms de parcelles
Les parcelles Mabounie utilisent une nomenclature spécifique qui est convertie:
- `mabou001` → `MAB01h`
- `mabou002` → `MAB02h`
- etc.

### Gestion des projections
Chaque site a sa propre zone UTM. Le script attribue automatiquement l'EPSG correct basé sur la localisation géographique.

### Fichier maître
`final_plot_name.csv` (racine du projet) contient la liste finale des parcelles retenues pour toutes les analyses.

## Dépendances R

```r
library(tidyverse)  # Manipulation de données
library(sf)         # Données spatiales
library(readxl)     # Lecture fichiers Excel
library(rio)        # Import/export flexible
```

## Utilisation

```r
# Depuis la racine du projet
source("0_Inventories_plot_preparation/plots_preparation.R")
```

**⚠️ Attention**: Les chemins sont codés en dur vers `E:/Arthur/OneDrive2/`. Ajuster selon votre configuration.

## Validation des données

Vérifier que:
- Tous les plots dans `final_plot_name.csv` ont des données dans `plots_inventories.csv`
- Les limites spatiales correspondent aux inventaires
- Les zones UTM sont correctement assignées
- Pas de valeurs manquantes critiques (plot_ref, site, dbh)
