# 3. Analyse des Trouées Forestières

## Vue d'ensemble

Ce module identifie et caractérise les trouées dans la canopée forestière à partir des CHM (Canopy Height Models). Les trouées sont définies comme des zones où la hauteur de végétation est inférieure à un seuil donné. Ce module est au cœur de l'analyse car les trouées reflètent la structure verticale des forêts.

## Objectif

Détecter automatiquement les trouées forestières à multiple seuils de hauteur (1-45m ou 1-60m), calculer leurs propriétés géométriques, et générer des profils verticaux qui caractérisent la structure forestière.

## Données d'entrée

### Rasters de hauteur
- `../2_ElevationData/CHM_final/{plot_ref}.tif` - Canopy Height Models

### Données spatiales
- `../0_Inventories_plot_preparation/final/plots_unique/{plot_ref}.gpkg` - Limites des parcelles
- `plots_info.csv` - Métadonnées des parcelles

## Scripts principaux

### `1_Gaps_identification.R`

**Fonction**: Identifie les trouées à multiple seuils de hauteur

**Algorithme principal - `get_gap_layer()`**:
```r
get_gap_layer(chm_layer, threshold, size)
```

**Étapes**:
1. **Reclassification binaire**: Pixels < `threshold` = trouée (1), sinon NA
2. **Détection des patches connectés**: Algorithme à 8 directions de connectivité
3. **Calcul de la surface** de chaque patch
4. **Filtrage par taille**: Conserve uniquement les trouées entre `size[1]` et `size[2]` m²
5. **Attribution d'ID unique** à chaque trouée

**Paramètres**:
- **Seuils de hauteur**: 1-45 m (ou 1-60 m selon les analyses)
- **Taille minimale**: 1 m² (résolution du raster)
- **Taille maximale**: 10⁸ m² (pas de limite supérieure pratique)

**Parallélisation**: Traitement simultané de toutes les parcelles

**Sortie**:
- `output/{plot_ref}/{plot_ref}_{height}m_gaps_height.tif`
- Un raster par seuil de hauteur (45 rasters par parcelle)

### `2_Gaps_data_frequencies.R`

**Fonction**: Calcule la distribution de taille des trouées

**Méthode - `get_gap_size_frequency_distribution()`**:
```r
get_gap_size_frequency_distribution(gaps_layer, type_of_return)
```

**Calcule**:
- **Fréquence** de chaque classe de taille de trouée
- **Distribution cumulative** des tailles
- **Statistiques** par parcelle, hauteur et buffer

**Buffers testés**:
- 0 m (parcelle exacte)
- 20 m (zone tampon)
- 50 m (zone tampon étendue)

**Sortie**:
- `output/gaps_size_frequency_data.csv`
- Distribution de taille pour chaque combinaison plot × hauteur × buffer

### `3_Gaps_metrics.R`

**Fonction**: Calcule la proportion de pixels en trouée

**Méthode - `calculate_gap_proportion()`**:
```r
calculate_gap_proportion(chm, height_threshold)
```

**Calcule**:
- **Proportion de trouées** = nombre de pixels < seuil / nombre total de pixels
- Pour chaque combinaison de:
  - Parcelle (n plots)
  - Hauteur (1-60 m)
  - Buffer (0, 20, 50, 100, 200 m)

**Sortie** principale:
- `output/results_gaps_metrics.csv`
  - Colonnes: `plot_name`, `height_aboveground`, `buffer`, `proportion`
  - Utilisé pour construire les **profils verticaux de trouées**

### `4_Gaps_graphics.R`

**Fonction**: Génère des visualisations des trouées

**Types de graphiques**:
- Distribution de taille des trouées
- Histogrammes et courbes cumulées
- Cartes de localisation des trouées
- Comparaisons entre parcelles

### `5_Gaps_advanced_analysis.R` (versions 1 & 2)

**Fonction**: Analyses avancées et métriques dérivées

**Calculs**:
- Métriques géométriques des trouées (périmètre, compacité)
- Indices de fragmentation
- Analyse spatiale des patterns de trouées

## Fichiers de sortie clés

### `output/results_gaps_metrics.csv` ⭐
**Profils verticaux de proportion de trouées**:
```
plot_name, height_aboveground, buffer, proportion
MAB01h, 1, 0, 0.05
MAB01h, 2, 0, 0.08
...
MAB01h, 60, 0, 0.95
```

**Usage**: Entrée principale pour le module 5 (analyses statistiques)

### `output/{plot_ref}/{plot_ref}_{height}m_gaps_height.tif`
**Rasters de trouées** par seuil:
- Valeur: ID de la trouée
- NA: Non-trouée ou hors parcelle
- Usage: Visualisation, analyses géométriques

### `output/gaps_size_frequency_data.csv`
**Distribution de taille des trouées**:
- Colonnes: `plot_name`, `height_aboveground`, `buffer`, `gap_area`, `gap_freq`
- Usage: Analyses de la dynamique forestière (théorie des trouées)

## Concept clé: Profil vertical

Le **profil vertical de trouées** est la courbe représentant la proportion de trouées en fonction de la hauteur:

```
Proportion = f(Hauteur)

Exemple conceptuel:
- À 1m: 5% de trouées (sous-bois ouvert)
- À 10m: 20% de trouées (strate intermédiaire)
- À 30m: 80% de trouées (canopée fermée)
- À 50m: 95% de trouées (rares émergents)
```

Ce profil caractérise la **stratification verticale** de la forêt.

## Paramètres importants

### Connectivité des trouées
```r
terra::patches(gaps, directions = 8)
```
- **8 directions**: Inclut les diagonales (connexion par coin)
- Définit ce qui constitue une "trouée unique"

### Résolution spatiale
- **1m × 1m**: Résolution des rasters CHM et gaps
- Compromise entre détail et temps de calcul

### Gamme de hauteurs
- **1-45 m**: Analyse standard
- **1-60 m**: Analyse étendue pour forêts hautes
- **Pas de 1 m**: Résolution fine de la stratification

## Dépendances R

```r
library(terra)        # Manipulation de rasters
library(tidyverse)    # Manipulation de données
library(sf)           # Données spatiales
library(rio)          # Import/export
library(foreach)      # Parallélisation
library(doParallel)   # Backend parallèle
library(progress)     # Barres de progression (script 3)
```

## Utilisation

### Workflow complet
```r
# 1. Identification des trouées (long! ~2-4h pour 50 plots)
source("3_Gaps/1_Gaps_identification.R")

# 2. Distribution de taille
source("3_Gaps/2_Gaps_data_frequencies.R")

# 3. Calcul des proportions (utilise script 3)
source("3_Gaps/3_Gaps_metrics.R")

# 4. Visualisations
source("3_Gaps/4_Gaps_graphics.R")
```

### Exécution par étapes
Il est recommandé d'exécuter les scripts séquentiellement car:
- Script 1 est très intensif (générer 45-60 rasters × n_plots)
- Scripts 2-3 dépendent des sorties du script 1

## Validation des données

### Vérifier les rasters de trouées
```r
library(terra)

# Charger un raster de trouées
gaps_10m <- rast("3_Gaps/output/MAB01h/MAB01h_10m_gaps_height.tif")

# Statistiques
freq(gaps_10m)  # Fréquence de chaque ID de trouée
plot(gaps_10m, main="Trouées à 10m - MAB01h")
```

### Vérifier les profils
```r
library(tidyverse)

# Charger les profils
profiles <- read_csv2("3_Gaps/output/results_gaps_metrics.csv")

# Plot pour une parcelle (buffer = 0)
profiles %>%
  filter(plot_name == "MAB01h", buffer == 0) %>%
  ggplot(aes(x = height_aboveground, y = proportion)) +
  geom_line() +
  labs(title = "Profil vertical - MAB01h",
       x = "Hauteur (m)", y = "Proportion de trouées")
```

### Valeurs attendues
- **Proportion à 1m**: 0-20% (sous-bois généralement fermé)
- **Proportion à 50m**: 80-99% (rares arbres dépassent 50m)
- **Forme de courbe**: Sigmoïde croissante
- **Point d'inflexion**: Souvent autour de 20-30m (transition sous-bois/canopée)

## Interprétation écologique

### Profil "abrupt" (pente forte)
- **Canopée homogène** et continue
- Strate dominante bien définie
- Typique des forêts matures à faible perturbation

### Profil "progressif" (pente douce)
- **Stratification complexe**
- Multiple strates de végétation
- Typique des forêts en régénération ou mixtes

### Position du point d'inflexion
- **Bas (< 15m)**: Canopée basse, forêt secondaire jeune
- **Haut (> 30m)**: Canopée haute, forêt mature

## Performance

### Temps de calcul (indicatif)
| Script | n=50 plots | Parallélisation |
|--------|------------|-----------------|
| Script 1 | 2-4h | Oui (75% cœurs) |
| Script 2 | 30-60 min | Non |
| Script 3 | 1-2h | Oui (75% cœurs) |

### Espace disque
- **Par parcelle**: ~200-500 MB (45-60 rasters)
- **50 parcelles**: ~10-25 GB
- Les rasters sont compressés (format GeoTIFF)

## Problèmes courants

### Script 1 très lent
**Solution**: C'est normal! Vérifier la parallélisation est active

### Profil plat ou aberrant
**Cause**: CHM de mauvaise qualité ou parcelle trop petite
**Solution**: Vérifier le CHM d'origine, exclure la parcelle si nécessaire

### Trouées fragmentées
**Cause**: Bruit dans le CHM ou résolution trop fine
**Solution**: Appliquer un filtre de lissage au CHM avant analyse

### Mémoire insuffisante
**Solution**: Réduire le nombre de cœurs dans la parallélisation
