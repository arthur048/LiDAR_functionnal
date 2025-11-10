# 2. Données d'Élévation

## Vue d'ensemble

Ce module génère les Modèles Numériques de Hauteur de Canopée (CHM - Canopy Height Model) à partir des nuages de points LiDAR. Les CHM représentent la hauteur de la végétation au-dessus du sol et sont essentiels pour l'identification des trouées forestières.

## Objectif

Transformer les nuages de points LiDAR 3D en rasters 2D de hauteur de canopée, normalisés et prêts pour l'analyse des trouées.

## Données d'entrée

### Nuages de points LiDAR
- `../1_LiDAR_extraction/output/{plot_ref}/{plot_ref}.las` - Fichiers LAS découpés (module 1)
- Points classifiés avec sol et végétation

### Données spatiales
- `plots_info_NASA2015.csv` - Métadonnées des parcelles
- `DRC01_lidar_transects.shp` / `DRC02_ferry_lines.shp` - Masques de référence

## Scripts principaux

### `ElevationData.R`

**Fonction**: Génère les CHM, DTM et normalise les nuages de points

**Workflow**:
1. **Lecture du nuage de points** LAS
2. **Filtrage du bruit** (points aberrants au-dessus du 95e percentile)
3. **Génération du DTM** (Digital Terrain Model - modèle numérique de terrain):
   - Algorithme: `tin()` - interpolation par triangulation
   - Résolution: adaptative selon les données
4. **Normalisation du nuage de points**:
   - Soustraction du DTM pour obtenir les hauteurs au-dessus du sol
5. **Génération du CHM**:
   - Algorithme: `p2r()` avec sous-cercle
   - Résolution: 1m (par défaut)
6. **Interpolation des valeurs manquantes** (pixels sans retour LiDAR)
7. **Masquage**:
   - Application du masque du transect/ferry line
   - Application du masque de la parcelle
8. **Export des rasters**

**Fonctions utilitaires**:
```r
filter_noise(las, sensitivity)   # Filtre les points aberrants
interpolate_na(rast)              # Interpole les pixels manquants
ensure_correct_crs(obj, epsg)    # Vérifie/corrige le CRS
get_reference_mask(...)           # Obtient le masque de référence
process_raster(raster, masks)     # Applique les masques
```

### `CHM_final_selection.R`

**Fonction**: Sélectionne les meilleurs CHM parmi plusieurs versions

**Workflow**:
1. **Lecture de** `selection_final_CHM.xlsx`
   - Colonne `keep`: 1 ou 2 pour les CHM à conserver
   - Colonne `origin_CHM`: dossier source du CHM
2. **Copie des CHM sélectionnés** vers `CHM_final/`
3. **Standardisation** des noms de fichiers

**Usage**: Permet de choisir entre plusieurs CHM générés avec différents paramètres.

### Scripts alternatifs

- `ElevationData_from-existing-CHM.R` - Réutilise des CHM existants
- `ElevationData_from-existing-CHM_AfriSAR.R` - Pour données AfriSAR
- `ElevationData_from-existing-CHM_NASA2015.R` - Pour données NASA 2015

## Configuration

### Parallélisation
```r
n_cores <- parallel::detectCores() / 2  # Utilise 50% des cœurs
cl <- makeCluster(floor(n_cores))
registerDoParallel(cl)
```

### Paramètres de génération

**DTM (Digital Terrain Model)**:
- Algorithme: Triangulation (TIN)
- Classification des points sol utilisée

**CHM (Canopy Height Model)**:
- Résolution: 1m × 1m
- Algorithme: Point-to-raster avec sous-cercle
- Filtrage du bruit: Valeurs > 95e percentile × sensibilité

## Structure de sortie

```
2_ElevationData/
├── CHM/               # CHM générés (intermédiaires)
├── DTM/               # Modèles numériques de terrain
├── LAS_normalized/    # Nuages de points normalisés
├── CHM_final/         # CHM finaux sélectionnés ⭐
└── selection_final_CHM.xlsx  # Table de sélection
```

## Fichiers de sortie

### `CHM_final/{plot_ref}.tif`
**Rasters de hauteur de canopée finaux**:
- **Résolution**: 1m
- **Valeurs**: Hauteur de végétation au-dessus du sol (m)
- **Format**: GeoTIFF
- **CRS**: UTM approprié à la zone
- **Utilisé par**: Module 3 (analyse des trouées)

### `DTM/{plot_ref}.tif`
**Modèles numériques de terrain**:
- Élévation du sol
- Interpolé à partir des points classifiés "sol"

### `LAS_normalized/{plot_ref}_normalized.las`
**Nuages de points normalisés**:
- Coordonnée Z = hauteur au-dessus du sol
- Utilisé pour validation et analyses complémentaires

## Traitement des données

### Gestion des valeurs manquantes
Les pixels sans retour LiDAR (sous couvert très dense) sont interpolés par:
1. **Focal mean** avec fenêtre 3×3
2. **Itération** jusqu'à ce qu'il ne reste plus de NA
3. **Préservation** des bordures avec masquage

### Filtrage du bruit
Points considérés comme bruit si:
- Z > p95 × sensibilité (défaut: 1.2)
- Retours isolés en haute altitude

### Double masquage
1. **Masque du transect/ferry**: Zone de couverture LiDAR
2. **Masque de la parcelle**: Limite exacte de l'inventaire
   → Garantit que seules les données fiables dans la parcelle sont conservées

## Dépendances R

```r
library(tidyverse)    # Manipulation de données
library(sf)           # Données spatiales
library(lidR)         # Traitement LiDAR
library(terra)        # Manipulation de rasters
library(foreach)      # Parallélisation
library(doParallel)   # Backend parallèle
```

## Utilisation

### Génération complète
```r
source("2_ElevationData/ElevationData.R")
```

### Sélection des CHM finaux
```r
# 1. Éditer selection_final_CHM.xlsx
#    - Marquer keep=1 ou keep=2 pour les CHM à conserver
#    - Indiquer l'origine (CHM, CHM_AfriSAR, etc.)

# 2. Exécuter la sélection
source("2_ElevationData/CHM_final_selection.R")
```

## Validation des données

### Vérifications essentielles
```r
library(terra)

# Charger un CHM
chm <- rast("2_ElevationData/CHM_final/{plot_ref}.tif")

# Statistiques de base
summary(chm)
global(chm, "mean", na.rm=TRUE)
global(chm, "max", na.rm=TRUE)

# Visualisation
plot(chm, main="CHM - {plot_ref}")
hist(chm, breaks=50, main="Distribution des hauteurs")

# Vérifier la résolution
res(chm)  # Doit être [1, 1]

# Vérifier le CRS
crs(chm)  # Doit être en UTM
```

### Valeurs typiques attendues
- **Hauteur moyenne**: 20-35 m (forêt mature)
- **Hauteur maximale**: 40-60 m (arbres émergents)
- **Proportion de NA**: < 5% après interpolation
- **Distribution**: Asymétrique vers les valeurs élevées

## Problèmes courants

### CHM avec trop de valeurs manquantes
**Cause**: Faible densité de points LiDAR
**Solution**: Vérifier la qualité du LAS d'origine, ajuster les paramètres d'interpolation

### Hauteurs aberrantes (> 80m)
**Cause**: Points de bruit non filtrés
**Solution**: Ajuster le paramètre `sensitivity` de `filter_noise()`

### Artéfacts en damier
**Cause**: Interpolation sur une zone trop large sans données
**Solution**: Resserrer le masquage ou améliorer le filtrage des points sol

### DTM mal généré
**Cause**: Classification du sol insuffisante
**Solution**: Revoir l'algorithme de classification ou utiliser un DTM externe
