# 4. Métriques de Parcelles

## Vue d'ensemble

Ce module calcule toutes les métriques forestières à l'échelle des parcelles: composition taxonomique, densité du bois, structure (surface terrière, biomasse), et tempérament écologique des espèces. Ces métriques seront corrélées avec la structure verticale dans le module 5.

## Objectif

Agréger les données d'inventaires individuels (arbres) en métriques synthétiques de parcelles et intégrer les traits fonctionnels des espèces (densité du bois, tempérament).

## Données d'entrée

### Inventaires de terrain
- `../0_Inventories_plot_preparation/final/plots_inventories.csv` - Données arbres (DBH, hauteur)
- `final_plot_name.csv` (racine) - Liste des parcelles retenues

### Bases de données de traits
- `../cofortraits.csv` (racine) - **Base de données locale** de densité du bois et traits fonctionnels
- **BIOMASS package databases** (automatique):
  - Global Wood Density Database (GWDD)
  - Base taxonomique (APG IV)

### Métriques des trouées
- `../3_Gaps/output/results_gaps_metrics.csv` - Profils verticaux

## Scripts et workflow

### `1_FieldDataInventories_computation.R`

**Fonction**: Traite les inventaires à l'échelle de l'arbre

**Pipeline de traitement**:

#### 1. Correction taxonomique
```r
correctTaxo(genus, species, useCache = TRUE)
```
- Standardise les noms d'espèces selon APG IV
- Corrige les synonymes et fautes de frappe
- Identifie les noms valides vs invalides

#### 2. Attribution de la densité du bois (WD)
**Hiérarchie des sources** (par priorité):
1. **CoForTraits** (base locale - données terrain)
2. **GWDD espèce** (Global Wood Density Database)
3. **GWDD genre** (moyenne du genre si espèce absente)
4. **Valeur par défaut** (0.55 g/cm³ si aucune donnée)

**Pourquoi prioriser CoForTraits?**
- Données mesurées localement en Afrique centrale
- Plus représentatives des conditions édaphiques et climatiques
- Réduction de l'incertitude

#### 3. Estimation de la hauteur
Pour les arbres sans mesure de hauteur:
```r
retrieveH(D = dbh, coord = coordinates, model = model)
```
- **Modèle local** si ≥ 50 arbres avec hauteur dans la parcelle
- **Modèle régional** (pantropical) sinon
- Allométries de Chave et al. (2014)

#### 4. Calcul des métriques individuelles
- **Surface terrière** (G): π × (DBH/2)²
- **Biomasse aérienne** (AGB): Allométrie de Chave et al. (2014)
  ```r
  computeAGB(D = dbh, WD = wood_density, H = height)
  ```

**Sortie**:
- `FieldDataInventories_computation.csv`
- Une ligne par arbre avec: taxonomie corrigée, WD, hauteur (mesurée/estimée), G, AGB

### `2_PlotsMetrics_FieldDataInventries.R`

**Fonction**: Agrège les métriques à l'échelle parcelle

**Métriques calculées**:
- **Densité de tiges**: Nombre d'arbres/ha
- **Surface terrière** (G): m²/ha
- **Biomasse aérienne** (AGB): Mg/ha
- **Densité du bois moyenne**: WD_mean
- **Densité du bois pondérée par G**: WD_BA (⭐ **variable cible principale**)
  ```r
  WD_BA = sum(WD × G) / sum(G)
  ```
- **Hauteur moyenne**, max, quantiles
- **Diamètre moyen**, max, quantiles
- **Richesse spécifique**

**Sortie**:
- `PlotsMetrics_FieldDataInventories.csv`
- Une ligne par parcelle

### `3_PlotsMetrics_Gaps.R`

**Fonction**: Agrège les métriques de trouées par parcelle

**À partir des profils verticaux**, calcule:
- **Proportion de trouées** à hauteurs fixes (10m, 15m, 20m, 30m)
- **Hauteurs caractéristiques** (h10, h25, h50, h75, h90):
  - h50 = hauteur où 50% de la parcelle est en trouée
- **Pente maximale** du profil (Smax) et hauteur associée (hSmax)
- **Aire sous la courbe** (AUC): Intégration numérique
- **Paramètres logistiques**: Ajustement d'une courbe sigmoïde
  ```r
  P(h) = a + (b - a) / (1 + exp(-k × (h - h0)))
  ```
  - a: asymptote inférieure
  - b: asymptote supérieure
  - h0: point d'inflexion
  - k: pente
- **Métriques dérivées**: IDV, ROS, IAS, ICT, etc.
- **Métriques zonales**: Découpe en 3 zones (basale, transition, supérieure)
- **Métriques différentielles**: Dérivées première et seconde
- **Métriques multi-échelles**: Rugosité à différentes échelles

**Sortie**:
- `plot_metrics_gaps.csv`
- Une ligne par parcelle avec ~50-100 métriques

### `4_FieldDataInventories_Analyse_Temperament.R`

**Fonction**: Attribue le tempérament écologique à chaque arbre

**Classification des espèces** (selon la littérature):
1. **Héliophile** (pioneer): Espèces pionnières strictes
   - Demandent pleine lumière
   - Régénération dans les grandes trouées
   - Ex: Musanga cecropioides, Macaranga spp.

2. **Héliophile non pionnière** (NPLD - Non-Pioneer Light Demander):
   - Tolèrent un peu d'ombre en jeunesse
   - Adultes nécessitent lumière
   - Ex: Terminalia superba, Triplochiton scleroxylon

3. **Tolérante à l'ombrage** (shade-tolerant):
   - Régénèrent sous couvert fermé
   - Croissance lente, bois dense
   - Ex: Staudtia kamerunensis, Guarea thompsonii

**Source**: Base de données d'experts + littérature (Hawthorne 1995, Poorter et al. 2006)

**Sortie**:
- `field_data_inventories_with_temperament.csv`
- Chaque arbre avec son tempérament

### `5_PlotsMetrics_Analyse_Temperament.R`

**Fonction**: Calcule les proportions de tempérament par parcelle

**Métriques calculées**:

#### Par surface terrière (G.rel):
- `prop_g_helio`: Proportion de G des héliophiles
- `prop_g_npld`: Proportion de G des NPLD
- `prop_g_shade`: Proportion de G des tolérants

#### Par nombre d'individus (ind.rel):
- `prop_ind_helio`: Proportion d'individus héliophiles
- `prop_ind_npld`: Proportion d'individus NPLD
- `prop_ind_shade`: Proportion d'individus tolérants

#### Par classes de diamètre:
- Proportions dans différentes classes de DBH
- Permet d'analyser la structure démographique

**Sortie**:
- `PlotsMetrics_Temperament_all.csv`
- Multiples lignes par parcelle (différentes subdivisions)

## Fichiers de sortie clés

### `PlotsMetrics_FieldDataInventories.csv` ⭐
**Métriques de composition et structure**:
- Variables d'intérêt: **WD_BA**, WD_mean, G, AGB, densité
- Usage: Variables cibles pour les corrélations

### `plot_metrics_gaps.csv` ⭐
**Métriques de structure verticale**:
- ~50-100 métriques dérivées des profils de trouées
- Usage: Variables explicatives pour les corrélations

### `PlotsMetrics_Temperament_all.csv` ⭐
**Proportions de tempérament**:
- Variables d'intérêt: **prop_g_helio**, prop_g_npld, **prop_g_shade**
- Usage: Variables cibles alternatives (liées à WD)

### `meanWD_by_species_CoFor.csv` (racine)
**Densités du bois moyennes par espèce** (CoForTraits):
- Utilisé comme référence de qualité

## Concepts clés

### Densité du bois pondérée (WD_BA)
**Variable cible principale** de l'étude:
```r
WD_BA = sum(WD_i × G_i) / sum(G_i)
```
- Donne plus de poids aux arbres dominants (gros diamètre)
- Meilleur indicateur de la biomasse et du carbone
- Corrélé avec le tempérament écologique

### Lien WD ↔ Tempérament
**Gradient écologique attendu**:
- **Héliophiles**: WD faible (~0.3-0.5)
  - Croissance rapide, bois léger
- **NPLD**: WD intermédiaire (~0.5-0.6)
- **Tolérants**: WD élevé (~0.6-0.9)
  - Croissance lente, bois dense

### Importance de la correction taxonomique
- **~15-30%** des noms d'espèces nécessitent correction
- Erreurs de frappe, synonymes, changements taxonomiques
- Impact direct sur l'attribution du WD

## Dépendances R

```r
library(tidyverse)    # Manipulation de données
library(BIOMASS)      # Correction taxonomique, allométries
library(rio)          # Import/export
library(sf)           # Données spatiales (pour coordonnées)
library(terra)        # Rasters (scripts gaps)
library(lidR)         # LiDAR (scripts gaps)
library(foreach)      # Parallélisation
library(doParallel)   # Backend parallèle
```

## Utilisation

### Workflow séquentiel complet
```r
# 1. Traitement des arbres individuels
source("4_Plots_metrics/1_FieldDataInventories_computation.R")

# 2. Agrégation des métriques de parcelles (inventaires)
source("4_Plots_metrics/2_PlotsMetrics_FieldDataInventries.R")

# 3. Agrégation des métriques de trouées
source("4_Plots_metrics/3_PlotsMetrics_Gaps.R")

# 4. Attribution du tempérament
source("4_Plots_metrics/4_FieldDataInventories_Analyse_Temperament.R")

# 5. Métriques de tempérament par parcelle
source("4_Plots_metrics/5_PlotsMetrics_Analyse_Temperament.R")
```

**⚠️ Ordre important**: Scripts 4-5 dépendent du script 1

## Validation des données

### Vérifier les densités du bois
```r
library(tidyverse)

# Charger les données
field_data <- read_csv2("4_Plots_metrics/FieldDataInventories_computation.csv")

# Distribution des WD
summary(field_data$WD)
hist(field_data$WD, breaks=50, main="Distribution des densités du bois")

# Proportion avec CoForTraits vs GWDD
table(field_data$WD_source)

# Vérifier les valeurs aberrantes
field_data %>% filter(WD < 0.2 | WD > 1.2)
```

### Vérifier les métriques de parcelles
```r
# Charger
plot_metrics <- read_csv2("4_Plots_metrics/PlotsMetrics_FieldDataInventories.csv")

# Gammes attendues
summary(plot_metrics$WD_BA)      # 0.45-0.75
summary(plot_metrics$G)          # 20-50 m²/ha
summary(plot_metrics$AGB)        # 200-600 Mg/ha
summary(plot_metrics$n_stems)    # 300-700 tiges/ha

# Corrélation WD_BA vs WD_mean
cor(plot_metrics$WD_BA, plot_metrics$WD_mean)  # ~0.8-0.95
```

### Vérifier les tempéraments
```r
# Charger
temperament <- read_csv2("4_Plots_metrics/PlotsMetrics_Temperament_all.csv") %>%
  filter(Type_subdivision == 3, Type_Proportion == "G.rel")

# Les proportions doivent sommer à 1
temperament %>%
  group_by(plot_ref) %>%
  summarise(total = sum(Proportion, na.rm=TRUE))

# Distribution des proportions
summary(temperament$Proportion[temperament$Temperament == "héliophile"])
summary(temperament$Proportion[temperament$Temperament == "tolérante à l'ombrage"])
```

## Problèmes courants

### WD non attribué (NA)
**Cause**: Espèce inconnue ou nom invalide
**Solution**: Vérifier la correction taxonomique, utiliser WD par défaut (0.55)

### Hauteurs estimées aberrantes
**Cause**: Diamètres extrêmes ou données manquantes
**Solution**: Vérifier les allométries, exclure les valeurs aberrantes

### Proportions de tempérament ne sommant pas à 1
**Cause**: Arbres sans tempérament attribué
**Solution**: Attribuer un tempérament par défaut ou exclure ces arbres

### Script 1 très lent
**Solution**: C'est normal (~30-60 min), la correction taxonomique interroge des bases de données en ligne avec cache

## Interprétation

### Parcelle à WD_BA élevé (> 0.65)
- **Dominance** d'espèces tolérantes à l'ombrage
- Forêt mature, peu perturbée
- **Hypothèse**: Structure verticale plus stratifiée?

### Parcelle à WD_BA faible (< 0.55)
- **Dominance** d'espèces héliophiles
- Forêt secondaire ou récemment perturbée
- **Hypothèse**: Structure verticale plus homogène?

### Lien attendu avec les trouées
**Hypothèse centrale du projet**:
- WD élevé ↔ Profils progressifs (stratification)
- WD faible ↔ Profils abrupts (canopée homogène)

→ Testé dans le module 5!
