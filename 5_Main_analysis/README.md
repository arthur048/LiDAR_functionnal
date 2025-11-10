# 5. Analyses Principales

## Vue d'ensemble

Ce module constitue le cœur analytique du projet. Il teste l'hypothèse centrale: **la structure verticale des forêts (profils de trouées) est-elle corrélée avec la composition fonctionnelle (densité du bois, tempérament des espèces)?**

## Objectif

Identifier et quantifier les relations entre:
- **Variables explicatives**: Métriques de structure verticale dérivées des profils de trouées
- **Variables cibles**: Densité du bois (WD_BA), proportions de tempérament (héliophiles, tolérants)

## Hypothèse scientifique

**H1**: Les forêts dominées par des espèces à bois dense (tolérantes à l'ombrage) présentent une structure verticale plus stratifiée (profils de trouées progressifs).

**H2**: Les forêts dominées par des espèces à bois léger (héliophiles) présentent une structure verticale plus homogène (profils de trouées abrupts).

## Données d'entrée

Toutes consolidées depuis les modules précédents:

### Métriques de parcelles (module 4)
- `../4_Plots_metrics/PlotsMetrics_FieldDataInventories.csv` - Composition, WD, structure
- `../4_Plots_metrics/plot_metrics_gaps.csv` - Métriques de trouées
- `../4_Plots_metrics/PlotsMetrics_Temperament_all.csv` - Proportions de tempérament
- `../4_Plots_metrics/field_data_inventories_with_temperament.csv` - Inventaires avec tempérament

### Profils bruts (module 3)
- `../3_Gaps/output/results_gaps_metrics.csv` - Profils verticaux complets

### Localisation
- `../plots_localisation.csv` (racine) - Coordonnées géographiques des parcelles

## Structure du module

```
5_Main_analysis/
├── functions/                  # Fonctions personnalisées
│   ├── graphics_profiles.R     # Visualisation des profils
│   ├── graphics_correlation.R  # Visualisation des corrélations
│   ├── metric_labels.R         # Noms lisibles des métriques
│   └── metric_labels_simplified.R
├── input/                      # Données consolidées
│   ├── plot_metrics.csv        # Métriques agrégées
│   ├── gaps_metrics.csv        # Profils verticaux
│   └── combined_data.csv       # Toutes les données jointes
├── output/
│   ├── data/                   # Résultats numériques
│   ├── figures/                # Graphiques
│   └── tables/                 # Tableaux de résultats
└── Scripts d'analyse (1-6)
```

## Scripts d'analyse

### `1_data_preparation.R`

**Fonction**: Prépare et consolide toutes les données

**Actions**:
1. **Charge** toutes les données des modules précédents
2. **Joint** les métriques de trouées avec les métriques d'inventaires
3. **Catégorise** les parcelles selon WD_BA:
   - Quartiles: 4 groupes équilibrés
   - Groupes uniformes: classes régulières
4. **Exporte** les données consolidées vers `input/`

**Sorties**:
- `input/plot_metrics.csv` - Une ligne par parcelle, toutes métriques
- `input/gaps_metrics.csv` - Profils verticaux
- `input/combined_data.csv` - Données complètes pour visualisations

### `2_profiles_visualization.R`

**Fonction**: Visualise les profils verticaux et leurs relations avec WD

**Analyses graphiques**:

#### A. Courbes de proportion de trouées
```r
plot_gap_curves_by_wd(data, wd_palette)
```
- Profils individuels par parcelle (lignes fines)
- Moyennes par catégorie de WD (lignes épaisses)
- Permet de visualiser si les catégories de WD se distinguent

#### B. Heatmap des profils
```r
create_profiles_heatmap(data)
```
- Lignes = hauteurs (0-60m)
- Colonnes = parcelles ordonnées par WD croissant
- Couleur = proportion de trouées
- Révèle des patterns visuels dans l'organisation verticale

#### C. Profils types
```r
create_typical_profiles(data, show_deviation = TRUE)
```
- Profil moyen par catégorie de WD
- Zones d'incertitude (écart-type ou erreur standard)
- Interpolation lisse pour comparaison

#### D. Signatures différentielles
```r
create_differential_signatures(data)
```
- Dérivée première (dP/dh): Taux de changement
- Identifie les zones de transition rapide
- Compare les pentes entre catégories de WD

**Sorties**: `output/figures/profiles/` (heatmaps, courbes, etc.)

### `3_correlation.R`

**Fonction**: Analyse les corrélations entre métriques

**Analyses**:

#### A. Corrélations univariées
Pour chaque **variable cible** (WD_BA, WD_mean, prop_g_helio, etc.):
1. **Calcul** des corrélations avec toutes les métriques de trouées
   - Pearson (linéaire)
   - Spearman (monotone)
2. **Test de significativité** (p-values)
3. **Classement** par force de corrélation
4. **Visualisation** des top N métriques

#### B. Matrices de corrélation
```r
corrplot(cor_matrix, method = "color")
```
- Toutes les métriques entre elles
- Identification des redondances
- Sélection des métriques indépendantes

#### C. Scatterplots
```r
ggplot(aes(x = metric, y = WD_BA)) +
  geom_point() +
  geom_smooth(method = "lm")
```
- Relations bivariées détaillées
- Identification des outliers
- Tests de linéarité

**Sorties**:
- `output/figures/correlation/{target_var}/` - Un dossier par variable cible
- `output/tables/correlation/` - Tableaux de corrélations

### `4_multivariate.R` et `5_multivariate_suite.R`

**Fonction**: Analyses multivariées

**Méthodes**:

#### A. Analyse en Composantes Principales (PCA)
```r
prcomp(metrics_data, scale = TRUE)
```
- **Réduction de dimensionnalité** des ~50-100 métriques de trouées
- **Identification** des axes principaux de variation
- **Projection** des parcelles dans l'espace réduit
- **Interprétation**: Quelles métriques expliquent la variabilité?

**Visualisations**:
- Biplot: Variables et individus
- Contribution des variables aux axes
- Pourcentage de variance expliquée

#### B. Clustering
```r
hclust(dist(data))
```
- **Regroupement** des parcelles similaires
- **Comparaison** avec les catégories de WD
- **Validation**: Les groupes naturels correspondent-ils au gradient de WD?

### `5_PCA_combined.R`

**Fonction**: PCA intégrative

- Combine métriques de trouées ET métriques d'inventaires
- Identifie les patterns globaux de co-variation
- Teste si WD structure l'espace multivarié

### `6_Modelisation.R`

**Fonction**: Modélisation prédictive

**Objectif**: Prédire WD_BA à partir des métriques de trouées

**Modèles testés**:
1. **Régression linéaire multiple**
   ```r
   lm(WD_BA ~ metric1 + metric2 + ...)
   ```
2. **Régression LASSO**
   - Sélection automatique de variables
   - Pénalisation pour éviter le sur-ajustement
3. **Random Forest** (si implémenté)
   - Non-linéaire
   - Importance des variables

**Validation**:
- Cross-validation
- R², RMSE, MAE
- Analyse des résidus

**Sorties**:
- Modèles sauvegardés
- Prédictions vs observations
- Importance des variables

## Fonctions utilitaires

### `functions/graphics_profiles.R`

Contient toutes les fonctions de visualisation des profils:
- `plot_gap_curves_by_wd()` - Courbes par catégorie
- `create_profiles_heatmap()` - Heatmap ordonnée
- `create_typical_profiles()` - Profils types avec incertitude
- `create_differential_signatures()` - Dérivées
- `plot_decompose_curve_segments()` - Décomposition en zones
- `plot_derivatives()` - Analyse différentielle détaillée
- `plot_multiscale_analysis()` - Lissage multi-échelle

### `functions/graphics_correlation.R`

Fonctions de visualisation des corrélations (à développer selon besoins).

### `functions/metric_labels.R`

**Fonction**: Traduction des codes de métriques en noms lisibles

```r
get_metric_labels()
```

Retourne un vecteur nommé:
```r
c(
  "prop_at_10m" = "Proportion de trouées à 10m",
  "h50" = "Hauteur à 50% de trouées",
  "Smax" = "Pente maximale",
  ...
)
```

**Usage**: Labels pour graphiques et tableaux

## Paramètres d'analyse

### Variables cibles principales
```r
target_vars <- c(
  "WD_BA",           # ⭐ Cible principale
  "WD_mean",
  "prop_g_helio",    # Proportion héliophiles (G)
  "prop_g_shade",    # Proportion tolérants (G)
  "prop_ind_helio",  # Proportion héliophiles (N)
  "prop_ind_shade"   # Proportion tolérants (N)
)
```

### Méthodes de corrélation
```r
corr_methods <- c("Pearson", "Spearman")
```

### Nombre de métriques à visualiser
```r
top_n_value <- 20  # Top 20 corrélations
```

### Résolution des graphiques
```r
dpi_value <- 600  # Publication quality
```

## Sorties principales

### Figures clés (exemples)
- `output/figures/profiles/gap_curves_by_wd_quartile.jpg`
- `output/figures/profiles/profiles_heatmap.jpg`
- `output/figures/correlation/WD_BA/top_correlations_pearson.jpg`
- `output/figures/correlation/WD_BA/scatterplot_h50_vs_WD_BA.jpg`

### Tableaux de résultats
- `output/tables/correlation/correlations_WD_BA.csv`
- `output/tables/correlation/top_metrics_all_targets.csv`

### Données intermédiaires
- `output/data/pca_results.RData`
- `output/data/model_predictions.csv`

## Dépendances R

```r
library(tidyverse)     # Manipulation, visualisation
library(ggpubr)        # Graphiques avancés
library(viridis)       # Palettes de couleurs
library(corrplot)      # Matrices de corrélation
library(patchwork)     # Assemblage de graphiques
library(ggridges)      # Density ridges
library(zoo)           # Séries temporelles (lissage)
library(fields)        # Interpolation spatiale
library(ggrepel)       # Labels non-superposés
library(cowplot)       # Assemblage de plots
library(foreach)       # Parallélisation (si nécessaire)
library(doParallel)    # Backend parallèle
```

## Utilisation

### Workflow séquentiel
```r
# Depuis 5_Main_analysis/

# 1. Préparation des données
source("1_data_preparation.R")

# 2. Visualisation des profils
source("2_profiles_visualization.R")

# 3. Analyses de corrélation
source("3_correlation.R")

# 4-5. Analyses multivariées
source("4_multivariate.R")
source("5_multivariate_suite.R")

# 6. Modélisation
source("6_Modelisation.R")
```

### Génération ciblée
Pour regénérer uniquement les graphiques de corrélation pour WD_BA:
```r
source("functions/metric_labels.R")
source("functions/graphics_correlation.R")

# Charger les données
plot_metrics <- read_csv2("input/plot_metrics.csv")

# Calculer et visualiser
# ... (code spécifique dans script 3)
```

## Interprétation des résultats

### Corrélations attendues (hypothèses)

**WD_BA élevé** devrait être associé avec:
- **h50 élevé**: Point d'inflexion plus haut (canopée haute)
- **Smax faible**: Pente plus douce (transition progressive)
- **AUC élevé**: Plus de volume sous la courbe (stratification)
- **Épaisseur de transition élevée**: Zone de transition large

**WD_BA faible** devrait être associé avec:
- **h50 bas**: Point d'inflexion plus bas
- **Smax élevé**: Transition abrupte
- **AUC faible**: Moins de stratification
- **Épaisseur de transition faible**: Transition courte

### Validation des hypothèses

**Si corrélations significatives (|r| > 0.5, p < 0.05)**:
- Confirme le lien structure verticale ↔ composition fonctionnelle
- Justifie l'utilisation du LiDAR pour inférer la composition

**Si corrélations faibles**:
- Autres facteurs dominent (sol, climat, historique)
- Échelle spatiale inadéquate (1 ha trop petit?)
- Métrique de WD incomplète (manque d'autres traits?)

## Problèmes courants

### Corrélations non significatives
**Causes possibles**:
- Taille d'échantillon insuffisante (< 30 parcelles)
- Hétérogénéité des sites (mélange de types de forêts)
- Qualité des données (erreurs de mesure)
**Solution**: Analyses stratifiées par site, augmenter l'échantillon

### Multicolinéarité
**Problème**: Métriques de trouées très corrélées entre elles
**Solution**: PCA pour réduire dimensionnalité, sélection de métriques indépendantes

### Outliers
**Identification**: Graphiques de résidus, distance de Cook
**Traitement**: Vérifier les données, analyses avec/sans outliers

### Graphiques illisibles
**Cause**: Trop de variables ou de parcelles
**Solution**: Ajuster `top_n_value`, filtrer les données

## Publications et communications

Les résultats de ce module constituent le matériel pour:
- **Figures principales** de l'article scientifique
- **Matériel supplémentaire** (heatmaps, tous les scatterplots)
- **Présentations** (courbes synthétiques, profils types)

### Recommandations graphiques
- **Pour article**: DPI 600, format JPEG/PNG
- **Pour présentation**: DPI 300, format PNG transparent
- **Palettes**: Viridis (colorblind-friendly)
- **Polices**: Arial ou Helvetica (universelles)
