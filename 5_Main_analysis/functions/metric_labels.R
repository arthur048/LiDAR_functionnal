#' Définit des noms lisibles pour les métriques utilisées dans les analyses
#' 
#' @return Un vecteur nommé contenant les correspondances entre codes et libellés
#' @export
get_metric_labels <- function() {
  metric_labels <- c(
    # Métriques basiques
    prop_at_10m = "Proportion de trouées à 10m",
    prop_at_15m = "Proportion de trouées à 15m",
    prop_at_20m = "Proportion de trouées à 20m",
    prop_at_30m = "Proportion de trouées à 30m",
    h10 = "Hauteur à 10% de trouées",
    h25 = "Hauteur à 25% de trouées",
    h50 = "Hauteur à 50% de trouées",
    h75 = "Hauteur à 75% de trouées",
    h90 = "Hauteur à 90% de trouées",
    Pmax = "Proportion maximale de trouées",
    Smax = "Pente maximale",
    hSmax = "Hauteur à pente maximale",
    auc = "Aire sous la courbe",
    
    # Métriques de logistique
    logistic_a = "Logistique - asymptote inf.",
    logistic_b = "Logistique - asymptote sup.",
    logistic_h0 = "Logistique - point d'inflexion",
    logistic_k = "Logistique - pente",
    logistic_R2 = "Logistique - R²",
    
    # Métriques combinées
    IDV = "Indice de Développement Vertical",
    ROS = "Ratio d'Occupation Stratifiée",
    IAS = "Indice d'Asymétrie Structurelle",
    ICT = "Indice de Contraste de Transition",
    IPN = "Indice de Progression Normalisée",
    ET = "Épaisseur de Transition",
    CSI = "Indice de Structure de Canopée",
    GVN = "Gradient Vertical Normalisé",
    RT = "Ratio de Transition",
    HPSI = "Harmonie Position-Pente",
    
    # Métriques zonales
    basal_thickness = "Épaisseur zone basale",
    transition_thickness = "Épaisseur zone de transition",
    upper_thickness = "Épaisseur zone supérieure",
    basal_norm_auc = "AUC normalisée zone basale",
    transition_norm_auc = "AUC normalisée zone transition",
    upper_norm_auc = "AUC normalisée zone supérieure",
    transition_slope = "Pente de la zone de transition",
    transition_r2 = "Linéarité de la transition (R²)",
    zone_heterogeneity_index = "Indice d'hétérogénéité zonale",
    
    # Métriques différentielles
    height_max_d1 = "Hauteur à pente maximale",
    max_d1 = "Pente maximale",
    height_max_d2 = "Hauteur à accélération max",
    max_d2 = "Accélération maximale",
    height_min_d2 = "Hauteur à décélération max",
    min_d2 = "Décélération maximale",
    asymmetry_ratio = "Ratio d'asymétrie",
    acceleration_phase = "Durée phase d'accélération",
    deceleration_phase = "Durée phase de décélération",
    phase_ratio = "Ratio des phases",
    
    # Métriques multi-échelles
    roughness_small = "Rugosité à petite échelle (2m)",
    roughness_medium = "Rugosité à échelle moyenne (10m)",
    roughness_large = "Rugosité à grande échelle (20m)",
    dominant_height_small = "Hauteur dominante à petite échelle",
    dominant_height_medium = "Hauteur dominante à échelle moyenne",
    dominant_height_large = "Hauteur dominante à grande échelle",
    energy_slope = "Pente spectrale d'énergie",
    energy_ratio_small_large = "Ratio d'énergie petite/grande échelle",
    
    # Métriques d'anisotropie
    anisotropy_mean = "Anisotropie moyenne",
    anisotropy_max = "Anisotropie maximale",
    height_max_anisotropy = "Hauteur d'anisotropie maximale",
    anisotropy_transition = "Anisotropie zone de transition",
    anisotropy_ratio = "Ratio d'anisotropie max/moyenne",
    
    # Indices optimisés
    ISA = "Indice de Stratification Asymétrique",
    ICV = "Indice de Complexité Verticale",
    IDT = "Indice de Différenciation Transitionnelle",
    
    # Variables de composition
    WD_BA = "Densité du bois pondérée par G",
    prop_g_helio = "Proportion G espèces héliophiles",
    prop_g_npld = "Proportion G espèces héliophiles non pionnières",
    prop_g_shade = "Proportion G espèces tolérantes à l'ombre",
    prop_ind_helio = "Proportion individus héliophiles",
    prop_ind_npld = "Proportion individus héliophiles non pionniers",
    prop_ind_shade = "Proportion individus tolérants à l'ombre",
    
    # Métriques transformées
    log_transition_thickness = "Log épaisseur zone transition",
    log_max_d1 = "Log pente maximale",
    ratio_h50_h90 = "Ratio h50/h90",
    ratio_basal_transition = "Ratio épaisseur basale/transition",
    h50_squared = "h50 au carré",
    anisotropy_squared = "Anisotropie au carré"
  )
  
  return(metric_labels)
}