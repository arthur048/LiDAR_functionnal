#' Définit des noms lisibles pour les métriques utilisées dans les analyses
#' 
#' @return Un vecteur nommé contenant les correspondances entre codes et libellés
#' @export
get_metric_labels <- function() {
  metric_labels <- c(
    # Métriques basiques
    prop_at_5m = "Proportion de trouées à 5m",
    prop_at_10m = "Proportion de trouées à 10m",
    prop_at_15m = "Proportion de trouées à 15m",
    prop_at_20m = "Proportion de trouées à 20m",
    prop_at_25m = "Proportion de trouées à 25m",
    prop_at_30m = "Proportion de trouées à 30m",
    prop_at_35m = "Proportion de trouées à 35m",
    prop_at_40m = "Proportion de trouées à 40m",
    h10 = "Hauteur à 10% de trouées",
    h25 = "Hauteur à 25% de trouées",
    h50 = "Hauteur à 50% de trouées",
    h75 = "Hauteur à 75% de trouées",
    h90 = "Hauteur à 90% de trouées",
    h100 = "Hauteur à 100% de trouées",
    
    # Aires sous la courbe
    auc_norm = "AUC normalisée (courbe entière)",
    auc = "AUC non normalisée (courbe entière)",
    auc_0_25_norm = "AUC normalisée (0% à 25% de trouées)",
    auc_0_25 = "AUC non normalisée (0% à 25% de trouées)",
    auc_25_75_norm = "AUC normalisée (25% à 75% de trouées)",
    auc_25_75 = "AUC non normalisée (25% à 75% de trouées)",
    auc_75_100_norm = "AUC normalisée (75% à 100% de trouées)",
    auc_75_100 = "AUC non normalisée (75% à 100% de trouées)",
    
    # Métriques zonales
    basal_thickness = "Épaisseur zone basale",
    transition_thickness = "Épaisseur zone de transition",
    upper_thickness = "Épaisseur zone supérieure",
    basal_auc_norm = "AUC normalisée zone basale",
    basal_auc = "AUC non normalisée zone basale",
    transition_auc_norm = "AUC normalisée zone transition",
    transition_auc = "AUC non normalisée zone transition",
    upper_auc_norm = "AUC normalisée zone supérieure",
    upper_auc = "AUC non normalisée zone supérieure",
    slope_basal = "Pente de la zone basale",
    slope_transition = "Pente de la zone de transition",
    slope_upper = "Pente de la zone supérieure",
    
    # Métriques différentielles
    height_max_d1 = "Hauteur à pente maximale (d1)",
    max_d1 = "Pente maximale (d1)",
    height_max_d2 = "Hauteur à accélération max (d2)",
    max_d2 = "Accélération maximale (d2)",
    height_min_d2 = "Hauteur à décélération max (d2)",
    min_d2 = "Décélération maximale (d2)",
    
    # Variables de composition
    WD_mean = "Densité moyenne du bois"
  )
  
  return(metric_labels)
}