rm(list=ls())
gc()

# Initialisation ----------------------------------------------------------

pkgs = c("rio","foreach", "doParallel","lidR", "terra", "sf", "tidyverse", "ggplot2", "BIOMASS", "here")

to_install = !pkgs %in% installed.packages() ; if(any(to_install)) {install.packages(pkgs[to_install])} ; inst = lapply(pkgs, library, character.only = TRUE) # load them

path0 = here("4_Plots_metrics")

path_output = file.path(path0, "output_temperament")
if(!dir.exists(path_output)){dir.create(path_output)}

setwd(path_output)

# Loading raw data -------------------------------------------------------------------------

plot_names <- read.csv2(here("final_plot_name.csv")) %>%
  pull(plot_name)

field_inventories = rio::import( # ATTENTION : bien utiliser la base de données où l'information relative aux tempéraments a été injectée, à partir du code n°4 'FieldDataInventories_Analyse_Temperament.R'
  paste0(path0, "field_data_inventories_with_temperament.csv"),
  sep = ";",
  dec = ","
) %>% 
  as_tibble() %>%
  filter(plot_ref %in% plot_names) # Charger les données d'inventaires pour les plots présents sous les survols LiDAR

PlotMetrics = rio::import(
  paste0(path0, "PlotsMetrics_FieldDataInventories.csv"),
  sep = ";",
  dec = ","
) %>% 
  as_tibble() %>%
  filter(plot_ref %in% plot_names) # Charger les métriques pour les plots présents sous les survols LiDAR

# Fonctions ---------------------------------------------------------------
# Fonction pour quantifier la proportion de tempérament par plot
# Possibilité de choisir le nombre de subdivision des tempéraments : 2 ou 3, et le type de pondération : surface terrière (G) ou individus (ind)
Temperament.by.plot_save = function(field_inventories,subdivision = 2, pond.rel_chr) { # Fonction pour quantifier le tempérament des espèces par parcelles
  
  #Le nombre de subdivision ({subdivision}) permet de choisir si on veut inclure les héliophiles non pionniers demandeurs de lumière dans la catégorie héliophile, ou la garder dans une catégorie séparées (3 subdivision)
  #La pondération pour le calcul de la proportion de tempérament ({pond.rel_chr}) permet de choisir si on souhaite calculer une proportion en fonction du nombre d'individus observés par plot (pond.rel_chr = "ind.rel") ou si on souhaite pondérer par la surface terrière (pond.rel_chr = "G.rel)
  # --> Il est nécessaire d'avoir une colonne qui représente la surface terrière proportionnelle au plot de chaque observation (G_arbre / sum(G_arbre du plot)), ou la proportion que représente un individus (1/ sum(observations / plots))
  # Ici il est choisit d'ignorer les arbres qui ont été catégorisé comme 'Inconnu' ; c'est à dire ceux qui sont absents de la base de données CoForTraits, mais également ceux dont on n'a pas pu déterminer le tempérament dans la base de données CoForTraits
  
  if (!pond.rel_chr %in% names(field_inventories)) {
    stop("pond.rel_chr doit être 'G.rel' ou 'ind.rel', une colonne qui donne une proportion relatif à un individu par rapport au reste de la parcelle")
  }
  
  pond.rel = sym(pond.rel_chr) #Transforme la chaîne de caractère en un symbole qui représente un objet de l'environnement
  
  if (subdivision == 2) { 
    
    Temp.By.Plot = field_inventories %>%
      group_by(plot_ref, Temperament) %>%
      summarise(Prop.Temp = sum(!!pond.rel, na.rm = TRUE) * 100, .groups = 'drop') %>% #On  utilise '!!' avant le symbole pour pouvoir le call dans les fonctions
      pivot_wider(names_from = Temperament, values_from = Prop.Temp, values_fill = list(count = 0)) %>%
      mutate_all(~ifelse(is.na(.), 0, .)) %>%
      mutate(
        helio = `héliophile` / (`héliophile` + `tolérante à l'ombrage`),
        ombre = `tolérante à l'ombrage` / (`héliophile` + `tolérante à l'ombrage`)
      ) %>%
      select(plot_ref, helio, ombre) %>%
      rename_with(.fn = ~paste0(., "_", pond.rel_chr), .cols = c("helio", "ombre"))
    
  } else if (subdivision == 3) {
    
    Temp.By.Plot = field_inventories %>%
      group_by(plot_ref, Temperament) %>%
      summarise(Prop.Temp = sum(!!pond.rel, na.rm = TRUE) * 100, .groups = 'drop') %>%
      pivot_wider(names_from = Temperament, values_from = Prop.Temp, values_fill = list(count = 0)) %>%
      mutate_all(~ifelse(is.na(.), 0, .)) %>%
      mutate(
        helio = `héliophile` / (`héliophile` + `tolérante à l'ombrage` + `héliophile non pionnière`),
        ombre = `tolérante à l'ombrage` / (`héliophile` + `tolérante à l'ombrage` + `héliophile non pionnière`),
        helio_np = `héliophile non pionnière` / (`héliophile` + `tolérante à l'ombrage` + `héliophile non pionnière`)
      ) %>%
      select(plot_ref, helio, ombre, helio_np) %>%
      rename_with(.fn = ~paste0(., "_", pond.rel_chr), .cols = c("helio", "ombre", "helio_np"))
    
  } else { Temp.By.Plot = "Subdivision doit être égal à 2 ou 3 (integer)" }
  
  return(Temp.By.Plot)
}
Temperament.by.plot = function(field_inventories,subdivision = 2, pond.rel_chr = "ind.rel") { # Fonction pour quantifier le tempérament des espèces par parcelles
  
  #Le nombre de subdivision ({subdivision}) permet de choisir si on veut inclure les héliophiles non pionniers demandeurs de lumière dans la catégorie héliophile, ou la garder dans une catégorie séparées (3 subdivision)
  #La pondération pour le calcul de la proportion de tempérament ({pond.rel_chr}) permet de choisir si on souhaite calculer une proportion en fonction du nombre d'individus observés par plot (pond.rel_chr = "ind.rel") ou si on souhaite pondérer par la surface terrière (pond.rel_chr = "G.rel)
  # --> Il est nécessaire d'avoir une colonne qui représente la surface terrière proportionnelle au plot de chaque observation (G_arbre / sum(G_arbre du plot)), ou la proportion que représente un individus (1/ sum(observations / plots))
  # Ici il est choisit d'ignorer les arbres qui ont été catégorisé comme 'Inconnu' ; c'est à dire ceux qui sont absents de la base de données CoForTraits, mais également ceux dont on n'a pas pu déterminer le tempérament dans la base de données CoForTraits
  
  if (!(pond.rel_chr %in% c('G.rel', 'ind.rel'))) {
    stop("pond.rel_chr doit être 'G.rel' ou 'ind.rel', selon le type de propotion : relative à un individu ou à la surface terrière (G) par rapport au reste de la parcelle")
  }
  
  if(pond.rel_chr == "G.rel"){
    
    G = "G"
    
    if (!G %in% names(field_inventories)) {
      stop("Une colonne nommée 'G' doit être présente et contenir la surface terrière individuel d'un arbre")
    }
  }
  
  pond.rel = sym(pond.rel_chr) #Transforme la chaîne de caractère en un symbole qui représente un objet de l'environnement
  
  if (subdivision == 2) { 
    
    Temp.By.Plot = field_inventories %>%
      filter(Temperament %in% c("héliophile" , "tolérante à l'ombrage")) %>%
      group_by(plot_ref) %>%
      mutate(
        G.rel = G / sum(G, na.rm = TRUE),
        ind.rel = 1/n()
      ) %>%
      ungroup() %>%
      group_by(plot_ref, Temperament) %>%
      summarise(Proportion = sum(!!pond.rel, na.rm = TRUE) * 100, .groups = 'drop') %>% #On  utilise '!!' avant le symbole pour pouvoir le call dans les fonctions
      mutate(
        Type_Proportion = pond.rel_chr,
        Type_subdivision = 2
      )
    
  } else if (subdivision == 3) {
    
    Temp.By.Plot = field_inventories %>%
      filter(Temperament %in% c("héliophile" , "tolérante à l'ombrage", "héliophile non pionnière")) %>%
      group_by(plot_ref) %>%
      mutate(
        G.rel = G / sum(G, na.rm = TRUE),
        ind.rel = 1/n()
      ) %>%
      ungroup() %>%
      group_by(plot_ref, Temperament) %>%
      summarise(Proportion = sum(!!pond.rel, na.rm = TRUE) * 100, .groups = 'drop') %>%
      dplyr::select(Temperament = Temperament, everything()) %>%
      mutate(
        Type_Proportion = pond.rel_chr,
        Type_subdivision = 3
      ) 
      
  } else { Temp.By.Plot = "Subdivision doit être égal à 2 ou 3 (integer)" }
  
  return(Temp.By.Plot)
}
# Fonction pour quantifier la proportion de tempérament par plot et par classe de valeur d'une colonne, ici avec objectif de faire des classes de DBH et de H mais la fonction est applicable à n'importe quelle métrique de terrain comme l'AGB, etc.
Temperament.by.plot_with.breaks = function(field_inventories,subdivision = 2, pond.rel_chr = "G.rel", breaks = FALSE, breaks_col = "dbh") { # Fonction pour quantifier le tempérament des espèces par parcelles
  
  #Le nombre de subdivision ({subdivision}) permet de choisir si on veut inclure les héliophiles non pionniers demandeurs de lumière dans la catégorie héliophile, ou la garder dans une catégorie séparées (3 subdivision)
  #La pondération pour le calcul de la proportion de tempérament ({pond.rel_chr}) permet de choisir si on souhaite calculer une proportion en fonction du nombre d'individus observés par plot (pond.rel_chr = "ind.rel") ou si on souhaite pondérer par la surface terrière (pond.rel_chr = "G.rel)
  # --> Il est nécessaire d'avoir une colonne qui représente la surface terrière proportionnelle au plot de chaque observation (G_arbre / sum(G_arbre du plot)), ou la proportion que représente un individus (1/ sum(observations / plots))
  # Ici il est choisit d'ignorer les arbres qui ont été catégorisé comme 'Inconnu' ; c'est à dire ceux qui sont absents de la base de données CoForTraits, mais également ceux dont on n'a pas pu déterminer le tempérament dans la base de données CoForTraits
  
  if (!pond.rel_chr %in% names(field_inventories)) {
    stop("pond.rel_chr doit être 'G.rel' ou 'ind.rel', une colonne qui donne une proportion relatif à un individu par rapport au reste de la parcelle")
  }
  
  pond.rel = sym(pond.rel_chr) #Transforme la chaîne de caractère en un symbole qui représente un objet de l'environnement
  
  if (!breaks_col %in% c("dbh", "h")) {
    stop("breaks_col must be either 'dbh' or 'h'")
  }
  
  # Define breaks and apply to the specified column
  breaks_def = if (breaks_col == "dbh") {
    cut(field_inventories$dbh, breaks = c(-Inf, breaks, Inf), include.lowest = TRUE)
  } else if (breaks_col == "h"){
    cut(field_inventories$h, breaks = c(-Inf, breaks, Inf), include.lowest = TRUE)
  } else { stop("La fonction n'est pas faites pour faire des intervalles autre que sur la colonne h ou dbh désolé, à creuser")}
  
  if (subdivision == 2) { 
    
    Temp.By.Plot = field_inventories %>%
      mutate(class_breaks = breaks_def) %>%
      filter(Temperament %in% c("héliophile" , "tolérante à l'ombrage")) %>%
      group_by(plot_ref, class_breaks) %>%
      mutate(
        G.rel = G / sum(G, na.rm = TRUE),
        ind.rel = 1/n()
      ) %>%
      ungroup() %>%
      group_by(plot_ref, Temperament, class_breaks) %>%
      summarise(Proportion = sum(!!pond.rel, na.rm = TRUE) * 100, .groups = 'drop') %>%
      mutate(
        Type_Proportion = pond.rel_chr,
        class_type = breaks_col,
        Type_subdivision = 2
      )
      
    
  } else if (subdivision == 3) {
    
    Temp.By.Plot = field_inventories %>%
      mutate(class_breaks = breaks_def) %>%
      filter(Temperament %in% c("héliophile" , "tolérante à l'ombrage", "héliophile non pionnière")) %>%
      group_by(plot_ref, class_breaks) %>%
      mutate(
        G.rel = G / sum(G, na.rm = TRUE),
        ind.rel = 1/n()
      ) %>%
      ungroup() %>%
      group_by(plot_ref, Temperament, class_breaks) %>%
      summarise(Proportion = sum(!!pond.rel, na.rm = TRUE) * 100, .groups = 'drop')  %>%
      select(Temperament = Temperament, everything()) %>%
      mutate(
        Type_Proportion = pond.rel_chr,
        class_type = breaks_col,
        Type_subdivision = 3
      )
    
  } else { stop("Subdivision doit être égal à 2 ou 3 (integer)") }
  
  return(Temp.By.Plot)
}

# Fonction pour générer un graphique de fréquence cumulée proportionnelle
generate_cumulative_proportion_plot <- function(data, column_name) {
  if (!column_name %in% names(data)) {
    stop("Le nom de la colonne spécifiée n'existe pas dans le dataframe.")
  }
  
  # Arrondir les valeurs de la colonne et calculer la fréquence cumulée proportionnelle
  data_for_plot <- data %>%
    mutate(RoundedValue = round(!!sym(column_name))) %>%
    group_by(RoundedValue) %>%
    summarise(Count = n()) %>%
    ungroup() %>%
    arrange(RoundedValue) %>%
    mutate(CumulativeCount = cumsum(Count),
           TotalCount = max(CumulativeCount),
           Proportion = CumulativeCount / TotalCount) %>%
    select(RoundedValue, Proportion)
  
  # Générer le graphique personnalisé
  ggplot(data_for_plot, aes(x = RoundedValue, y = Proportion)) +
    geom_line() + # Ne montre que la ligne
    scale_x_continuous(breaks = function(x) pretty(x, n = 10)) + # Ajuste les breaks de l'axe des x
    scale_y_continuous(labels = scales::percent_format(), breaks = seq(0, 1, by = 0.1)) + # Définit les breaks de l'axe des y tous les 10%
    labs(
      title = paste("Fréquence cumulée proportionnelle pour", column_name),
      x = paste("Valeur de", column_name),
      y = "Fréquence Cumulée Proportionnelle (%)"
    ) +
    theme_minimal()
}

# Quantification du tempérament par parcelle sans distinction des catégories de hauteur ou de diamètre ----

# Combinaison des résultats dans un seul tibble
PlotsMetrics_Temperament_all <- bind_rows(Temperament.by.plot(field_inventories, subdivision = 3, pond.rel_chr = 'G.rel')) %>%
  bind_rows(Temperament.by.plot(field_inventories, subdivision = 3, pond.rel_chr = 'ind.rel'))
  

rio::export(
  PlotsMetrics_Temperament_all,
  paste0(path0, "PlotsMetrics_Temperament_all.csv"),
  sep = ";",
  dec = ","
)

# Graphique pour comparer l'influence du choix de la métrique relative : G ou ind ----

data_graph <- PlotsMetrics_Temperament_all %>%
  mutate(
    measure = paste(Temperament)
  ) 
  

p1 = ggplot(data_graph, aes(x = Type_Proportion, y = Proportion, fill = Type_Proportion)) +
  geom_boxplot() +
  facet_wrap(~measure, scales = "free") +
  theme_light() +
  labs(
    title = "Comparaison de l'influence de la pondération par plot :relative au nombre d'individus \n(ind.rel) ou à la surface terrière (G.rel) sur la proportion de tempérament par plot",
    y = "Proportion tempérament",
    x = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )


ggsave(paste0(path0, "Temperament_influence_pondérationPlots_boxplots.svg"), plot = p1, width = 15, height = 10, device = "svg")

