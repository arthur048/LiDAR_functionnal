rm(list=ls())
gc()

#installr::updateR()
#update.packages()

#library(help = "lidR")

# Open the file in RStudio to edit it : ---> usethis::edit_file("~/AppData/Roaming/RStudio/templates/default.R") <---

# Initialisation ----------------------------------------------------------

pkgs.a = c("rio", "lubridate", "doBy", "rstatix", "randomForest", "blockCV", "Metrics", "Boruta",
           "caret", "PMCMRplus","ggpubr", "multcomp", "emmeans","ggpol", "rstatix", "ggpmisc", "corrplot", "FactoMineR", "factoextra",
           "ggthemes", "svglite", "viridis", "usethis", "broom", "sysfonts", "showtext", "broom", "lidR", "rgdal", "parallel", "remotes",
           "raster", "sf", "stars", "terra", "units", "tidyverse", "devtools", "BIOMASS", "igraph", "poweRlaw", "gridExtra")
pkgs.priority = c("raster", "stars", "lidR", "terra", "sf", "tidyverse")
pkgs = c(pkgs.a, pkgs.priority)
to_install = !pkgs %in% installed.packages()
if(any(to_install)) {
  install.packages(pkgs[to_install])
}
inst = lapply(pkgs, library, character.only = TRUE) # load them


path0 = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/6_Temperament_exploration/"

path_output = paste0(path0, "output_exploratory_analysis/")
if(!dir.exists(path_output)){dir.create(path_output)}

setwd(path0)

path_grid = paste0(path0, "grid_gapsMetrics/")
path_metrics = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/4_Plots_metrics/"

# Fonctions ---------------------------------------------------------------

create_files_path_list = function(plot_name, params) {
  library(purrr)
  library(stringr)
  library(dplyr)
  
  files_paths = map2(params$base_path, params$to_grepl, ~list.dirs(path = .x, full.names = TRUE, recursive = TRUE)) %>%
    map2(., params$to_grepl, ~.[grepl(.y, .)]) %>%
    imap(~map(.x, function(folder) {
      sapply(plot_name, function(name) {
        list.files(path = folder, pattern = paste0(name, params$file_extension[.y]), full.names = TRUE)
      }, simplify = FALSE)
    })) %>%
    set_names(params$list_name) %>%
    map(function(files_list) {
      map(plot_name, function(name) {
        files = unlist(files_list) # Convertit la liste en vecteur de caractères
        matched_files = files[str_detect(files, name)]
        if (length(matched_files) > 0) matched_files else NA
      })
    }) %>%
    map(function(list) Filter(function(x) !is.na(x), unlist(list, recursive = FALSE))) %>%
    map(unique) # Supprime les doublons
  
  # Réordonner chaque sous-liste selon l'ordre de plot_name
  files_paths_ordered = map(files_paths, function(file_list) {
    file_order = sapply(plot_name, function(name) {
      grep(name, file_list)
    }) %>% unlist()
    
    # Utiliser unique pour supprimer les index en double causés par des correspondances multiples
    file_list[unique(file_order)]
  })
  
  return(files_paths_ordered)
}

# Loading raw data -------------------------------------------------------------------------

plot_name = "E:/Arthur/OneDrive - Universite de Liege/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Plots/plots_unique" %>%
  list.files(pattern = "\\.shp$", full.names = TRUE) %>%
  unlist() %>%
  unique() %>%
  map_chr(~sub("^.+/(plot_[^/]+)\\.shp$", "\\1", .))

plots_toExclude = c("plot_JEU2","plot_JEU3","plot_JEU4","plot_JEU5")

if(T){
  params_files = list(
    base_path = paste0(
      "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/",
      c("2_ElevationData/","2_ElevationData/","2_ElevationData/", "3_ITD_ITS_Gaps/", "3_ITD_ITS_Gaps/", "3_ITD_ITS_Gaps/")),
    to_grepl = c("CHM_", "DTM_", "nLAS_filtered_", "Crowns_", "Gaps_", "Trees_"),
    file_extension = c("\\.tif$", "\\.tif$", "\\.las$", "\\.shp$", "\\.shp$", "\\.shp$"),
    list_name = c("list.chm", "list.dtm", "list.nlasFiltered", "list.crowns", "list.gaps", "list.trees")
  )
  
  files_paths = create_files_path_list(plot_name, params_files)
}

plot_name = unique(sub("^.+/(plot_[^/]+)\\.tif$", "\\1", files_paths$list.chm))  %>% .[!. %in% plots_toExclude]

shp_files_paths = "E:/Arthur/OneDrive - Universite de Liege/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Plots/plots_unique" %>%
  list.files(pattern = "\\.shp$", full.names = TRUE) %>%
  unlist()%>%
  .[grepl(paste(plot_name, collapse="|"), .)]

FieldInventories = rio::import( # ATTENTION : bien utiliser la base de données où l'information relative aux tempéraments a été injectée, à partir du code n°4 'FieldDataInventories_Analyse_Temperament.R'
  paste0(path_metrics, "FieldDataInventories_computation_with_temperament.csv"),
  sep = ";",
  dec = ","
) %>% 
  as_tibble() %>%
  filter(plot %in% plot_name & !plot %in% plots_toExclude) # Charger les données d'inventaires pour les plots présents sous les survols LiDAR

PlotMetrics = rio::import(
  paste0(path_metrics, "PlotsMetrics_FieldDataInventories.csv"),
  sep = ";",
  dec = ","
) %>% 
  as_tibble() %>%
  filter(plot %in% plot_name & !plot %in% plots_toExclude) # Charger les métriques pour les plots présents sous les survols LiDAR

PlotMetricsLiDAR = rio::import(
  paste0(path_metrics, "PlotsMetrics_LiDAR.csv"),
  sep = ";",
  dec = ","
) %>% 
  as_tibble() %>%
  filter(plot %in% plot_name & !plot %in% plots_toExclude) # Charger les métriques pour les plots présents sous les survols LiDAR

PlotTemperament_all = rio::import("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/4_Plots_metrics/PlotsMetrics_Temperament_all.csv", sep = ';', dec = ',') %>%
  as_tibble() %>%
  filter(plot %in% plot_name & !plot %in% plots_toExclude)

PlotTemperament_breaks = rio::import("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/4_Plots_metrics/PlotsMetrics_Temperament_all_withBreaks.csv", sep = ';', dec = ',') %>%
  as_tibble() %>%
  filter(plot %in% plot_name & !plot %in% plots_toExclude)

#save.image("")

### Répartition tempérament par classe -------------------------------

PlotTemperament_breaks %>%
  filter(plot == "plot_123") %>%
  mutate(class_breaks = factor(class_breaks, levels = c("[-Inf,10]", "[-Inf,20]", "(10,15]", "(15,25]", "(20,40]", "(25, Inf]", "(40, Inf]"))) %>%
  ggplot(aes(x = class_breaks, y = Proportion, fill = Temperament)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Type_subdivision + class_type + Type_Proportion  , scales = "free", ncol = 4) +
  theme_minimal() +
  labs(title = "Cumulative Proportion by Temperament across class_breaks", 
       y = "Total Proportion", 
       x = "Class Breaks")

PlotTemperament_breaks %>%
  mutate(class_breaks = factor(class_breaks, levels = c("[-Inf,10]", "[-Inf,20]", "(10,15]", "(15,25]", "(20,40]", "(25, Inf]", "(40, Inf]"))) %>%
  ggplot(aes(x = class_breaks, y = Proportion, fill = Temperament)) +
  geom_boxplot(outlier.shape = NA) + # Boxplot for variability
  facet_wrap(~Type_subdivision + class_type + Type_Proportion, scales = "free", ncol = 4) +
  theme_minimal() +
  labs(title = "Variability of Proportion by Temperament across class_breaks", 
       y = "Proportion", 
       x = "Class Breaks") +
  coord_cartesian(ylim = c(0, 100)) # Ensuring proportions are between 0 and 100

plots_list = lapply(unique(PlotTemperament_breaks$plot), function(plot_name) {
  data_plot = PlotTemperament_breaks %>%
    filter(plot == plot_name)
  
  data_plot %>%
    mutate(class_breaks = factor(class_breaks, levels = c("[-Inf,10]", "[-Inf,20]", "(10,15]", "(15,25]", "(20,40]", "(25, Inf]", "(40, Inf]"))) %>%
    ggplot(aes(x = class_breaks, y = Proportion, fill = Temperament)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~Type_subdivision + class_type + Type_Proportion, scales = "free", ncol = 4) +
    theme_minimal() +
    labs(title = paste("Cumulative Proportion by Temperament across class_breaks for", plot_name), 
         y = "Total Proportion", 
         x = "Class Breaks")
}) %>% map(~print(.))

### Analyse multivariée position des plots par proportion de tempérament  --------

data_pca = PlotTemperament_breaks %>%
  filter(Type_subdivision == 3, class_type == "dbh", Type_Proportion == "ind.rel") %>%
  mutate(Temperament_classDBH = paste(Temperament, class_breaks)) %>%
  select(-Type_Proportion, -class_type, -Type_subdivision, -Temperament, -class_breaks) %>%
  pivot_wider(names_from = Temperament_classDBH, values_from = Proportion, values_fill = list(Proportion = 0)) %>%
  left_join(PlotMetrics %>% select(plot, WD.BA), by = "plot")

# PCA analysis
pca_result = dudi.pca(data_pca %>% select(-plot) %>% mutate(across(everything(), as.numeric)), scale = TRUE, scannf = FALSE, nf = 5)

fviz_eig(pca_result)
fviz_pca_ind(pca_result,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(pca_result,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_biplot(pca_result, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val

# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

### WD.BA en fonction des pionniers par breaks ---------------------------

data = PlotTemperament_breaks %>%
  filter(Type_Proportion == "G.rel",
         class_type == "dbh",
         Type_subdivision == 3,
         Temperament == "héliophile") %>%
  select(plot, class_breaks, Proportion) %>%
  pivot_wider(names_from = class_breaks, values_from = Proportion, values_fill = 0) %>%
  left_join(PlotMetrics %>% select(plot, WD.BA))

plot_h = ggplot(data) +
  geom_point(aes(x = `[-Inf,20]`, y = WD.BA, color = "`[-Inf,20]`")) +
  geom_smooth(aes(x = `[-Inf,20]`, y = WD.BA, color = "`[-Inf,20]`"), method = "lm", se = FALSE) +
  geom_point(aes(x = `(20,40]`, y = WD.BA, color = "`(20,40]`")) +
  geom_smooth(aes(x = `(20,40]`, y = WD.BA, color = "`(20,40]`"), method = "lm", se = FALSE) +
  geom_point(aes(x = `(40, Inf]`, y = WD.BA, color = "`(40, Inf]`")) +
  geom_smooth(aes(x = `(40, Inf]`, y = WD.BA, color = "`(40, Inf]`"), method = "lm", se = FALSE) +
  scale_color_manual(values = c("`[-Inf,20]`" = "blue", "`(20,40]`" = "red", "`(40, Inf]`" = "green"), name = "Pioneer Percent Range") +
  scale_x_continuous(name = "Pioneer Percent", limits = c(min(data$`[-Inf,20]`, data$`(20,40]`, data$`(40, Inf]`, na.rm = TRUE), max(data$`[-Inf,20]`, data$`(20,40]`, data$`(40, Inf]`, na.rm = TRUE))) +
  labs(y = "WD.BA", title = "Scatter plots with regression lines for different Pioneer Percent ranges") +
  theme_minimal()

plot_dbh = ggplot(data) +
  geom_point(aes(x = `[-Inf,20]`, y = WD.BA, color = "`[-Inf,20]`")) +
  geom_smooth(aes(x = `[-Inf,20]`, y = WD.BA, color = "`[-Inf,20]`"), method = "lm", se = FALSE) +
  geom_point(aes(x = `(20,40]`, y = WD.BA, color = "`(20,40]`")) +
  geom_smooth(aes(x = `(20,40]`, y = WD.BA, color = "`(20,40]`"), method = "lm", se = FALSE) +
  geom_point(aes(x = `(40, Inf]`, y = WD.BA, color = "`(40, Inf]`")) +
  geom_smooth(aes(x = `(40, Inf]`, y = WD.BA, color = "`(40, Inf]`"), method = "lm", se = FALSE) +
  scale_color_manual(values = c("`[-Inf,20]`" = "blue", "`(20,40]`" = "red", "`(40, Inf]`" = "green"), name = "Pioneer Percent Range") +
  scale_x_continuous(name = "Pioneer Percent", limits = c(min(data$`[-Inf,20]`, data$`(20,40]`, data$`(40, Inf]`, na.rm = TRUE), max(data$`[-Inf,20]`, data$`(20,40]`, data$`(40, Inf]`, na.rm = TRUE))) +
  labs(y = "WD.BA", title = "Scatter plots with regression lines for different Pioneer Percent ranges") +
  theme_minimal()
