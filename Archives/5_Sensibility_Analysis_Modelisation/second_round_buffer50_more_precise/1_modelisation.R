rm(list=ls())
gc()

create_dir = function(path) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}
create_empty_list = function(niveaux, index) {
  if (index > length(niveaux)) {
    return(list())
  } else {
    currentLevel = niveaux[[index]]
    return(lapply(currentLevel, function(level) create_empty_list(niveaux, index + 1)) %>% setNames(currentLevel))
  }
}

# Initialisation et path preparation ----------------------------------------------------------

pkgs = c("terra", "tidyverse", "sf", "rio", "caret", "progress")

to_install = !pkgs %in% installed.packages() ; if(any(to_install)) {install.packages(pkgs[to_install])} ; inst = lapply(pkgs, library, character.only = TRUE) # load them

path0 = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/5_Sensibility_Analysis_Modelisation/second_round_buffer50_more_precise/" ; setwd(path0)
path_elevation = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/2_ElevationData/"
path_gaps = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/3_Gaps/output/"
path_metrics = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/4_Plots_metrics/"

path_output = paste0(path0,"output/") ; lapply(path_output, create_dir)

study_area = c("Malebo", "Yangambi", "YangambiNestor")
study_area_epsg = c("32733", "32635", "32635")

path_chm_folder = paste0(path_elevation, "CHM_", study_area, "/")
path_chm = lapply(path_chm_folder, function(folder) {
  list.files(path = folder, pattern = "\\.tif$", full.names = TRUE)
}) %>% unlist()

plot_to_exclude = c("plot_JEU2","plot_JEU3","plot_JEU4","plot_JEU5")

plot_name_short = path_chm %>% gsub(pattern = ".*plot_|\\.tif$", replacement = "")

plot_name = paste0("plot_", plot_name_short)

plot_to_consider = plot_name[!plot_name %in% plot_to_exclude]

path_plot = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Plots/plots_unique/" %>%
  list.files(pattern = "\\.shp$", full.names = TRUE) %>%
  .[grep(paste(plot_name, collapse = "|"), .)]

# Data reading -------------------------------------------------------------------------

gaps_freq = rio::import(paste0(path_gaps, "results_gaps_freq.csv"), sep = ";", dec = ",") %>% as_tibble()
gaps_metrics = rio::import(paste0(path_gaps, "results_gaps_metrics.csv"), sep = ";", dec = ",") %>% as_tibble()

FieldInventories = rio::import( # ATTENTION : bien utiliser la base de données où l'information relative aux tempéraments a été injectée, à partir du code n°4 'FieldDataInventories_Analyse_Temperament.R'
  paste0(path_metrics, "FieldDataInventories_computation_with_temperament.csv"),
  sep = ";",
  dec = ","
) %>% as_tibble()

PlotMetrics = rio::import(paste0(path_metrics, "PlotsMetrics_FieldDataInventories.csv"),sep = ";",dec = ",") %>% as_tibble()
PlotTemperament = rio::import(paste0(path_metrics, "PlotsMetrics_Temperament_all.csv"), sep = ';', dec = ',') %>% as_tibble() 

# Data preparation for modelisation ---------------------------------------

type_subdivision = 3
type_proportion = "G.rel"

data_temperament = PlotTemperament %>%
  mutate(Temperament = case_when(
    Temperament == "héliophile" ~ "helio",
    Temperament == "tolérante à l'ombrage" ~ "ombre",
    Temperament == "héliophile non pionnière" ~ "helio.np",
    TRUE ~ as.character(Temperament)  # Fallback, au cas où il y aurait d'autres valeurs
  )) %>%
  filter(Type_Proportion == type_proportion, Type_subdivision == type_subdivision) %>%
  pivot_wider(names_from = Temperament, values_from = Proportion)

heights = seq(4,22)
bufs = 50
min_gap_size = seq(1,120)
max_gap_size = 2.5 * 10^3

predictors = c("proportion", "lambda")
respons = c("WD.BA", "helio", "ombre", "helio.np")

data_modelisation = gaps_metrics %>%
  filter(plot_name %in% plot_to_consider, height_aboveground %in% heights, buffer %in% bufs, xmin %in% min_gap_size) %>%
  mutate(has_lambda = !is.na(lambda)) %>%
  left_join(PlotMetrics %>% select(plot_name = plot, any_of(respons))) %>%
  left_join(data_temperament %>% select(plot_name = plot, any_of(respons))) %>%
  mutate(proportion = if_else(is.na(proportion), 0, proportion))

niveaux_model = list(
  niveau1 = respons,
  niveau2 = predictors,
  niveau3 = bufs %>% as.character(),
  niveau4 = heights %>% as.character(),
  niveau5 = min_gap_size %>% as.character()
)

model_list = create_empty_list(niveaux_model, 1)

# Loop for linear model creation ------------------------------------------

z = 0 ; a=1 ; b=2 ; c=3 ; d=20; e = 10 

maxz = length(niveaux_model[[1]])*length(niveaux_model[[3]])*length(niveaux_model[[4]])*length(niveaux_model[[5]])
pb = progress_bar$new(format = "  Progress: [:bar] :percent :elapsed/:eta", total = maxz, clear = FALSE)

for(a in 1:length(niveaux_model[[1]])){
  for (c in 1:length(niveaux_model[[3]])){
    for (d in 1:length(niveaux_model[[4]])){
      for (e in 1:length(niveaux_model[[5]])){
        
        z = z + 1 
        
        data_toUse = data_modelisation %>%
          filter(height_aboveground == niveaux_model[["niveau4"]][[d]], xmin == niveaux_model[["niveau5"]][[e]], buffer == niveaux_model[["niveau3"]][[c]]) %>%
          dplyr::select(respons = niveaux_model[["niveau1"]][a], all_of(predictors))
        
        test_NA = data_toUse %>% 
          dplyr::select(lambda) %>%
          is.na() %>%
          any()
        
        if(test_NA){
          model_list[[a]][["lambda"]][[c]][[d]][[e]] = NA
        } else { 
          model_list[[a]][["lambda"]][[c]][[d]][[e]] = train(respons ~ lambda + proportion, data = data_toUse, method = 'lm', trControl = trainControl(method = "LOOCV"), preProcess = NULL)
        }
        
        pb$tick()  # Update progress bar
      }
    }
  }
}

saveRDS(model_list, file = paste0(path_output, "model.list.rds"))
