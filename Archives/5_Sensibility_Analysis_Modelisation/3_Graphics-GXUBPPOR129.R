rm(list=ls())
gc()

create_dir = function(path) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}

# Initialisation et path preparation ----------------------------------------------------------

pkgs = c("terra", "tidyverse", "sf", "rio", "caret", "progress", "viridis")

to_install = !pkgs %in% installed.packages() ; if(any(to_install)) {install.packages(pkgs[to_install])} ; inst = lapply(pkgs, library, character.only = TRUE) # load them

path0 = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/5_Sensibility_Analysis_Modelisation/" ; setwd(path0)

path_figures = paste0("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/Figures/Sensibility_Analysis_Modelisation/") ; lapply(path_figures, create_dir)

path_results = paste0(path0,"output/")
path_gaps = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/3_Gaps/output/"

gaps_freq = rio::import(paste0(path_gaps, "results_gaps_freq.csv"), sep = ";", dec = ",") %>% as_tibble()
results_sensibility = rio::import(paste0(path_results, "results_sensibility_analysis_modelisation.csv"), sep = ";", dec = ",") %>% as_tibble()

plot_to_exclude = c("plot_JEU2","plot_JEU3","plot_JEU4","plot_JEU5")
plot_name = gaps_freq %>% filter(!plot_name %in% plot_to_exclude) %>% pull(plot_name) %>% unique()


# Exploring sensibility analysis ------------------------------------------

wd_ba = results_sensibility %>% filter(respons == "WD.BA", !is.na(RMSE), predictors != "proportion", buffer == 50)
helio = results_sensibility %>% filter(respons == "helio")
ombre = results_sensibility %>% filter(respons == "ombre")
helio_np = results_sensibility %>% filter(respons == "helio.np")

# Graphics sensibility analysis -------------------------------------------------------------------------

resp = "WD.BA"  # "helio"    "ombre"    "helio.np"
#resp = "helio"

pred = "proportion + lambda" 
# pred = "proportion"   

buf = 50 #  0  20  50 100 200

data = results_sensibility %>%
  filter(respons == resp, predictors == pred, buffer == buf)

ggplot(data, aes(x = xmin, y = height_aboveground, size = RMSE)) +
  geom_point(aes(color = RMSE), alpha = 0.7) +
  scale_color_viridis(limits = c(min(data$RMSE), max(data$RMSE))) +
  scale_size_continuous(range = c(2, 10), limits = c(min(data$RMSE), max(data$RMSE))) +
  theme_minimal() +
  labs(x = "Minimum Gap Size (m²)", y = "Height Threshold (m)")
