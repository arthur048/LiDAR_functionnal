rm(list=ls())
gc()

create_dir = function(path) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}

# Initialisation et path preparation ----------------------------------------------------------

pkgs = c("lidR", "terra", "tidyverse", "sf", "rio",  "sysfonts", "showtext", "ggpubr", "viridis",
         "tmap", "tidyterra", "ggspatial", "ggpubr", "cowplot", "gridExtra")

to_install = !pkgs %in% installed.packages() ; if(any(to_install)) {install.packages(pkgs[to_install])} ; inst = lapply(pkgs, library, character.only = TRUE) # load them

font_paths(new = "E:/Arthur/OneDrive2/R"); font_add("LM12", regular = "lmroman12-regular.otf", bold = "lmroman12-bold.otf") ; showtext_auto()

path0 = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/Figures/" ; setwd(path0)
path_out = paste0(path0, "Autres_Figures/")
path_data = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/"

path_plots = paste0(path_data, "0_Plots/plots_unique/")
path_lidar = paste0(path_data, "1_LiDAR_extraction/")
path_rgb = paste0(path_data, "1a_RGB_extraction/RGB/")
path_elevation = paste0(path_data, "2_ElevationData/")
path_gaps = paste0(path_data, "3_Gaps/output/")
path_metrics = paste0(path_data, "4_Plots_metrics/")
path_results_sensibility_analysis = paste0(path_data, "5_Sensibility_Analysis_Modelisation/")

gaps_freq = rio::import(paste0(path_gaps, "results_gaps_freq.csv"), sep = ";", dec = ",") %>% as_tibble()
gaps_metrics = rio::import(paste0(path_gaps, "results_gaps_metrics.csv"), sep = ";", dec = ",") %>% as_tibble()
results_sensibility = rio::import(paste0(path_results_sensibility_analysis, "output/results_sensibility_analysis_modelisation.csv"), sep = ";", dec = ",") %>% as_tibble()
#model_list = readRDS(paste0(path_results_sensibility_analysis, "model.list.rds"))

plot_to_exclude = c("plot_JEU2","plot_JEU3","plot_JEU4","plot_JEU5")
plot_name_full = gaps_freq %>% pull(plot_name) %>% unique()
plot_name = plot_name_full[!plot_name_full %in% plot_to_exclude]
plot_nestor = c("plot_JEU-L2",  "plot_JEU-L6B", "plot_MIX-L6A", "plot_MIX-L7A")

FieldInventories = rio::import( # ATTENTION : bien utiliser la base de données où l'information relative aux tempéraments a été injectée, à partir du code n°4 'FieldDataInventories_Analyse_Temperament.R'
  paste0(path_metrics, "FieldDataInventories_computation_with_temperament.csv"),
  sep = ";",
  dec = ","
) %>% as_tibble()

plot_location = FieldInventories %>%
  select(plot, site) %>%
  unique()

PlotMetrics = rio::import(paste0(path_metrics, "PlotsMetrics_FieldDataInventories.csv"),sep = ";",dec = ",") %>% as_tibble()
PlotLiDARMetrics = rio::import(paste0(path_metrics, "PlotsMetrics_LiDAR.csv"),sep = ";",dec = ",") %>% as_tibble()
PlotTemperament = rio::import(paste0(path_metrics, "PlotsMetrics_Temperament_all.csv"), sep = ';', dec = ',') %>% as_tibble() 

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

data_temperament_ind.rel = PlotTemperament %>%
  mutate(Temperament = case_when(
    Temperament == "héliophile" ~ "helio",
    Temperament == "tolérante à l'ombrage" ~ "ombre",
    Temperament == "héliophile non pionnière" ~ "helio.np",
    TRUE ~ as.character(Temperament)  # Fallback, au cas où il y aurait d'autres valeurs
  )) %>%
  filter(Type_Proportion == "ind.rel", Type_subdivision == type_subdivision) %>%
  pivot_wider(names_from = Temperament, values_from = Proportion)

# Figures relation height - temperament by DBH breaks-------------------------------------

color = c("#008fba", "#6ACC64")

temp = c("tolérante à l'ombrage", "héliophile", "héliophile non pionnière")
level_temp = c("species", "genus", "family")

plot_to_remove_height_allom = paste0("plot_", c("GIL1", "GIL2", "JEU-L2", "JEU-L6B", "JEU4", "MIX-L6A", "MIX-L7A", "MIX2", "MIX6", "JEU2"))

data_graph = FieldInventories %>%
  filter(!plot %in% plot_to_remove_height_allom, level_Temperament.3subdiv %in% level_temp, Temperament.3subdiv %in% temp, dbh > 10) %>%
  select(dbh, h, Temperament.3subdiv) %>%
  mutate(DBH = factor(case_when(
    dbh <= 20 ~ "DBH : 10-20 cm",
    dbh > 20 & dbh <= 50 ~ "DBH : 20-50 cm",
    dbh > 50 ~ "DBH >50 cm"
  ), levels = c("DBH : 10-20 cm", "DBH : 20-50 cm", "DBH >50 cm"))) %>%
  mutate(Temperament.3subdiv = factor(case_when(
    Temperament.3subdiv == "héliophile" ~ "Pioneer Species",
    TRUE ~ "Other Species"
  ), levels = c("Pioneer Species", "Other Species")))


# Creating the boxplots with the adjusted DBH ordering
plot_1 = data_graph %>%
  ggplot(aes(x = Temperament.3subdiv, y = h)) +
  geom_boxplot() +
  facet_wrap(~DBH, scales = "free_x") +
  scale_y_continuous(breaks = seq(0, ceiling(max(data_graph$h) / 5) * 5, by = 5)) +
  #scale_fill_manual(values = color) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, size = 13),
    legend.title = element_text(size = 13, hjust = 0.5),
    legend.text = element_text(size = 11),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) + 
  labs(x = "none", y = "Height Aboveground (m)", title = "")

print(plot_1)

file_path = paste0(path_out, "FigSup_Distribution of Height by Temperament and DBH Category.svg")
ggsave(file_path, 
       plot = plot_1, 
       device = "svg", 
       width = 26, 
       height = 16, 
       units = "cm")

# Figures relation height - temperament by DBH continuous ; -------------------------------------

temp = c("tolérante à l'ombrage", "héliophile", "héliophile non pionnière")
level_temp = c("species", "genus", "family")

if(F){
  #Code pour inspection visuelle et identification des plots dont la hauteur a été dérivée d'une équation allométrique et non mesurée sur le terrain
  ggplot(data_graph, aes(x = dbh, y = h)) +
    geom_point() +
    facet_wrap(~ plot) +
    theme_minimal() +
    labs(title = "Relation entre la hauteur (h) et le diamètre à hauteur de poitrine (dbh) par parcelle",
         x = "DBH (cm)",
         y = "Height (m)")
  
  
  i = 1
  ggplot(FieldInventories %>% filter(plot == unique(FieldInventories$plot)[i]), aes(x = dbh, y = h)) +
    geom_point() +
    theme_minimal() +
    labs(title = unique(FieldInventories$plot)[i],
         x = "DBH (cm)",
         y = "Height (m)")
  i = i + 1
  
}

plot_to_remove_height_allom = paste0("plot_", c("GIL1", "GIL2", "JEU-L2", "JEU-L6B", "JEU4", "MIX-L6A", "MIX-L7A", "MIX2", "MIX6", "JEU2"))

data_graph = FieldInventories %>%
  filter(!plot %in% plot_to_remove_height_allom,level_Temperament.3subdiv %in% level_temp, Temperament.3subdiv %in% temp, dbh > 10) %>%
  select(plot, dbh, h, Temperament.3subdiv)  %>%
  mutate(new_temperament = factor(case_when(
    Temperament.3subdiv == "héliophile" ~ "Pioneer Species",
    TRUE ~ "Other Species"
  ), levels = c("Pioneer Species", "Other Species")))

model1_exp <- lm(log(h) ~ log(dbh), data = data_graph %>% filter(new_temperament == 'Pioneer Species'))
model2_exp <- lm(log(h) ~ log(dbh), data = data_graph %>% filter(new_temperament == 'Other Species'))

model_params_exp = tibble(
  model = c("model1", "model2"),
  intercept = c(exp(coef(model1_exp)[1]), exp(coef(model2_exp)[1])),
  exponent = c(coef(model1_exp)[2], coef(model2_exp)[2]),
  p_value = c(summary(model1_exp)$coefficients[2,4], summary(model2_exp)$coefficients[2,4]),
  adjusted_r_squared = c(summary(model1_exp)$adj.r.squared, summary(model2_exp)$adj.r.squared)
)  %>%
  mutate(
    p_value_text = case_when(
      p_value > 0.05 ~ "> 0.05",
      p_value <= 0.05 & p_value > 0.01 ~ "< 0.05",
      p_value <= 0.01 & p_value > 0.001 ~ "< 0.01",
      p_value <= 0.001 ~ "< 0.001"
    )
  )

min_simulation = 10
max_simulation = 100

data_model_exp = tibble(
  dbh = seq(min_simulation,max_simulation, by = 0.01),
  h_pioneer = model_params_exp[1,]$intercept * (dbh^model_params_exp[1,]$exponent),
  h_others = model_params_exp[2,]$intercept * (dbh^model_params_exp[2,]$exponent),
) %>%
  pivot_longer(
    cols = c(h_pioneer, h_others),
    names_to = "Temperament",
    values_to = "h",
    names_prefix = "h_",
    names_transform = list(Temperament = ~ case_when(
      .x == "pioneer" ~ "Pioneer Species",
      .x == "others" ~ "Other Species"
    ))
  )

plot_2 = ggplot() +
  #geom_point(data = data_graph %>% filter(dbh >= min_simulation, dbh <= max_simulation, h < max(data_model_exp$h)), aes(x = dbh, y = h, shape = new_temperament), size = 1.5, alpha = 0.5) +
  geom_line(data = data_model_exp, aes(x = dbh, y = h, color = Temperament, group = Temperament), linewidth = 1) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, ceiling(max(data_model_exp$dbh) / 10) * 10, by = 10)) +
  scale_y_continuous(breaks = seq(0, ceiling(max(data_model_exp$h) / 5) * 5, by = 5), limit = range(seq(0, ceiling(max(data_model_exp$h) / 5) * 5, by = 5))) +
  scale_shape_manual(values = c(16, 17)) +  # Assigner des formes différentes aux catégories
  scale_color_manual(values = c("Pioneer Species" = "#5ec962", "Other Species" = "#3b528b")) +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0, size = 11),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    legend.position = c(0.75, 0.10),
    legend.justification = c("left", "bottom")
  ) + 
  labs(
    x = "DBH (cm)",
    y = "Height Aboveground (m)",
    title = "Comparison of Height ~ DBH Relationship for Pioneer and Other Species",
    color = "Species Light Requirements"
    #subtitle = bquote(atop(
     # "Pioneer Species: Height = " ~ .(round(model_params_exp$intercept[1], 2)) ~ "*" ~ "DBH"^.(round(model_params_exp$exponent[1], 2)) ~ "; Adj R"^2 ~ "=" ~ .(round(model_params_exp$adjusted_r_squared[1], 2)) ~ ", P-value =" ~ .(model_params_exp$p_value_text[1]),
      #"Other Species: Height = " ~ .(round(model_params_exp$intercept[2], 2)) ~ "*" ~ "DBH"^.(round(model_params_exp$exponent[2], 2)) ~ "; Adj R"^2 ~ "=" ~ .(round(model_params_exp$adjusted_r_squared[2], 2)) ~ ", P-value =" ~ .(model_params_exp$p_value_text[2])
    #))
  )

plot(plot_2)


file_path = paste0(path_out, "FigSup_Distribution of Height and DBH by Temperament Continuous_IUFRO.svg")
ggsave(file_path, 
       plot = plot_2, 
       device = "svg", 
       width = 26, 
       height = 13, 
       units = "cm")


# Define new lighter colors for points
data_graph = data_graph %>%
  mutate(point_color = ifelse(new_temperament == "Pioneer Species", "#8ddc91", "#6a75a6"))

data_model_exp = data_model_exp %>%
  mutate(line_color = ifelse(Temperament == "Pioneer Species", "#5ec962", "#3b528b"))

plot_3 = ggplot() +
  geom_point(data = data_graph %>% filter(dbh >= min_simulation, dbh <= max_simulation, h < max(data_model_exp$h)), 
             aes(x = dbh, y = h, shape = new_temperament, color = point_color), size = 1.5, alpha = 0.8) +
  geom_line(data = data_model_exp, 
            aes(x = dbh, y = h, color = line_color, group = Temperament), linewidth = 1) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(10, ceiling(max(data_model_exp$dbh) / 10) * 10, by = 10)) +
  scale_y_continuous(breaks = seq(0, ceiling(max(data_model_exp$h) / 5) * 5, by = 5), limits = range(seq(0, ceiling(max(data_model_exp$h) / 5) * 5, by = 5))) +
  scale_color_identity() +
  scale_shape_manual(values = c("Pioneer Species" = 17, "Other Species" = 19)) +
  theme(
    text = element_text(color = "black",size = 12),
    axis.text.x = element_text(color = "black",size = 11),
    axis.text.y = element_text(color = "black",size = 11),
    plot.title = element_text(color = "black",hjust = 0.5, size = 13),
    plot.subtitle = element_text(color = "black",hjust = 0, size = 11),
    legend.title = element_text(color = "black",size = 13),
    legend.text = element_text(color = "black",size = 11),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    legend.position = c(0.75, 0.05),
    legend.justification = c("left", "bottom")
  ) +
  labs(
    x = "DBH (cm)",
    y = "Height Aboveground (m)",
    title = "Comparison of Height ~ DBH Relationship for Pioneer and Other Species",
    color = "Species Light Requirements",
    shape = "Species Light Requirements"
  ) +
  guides(
    color = guide_legend(override.aes = list(color = c("#5ec962", "#3b528b"))),
    shape = guide_legend(override.aes = list(shape = c(17, 19), color = c("#5ec962", "#3b528b")))
  )

plot_3


plot(plot_3)

file_path2 = paste0(path_out, "FigSup_Distribution of Height and DBH by Temperament Continuous_withPoints_IUFRO.svg")
ggsave(file_path2, 
       plot = plot_3, 
       device = "svg", 
       width = 26, 
       height = 13, 
       units = "cm")

# Figure relation WD mean ~ WD ba -----------------------------------------

data_graph = PlotMetrics %>%
  filter(plot %in% plot_name) %>%
  select(plot, site, BA = G, WD.BA, WD.mean)

model <- lm(WD.BA ~ WD.mean, data = data_graph)
summary(model)

model_params = tibble(
  intercept = coef(model)[1],
  slope = coef(model)[2],
  p_value = summary(model)$coefficients[2,4],
  adjusted_r_squared = summary(model)$adj.r.squared
)  %>%
  mutate(
    p_value_text = case_when(
      p_value > 0.05 ~ "> 0.05",
      p_value <= 0.05 & p_value > 0.01 ~ "< 0.05",
      p_value <= 0.01 & p_value > 0.001 ~ "< 0.01",
      p_value <= 0.001 ~ "< 0.001"
    )
  )

plot_3 = ggplot(data_graph, aes(x = WD.mean, y = WD.BA)) + 
  geom_point(size = 1.5, alpha = 1, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  scale_x_continuous(breaks = seq(0.50, 0.75, by = 0.05), limit = c(0.53, 0.73)) +
  scale_y_continuous(breaks = seq(0.50, 0.75, by = 0.05), limit = c(0.53, 0.78)) +
  theme_minimal() + 
  theme(
    text = element_text(size = 12),
    plot.subtitle = element_text(hjust = 0, size = 11),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5)
  ) + 
  labs(x = bquote(WD[mean] ~ "(g" ~ cm^-3 ~ ")"),
       y = bquote(WD[BA] ~ "(g" ~ cm^-3 ~ ")"),
       subtitle = bquote(atop(
         WD[BA] ~ "=" ~ .(round(model_params$slope[1], 2)) ~ "*" ~ WD[mean] ~ "; Adj R"^2 ~ "=" ~ .(round(model_params$adjusted_r_squared[1], 2)) ~ ", P-value =" ~ .(model_params$p_value_text[1])
       ))
  )

plot_3

file_path = paste0(path_out, "FigSup_Relation_WDba_WDmean.svg")
ggsave(file_path, 
       plot = plot_3, 
       device = "svg", 
       width = 26, 
       height = 16, 
       units = "cm")

# Figures relation BA ~ Wood Density ------------------------------------

data_graph = PlotMetrics %>%
  filter(plot %in% plot_name) %>%
  select(plot, site, BA = G, WD.BA, WD.mean)

plots_list = list(
  ggplot(data_graph, aes(x = WD.BA, y = BA)) + 
    geom_point() + 
    theme_minimal() + 
    theme(
      text = element_text(size = 12),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks.y = element_line(color = "black", size = 0.5),
      axis.ticks.x = element_line(color = "black", size = 0.5)
    ) + 
    labs(x = bquote(WD[BA] ~ "(g" ~ cm^-3 ~ ")"),
         y = "BA (m²)"
    )
  ,
  
  ggplot(data_graph, aes(x = WD.mean, y = BA)) + 
    geom_point() + 
    theme_minimal() + 
    theme(
      text = element_text(size = 12),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks.y = element_line(color = "black", size = 0.5),
      axis.ticks.x = element_line(color = "black", size = 0.5)
    ) + 
    labs(x = bquote(WD[mean] ~ "(g" ~ cm^-3 ~ ")"),
         y = ""
    )
)

plot4 = cowplot::plot_grid(plotlist = plots_list,
                           ncol = 2, 
                           align = "hv", 
                           axis = "rlbt")

plot4

ggsave(paste0(path_out, "FigSup_BA_WD.svg"), 
       plot = plot4, 
       device = "svg", 
       width = 26, 
       height = 12, 
       units = "cm")




# Figures relation WD / Temperament (plot 1 = BA weighted ; plot 2 = simple proportion et WDmean---------------------------------------

data_graph = PlotMetrics %>%
  filter(plot %in% plot_name)  %>%
  left_join(data_temperament) %>%
  select(plot, site, WD.BA, WD.mean, helio, helio.np, ombre) %>%
  mutate(helio = helio / 100,
         helio.np = helio.np / 100,
         ombre = ombre / 100)

data_graph2 = PlotMetrics %>%
  filter(plot %in% plot_name)  %>%
  left_join(data_temperament_ind.rel) %>%
  select(plot, site, WD.BA, WD.mean, helio, helio.np, ombre) %>%
  mutate(helio = helio / 100,
         helio.np = helio.np / 100,
         ombre = ombre / 100)

model1 <- lm(WD.BA ~ helio, data = data_graph) ; summary(model1)
model2 <- lm(WD.BA ~ helio.np, data = data_graph) ; summary(model2)
model3 <- lm(WD.BA ~ ombre, data = data_graph) ; summary(model3)
model4 <- lm(WD.mean ~ helio, data = data_graph2) ; summary(model4)
model5 <- lm(WD.mean ~ helio.np, data = data_graph2) ; summary(model5)
model6 <- lm(WD.mean ~ ombre, data = data_graph2) ; summary(model6)


model_params = tibble(
  model = c("WD.BA_helio", "WD.BA_helio.np", "WD.BA_ombre", "WD.mean_helio", "WD.mean_helio.np", "WD.mean_ombre"),
  intercept = c(coef(model1)[1], coef(model2)[1], coef(model3)[1], coef(model4)[1], coef(model5)[1], coef(model6)[1]),
  slope = c(coef(model1)[2], coef(model2)[2], coef(model3)[2], coef(model4)[2], coef(model5)[2], coef(model6)[2]),
  p_value = c(summary(model1)$coefficients[2,4], summary(model2)$coefficients[2,4], summary(model3)$coefficients[2,4], summary(model4)$coefficients[2,4], summary(model5)$coefficients[2,4], summary(model6)$coefficients[2,4]),
  adjusted_r_squared = c(summary(model1)$adj.r.squared, summary(model2)$adj.r.squared, summary(model3)$adj.r.squared, summary(model4)$adj.r.squared, summary(model5)$adj.r.squared, summary(model6)$adj.r.squared)
)  %>%
  mutate(
    p_value_text = case_when(
      p_value > 0.05 ~ "> 0.05",
      p_value <= 0.05 & p_value > 0.01 ~ "< 0.05",
      p_value <= 0.01 & p_value > 0.001 ~ "< 0.01",
      p_value <= 0.001 ~ "< 0.001"
    )
  )

plots_list1 = list(
  ggplot(data_graph, aes(x = helio, y = WD.BA)) + 
    geom_point(size = 1.5, alpha = 1, color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
    theme_minimal() + 
    theme(
      text = element_text(size = 12),
      plot.subtitle = element_text(hjust = 0, size = 11),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5)
    ) + 
    labs(x = bquote(Pioneer[BA]),
         y = bquote(WD[BA] ~ "(g" ~ cm^-3 ~ ")"),
         subtitle = bquote(atop(
           WD[BA] ~ "=" ~ .(round(model_params$intercept[1], 2))  ~ .(round(model_params$slope[1], 2)) ~ "*" ~ Pioneer[BA],
           Adj ~ R^2 ~ "=" ~ .(round(model_params$adjusted_r_squared[1], 2)) ~ ", P-value" ~ .(model_params$p_value_text[1])
         ))
    )
  ,
  
  ggplot(data_graph, aes(x = helio.np, y = WD.BA)) + 
    geom_point(size = 1.5, alpha = 1, color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
    theme_minimal() + 
    theme(
      text = element_text(size = 12),
      plot.subtitle = element_text(hjust = 0, size = 11),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5)
    ) + 
    labs(x = bquote(NPLD[BA]),
         y = "",
         subtitle = bquote(atop(
           WD[BA] ~ "=" ~ .(round(model_params$intercept[2], 2)) ~ .(round(model_params$slope[2], 2)) ~ "*" ~ NPLD[BA],
           Adj ~ R^2 ~ "=" ~ .(round(model_params$adjusted_r_squared[2], 2)) ~ ", P-value" ~ .(model_params$p_value_text[2])
         ))
    )
  ,
  
  ggplot(data_graph, aes(x = ombre, y = WD.BA)) + 
    geom_point(size = 1.5, alpha = 1, color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
    theme_minimal() + 
    theme(
      text = element_text(size = 12),
      plot.subtitle = element_text(hjust = 0, size = 11),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5)
    ) + 
    labs(x = bquote(Shade-Tolerant[BA]),
         y = "",
         subtitle = bquote(atop(
           WD[BA] ~ "=" ~ .(round(model_params$intercept[3], 2)) ~ " + " ~ .(round(model_params$slope[3], 2)) ~ "*" ~ Shade-Tolerant[BA],
           Adj ~ R^2 ~ "=" ~ .(round(model_params$adjusted_r_squared[3], 2)) ~ ", P-value" ~ .(model_params$p_value_text[3])
         ))
    )
)

plots_list2 = list(
  ggplot(data_graph2, aes(x = helio, y = WD.mean)) + 
    geom_point(size = 1.5, alpha = 1, color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
    theme_minimal() + 
    theme(
      text = element_text(size = 12),
      plot.subtitle = element_text(hjust = 0, size = 11),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5)
    ) + 
    labs(x = bquote(Pioneer),
         y = bquote(WD[mean] ~ "(g" ~ cm^-3 ~ ")"),
         subtitle = bquote(atop(
           WD[mean] ~ "=" ~ .(round(model_params$intercept[4], 2)) ~ .(round(model_params$slope[4], 2)) ~ "*" ~ Pioneer,
           Adj ~ R^2 ~ "=" ~ .(round(model_params$adjusted_r_squared[4], 2)) ~ ", P-value" ~ .(model_params$p_value_text[4])
         ))
    )
  ,
  
  ggplot(data_graph2, aes(x = helio.np, y = WD.mean)) + 
    geom_point(size = 1.5, alpha = 1, color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
    theme_minimal() + 
    theme(
      text = element_text(size = 12),
      plot.subtitle = element_text(hjust = 0, size = 11),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5)
    ) + 
    labs(x = bquote(NPLD),
         y = "",
         subtitle = bquote(atop(
           WD[mean] ~ "=" ~ .(round(model_params$intercept[5], 2)) ~ .(round(model_params$slope[5], 2)) ~ "*" ~ NPLD,
           Adj ~ R^2 ~ "=" ~ .(round(model_params$adjusted_r_squared[5], 2)) ~ ", P-value" ~ .(model_params$p_value_text[5])
         ))
    )
  ,
  
  ggplot(data_graph2, aes(x = ombre, y = WD.mean)) + 
    geom_point(size = 1.5, alpha = 1, color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
    theme_minimal() + 
    theme(
      text = element_text(size = 12),
      plot.subtitle = element_text(hjust = 0, size = 11),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5)
    ) + 
    labs(x = bquote(Shade-Tolerant),
         y = "",
         subtitle = bquote(atop(
           WD[mean] ~ "=" ~ .(round(model_params$intercept[6], 2)) ~ "+" ~ .(round(model_params$slope[6], 2)) ~ "*" ~ Shade-Tolerant,
           Adj ~ R^2 ~ "=" ~ .(round(model_params$adjusted_r_squared[6], 2)) ~ ", P-value" ~ .(model_params$p_value_text[6])
         ))
    )
)


plot5 = cowplot::plot_grid(plotlist = plots_list1,
                           ncol = 3, 
                           align = "hv", 
                           axis = "rlbt")

plot5

plot6 = cowplot::plot_grid(plotlist = plots_list2,
                           ncol = 3, 
                           align = "hv", 
                           axis = "rlbt")

plot6

ggsave(paste0(path_out, "FigSup_WDBA_Temperament.svg"), 
       plot = plot5, 
       device = "svg", 
       width = 27, 
       height = 10, 
       units = "cm")

ggsave(paste0(path_out, "FigSup_WDmean_Temperament.svg"), 
       plot = plot6, 
       device = "svg", 
       width = 27, 
       height = 10, 
       units = "cm")

# Figures gaps succession profil LiDAR ------------------------------------

if(F){
  
  # Sélection des plots à illustrer M2P14, 151, M2P14, aussi via QGIS inspection visuelle
  
  plot_to_consider = plot_name 
  
  exploration_data = PlotMetrics %>%
    filter(plot %in% plot_to_consider) %>%
    left_join(PlotLiDARMetrics) %>%
    left_join(
      gaps_metrics %>% filter(buffer == 50, height_aboveground == 16, plot_name %in% plot_to_consider, xmin == 10), by = c("plot" = "plot_name")
    ) %>%
    left_join(data_temperament) %>%
    select(plot, H.mean, H.Lor, QMD, WD.BA, WD.mean, AGB_Chave, TCH, proportion, lambda, frequency, mean_area, helio, helio.np, ombre)
  
  
  ggplot(exploration_data, aes(x = WD.BA, y = proportion)) +
    geom_point() +  # Adding points
    geom_text(aes(label = plot), nudge_y = 0.02, check_overlap = TRUE) +  # Adding labels with a slight nudge
    labs(x = "WD.BA", y = "Proportion", title = "Scatter plot of Proportion vs WD.BA") +
    theme_minimal()
  
  ggplot(exploration_data, aes(x = helio, y = proportion)) +
    geom_point() +  # Adding points
    geom_text(aes(label = plot), nudge_y = 0.02, check_overlap = TRUE) +  # Adding labels with a slight nudge
    labs(x = "helio", y = "Proportion", title = "Scatter plot of Proportion vs helio") +
    theme_minimal()
  
  plot(exploration_data$lambda ~ exploration_data$WD.BA)
  lm(exploration_data$WD.BA ~ exploration_data$proportion + exploration_data$lambda) %>% summary()
} 

custom_colors = colorRampPalette(colors = c("#008fba", "#6ACC64", "yellow", "#EE854A", "#720d00"))(64)
colors_0_16 = colorRampPalette(colors = c("#e3abdb", "#4878D0"))(32)
combined_colors = c(colors_0_16, custom_colors)

plots_to_map = c("plot_M2P14", "plot_M2P11", "plot_M2P10")
plot_location_corresponding = plot_location %>%
  filter(plot %in% plots_to_map)

data_annotation = PlotMetrics %>%
  filter(plot %in% plots_to_map) %>%
  left_join(PlotLiDARMetrics) %>%
  left_join(
    gaps_metrics %>% filter(buffer == 50, height_aboveground == 16, plot_name %in% plots_to_map, xmin == 10), by = c("plot" = "plot_name")
  ) %>%
  left_join(data_temperament) %>%
  select(plot, H.mean, H.Lor, QMD, WD.BA, WD.mean, AGB_Chave, TCH, proportion, lambda, frequency, mean_area, helio, helio.np, ombre)

height = 16
min_gap_size = 10

path_chm = paste0(path_elevation, "CHM_", plot_location_corresponding$site, "/", plot_location_corresponding$plot, ".tif")
path_lidar = paste0(path_elevation, "nLAS_filtered_", plot_location_corresponding$site, "/", plot_location_corresponding$plot, ".las")
path_images = paste0(path_rgb, plot_location_corresponding$plot, ".tif")
path_ouverture = paste0(path_gaps, plot_location_corresponding$plot, "/", plot_location_corresponding$plot, "_gaps_height_", height, "_m.tif")
path_shp = paste0(path_plots, plot_location_corresponding$plot, ".shp")

plot_raster_list = list()
plot_lidar_list = list()

for(i in 1:length(plots_to_map)){
  plot_to_consider = plots_to_map[i]
  annotation = data_annotation %>% filter(plot == plot_to_consider)
  
  shp = path_shp[str_detect(path_shp, pattern = plot_to_consider)] %>% sf::st_read() %>% select(geometry)
  shp50 = shp %>% sf::st_buffer(dist = 50)
  
  chm_full = path_chm[str_detect(path_chm, pattern = plot_to_consider)] %>% terra::rast() %>% terra::crop(shp50, mask = TRUE)
  chm_full = ifel(chm_full < 0, 0, chm_full)

  plot_raster = ggplot() + 
    geom_spatraster(data = chm_full) +
    scale_fill_gradientn(colors = custom_colors, breaks = c(16, seq(30,50, by = 10)), limits = c(16, 50), na.value = "white") + 
    geom_sf(data = shp, fill = NA, color = "black", linewidth = 0.5, linetype = "solid", alpha = 1) + 
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.title = element_text(size = 13, hjust = 0.5),
      legend.text = element_text(size = 11),
      legend.title.position = "top"
    ) +
    labs(
      fill = "Height Aboveground (m)"
    )
  
  plot_raster_list[[i]] = plot_raster ; names(plot_raster_list)[[i]] = plot_to_consider
  plot(plot_raster)
  
  lidar = path_lidar[str_detect(path_lidar, pattern = plot_to_consider)] %>% lidR::readLAS() %>% lidR::clip_roi(shp50)
  p1= c(lidar@header$`Min X`, lidar@header$`Min Y`)
  p2= c(lidar@header$`Max X`, lidar@header$`Max Y`)
  transect = clip_transect(lidar, p1, p2, width = 8, xz = TRUE)
  data_plot_lidar = transect@data
  
  max_height_points = data_plot_lidar %>%
    mutate(X_interval = floor(X / 1) * 1) %>%
    group_by(X_interval) %>%
    summarize(Z_max = max(Z)) %>%
    mutate(Z_max = Z_max + 0.2) %>%  # Augmenter les hauteurs de 0,2
    ungroup() %>%
    complete(X_interval = seq(min(X_interval), max(X_interval), by = 0.1)) %>%
    fill(Z_max, .direction = "downup") %>%
    mutate(Z_smooth = zoo::na.approx(Z_max, X_interval, rule = 2))
  
  spline_interpolation = as_tibble(spline(max_height_points$X_interval, max_height_points$Z_smooth, n = 500))
  
  plot_lidar = ggplot(data_plot_lidar, aes(X, Z, color = Z)) + 
    geom_point(size = 0.3) + 
    geom_line(data = spline_interpolation, aes(x = x, y = y), color = "black", linewidth = 0.5, alpha = 1) + 
    geom_hline(yintercept = height, color = "black", linewidth = 0.5, alpha = 1) +
    coord_equal() + 
    theme_minimal() +
    scale_color_gradientn(colors = combined_colors, 
                          name = "Height Aboveground (m)", 
                          breaks = seq(0,50, by = 10), 
                          limits = c(0, 50)) +
    scale_y_continuous(expand = c(0, 0), 
                       breaks = seq(0,50, by = 10), 
                       limits = c(0, 50), 
                       name = "Height Aboveground (m)") +
    scale_x_continuous(expand = c(0,0), breaks = seq(0, 250, by = 50), limits = c(0, 250), name = "Distance (m)") +
    theme(
      text = element_text(size = 12),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, size = 13),
      legend.position = "bottom",
      legend.justification = "center",
      legend.direction = "horizontal",
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks.y = element_line(color = "black", size = 0.5),
      axis.ticks.x = element_line(color = "black", size = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),  # Supprimer les lignes de fond sur l'axe des x
      panel.grid.minor.x = element_blank(),  # Supprimer les lignes de fond sur l'axe des x
      panel.grid.minor.y = element_blank(),
      legend.title = element_text(size = 13, hjust = 0.5),
      legend.text = element_text(size = 11),
      legend.title.position = "top"  
    )  +
    labs(
      fill = "Height Aboveground (m) :",
      title = bquote(atop(
        "Pioneer species" == .(round(annotation$helio, 2)) ~ "%" ~ ", " ~ 
          "WD"[BA] == .(round(annotation$WD.BA, 2)) ~ " g.cm"^"-3" ~ ", " ~ 
          "WD"[mean] == .(round(annotation$WD.mean, 2)) ~ " g.cm"^"-3",
        "Canopy Openings at 16 m" == .(round(annotation$proportion, 2) * 100) ~ "%, " ~ 
          lambda == .(round(annotation$lambda, 2))
      ))
    )
  
  plot_lidar_list[[i]] = plot_lidar ; names(plot_lidar_list)[[i]] = plot_to_consider
  plot(plot_lidar)

}

# Modify plot lists to remove legends, leaving legends in one plot for extraction

plot_a <- plot_lidar_list[[1]] + theme(legend.position = "none", plot.margin = margin(t = 0.1, b = 0.1, unit = "cm"))
plot_b <- plot_lidar_list[[2]] + theme(legend.position = "none", plot.margin = margin(t = 0.1, b = 0.1, unit = "cm"))
plot_c <- plot_lidar_list[[3]]  + theme(legend.position = "none", plot.margin = margin(t = 0.1, b = 0.1, unit = "cm"))
plot_d <- plot_raster_list[[1]] + theme(legend.position = "none", plot.margin = margin(t = 0.1, b = 0.1, unit = "cm"))
plot_e <- plot_raster_list[[2]] + theme(legend.position = "none", plot.margin = margin(t = 0.1, b = 0.1, unit = "cm"))
plot_f <- plot_raster_list[[3]]  + theme(legend.position = "none", plot.margin = margin(t = 0.1, b = 0.1, unit = "cm"))

legend_left <- cowplot::get_plot_component(plot_lidar_list[[3]], 'guide-box-bottom', return_all = TRUE) 
legend_right <- cowplot::get_plot_component(plot_raster_list[[3]], 'guide-box-bottom', return_all = TRUE)

left_column_no_legend <- cowplot::plot_grid(
  plot_a, plot_b, plot_c, 
  ncol = 1,
  align = 'hv', 
  axis = 'rlbt'
)

right_column_no_legend <- cowplot::plot_grid(
  plot_d, plot_e, plot_f, 
  ncol = 1,
  align = 'hv', 
  axis = 'rlbt'
)

combined_plots <- cowplot::plot_grid(
  left_column_no_legend,
  right_column_no_legend,
  ncol = 2,
  rel_widths = c(10, 3),
  align = "hv", 
  axis = "rlbt"
)
  
final_plots = cowplot::plot_grid(
  tmp,
  legend_left,
  ncol = 1,
  nrow = 2,
  rel_heights = c(15,1),
  align = "hv", 
  axis = "rlbt"
)

# Afficher le plot final
print(final_plots)

ggsave(paste0(path_out, "Figure3_Gaps_and_LiDAR_profil_by_WD_BA_and_Forest_Succession.svg"), 
       plot = final_plots, 
       device = "svg", 
       width = 26, 
       height = 26, 
       units = "cm")


