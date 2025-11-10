rm(list=ls())
gc()

create_dir = function(path) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}

# Initialisation et path preparation ----------------------------------------------------------

pkgs = c("terra", "tidyverse", "sf", "rio", "caret", "progress", "viridis", "sysfonts", "showtext")

to_install = !pkgs %in% installed.packages() ; if(any(to_install)) {install.packages(pkgs[to_install])} ; inst = lapply(pkgs, library, character.only = TRUE) # load them

font_paths(new = "E:/Arthur/OneDrive2/R"); font_add("LM12", regular = "lmroman12-regular.otf", bold = "lmroman12-bold.otf") ; showtext_auto()

path0 = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/5_Sensibility_Analysis_Modelisation/second_round_buffer50_more_precise/" ; setwd(path0)

path_figures = paste0("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/Figures/5_Sensibility_Analysis_Modelisation/") ; lapply(path_figures, create_dir)

path_results = paste0(path0,"output/")
path_gaps = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/3_Gaps/output/"

gaps_freq = rio::import(paste0(path_gaps, "results_gaps_freq.csv"), sep = ";", dec = ",") %>% as_tibble()
results_sensibility = rio::import(paste0(path_results, "results_sensibility_analysis_modelisation.csv"), sep = ";", dec = ",") %>%
  as_tibble() %>%
  mutate(
    RMSE = case_when(
      respons == "WD.BA" ~ RMSE,
      TRUE ~ RMSE/100 
    )
  )

plot_to_exclude = c("plot_JEU2","plot_JEU3","plot_JEU4","plot_JEU5")
plot_name = gaps_freq %>% filter(!plot_name %in% plot_to_exclude) %>% pull(plot_name) %>% unique()


# Exploring sensibility analysis ------------------------------------------

wd_ba = results_sensibility %>% filter(respons == "WD.BA", !is.na(RMSE), predictors != "proportion", buffer == 50)
helio = results_sensibility %>% filter(respons == "helio")
ombre = results_sensibility %>% filter(respons == "ombre")
helio_np = results_sensibility %>% filter(respons == "helio.np")

# Graphics sensibility analysis -------------------------------------------------------------------------

#resp = "WD.BA"  # "helio"    "ombre"    "helio.np"
resp = "helio"

pred = "proportion + lambda" 
buf = 50
max_min_gap_size = 50
resp_letters <- setNames(c("(A)", "(B)", "(C)", "(D)"), c("WD.BA", "helio", "helio.np", "ombre"))

data <- results_sensibility %>%
  filter(respons == resp, predictors == pred, buffer == buf, xmin <= max_min_gap_size) %>%
  mutate(
    `1 - RMSE` = 1 - RMSE,
    `1 - rRMSE` = 1 - rRMSE,
    resp_letter = resp_letters[respons]
  )

ggplot(data, aes(x = xmin, y = height_aboveground, size = `1 - RMSE`, color = `1 - RMSE`)) +
  geom_point(data = data %>% filter(P_Value > 0.05), aes(shape = "significant"), color = "red", size = 2, alpha = 0.6, stroke = 1) +
  geom_point(alpha = 0.7) +  # Base layer of points with original formatting
  geom_text(aes(label = unique(resp_letter), family = "LM12"), x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 6, color = "black", check_overlap = TRUE) +
  scale_shape_manual(values = c(significant = 4), labels = c("P-value > 0.05")) +
  guides(shape = guide_legend(title = NULL)) +
  scale_color_viridis(limits = c(min(data$`1 - RMSE`), max(data$`1 - RMSE`))) +
  scale_size_continuous(
    range = c(1, 4), 
    limits = c(min(data$`1 - RMSE`), max(data$`1 - RMSE`)),
    guide = guide_legend(reverse = TRUE)
  ) +
  theme_classic(base_family = "LM12") +
  labs(x = "Minimum Gap Size (m²)", y = "Height Threshold (m)") +
  scale_x_continuous(breaks = c(1, seq(5, max_min_gap_size, by = 5)), limits = c(1, max_min_gap_size)) +
  scale_y_continuous(breaks = seq(1, 25, by = 1), limits = c(1, 25)) +
  theme(
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    axis.line = element_line(color = "black", linewidth = 0.4)
  )


# Save the plot
ggsave(paste0(path_figures, resp, "_sensibility_analysis_gaps_definition.svg"),
       plot = last_plot(),
       device = "svg",
       width = 24,
       height = 12,
       units = "cm")


# Plots for IUFRO Poster --------------------------------------------------

resp = "WD.BA"
pred = "proportion + lambda" 
buf = 50
max_min_gap_size = 50

data <- results_sensibility %>%
  filter(respons == resp, predictors == pred, buffer == buf, xmin <= max_min_gap_size) %>%
  mutate(RMSE = -RMSE)

ggplot(data, aes(x = xmin, y = height_aboveground, size = RMSE, color = RMSE)) +
  geom_point(alpha = 1) +  # Base layer of points with original formatting
  geom_point(data = data %>% filter(P_Value > 0.05), aes(shape = "significant"), color = "red", size = 2, alpha = 1, stroke = 1) +
  scale_shape_manual(values = c(significant = 4), labels = c("P-value > 0.05")) +
  guides(shape = guide_legend(title = NULL)) +
  scale_color_viridis(
    limits = c(-max(data$RMSE), -min(data$RMSE)),
    labels = function(x) -x,
    guide = guide_colorbar(title = "Model RMSE (g/cm³)")
  ) +
  scale_size_continuous(
    range = c(1, 4), 
    limits = c(-max(data$RMSE), -min(data$RMSE)),
    guide = guide_legend(title = "Model RMSE (g/cm³)", reverse = TRUE),
    labels = function(x) -x
  ) +
  theme_minimal() +
  labs(x = "Minimum Gap Size (m²)",
       y = "Height Aboveground (m)",
       title = "Canopy Gaps Definition for Wood Density Estimations ") +
  scale_x_continuous(breaks = c(1, seq(5, max_min_gap_size, by = 5)), limits = c(1, max_min_gap_size)) +
  scale_y_continuous(breaks = seq(4, 20, by = 1), limits = c(4, 20)) +
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
  )


# Save the plot
ggsave(paste0(path_figures, resp, "_sensibility_analysis_gaps_definition_IUFRO.svg"),
       plot = last_plot(),
       device = "svg",
       width = 24,
       height = 12,
       units = "cm")
ggsave(paste0(path_figures, resp, "_sensibility_analysis_gaps_definition_IUFRO.jpeg"),
       plot = last_plot(),
       device = "jpeg",
       width = 24,
       height = 12,
       units = "cm")


