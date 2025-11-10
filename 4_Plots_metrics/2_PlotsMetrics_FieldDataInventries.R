rm(list=ls())
gc()

# Initialisation ----------------------------------------------------------
pkgs <- c("tidyverse", "rio", "BIOMASS", "foreach", "doParallel", "sf", "terra", "lidR")

to_install <- !pkgs %in% installed.packages()
if(any(to_install)) install.packages(pkgs[to_install])
inst <- lapply(pkgs, library, character.only = TRUE)

# Définition du chemin de base
path0 <- "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/4_Plots_metrics/"

# Loading raw data -------------------------------------------------------------------------
plot_name = read.csv2("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/final_plot_name.csv") %>% pull(plot_name)

field_inventories = rio::import(
  paste0(path0, "FieldDataInventories_computation.csv"),
  sep = ";",
  dec = ","
) %>%
  as_tibble()

plots_info <- rio::import(
  "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/final/plots_info.csv",
  dec = ",", sep = ";") %>%
  as_tibble() %>%
  filter(plot_ref %in% plot_name)

plots_afrisar_wd <- rio::import("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/final/plots_AfriSAR_WD.csv", 
                                dec = ",", sep = ";") %>%
  as_tibble() %>%
  dplyr::filter(plot_ref %in% plot_name)

# Plots field metrics -------------------------------------------------
#  --> attention aux outliers qui ne sont pas filtrés

plots_metrics = field_inventories %>%
  group_by(plot_ref, site) %>%
  summarise(
    H_mean = mean(h),
    H_95 = quantile(h, 0.95),
    H_Lor = sum(G_rel * h) / sum(G_rel),
    n_tree = n(),
    DBH_mean = mean(dbh),
    QMD = sqrt(sum(dbh^2) / n_tree),
    G = sum(G),
    WD_BA = sum(WD_pond),
    WD_mean = mean(meanWD),
    AGB_Chave = sum(AGB_Chave),
    .groups = "drop"
  ) %>%
  bind_rows(
    # Filter plots_afrisar_wd to include only plot_ref values not in field_inventories
    plots_afrisar_wd %>%
      rename(WD_mean = WD) %>%
      # This anti_join keeps only rows with plot_ref values not in field_inventories
      anti_join(field_inventories, by = "plot_ref") %>%
      mutate(
        site = NA_character_,
        H_mean = NA_real_,
        H_95 = NA_real_,
        H_Lor = NA_real_,
        n_tree = NA_integer_,
        DBH_mean = NA_real_,
        QMD = NA_real_,
        G = NA_real_,
        WD_BA = NA_real_,
        AGB_Chave = NA_real_
      )
  )

# Plots top 20 field metrics ------------------------------------------------

plots_metrics_20 = field_inventories %>%
  as_tibble() %>%
  filter(h/dbh < 1.2) %>%
  group_by(plot_ref) %>%
  arrange(desc(dbh)) %>%
  slice(1:20) %>%
  summarise(
    H_mean_20 = mean(h),
    WD_mean_20 = mean(meanWD),
    QMD_20 = sqrt(
      sum(dbh ^ 2) / n()
    ),
    AGB_Bastin_20 = 0.0735 * (QMD_20 * H_mean_20 * WD_mean_20) ^ 1.1332
  ) %>%
  ungroup()

plots_metrics_final = plots_metrics %>% 
  left_join(plots_metrics_20)

rio::export(
  plots_metrics_final,
  paste0(path0, "PlotsMetrics_FieldDataInventories.csv"),
  sep = ";",
  dec = ","
)

# comparaison du WD fournit par AfriSAR et celui calculé ici ------------------------------------------------

wd_comparison = plots_afrisar_wd %>%
  select(plot_ref, WD) %>%
  semi_join(field_inventories, by = "plot_ref") %>%
  left_join(plots_metrics_final %>% select(plot_ref, WD_mean, WD_BA)) %>%
  mutate(
    WD_mean_diff = WD_mean - WD,
    WD_BA_diff = WD_BA - WD
  )

if(T){
  library(ggplot2)
  library(patchwork)
  library(gridExtra)
  library(broom)
  library(dplyr)
  
  # Calculate regression statistics for WD vs WD_mean
  lm_mean <- lm(WD_mean ~ WD, data = wd_comparison)
  lm_mean_stats <- glance(lm_mean)
  lm_mean_coef <- tidy(lm_mean)
  mae_mean <- mean(abs(wd_comparison$WD_mean_diff))
  
  # Calculate regression statistics for WD vs WD_BA
  lm_ba <- lm(WD_BA ~ WD, data = wd_comparison)
  lm_ba_stats <- glance(lm_ba)
  lm_ba_coef <- tidy(lm_ba)
  mae_ba <- mean(abs(wd_comparison$WD_BA_diff))
  
  # Determine axis limits
  min_wd <- min(c(wd_comparison$WD, wd_comparison$WD_mean, wd_comparison$WD_BA))
  max_wd <- max(c(wd_comparison$WD, wd_comparison$WD_mean, wd_comparison$WD_BA))
  # Round to nearest 0.05 below/above
  min_limit <- floor(min_wd*200)/200
  max_limit <- ceiling(max_wd*200)/200
  # Create sequence by 0.05
  breaks_seq <- seq(min_limit, max_limit, by = 0.02)
  
  # Create annotation text for the plots
  mean_annot <- paste0(
    "y = ", round(lm_mean_coef$estimate[1], 3), 
    " + ", round(lm_mean_coef$estimate[2], 3), "x\n",
    "R² = ", round(lm_mean_stats$r.squared, 3), "\n",
    "MAE = ", round(mae_mean, 3)
  )
  
  ba_annot <- paste0(
    "y = ", round(lm_ba_coef$estimate[1], 3), 
    " + ", round(lm_ba_coef$estimate[2], 3), "x\n",
    "R² = ", round(lm_ba_stats$r.squared, 3), "\n",
    "MAE = ", round(mae_ba, 3)
  )
  
  # Create the first plot: WD vs WD_mean
  p1 <- ggplot(wd_comparison, aes(x = WD, y = WD_mean)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    annotate("text", x = min_limit + 0.05, y = max_limit - 0.05, 
             label = mean_annot, hjust = 0, vjust = 1, size = 4) +
    labs(
      title = "WD vs WD_mean",
      x = "WD (plots_afrisar_wd)",
      y = "WD_mean (plots_metrics_final)"
    ) +
    scale_x_continuous(limits = c(min_limit, max_limit), breaks = breaks_seq) +
    scale_y_continuous(limits = c(min_limit, max_limit), breaks = breaks_seq) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Create the second plot: WD vs WD_BA
  p2 <- ggplot(wd_comparison, aes(x = WD, y = WD_BA)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    annotate("text", x = min_limit + 0.05, y = max_limit - 0.05, 
             label = ba_annot, hjust = 0, vjust = 1, size = 4) +
    labs(
      title = "WD vs WD_BA",
      x = "WD (plots_afrisar_wd)",
      y = "WD_BA (plots_metrics_final)"
    ) +
    scale_x_continuous(limits = c(min_limit, max_limit), breaks = breaks_seq) +
    scale_y_continuous(limits = c(min_limit, max_limit), breaks = breaks_seq) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Combine the plots side by side
  combined_plot <- p1 + p2 + 
    plot_layout(ncol = 2, widths = c(1, 1)) +
    plot_annotation(
      title = "Comparison of Wood Density Measurements",
      subtitle = "Blue line = linear regression, Red dashed line = 1:1 relationship",
      theme = theme(plot.title = element_text(hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5))
    )
}

combined_plot

ggsave(
  filename = "wood_density_comparison.png",
  plot = combined_plot,
  width = 10,  # Width in inches
  height = 5,  # Height in inches
  dpi = 300    # Resolution in dots per inch
)
