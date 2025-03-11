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
    plots_afrisar_wd %>%
      rename(WD_mean = WD) %>%
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

rio::export(
  plots_metrics %>% left_join(plots_metrics_20),
  paste0(path0, "PlotsMetrics_FieldDataInventories.csv"),
  sep = ";",
  dec = ","
)
