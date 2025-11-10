#### 1️⃣ INITIALIZATION ####
# Clean environment and memory
rm(list = ls())
gc()

# Définition des packages nécessaires
pkgs <- c("tidyverse","FactoMineR", "factoextra", "corrplot", "patchwork", "ade4", "vegan", "broom", "here", "caret")

# Automatic installation and loading of missing packages
to_install <- !pkgs %in% installed.packages()
if(any(to_install)) {install.packages(pkgs[to_install])}
inst <- lapply(pkgs, library, character.only = TRUE)

#### 2️⃣ PATH DEFINITION ####
# Helper function to build paths (adapt to your environment)
project_path <- function(...) {
  here(...)
}

# Base directory definition
project_dir <- here()

# Essential paths
path_main <- project_path("5_Main_analysis")
path_input <- file.path(path_main, "input")
path_output <- file.path(path_main, "output")

# Paths for this specific analysis
path_data <- file.path(path_output, "data/modelisation")
path_figures <- file.path(path_output, "figures/modelisation")
path_table <- file.path(path_output, "tables/modelisation")

# Folder creation
dir.create(path_data, recursive = TRUE, showWarnings = FALSE)
dir.create(path_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(path_table, recursive = TRUE, showWarnings = FALSE)

#### 3️⃣ DATA LOADING ####
# Load main data
plot_metrics <- read_csv2(file.path(path_input, "plot_metrics.csv"))
inventory_temperament <- read_csv2(file.path(path_input, "inventory_with_temperament.csv"))

#### 4️⃣ DATA PREPARATION ####
# Filter plots with temperament data
plots_with_temperament <- inventory_temperament %>%
  pull(plot_ref) %>%
  unique()

# Filter plot_metrics to keep only these plots
plot_data_complete <- plot_metrics %>%
  filter(plot_name %in% plots_with_temperament)

# Define columns for analyses
metadata_cols <- c("plot_name", "WD_BA", "prop_g_helio", "prop_g_npld", "prop_g_shade", 
                   "prop_ind_helio", "prop_ind_npld", "prop_ind_shade")

# Canopy structural metrics
gap_metrics_to_remove <- c("basal_thickness", "transition_thickness", "upper_thickness", 
                           "auc_0_25_norm", "auc_25_75_norm", "auc_75_100_norm",
                           "basal_auc_norm", "upper_auc_norm", "transition_auc_norm",
                           "basal_auc", "upper_auc", "transition_auc", 	
                           "height_min_d2", "min_d2", "height_max_d2", "max_d2")

gap_metrics_cols <- setdiff(names(plot_data_complete), c(metadata_cols, gap_metrics_to_remove))

# Metadata with wood density and succession
metadata <- plot_data_complete %>%
  select(all_of(metadata_cols))

plots_name <- metadata %>%
  pull(plot_name)

# 1. Clean taxonomic data
inventory_clean <- inventory_temperament %>%
  filter(!is.na(species_full) &
           species_full != "" &
           species_full != "NA NA" &
           species_full != "sp.2 NA" &
           species_full != "indet indet " &
           species_full != "bois mort" &
           species_full != "Unknown NA"
  )

# 2. Build species × plots contingency table
species_plot_table <- inventory_clean %>%
  count(plot_ref, species_full) %>%
  pivot_wider(names_from = species_full, 
              values_from = n,
              values_fill = 0) %>%
  column_to_rownames("plot_ref")

# 2.a Table espèces × parcelles pondérée par DBH
if(T){
  species_dbh_plot_table <- inventory_clean %>%
    # Somme des sections par espèce et parcelle (DBH en cm)
    group_by(plot_ref, species_full) %>%
    summarise(basal_area = sum(G), .groups = "drop") %>% # Surface terrière en m²
    # Convertir au format large (contingence)
    pivot_wider(names_from = species_full, 
                values_from = basal_area,
                values_fill = 0) %>%
    # Conserver plot_ref comme identifiant de ligne
    column_to_rownames("plot_ref")
  
  genus_plot_table <- inventory_clean %>%
    count(plot_ref, genus) %>%
    pivot_wider(names_from = genus, 
                values_from = n,
                values_fill = 0) %>%
    column_to_rownames("plot_ref")
  
  #Test table avec uniquement le genre
  genus_wd_table <- inventory_clean %>%
    group_by(genus) %>%
    reframe(wd = min(meanWD, na.rm = TRUE)) %>%
    distinct(genus, .keep_all = TRUE) %>%
    column_to_rownames("genus")
}

# 3. Build species x wood density table
species_wd_table <- inventory_clean %>%
  group_by(species_full) %>%
  reframe(wd = min(meanWD, na.rm = TRUE)) %>%
  distinct(species_full, .keep_all = TRUE) %>%
  column_to_rownames("species_full")

# 4. Build structure x plots table

structure_plots_table <- plot_data_complete %>%
  mutate(h100 = ifelse(is.na(h100), 1, h100)) %>%
  select(plot_name , all_of(gap_metrics_cols)) %>%
  column_to_rownames("plot_name")

#### 5️⃣ MULTIVARIATE ANALYSIS (NSCA, CA, PCA) ####

# =============================================================================
# A. NSCA - Analyse de Correspondance Non-Symétrique
# =============================================================================

# Effectuer la NSCA sur les données floristiques
nsca_result <- dudi.nsc(species_plot_table, scannf = FALSE, nf = min(nrow(species_plot_table) - 1, ncol(species_plot_table) - 1))
nsca_dbh_result <- dudi.nsc(species_dbh_plot_table, scannf = FALSE, nf = min(nrow(species_dbh_plot_table) - 1, ncol(species_dbh_plot_table) - 1))
nsca_genus_result <- dudi.nsc(genus_plot_table, scannf = FALSE, nf = min(nrow(genus_plot_table) - 1, ncol(genus_plot_table) - 1))

# Scores des parcelles sur les axes NSCA
nsca_site_scores <- nsca_result$li
nsca_site_scores_dbh <- nsca_dbh_result$li
nsca_site_scores_genus <- nsca_genus_result$li

colnames(nsca_site_scores) <- paste0("NSCA", 1:ncol(nsca_site_scores))
colnames(nsca_site_scores_dbh) <- paste0("NSCA_DBH", 1:ncol(nsca_site_scores_dbh))
colnames(nsca_site_scores_genus) <- paste0("NSCA_GENUS", 1:ncol(nsca_site_scores_genus))

print(summary(nsca_result))
print(summary(nsca_dbh_result))
print(summary(nsca_genus_result))

# =============================================================================
# B. CA - Correspondance Analysis
# =============================================================================

# Effectuer la CA sur les données floristiques
ca_result <- dudi.coa(species_plot_table, scannf = F, nf = min(nrow(species_plot_table) - 1, ncol(species_plot_table) - 1))
ca_genus_result <- dudi.coa(genus_plot_table, scannf = F, nf = min(nrow(genus_plot_table) - 1, ncol(genus_plot_table) - 1))

# Scores des parcelles sur les axes CA
ca_site_scores <- ca_result$li
ca_site_scores_genus <- ca_genus_result$li

colnames(ca_site_scores) <- paste0("CA", 1:ncol(ca_site_scores))
colnames(ca_site_scores_genus) <- paste0("CA_GENUS", 1:ncol(ca_site_scores_genus))

print(summary(ca_result))
print(summary(ca_genus_result))

# =============================================================================
# C. PCA - Principal Component Analysis
# =============================================================================

# Réaliser la PCA sur les données de structure de canopée
pca_result <- dudi.pca(structure_plots_table, scannf = F, nf = min(nrow(structure_plots_table) - 1, ncol(structure_plots_table) - 1))

# Scores des parcelles sur les axes CA
pca_site_scores <- pca_result$li

colnames(pca_site_scores) <- paste0("PCA", 1:ncol(pca_site_scores))

print(summary(pca_result))

fviz_eig(pca_result)

fviz_pca_ind(pca_result,
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(pca_result,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


#### EXPLORATION MODELISATION OF WD ####

# =============================================================================
# A. structure = f(composition) --> PCA = f(NSCA, CA)
# =============================================================================

data_model1 = tibble(
  plot_name = plot_data_complete$plot_name,
  WD = plot_data_complete$WD_BA,
  PCA1 = pca_site_scores$PCA1,
  PCA2 = pca_site_scores$PCA2,
  CA1 = ca_site_scores$CA1,
  CA2 = ca_site_scores$CA2,
)

# =============================================================================
# EXPLORATION
# =============================================================================


hist(data_model1$PCA1, main = "Distribution de PCA1", xlab = "PCA1")
qqnorm(data_model1$PCA1, main = "Q-Q plot PCA1")
qqline(data_model1$PCA1, col = "red")

hist(data_model1$PCA2, main = "Distribution de PCA2", xlab = "PCA2")
qqnorm(data_model1$PCA2, main = "Q-Q plot PCA2")
qqline(data_model1$PCA2, col = "red")

hist(data_model1$CA1, main = "Distribution de CA1", xlab = "CA1")
qqnorm(data_model1$CA1, main = "Q-Q plot CA1")
qqline(data_model1$CA1, col = "red")

hist(data_model1$CA2, main = "Distribution de CA2", xlab = "CA2")
qqnorm(data_model1$CA2, main = "Q-Q plot CA2")
qqline(data_model1$CA2, col = "red")

plot(data_model1$PCA1 ~ data_model1$CA1)
plot(data_model1$PCA1 ~ data_model1$CA1)

plot(data_model1$PCA2 ~ data_model1$CA1)
plot(data_model1$PCA2 ~ data_model1$CA2)

lm = lm(WD ~ CA1 + CA2, data = data_model1)
summary(lm)
plot(WD ~ CA1, data = data_model1)
plot(WD ~ CA2, data = data_model1)

lm = lm(WD ~ PCA1, data = data_model1)
summary(lm)
plot(WD ~ PCA1, data = data_model1)
plot(WD ~ CA2, data = data_model1)

lm = lm(WD_BA ~ prop_at_20m, data = plot_metrics)
summary(lm)

lm = lm(WD ~ PCA1 + PCA2 + CA1 + CA2, data = data_model1)
summary(lm)
plot(lm)

