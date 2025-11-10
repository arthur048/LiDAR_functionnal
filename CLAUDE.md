# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a doctoral research project analyzing forest structure using LiDAR data and field inventories in Central African forests. The workflow processes LiDAR point clouds to identify canopy gaps and correlates them with tree species traits (particularly wood density and species temperament/light requirements) to understand relationships between forest vertical structure and functional composition.

## Environment Setup

### R Environment
- R version: 4.4.0
- Package management: `renv` (lockfile present)
- To restore the R environment: `Rscript -e "renv::restore()"`
- Project uses `.Rproj.user` for RStudio configuration

### Key R Packages
Core dependencies (see renv.lock for complete list):
- Spatial: `sf`, `terra`, `lidR`
- Data manipulation: `tidyverse`, `rio`
- Forest analysis: `BIOMASS` (taxonomic correction, height models, AGB estimation)
- Performance: `foreach`, `doParallel`

## Project Workflow Architecture

The analysis follows a numbered pipeline structure (0-5), where each step depends on outputs from previous steps:

### 0_Inventories_plot_preparation/
Prepares field plot data and spatial boundaries.
- Input: Shapefiles from AfriSAR campaigns, field inventory Excel/CSV files
- Output: `final/plots_inventories.csv`, `final/plots_info.csv`, `final/plots_limites.gpkg`
- Key script: `plots_preparation.R`
- Handles multiple data sources (Rabi, Mondah, Mabounie, Lope, Yangambi, Malebo)
- Assigns UTM zones and standardizes plot references

### 1_LiDAR_extraction/
Extracts LiDAR point clouds for each plot from raw LiDAR data.
- Input: Raw LAS files from NASA 2015 campaign, plot boundaries from step 0
- Output: Clipped LAS files per plot in `output/`
- Key script: `LiDAR_extraction.R`
- Uses parallel processing (`set_lidr_threads()`)
- Distinguishes between transect plots and ferry line plots

### 2_ElevationData/
Generates Canopy Height Models (CHM) from LiDAR point clouds.
- Input: Clipped LAS files from step 1
- Output: CHM rasters (.tif) in `CHM_final/`
- Multiple scripts for different data sources (NASA2015, AfriSAR)
- Key script: `CHM_final_selection.R` selects best CHM per plot

### 3_Gaps/
Identifies and characterizes canopy gaps from CHM.
- Input: CHM rasters from step 2
- Output: Gap rasters and metrics in `output/`
- Key scripts:
  - `1_Gaps_identification.R`: Detects gaps at multiple height thresholds (1-45m) using custom `get_gap_layer()` function
  - `2_Gaps_data_frequencies.R`: Gap size-frequency distributions
  - `3_Gaps_metrics.R`: Calculates gap metrics (area, perimeter, shape indices)
- Uses parallel processing for multiple plots/thresholds

### 4_Plots_metrics/
Computes plot-level metrics from field inventories and gaps.
- Input: Field inventories from step 0, gap metrics from step 3, CoForTraits database
- Output: Multiple CSV files with plot-level aggregations
- Key scripts:
  - `1_FieldDataInventories_computation.R`:
    - Taxonomic correction via `BIOMASS::correctTaxo()`
    - Wood density assignment (prioritizes CoForTraits local data over GWDD)
    - Height estimation for missing values
    - Calculates basal area, AGB per tree
  - `2_PlotsMetrics_FieldDataInventries.R`: Aggregates to plot level
  - `3_PlotsMetrics_Gaps.R`: Aggregates gap metrics per plot
  - `4_FieldDataInventories_Analyse_Temperament.R`: Assigns species temperament (shade-tolerant/NPLD/pioneer)
  - `5_PlotsMetrics_Analyse_Temperament.R`: Calculates temperament-based metrics per plot
- Critical data files:
  - `cofortraits.csv`: Local wood density database (root directory)
  - `meanWD_by_species_CoFor.csv`: Species-averaged wood densities

### 5_Main_analysis/
Main statistical analyses and manuscript figures.
- Input: All plot metrics from step 4
- Key script: `1_data_preparation.R` merges all datasets
- Structure:
  - `functions/`: Custom analysis functions
  - `input/`: Consolidated input data
  - `output/data/`: Analysis results
  - `output/figures/`: Publication-ready figures
  - `output/tables/`: Results tables

## Important Data Handling Patterns

### File Paths
Scripts use absolute paths pointing to `E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/`. When modifying paths, maintain consistency across dependent scripts.

### CSV Encoding
- Delimiter: `;` (semicolon)
- Decimal separator: `,` (comma - European convention)
- Use `read_csv2()` or `rio::import(sep=";", dec=",")` for reading
- Encoding: UTF-8 for species names with special characters

### Spatial Data
- UTM projections vary by site (UTM 32-36, North/South)
- EPSG codes stored in named vectors (e.g., `utm_epsg` in scripts)
- Always check CRS consistency when combining spatial datasets
- Use `st_make_valid()` for geometry issues

### Parallel Processing
Most computationally intensive scripts use parallel processing:
```r
n_cores <- parallel::detectCores() - 2  # Leave cores for system
cl <- makeCluster(n_cores)
registerDoParallel(cl)
# ... foreach loop
stopCluster(cl)  # Always close clusters
```

## Custom Functions

### Gap Analysis (`3_Gaps/1_Gaps_identification.R`)
- `get_gap_layer(chm_layer, threshold, size)`: Core gap detection function
  - Reclassifies CHM below threshold as gaps
  - Uses 8-direction connectivity (`terra::patches(directions=8)`)
  - Filters by minimum/maximum area
  - Returns SpatRaster with unique gap IDs

### Species Temperament
Species classified into three groups based on light requirements:
- "héliophile" (pioneer): Light-demanding pioneers
- "héliophile non pionnière" (NPLD): Non-pioneer light-demanders
- "tolérante à l'ombrage" (shade-tolerant): Shade-tolerant species

Proportions calculated by basal area (`G.rel`) and stem count (`ind.rel`)

## Common Development Commands

### Running Individual Scripts
```r
# From R console
source("path/to/script.R")
```

### Updating Dependencies
```r
# Snapshot current package versions
renv::snapshot()

# Update packages
renv::update()
```

### Spatial Data Inspection
```r
# Quick plot inspection
library(terra)
r <- rast("path/to/file.tif")
plot(r)

library(sf)
s <- st_read("path/to/file.gpkg")
plot(st_geometry(s))
```

## Critical Data Dependencies

### External Databases
- **CoForTraits** (`cofortraits.csv`): Local trait database for Central African species - ESSENTIAL for accurate wood density values
- **BIOMASS package databases**: Taxonomic reference and pantropical allometric equations

### Plot Reference System
- Plot names follow site-specific conventions (e.g., MAB01h, Betamba_1)
- Master list in `final_plot_name.csv` at root
- Always filter datasets by this list for consistency across analyses

## Known Constraints

1. **Memory intensive**: LiDAR processing requires substantial RAM (45 height thresholds × n plots × gap detection)
2. **Processing time**: Full gap detection pipeline can take hours even with parallelization
3. **Data size**: LAS files and rasters excluded from git (see .gitignore)
4. **Path dependencies**: Scripts assume specific directory structure - moving files requires updating multiple paths
5. **Language mixing**: Comments and variable names mix French and English - species temperament categories are in French

## Workflow Execution Order

To reproduce the full analysis from scratch:
1. Run `0_Inventories_plot_preparation/plots_preparation.R`
2. Run `1_LiDAR_extraction/LiDAR_extraction.R` (requires raw LAS files)
3. Run `2_ElevationData/CHM_final_selection.R`
4. Run `3_Gaps/1_Gaps_identification.R` through `3_Gaps/3_Gaps_metrics.R`
5. Run all scripts in `4_Plots_metrics/` sequentially (1-5)
6. Run analyses in `5_Main_analysis/`

Each step generates outputs required by subsequent steps - do not skip stages.
