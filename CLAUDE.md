# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a doctoral research project analyzing forest structure using LiDAR data and field inventories in Central African forests. The workflow processes LiDAR point clouds to identify canopy gaps and correlates them with tree species traits (particularly wood density and species temperament/light requirements) to understand relationships between forest vertical structure and functional composition.

## Environment Setup

### R Environment
- R version: 4.4.0
- Package management: `renv` (lockfile present)
- To restore the R environment: `Rscript -e "renv::restore()"`
- **R Project**: `LiDAR_functionnal.Rproj` at project root
- Project uses `.Rproj.user` for RStudio configuration

### Key R Packages
Core dependencies (see renv.lock for complete list):
- Spatial: `sf`, `terra`, `lidR`
- Data manipulation: `tidyverse`, `rio`
- Forest analysis: `BIOMASS` (taxonomic correction, height models, AGB estimation)
- Performance: `foreach`, `doParallel`
- **Path management**: `here` ⭐ (ESSENTIAL for relative paths)

### Path Management with here()
This project uses `here::here()` for all internal paths to ensure portability:
- ✅ **Always use**: `here("module", "subfolder", "file.ext")` for project files
- ❌ **Never hardcode**: Absolute paths like `E:/Arthur/...`
- ⚠️ **Exception**: External data paths (raw LiDAR, original datasets) remain absolute but are clearly marked

**Example**:
```r
# ✅ GOOD - Relative path with here()
plots <- st_read(here("0_Inventories_plot_preparation", "final", "plots_limites.gpkg"))

# ❌ BAD - Hardcoded path
plots <- st_read("E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/0_Inventories_plot_preparation/final/plots_limites.gpkg")

# ⚠️ EXTERNAL DATA - Keep absolute but document clearly
# EXTERNAL DATA PATH: Raw LiDAR data - adjust to your local setup
lidar_base <- "E:/Arthur/Doctorat_DataAnnexe/LiDAR_RDC_2015/00-raw-data"
```

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

## Coding Conventions and Style Guide

### Script Structure (ALWAYS Follow)
Every R script follows this standard structure:

```r
# 1. Environment reset
rm(list=ls())
gc()

# 2. Package initialization
pkgs <- c("tidyverse", "sf", ...)
to_install <- !pkgs %in% installed.packages()
if(any(to_install)) {install.packages(pkgs[to_install])}
inst <- lapply(pkgs, library, character.only = TRUE)

# 3. Path definitions
path_base <- "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles"
path0 <- "..."  # Working directory

# 4. Custom functions (if needed)
my_function <- function(x, y) {
  # Function body
}

# 5. Main code with clear sections
```

### Section Headers

Use consistent hierarchical section markers:

#### Top-level sections (major parts):
```r
# =================================================================
# TITLE IN ALL CAPS
# Description of what this section does
# =================================================================
```

#### Secondary sections with emojis:
```r
#### 1️⃣ INITIALIZATION ####
#### 2️⃣ DATA LOADING ####
#### 3️⃣ ANALYSIS ####
```

#### Subsections:
```r
### 3A. Sub-analysis name ###
```

#### Minor sections:
```r
# Simple comment ----
```

### Naming Conventions

**Variables and functions**: `snake_case`
```r
plot_name <- "MAB01h"
gap_metrics <- calculate_gap_proportion(chm, threshold)
height_aboveground <- 10
```

**Function names**: Descriptive verbs
```r
get_gap_layer()            # Retrieval
calculate_gap_proportion()  # Computation
create_profiles_heatmap()   # Creation
filter_noise()              # Transformation
```

**File names**: Descriptive with underscores
```r
plots_inventories.csv
results_gaps_metrics.csv
PlotsMetrics_FieldDataInventories.csv  # Mixed case for emphasis OK
```

### Code Style

#### Spacing
```r
# YES - spaces around operators
x <- 5 + 3
result <- mean(data, na.rm = TRUE)

# NO - no spaces
x<-5+3
result<-mean(data,na.rm=TRUE)
```

#### Assignment
Use `<-` for assignment (not `=`):
```r
# YES
data <- read_csv2("file.csv")

# NO
data = read_csv2("file.csv")
```

#### Pipe operator
Use `%>%` extensively for data manipulation:
```r
data %>%
  filter(dbh >= 10) %>%
  mutate(ba = pi * (dbh/2)^2) %>%
  group_by(plot_ref) %>%
  summarise(total_ba = sum(ba))
```

#### Line length
Keep lines reasonably short. Break long pipe chains:
```r
# Good
result <- data %>%
  filter(condition1) %>%
  mutate(new_var = calculation) %>%
  group_by(group_var) %>%
  summarise(stat = mean(value))

# Avoid very long single lines
```

### Documentation Style

#### Function documentation (Roxygen-style)
```r
#' Brief description of function
#'
#' Longer description with details about what the function does
#' and when to use it.
#'
#' @param param1 Description of first parameter
#' @param param2 Description of second parameter
#' @return Description of return value
#' @export (if applicable)
my_function <- function(param1, param2) {
  # Implementation
}
```

#### Inline comments
```r
# Descriptive comment explaining WHY, not just WHAT
chm <- classify(chm, rcl = matrix(...))  # Reclassify based on threshold

# For complex logic, explain the reasoning
# We use 8-direction connectivity because corner pixels should be
# considered part of the same gap (follows standard forest gap definition)
gaps <- patches(chm, directions = 8)
```

#### Section comments
Always in French, with clear hierarchy:
```r
# === LECTURE ET PRÉPARATION DES DONNÉES ===
# Cette section charge les données brutes et effectue...

#### 3️⃣ CHARGEMENT DES DONNÉES ####
# Charger les métriques de trouées calculées précédemment
```

### Language Usage

**Comments**: Mix of French and English - French for major sections and explanations, English for code-level comments is acceptable

**Variable names**: Mostly English (`height_aboveground`, `plot_name`) but French is OK for domain-specific terms (`trouées` might appear in comments)

**Function parameter names**: English

**Output messages**: French
```r
cat("Répertoire créé:", dir, "\n")
cat("Traitement de la parcelle:", plot_name, "\n")
```

### Error Handling and Validation

Always validate inputs:
```r
# Check required columns
required_cols <- c("height_aboveground", "proportion", "plot_name")
if (!all(required_cols %in% colnames(data))) {
  stop("Les données doivent contenir les colonnes: ", paste(required_cols, collapse = ", "))
}

# Check for sufficient data
if (sum(valid_data) < 3) {
  warning("Pas assez de données pour le plot ", plot_name)
  return(NULL)
}
```

Use `tryCatch()` for robust parallel processing:
```r
results <- foreach(...) %dopar% {
  tryCatch({
    # Main processing
  }, error = function(e) {
    message("Erreur pour ", plot_name, ": ", e$message)
    return(NULL)
  })
}
```

### Data I/O Patterns

**Reading CSV**:
```r
# ALWAYS use read_csv2() or specify delimiters explicitly
data <- read_csv2("file.csv")  # Assumes sep=";", dec=","

# Or with rio
data <- rio::import("file.csv", sep = ";", dec = ",")
```

**Writing CSV**:
```r
# Use write_csv2() for consistency
write_csv2(data, "output.csv")

# Or with rio
rio::export(data, "output.csv", sep = ";", dec = ",")
```

**Reading spatial data**:
```r
# Use sf for vectors
plots <- st_read("plots.gpkg", quiet = TRUE)  # quiet=TRUE to reduce console output

# Use terra for rasters
chm <- rast("chm.tif")
```

### Parallel Processing Pattern

Standard pattern used throughout:
```r
# Configuration
n_cores <- parallel::detectCores()
num_workers <- max(1, floor(n_cores * 0.75))  # Use 75% of cores

# Initialization
cl <- makeCluster(num_workers)
registerDoParallel(cl)

# Parallel loop
results <- foreach(
  i = seq_along(items),
  .packages = c("terra", "tidyverse", "sf"),  # Specify all needed packages
  .errorhandling = "pass"  # Continue on errors
) %dopar% {
  tryCatch({
    # Processing code
  }, error = function(e) {
    message("Erreur: ", e$message)
    return(NULL)
  })
}

# ALWAYS close cluster
stopCluster(cl)
```

### Progress Reporting

For long-running operations, provide feedback:
```r
cat(glue::glue("Traitement de {length(plot_names)} parcelles...\n"))
cat(glue::glue("Temps écoulé: {round(difftime(Sys.time(), time_start, units='mins'), 1)} min\n"))
```

For iterative processes:
```r
for (i in seq_along(items)) {
  if (i %% 10 == 0) {
    cat(sprintf("Progression: %d/%d (%.1f%%)\n", i, length(items), i/length(items)*100))
  }
  # Processing
}
```

### Graphics and Visualization

#### ggplot2 style
```r
# Standard pattern
p <- ggplot(data, aes(x = var1, y = var2, color = category)) +
  geom_point() +
  scale_color_manual(values = palette) +
  labs(
    title = "Titre descriptif",
    x = "Label axe X",
    y = "Label axe Y",
    color = "Légende"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold")
  )
```

#### Color palettes
Prefer `viridis` for continuous scales (colorblind-friendly):
```r
scale_color_viridis_c(option = "viridis")  # or "plasma", "turbo", etc.
scale_fill_viridis_d(option = "viridis")
```

#### Saving plots
```r
ggsave(
  filename = file.path(output_path, "figure_name.jpg"),
  plot = p,
  width = 12,
  height = 8,
  dpi = 600  # Publication quality
)
```

### Utility Functions

Create helper functions for repeated operations:
```r
# Directory creation
create_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

# Project path helper
project_path <- function(...) {
  file.path(project_dir, ...)
}
```

### Configuration and Parameters

Group parameters at the top of scripts:
```r
#### 4️⃣ PRÉTRAITEMENT ####
# Définition des paramètres
height_aboveground <- seq.int(1, 60)  # Seuils de hauteur à tester
buffers <- c(0, 20, 50, 100, 200)    # Tailles de buffer à appliquer
dpi_value <- 600                      # Résolution des graphiques
top_n_value <- 20                     # Nombre de métriques à visualiser
```

### Best Practices Specific to This Project

1. **Always filter by `final_plot_name.csv`**: Ensure consistency across all analyses
   ```r
   plot_names <- read_csv2("final_plot_name.csv") %>% pull(plot_name)
   data <- data %>% filter(plot_ref %in% plot_names)
   ```

2. **Preserve taxonomic corrections**: Use `useCache = TRUE` with BIOMASS functions
   ```r
   taxo_info <- correctTaxo(genus, species, useCache = TRUE, verbose = FALSE)
   ```

3. **Prioritize CoForTraits for wood density**: Always check local database first

4. **Document buffer choices**: When using spatial buffers, explain why (e.g., edge effects)

5. **Validate CRS consistency**: Before spatial operations
   ```r
   if(st_crs(obj1) != st_crs(obj2)) {
     obj2 <- st_transform(obj2, st_crs(obj1))
   }
   ```

6. **Close parallel clusters**: Memory leaks if forgotten
   ```r
   stopCluster(cl)  # NEVER forget this!
   ```

### Anti-patterns to Avoid

❌ **Don't**: Hardcode plot names directly
```r
plots <- c("MAB01h", "MAB02h", ...)  # BAD
```
✅ **Do**: Read from master file
```r
plots <- read_csv2("final_plot_name.csv") %>% pull(plot_name)  # GOOD
```

❌ **Don't**: Use `read.csv()` with defaults
```r
data <- read.csv("file.csv")  # Will fail with European format
```
✅ **Do**: Use `read_csv2()` or specify delimiters
```r
data <- read_csv2("file.csv")  # Correct for European format
```

❌ **Don't**: Ignore failed iterations silently
```r
for (i in plots) {
  process(i)  # Silent failure
}
```
✅ **Do**: Report errors
```r
for (i in plots) {
  tryCatch({
    process(i)
  }, error = function(e) {
    message("Erreur pour ", i, ": ", e$message)
  })
}
```

## Development Workflow

When adding new analyses:
1. **Start from template**: Use existing scripts as templates for structure
2. **Test on subset**: Use 2-3 plots before running on full dataset
3. **Document parameters**: Add clear comments for all tunable parameters
4. **Save intermediate results**: Don't recompute expensive operations
5. **Version outputs**: Add dates or versions to output files if iterating

## Module-Specific README Files

Each numbered module (0-5) now has a detailed README.md explaining:
- Purpose and objectives
- Input/output data
- Script descriptions and workflows
- Usage examples
- Validation procedures
- Common problems and solutions

**Refer to these READMEs for user-facing documentation**. The CLAUDE.md focuses on coding conventions and developer guidance.

## Git/GitHub Workflow (MANDATORY)

### Default Branch Strategy
**ALWAYS work on the `claude` branch** when making changes via Claude Code:
```bash
# Check current branch
git branch

# Create claude branch (first time only)
git checkout -b claude

# Switch to claude branch
git checkout claude
```

### Commit Workflow

#### When to Commit
Commit changes when:
1. ✅ Completing a logical unit of work (e.g., updating one module's paths)
2. ✅ After successfully testing a new feature
3. ✅ Before switching to a different task
4. ✅ At the end of a working session

#### Commit Message Convention
Follow this structure:
```
<type>(<scope>): <short summary>

<optional detailed description>

<optional footer>
```

**Types**:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `refactor`: Code refactoring (no functionality change)
- `style`: Formatting changes (no code logic change)
- `chore`: Maintenance tasks
- `test`: Adding or updating tests

**Examples**:
```bash
# Good commit messages
git commit -m "refactor(paths): migrate module 0 to use here() for relative paths"
git commit -m "docs(readme): update module 3 README with workflow details"
git commit -m "feat(gaps): add multi-scale roughness metrics"
git commit -m "fix(plots): correct EPSG assignment for Mabounie plots"

# Bad commit messages (avoid)
git commit -m "update"
git commit -m "fix bug"
git commit -m "changes"
```

### Standard Git Commands

#### Check Status
```bash
# See what files have changed
git status

# See detailed changes
git diff

# See changes for a specific file
git diff path/to/file.R
```

#### Stage and Commit
```bash
# Stage specific files
git add path/to/file1.R path/to/file2.R

# Stage all changes in a module
git add 3_Gaps/

# Stage all R scripts
git add "*.R"

# Commit with message
git commit -m "refactor(module3): update paths to use here()"

# Amend last commit (if needed)
git commit --amend -m "refactor(module3): update paths and add documentation"
```

#### Push to Remote
```bash
# Push claude branch to remote (first time)
git push -u origin claude

# Subsequent pushes
git push

# Force push (USE WITH CAUTION, only on claude branch)
git push --force  # Only if you've amended commits
```

#### View History
```bash
# See commit history
git log --oneline -10

# See files changed in last commit
git show --name-only

# See detailed changes in last commit
git show
```

### Working with Remote

#### Before Starting Work
```bash
# Fetch latest changes from remote
git fetch origin

# Check if claude branch is up to date
git status

# Pull changes if behind
git pull origin claude
```

#### Regular Sync Pattern
```bash
# 1. Check status
git status

# 2. Stage changes
git add .

# 3. Commit
git commit -m "type(scope): description"

# 4. Push
git push
```

### Merge Workflow (For User)
When `claude` branch is ready to merge to `main`:
```bash
# Switch to main
git checkout main

# Pull latest main
git pull origin main

# Merge claude branch
git merge claude

# Resolve conflicts if any
# ... edit conflicted files ...
git add <resolved-files>
git commit

# Push merged main
git push origin main
```

### Conflict Resolution
If conflicts occur:
1. Git will mark conflicts in files with `<<<<<<<`, `=======`, `>>>>>>>`
2. Edit files to resolve conflicts
3. Remove conflict markers
4. Test that code still works
5. Stage resolved files: `git add <file>`
6. Complete merge: `git commit`

### .gitignore Best Practices
The project `.gitignore` excludes:
- Large data files (`*.las`, `*.tif`, `*.lax`)
- R temporary files (`.Rhistory`, `.RData`, `.Rproj.user/`)
- Output directories with large results (should be documented in README)

**Never commit**:
- Raw LiDAR files (too large)
- Intermediate rasters (regenerate-able)
- Personal configuration files
- Sensitive data

**Always commit**:
- R scripts
- Documentation (`.md` files)
- Small reference data (`*.csv` < 10MB)
- Configuration files (`.Rproj`, `renv.lock`)

### Branch Protection
- `main` branch: Production-ready code only
- `claude` branch: Active development by Claude Code
- Feature branches: For specific experimental features (optional)

### Emergency Procedures

#### Undo Last Commit (not pushed)
```bash
# Keep changes in working directory
git reset --soft HEAD~1

# Discard changes completely
git reset --hard HEAD~1
```

#### Undo Last Push (dangerous!)
```bash
# Only on claude branch, never on main!
git reset --hard HEAD~1
git push --force origin claude
```

#### Discard All Local Changes
```bash
# Discard changes in specific file
git checkout -- path/to/file.R

# Discard all changes
git reset --hard HEAD
```

#### Recover Lost Commits
```bash
# Show reflog (history of HEAD movements)
git reflog

# Recover to a specific commit
git reset --hard <commit-hash>
```

### Best Practices Summary

1. **Always work on `claude` branch**
2. **Commit frequently** with clear messages
3. **Push regularly** to backup work
4. **Never force push to `main`**
5. **Test before committing**
6. **Use descriptive commit messages**
7. **Keep commits atomic** (one logical change per commit)
8. **Don't commit large binary files**
9. **Pull before starting new work**
10. **Review changes before committing** (`git diff`)

### Daily Workflow Example
```bash
# Morning - start work
git checkout claude
git pull origin claude

# Make changes, test
# ... edit files ...

# Commit progress
git add 3_Gaps/1_Gaps_identification.R
git commit -m "refactor(gaps): migrate to here() paths"

# Continue working
# ... more edits ...

# Commit again
git add 3_Gaps/2_Gaps_data_frequencies.R 3_Gaps/3_Gaps_metrics.R
git commit -m "refactor(gaps): complete path migration for core scripts"

# End of day - push all work
git push

# Verify push succeeded
git status
```

This workflow ensures:
- ✅ All changes are tracked
- ✅ Work is backed up remotely
- ✅ Easy to review what changed
- ✅ Main branch stays clean
- ✅ Easy rollback if needed
