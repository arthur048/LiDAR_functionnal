rm(list=ls())
gc()

# Initialisation ----------------------------------------------------------
pkgs = c("terra", "tidyverse", "sf", "rio", "foreach", "doParallel")

to_install = !pkgs %in% installed.packages() ; if(any(to_install)) {install.packages(pkgs[to_install])} ; inst = lapply(pkgs, library, character.only = TRUE) # load them

path0 = "here()3_Gaps/" ; setwd(path0)

# Fonctions d'analyses des trouées (équivalent du package ForestGa --------

create_dir = function(path) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}

# Mes propres fonctions pour l'analyse des gaps
get_gap_layer = function(chm_layer, threshold = 10, size = c(1, 10^8)) {
  
  ##### Cette fonction prend en entrée un CHM de type SpatRaster avec des NA dans les zones à ne pas analyser,
  ##### Et renvoie en sortie une couche avec les gaps identifié sous format SpatRaster avec un numéro associé à chaque gap
  ##### La fonction demande en entrée a minimum height abovegrount to define gap, and the minimum and maximum area of gaps (m²)
  ##### La couche d'entrée doit avoir une projection en format métrique et non degrés décimaux 
  
  # Reclassifying based on threshold
  chm_layer = terra::classify(chm_layer,
                              rcl = c(-Inf, threshold, 1,
                                      threshold, Inf, NA) %>% matrix(ncol = 3, byrow = TRUE),
                              right = TRUE # > threshold = NA, <= threshold = 1
  )
  
  # Clumping connected areas
  gaps = terra::patches(chm_layer, directions = 8)
  
  # Frequency and size calculation
  rcl = terra::freq(gaps)
  rcl[, 3] <- rcl[, 3] * terra::res(chm_layer)[1]^2
  
  # Reclassify gaps based on size
  z = terra::classify(gaps,
                      rcl = rcl[,2:3],
                      right = NA
  )
  
  # Remove cells based on min and max size
  gaps[values(z) > size[2]] = NA
  gaps[values(z) < size[1]] = NA
  
  gaps = terra::patches(gaps, directions = 8)
  names(gaps) = "gaps"
  
  return(gaps)
}

# Reading data -------------------------------------------------------------------------

path_output = paste0(path0,"output/") ; lapply(path_output, create_dir)
path_elevation = "here()2_ElevationData/"

path_chm_folder = paste0(path_elevation, "CHM_final/")
path_chm = lapply(path_chm_folder, function(folder) {
  list.files(path = folder, pattern = "\\.tif$", full.names = TRUE)
}) %>% unlist()

# Lecture du fichier plots_info
plots_info <- read_csv2("here()0_Inventories_plot_preparation/final/plots_info.csv")

plot_name = basename(gsub("\\.tif$", "", path_chm))

path_outputs = paste0(path_output, plot_name) ; lapply(path_outputs, create_dir)

# Gaps detection ----------------------------------------------------------

height_aboveground = seq.int(1:45)

time_start_all = Sys.time()
i = 1 ; height = 10 
# Configuration de la parallélisation
n_cores <- parallel::detectCores() - 2
cl <- makeCluster(floor(n_cores))
registerDoParallel(cl)

results <- foreach(
  i = 1:length(plot_name),
  .packages = c("terra", "tidyverse", "sf"),
  .errorhandling = "pass"
) %dopar% {
  tryCatch({
    plot_name_part = plot_name[[i]]
    path_output_part = path_outputs[basename(path_outputs) == plot_name_part]
    path_chm_part = path_chm[basename(path_chm) == paste0(plot_name_part, ".tif")]
    
    chm = rast(path_chm_part)
    for (height in height_aboveground){
      gaps = get_gap_layer(chm_layer = chm, threshold = height)
      writeRaster(gaps, paste0(path_output_part, "/", plot_name_part, "_", height,"m_gaps_height.tif"), overwrite = TRUE)
    }
    
    return(list(
      status = "success",
      plot_name = plot_name_part
    ))
    
  }, error = function(e) {
    return(list(
      status = "error",
      plot_name = plot_name_part,
      message = as.character(e)
    ))
  })
}

stopCluster(cl)

results_df <- bind_rows(lapply(results, function(x) {
  if(is.list(x)) {
    as_tibble(x)
  } else {
    tibble(status = "error", message = "Unknown error")
  }
}))

