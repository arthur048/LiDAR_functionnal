# Script pour mettre à jour tous les chemins hardcodés vers here()
# Ce script est à exécuter UNE SEULE FOIS pour la migration

library(here)
library(stringr)

# Pattern à rechercher
old_pattern <- "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/"

# Fonction pour remplacer les chemins dans un fichier
update_file_paths <- function(file_path) {
  cat("Traitement de:", file_path, "\n")

  # Lire le fichier
  content <- readLines(file_path, warn = FALSE)

  # Ajouter library(here) si nécessaire
  has_here <- any(str_detect(content, "library\\(here\\)"))
  if (!has_here && any(str_detect(content, old_pattern))) {
    # Trouver où insérer library(here) - après les autres library()
    lib_lines <- which(str_detect(content, "^library\\("))
    if (length(lib_lines) > 0) {
      insert_pos <- max(lib_lines)
      content <- c(content[1:insert_pos],
                   "library(here)",
                   content[(insert_pos+1):length(content)])
    }
  }

  # Remplacer les chemins
  # Pattern 1: st_read("E:/Arthur/.../file.shp") -> st_read(here("path", "file.shp"))
  # Pattern 2: read_csv2("E:/Arthur/.../file.csv") -> read_csv2(here("path", "file.csv"))
  # etc.

  for (i in seq_along(content)) {
    line <- content[i]

    # Si la ligne contient le pattern ancien
    if (str_detect(line, fixed(old_pattern))) {
      # Ne pas toucher aux commentaires d'explication
      if (str_detect(line, "^#")) next

      # Ne pas toucher aux chemins externes (Doctorat_DataAnnexe, OriginalDataFiles)
      if (str_detect(line, "Doctorat_DataAnnexe|OriginalDataFiles")) next

      # Extraire les parties du chemin après LiDAR_functionnal/
      matches <- str_match_all(line, paste0(fixed(old_pattern), "([^\"']+)"))[[1]]

      if (nrow(matches) > 0) {
        for (j in 1:nrow(matches)) {
          old_path <- matches[j, 1]
          relative_path <- matches[j, 2]

          # Convertir le chemin en arguments here()
          path_parts <- str_split(relative_path, "/")[[1]]
          here_args <- paste0('"', path_parts, '"', collapse = ", ")
          new_path <- paste0('here(', here_args, ')')

          # Remplacer dans la ligne
          line <- str_replace(line, fixed(old_path), new_path)
        }
        content[i] <- line
      }
    }
  }

  # Écrire le fichier modifié
  writeLines(content, file_path)
  cat("  -> Mise à jour terminée\n")
}

# Liste des scripts R à traiter (excluant les archives et renv)
scripts_to_update <- list.files(
  path = here(),
  pattern = "\\.R$",
  recursive = TRUE,
  full.names = TRUE
) %>%
  .[!str_detect(., "renv/|Archives/|archive|update_paths\\.R")]

# Afficher le nombre de fichiers
cat("Nombre de fichiers à traiter:", length(scripts_to_update), "\n\n")

# Traiter chaque fichier
walk(scripts_to_update, update_file_paths)

cat("\nMise à jour terminée!\n")
