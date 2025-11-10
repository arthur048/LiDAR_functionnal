rm(list=ls())
gc()


create_dir = function(path) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}

# Initialisation et path preparation ----------------------------------------------------------

pkgs = c("terra", "tidyverse", "sf", "rio", "caret", "progress")

to_install = !pkgs %in% installed.packages() ; if(any(to_install)) {install.packages(pkgs[to_install])} ; inst = lapply(pkgs, library, character.only = TRUE) # load them

path0 = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/5_Sensibility_Analysis_Modelisation/" ; setwd(path0)
path_output = paste0(path0,"output/") ; lapply(path_output, create_dir)

path_gaps = "E:/Arthur/OneDrive2/R/DoctoratGIS/WorkingFiles/LiDAR_functionnal/3_Gaps/output/"
gaps_freq = rio::import(paste0(path_gaps, "results_gaps_freq.csv"), sep = ";", dec = ",") %>% as_tibble()
gaps_metrics = rio::import(paste0(path_gaps, "results_gaps_metrics.csv"), sep = ";", dec = ",") %>% as_tibble()

plot_to_exclude = c("plot_JEU2","plot_JEU3","plot_JEU4","plot_JEU5")
plot_name = gaps_freq %>% filter(!plot_name %in% plot_to_exclude) %>% pull(plot_name) %>% unique()

model_list = readRDS(file = paste0(path0, "model.list.rds"))

# Data preparation for results extraction ---------------------------------------

create_empty_structure <- function(lst, current_depth = 1, max_depth = Inf) {
  # Base case: if current depth exceeds max_depth or if the item is not a list, stop recursion
  if (current_depth > max_depth || !is.list(lst)) {
    return(NULL)
  }
  # Recursive case: create a new list with the same structure but empty
  setNames(lapply(lst, create_empty_structure, current_depth = current_depth + 1, max_depth = max_depth), names(lst))
}

result_list <- create_empty_structure(model_list, max_depth = 5)


# Extraction des résultats ------------------------------------------------

maxz = length(model_list)*length(model_list[[1]])*length(model_list[[1]][[1]])*length(model_list[[1]][[1]][[1]])*length(model_list[[1]][[1]][[1]][[1]])
pb = progress_bar$new(format = "  Progress: [:bar] :percent :elapsed/:eta", total = maxz, clear = FALSE)

results = tibble()

z = 0 ; a=1 ; b=2 ; c=3 ; d=16; e = 1 
for(a in 1:length(model_list)){
  for(b in 1:length(model_list[[1]])){
    for (c in 1:length(model_list[[1]][[1]])){
      for (d in 1:length(model_list[[1]][[1]][[1]])){
        for (e in 1:length(model_list[[1]][[1]][[1]][[1]])){
          
          z = z + 1 
          
          model = model_list[[a]][[b]][[c]][[d]][[e]]
          
          test_model = inherits(model, "train")
          
          if(test_model){
            
            coefficients_str <- paste(names(model[["finalModel"]][["coefficients"]]), 
                                      model[["finalModel"]][["coefficients"]], 
                                      sep = "=", collapse = " + ")
            
            model_summary <- summary(model)
            
            overall_p_value <- pf(model_summary$fstatistic["value"], 
                                  model_summary$fstatistic["numdf"], 
                                  model_summary$fstatistic["dendf"], 
                                  lower.tail = FALSE) %>% as.numeric()
            
            
            result_list[[a]][[b]][[c]][[d]][[e]] = tibble(
              respons = names(result_list)[[a]],
              predictors = if (names(result_list[[a]])[[b]] == "lambda") "proportion + lambda" else names(result_list[[a]])[[b]],
              buffer = names(result_list[[a]][[b]])[[c]],
              height_aboveground = names(result_list[[a]][[b]][[c]])[[d]],
              xmin = names(result_list[[a]][[b]][[c]][[d]])[[e]],
              RMSE = model$results$RMSE,
              rRMSE = model$results$RMSE / mean(model$pred$obs),
              Rsquared = model$results$Rsquared,
              MAE = model$results$MAE,
              P_Value = overall_p_value,
              model = coefficients_str
            )
            
          } else { 
            
            result_list[[a]][[b]][[c]][[d]][[e]] = tibble(
              respons = names(result_list)[[a]],
              predictors = if (names(result_list[[a]])[[b]] == "lambda") "proportion + lambda" else names(result_list[[a]])[[b]],
              buffer = names(result_list[[a]][[b]])[[c]],
              height_aboveground = names(result_list[[a]][[b]][[c]])[[d]],
              xmin = names(result_list[[a]][[b]][[c]][[d]])[[e]],
              RMSE = NA,
              rRMSE = NA,
              Rsquared = NA,
              MAE = NA,
              P_Value = NA,
              model = NA
            )
          }
          
          results = results %>% bind_rows(result_list[[a]][[b]][[c]][[d]][[e]])
          
          pb$tick()  # Update progress bar
          
        }
      }
    }
  }
}

saveRDS(result_list, file = paste0(path_output, "result_list.rds"))
rio::export(results, paste0(path_output, "results_sensibility_analysis_modelisation.csv"), dec = ",", sep = ";", append = FALSE)
