# Using R version 3.5.2

#' Return a dataframe that denotes the Red List Status of each
#' "species" over a specified interval of time based on change in abundance

#' @param indicators_project file path where all project data is kept
#' @param scenario string that denotes which impact directory you want to look in
#' @param simulation_number string that denotes which simulation directory you 
#' want to prepare the data from
#' @param format string that specifies whether you want long or wide format
#' output data. Defaults to long.
#' 
#' TODO: Not sure method for calculating percentage change is correct, need
#' to read RL guidelines more carefully (might need to add scenario as an input?)
#' 
#' TODO: Exponentiate abundance figures?
#' 
#' TODO: Check RL status thresholds
#' 
#' TODO: This is very slow, reading in the abundance data takes a long time. 
#' See if there's a way to optimise it. 

prepare_rli_inputs <- function(indicators_project, scenario, simulation_number, format = "long"){
  
  require(stringr)
  require(tidyverse)
  require(reshape2)
  
  #source("C:/Users/ssteven/Dropbox/Deakin/Serengeti-analysis/BiodiversityIndicators-project-code/madingley_code/prepare_data_functions.R")
  
  # Find data
  
    all_scenario_folders <- list.dirs(file.path(indicators_project, "Serengeti",
                                              "Outputs_from_adaptor_code",
                                              "map_of_life"), recursive = FALSE)
  
    target_scenario_folder <- str_subset(all_scenario_folders, scenario)
    
    ## Need to fix this
    
    scenario_simulation_folders <- list.dirs(target_scenario_folder, recursive = FALSE)
    
    processed_simulation_outputs <- str_subset(scenario_simulation_folders, simulation_number )
    
    files <- list.files(processed_simulation_outputs)

    generation_files <- str_subset(files, "generation_lengths")
    generation_lengths <- generation_files[!str_detect(generation_files, ".csv")]
    generation <- readRDS(file.path(processed_simulation_outputs, generation_lengths))
    
    abundance_files <- str_subset(files, "adult_abundance")
    abundance <- abundance_files[!str_detect(abundance_files, ".csv")]
    data <- readRDS(file.path(processed_simulation_outputs, abundance))
    
    print("testing function with adult abundance files when this is fixed")
  
  # Create or set output folder
    
    output_folder <- file.path(indicators_project, 
                               "/Serengeti/Outputs_from_indicator_code",
                               paste(scenario, "indicator_outputs" , sep ="_"))
    
    if( !dir.exists( file.path(output_folder) ) ) {
      dir.create( file.path(output_folder), recursive = TRUE )
   
    }
  
  # Tidy adaptor output data (remove burn in, convert to annual timesteps)
  
    tmp <- data %>%
    remove_burn_in(12*1000) %>%
    convert_timesteps(12, mean) %>%
    as.data.frame() %>%
    dplyr::mutate(functional_group_index = row.names(.))
  
  # Add timeframe over which to assess decline, using generation lengths
    
    annual_data <- generation %>%
                   dplyr::mutate(generation_by_three = gen_length_mean * 3) %>%
                   dplyr::mutate(timeframe = ifelse(generation_by_three > 10 , 
                                                    round(generation_by_three), 10)) %>%
                   dplyr::select(functional_group_index, timeframe) %>%
                   merge(tmp, by = "functional_group_index")
    
    names(annual_data) <- c("functional_group_index","timeframe",1:300)

    rm(tmp, data, generation)
    
 # exponentiate?
    
 # Calculate percentage change over timeframe for each row
    
    # Convert to long format and split into individual dataframes per species
    # (allows us to apply a different assessment timeframe per species)
    
    annual_data_long <- melt(annual_data, id = c("functional_group_index","timeframe"))
    
    annual_data_long$timeframe <- as.integer(annual_data_long$timeframe)
    
    annual_data_long <- annual_data_long %>%
                        dplyr::mutate(functional_group = case_when(grepl("10.", functional_group_index) ~ "herbivorous endotherm",
                                                   grepl("11.", functional_group_index) ~ "carnivorous endotherm",
                                                   grepl("12.", functional_group_index) ~ "omnivorous endotherm",
                                                   grepl("13.16.", functional_group_index) ~ "herbivorous ectotherm",
                                                   grepl("14.17.", functional_group_index) ~ "carnivorous ectotherm",
                                                   grepl("15.18.", functional_group_index) ~ "omnivorous ectotherm"))
    
    
    
    annual_by_species <- split(annual_data_long, annual_data_long$functional_group_index)
    
    rm(annual_data, annual_data_long)
    
    # Calculate percentage change and assign Red List Status according to magnitude
    # of change
    
    rl_status_list <- list()
    
    for (i in seq_along(annual_by_species)) {
      
    t <- annual_by_species[[i]]$timeframe[1]
      
    rl_status_list[[i]] <- annual_by_species[[i]]  %>%
      dplyr::mutate(Diff = 1 - (value/lag(value, 10))) %>% ## TODO: fix timeframe
      dplyr::mutate(RL_status = ifelse(Diff >= -0.40, "LC",
                                ifelse(Diff <= -0.40 & Diff >= -0.50, "NT",
                                ifelse(Diff < -0.50 & Diff >= -0.70, "VU",
                                ifelse(Diff < -0.70 & Diff >= -0.90, "EN",
                                ifelse(Diff < -0.90 & Diff >= -0.99, "CR",
                                ifelse(Diff <= -0.99, "EX",
                                "NA"))))))) 
    
    #df.new <- df[seq(1, nrow(df), t), ] # Sample by timeframe?
    
    }  
    
    rm(annual_by_species)
    
    if (format == "wide") {
    
    # Convert each species dataframe back into wide format (row is species,
    # columns are timesteps)
    
    rl_status_wide_list <- list()
    
    for (i in seq_along(rl_status_list)){
    
    rl_status_wide_list[[i]] <- rl_status_list[[i]] %>%
                                dplyr::select(functional_group_index,variable, 
                                              RL_status) %>%
                                spread(variable, RL_status, fill = NA, 
                                       convert = FALSE)
  }
    
    # Combine each species back into one data frame
    
    rl_status_wide <- do.call(rbind, rl_status_wide_list)
    
    rm(rl_status_wide_list)
    
    return(rl_status_wide)
    
    # Save in output folder
    
    saveRDS( rl_status_wide, file = file.path(output_folder,paste(scenario, simulation_number, "RLI_inputs_wide", sep = "_" )))
    write.csv( rl_status_wide, file = file.path(output_folder,paste(scenario, simulation_number, "RLI_inputs_wide.csv", sep = "_" )))

    } else {
      
      rl_status_long <- do.call(rbind, rl_status_list) %>%
                        dplyr::select(functional_group_index, functional_group, variable, RL_status)
      
      names(rl_status_long) <- c("functional_group_index", "functional_group", "timestep", "RL_status")
      
      return(rl_status_long)
      
      # Save in output folder
      
      saveRDS( rl_status_long, file = file.path(output_folder,paste(scenario, simulation_number, "RLI_inputs_long", sep = "_" )))
      write.csv( rl_status_long, file = file.path(output_folder,paste(scenario, simulation_number, "RLI_inputs_long.csv", sep = "_" )))
      
  
}
  
}
          


