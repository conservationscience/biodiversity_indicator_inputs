# Using R version 3.5.2

#' Return a dataframe that ...

#' @param indicators_project file path where all project data is kept
#' @param scenario string that denotes which impact directory you want to look in
#' @param simulation_number string that denotes which simulation directory you 
#' want to prepare the data from
#' @param variable string that denotes whether you want to use abundance or biomass
#' @param burnin integer, specifies how many burnin timesteps to remove
#' @param interval integer, specifies how many monthly timesteps you want to 
#' aggregate by.  Should usually be 12 if you want to convert monthly to annual
#' @param func function to use when converting monthly timesteps (eg mean, min,
#' max, sample etc)
#' @return a dataframe (or list of dataframes, see TODO section) where the columns
#' are mean monthly biomass for that year (or whatever specified interval is)


#' TODO: At this stage the inputs and outputs of this function are for a single
#' replicate.  But I have set it up as a list in case later on it makes more
#' sense to average over the replicates before calculating the indicator
#' 
#' TODO: Add default values to variable, interval and function (adult_biomass,
#' 12 and mean respectively) unless specified.
#' 

# indicators_project <- "N:/Quantitative-Ecology/Indicators-Project/" # File path for entire project directory
# location <- 'Serengeti' # Modelled location you want to process
# scenario <- 'Test_runs' # Scenario you want to process
# simulation <- 'aa_BuildModel/' # Model simulation number aka 'BuildModel' directory you want to process
# burnin <- 0.5*12 # number of years burn in * 12 (to convert to monthly)
# simulation_number <- "aa" # Number of your buildmodel file
# variable <- "adult_biomass"
# interval <- 3
# func <- mean

prepare_proportion_total_biomass_inputs <- function(indicators_project, scenario, simulation_number, variable, burnin, interval, func){
  
  require(stringr)
  require(tidyverse)
  require(reshape2)
  
  
  # Find data
  
  all_scenario_folders <- list.dirs(file.path(indicators_project, "Serengeti",
                                              "Outputs_from_adaptor_code",
                                              "map_of_life"), recursive = FALSE)
  
  target_scenario_folder <- str_subset(all_scenario_folders, scenario)
  
  scenario_simulation_folders <- list.dirs(target_scenario_folder, recursive = FALSE)
  
  processed_simulation_outputs <- str_subset(scenario_simulation_folders, simulation_number )
  
  all_files <- list.files(processed_simulation_outputs)
  
  target_files <- str_subset(all_files, variable)
  target_rds_files <- target_files[!str_detect(target_files, ".csv")]
  
  replicates <- list()
  
  for (i in seq_along(target_rds_files)) {
    
  replicates[[i]] <- readRDS(file.path(processed_simulation_outputs, target_rds_files[[i]]))
  
  }
  
  # Create or set output folder
  
  output_folder <- file.path(indicators_project, 
                             "/Serengeti/Outputs_from_indicator_code",
                             paste(scenario, "indicator_outputs" , sep ="_"))
  
  if( !dir.exists( file.path(output_folder) ) ) {
    dir.create( file.path(output_folder), recursive = TRUE )
    
  }

  # Function to remove burnin
  
  remove_burn_in <- function(data, burnin) {
      
      data[,(burnin + 1):ncol(data)]
      
  }
    
  replicates_no_burnin <- lapply(replicates, remove_burn_in, burnin )

  # Function to convert monthly timesteps to yearly by taking the mean of a specified interval (12 to convert monthly to yearly)
    
  convert_timesteps <- function(dataframe, interval, func){
      
      monthly_matrix <- t(dataframe)
      
      n <- interval
      
      time_converted_matrix <- t(aggregate(monthly_matrix,list(rep(1:(nrow(monthly_matrix) %/% 
                                                                        n + 1), each = n, len = nrow(monthly_matrix))), func, na.rm = TRUE))
      
      time_converted_matrix <- time_converted_matrix[-1,]
      
      time_converted_matrix[is.nan(time_converted_matrix)] = NA
      
      return(time_converted_matrix)
      
    }
    
    # Loop through each replicate and convert biomass per month to mean annual biomass
    
    if (interval > 1) {
      
      time_converted_replicates <- list()
      
      for (i in seq_along(replicates_no_burnin)) {
        
        time_converted_replicates[[i]] <- 
          convert_timesteps(replicates_no_burnin[[i]], interval, func)
      }
      
      proportion_total_biomass_inputs <- time_converted_replicates
      
      saveRDS( proportion_total_biomass_inputs, file = file.path(output_folder,
      paste(scenario, simulation_number, "proportion_total_biomass_inputs", sep = "_" )))
      
      return(proportion_total_biomass_inputs)
     
    } else if (interval == 1 ) {
      
      proportion_total_biomass_inputs <- replicates_no_burnin
      
      saveRDS( proportion_total_biomass_inputs, 
               file = file.path(output_folder,paste(scenario, simulation_number, 
               "proportion_total_biomass_inputs", sep = "_" )))
     
      return(proportion_total_biomass_inputs)
      
    }
}
    
 
  
  # Take the mean of all replicates to create one dataframe of proportion biomass
  # at each annual timestep
  
  # proportion_total_biomass <-  proportion_total_biomass_replicate_df %>%
  #                              mutate(relative_proportion_biomass = 
  #                              rowMeans(proportion_total_biomass_replicate_df,
  #                                             na.rm = TRUE)) %>%
  #                              mutate(year = c(1:nrow(proportion_total_biomass_replicate_df))) %>%
  #                              select(relative_proportion_biomass, year)
  
  


  

  