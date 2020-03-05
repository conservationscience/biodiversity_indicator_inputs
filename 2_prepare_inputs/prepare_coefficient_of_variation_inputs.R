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
#' replicate.  But may need to be altered for replicates later
#' 
#' TODO: Could all be done better.  Make the way we split groups more flexible via
#' some arguments, currently it is hard coded so you can only split by model
#' functional groups. Might need to wait until after species matching
#' 
#' TODO: Add default values to variable, interval and function (adult_biomass,
#' 12 and mean respectively) unless specified.
#' 

prepare_coefficient_of_variation_inputs <- function(indicators_project, scenario, 
                          simulation_number, variable, burnin, interval, func){
  
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
  
  coefficient_of_variation_inputs <- replicates_no_burnin[[1]]
  
  saveRDS( coefficient_of_variation_inputs, 
           file = file.path(output_folder,paste(scenario, simulation_number, 
                  "coefficient_of_variation_inputs", sep = "_" )))
  
  return(coefficient_of_variation_inputs)

 
}

# Test function

# 
# indicators_project <- "N:/Quantitative-Ecology/Indicators-Project/" # File path for entire project directory
# location <- 'Serengeti' # Modelled location you want to process
# scenario <- 'Test_runs' # Scenario you want to process
# simulation <- 'aa_BuildModel/' # Model simulation number aka 'BuildModel' directory you want to process
# burnin <- 0.5*12 # number of years burn in * 12 (to convert to monthly)
# simulation_number <- "aa" # Number of your buildmodel file
# variable <- "adult_biomass"
# interval <- 3
# func <- mean
# pre_exploitation_period <- 2
# colour_scheme <- c("black", "red", "blue","turquoise3", "purple","green3", "pink")
# 
# test <- prepare_coefficient_of_variation_inputs(indicators_project, scenario,
#         simulation_number, variable, burnin, interval, func)
#  
  
