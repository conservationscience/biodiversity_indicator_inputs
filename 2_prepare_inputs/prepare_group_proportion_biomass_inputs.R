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


prepare_group_proportion_biomass_inputs <- function(indicators_project, scenario, 
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
  
  # Read in adult biomass over time for all groups (each simulation should be
  # a matrix where rows = 'species' and columns = 'monthly timesteps for entire
  # model duration including burnin period)
  
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
  
  # Remove the burnin time steps from the matrix
  
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
    
  #  if (interval > 1) {
      
      time_converted_replicates <- list()
      
      for (i in seq_along(replicates_no_burnin)) {
        
        time_converted_replicates[[i]] <- 
          convert_timesteps(replicates_no_burnin[[i]], interval, func)
      }
      
   # }
      
      # Subset dataframes by functional group
      
      time_converted_replicates <- time_converted_replicates[[1]]
      
      herbivore_endotherm_replicates <- time_converted_replicates[grepl("10.",
                                        rownames(time_converted_replicates)),]
      
      herbivore_ectotherm_replicates <- time_converted_replicates[grepl("13.16.",
                                        rownames(time_converted_replicates)),]
      
      omnivore_endotherm_replicates <- time_converted_replicates[grepl("12.",
                                        rownames(time_converted_replicates)),]
      
      omnivore_ectotherm_replicates <- time_converted_replicates[grepl("15.18.",
                                        rownames(time_converted_replicates)),]
      
      carnivore_endotherm_replicates <- time_converted_replicates[grepl("11.",
                                        rownames(time_converted_replicates)),]
      
      carnivore_ectotherm_replicates <- time_converted_replicates[grepl("14.17.",
                                        rownames(time_converted_replicates)),]

      
      group_proportion_biomass_inputs <- list(herbivore_endotherm_replicates,
                                                 herbivore_ectotherm_replicates,
                                                 omnivore_endotherm_replicates,
                                                 omnivore_ectotherm_replicates,
                                                 carnivore_endotherm_replicates,
                                                 carnivore_ectotherm_replicates)
      
      names(group_proportion_biomass_inputs) <- c("herbivorous endotherms",
                                                  "herbivorous ectotherms",
                                                  "omnivorous endotherms",
                                                  "omnivorous ectotherms",
                                                  "carnivorous endotherms",
                                                  "carnivorous ectotherms")
      
      saveRDS( group_proportion_biomass_inputs, 
               file = file.path(output_folder,paste(scenario, simulation_number, 
               "group_proportion_biomass_inputs", sep = "_" )))
      
      return(group_proportion_biomass_inputs) 
    
}

# Test function


# indicators_project <- "N:/Quantitative-Ecology/Indicators-Project/" # File path for entire project directory
# location <- 'Serengeti' # Modelled location you want to process
# scenario <- 'Test_runs' # Scenario you want to process
# simulation <- 'aa_BuildModel/' # Model simulation number aka 'BuildModel' directory you want to process
# burnin <- 0.5*12 # number of years burn in * 12 (to convert to monthly)
# simulation_number <- "aa" # Number of your buildmodel file
# variable <- "adult_biomass"
# interval <- 3
# func <- mean
# 
# test <- prepare_group_proportion_biomass_inputs(indicators_project, scenario,
#                               simulation_number, variable, burnin, interval, func)
 
  
