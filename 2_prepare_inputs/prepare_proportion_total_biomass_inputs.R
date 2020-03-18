# Using R version 3.5.2

#' Return a dataframe that ...

#' @param inputs file path for an individual, processed dataframe
#' @param outputs file path where you would like to store your outputs
#' @param simulation_number string that denotes the simulation number the data came from (eg "001")
#' @param replicate_number string that denotes what number replicate within the simulation the data is
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

# indicator <- "proportion_total_biomass"
# variable <- "biomass"
# interval <- 3
# func <- mean
# simulation_number <- "ae"
# replicate_numbers <- 0: (length(test_input) - 1)
# replicate_number <- as.character(replicate_numbers[1])
# test_input <- "N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_adaptor_code/map_of_life/Test_runs/ae_BuildModel/MassBinsOutputs_NI_0_Cell0_biomass.rds"
# test_output <- file.path(IndicatorsProject, location, "Outputs_from_indicator_code/inputs", indicator,
#                          scenarios[1])
# x <- prepare_proportion_total_biomass_inputs(test_input, test_output,simulation_number, 
#                                              cell_number, 
#                                              burnin, interval, func )

prepare_proportion_total_biomass_inputs <- function(inputs, outputs, simulation_number, 
                                                    replicate_number, burnin, interval, func){
  
  require(stringr)
  require(tidyverse)
  require(reshape2)
  
  scenario <- basename(outputs)
  
  replicate <- readRDS(inputs)
  
  # Create or set output folder
  
  output_folder <- outputs
  
  if( !dir.exists( file.path(output_folder) ) ) {
    dir.create( file.path(output_folder), recursive = TRUE )
    
  }

  # Function to remove burnin
  
  remove_burn_in <- function(data, burnin) {
      
      data[,(burnin + 1):ncol(data)]
      
  }
    
  #replicates_no_burnin <- lapply(replicates, remove_burn_in, burnin )
  
  replicate_no_burnin <- remove_burn_in(replicate, burnin)

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
      
      proportion_total_biomass_inputs <- convert_timesteps(replicate_no_burnin, interval, func)
      
      saveRDS( proportion_total_biomass_inputs, file = file.path(output_folder,
      paste(scenario, simulation_number, replicate_number, "proportion_total_biomass_inputs", sep = "_" ))) 
      
      return(proportion_total_biomass_inputs)
     
    } else if (interval == 1 ) {
      
      proportion_total_biomass_inputs <- replicate_no_burnin
      
      saveRDS( proportion_total_biomass_inputs, 
               file = file.path(output_folder,paste(scenario, simulation_number, cell_number,  
               "proportion_total_biomass_inputs", sep = "_" )))
     
      return(proportion_total_biomass_inputs)
      
    }
}
    

  
  


  

  