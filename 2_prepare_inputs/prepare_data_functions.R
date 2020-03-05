
## TODO: Look at split groups code for rowname matrix subsetting
## Also use which function in exploratory functions code - can use it to
## subset the matrix, but also add a new column with NAs first then use which
## to add information for selected rows only
## Also look at these examples to document what the functions do (set parameters)



ListMassBinsFiles <- function(resultsDir){
  
  files<-dir(resultsDir)
  files<-files[grep("MassBinsOutputs",files)]
  files<-files[grep("Cell",files)]
  files<-files[grep("rds",files)]
  
  return(files)
  
}

# TO DO: Fix this function so it renames the elements in the list

get_rds_files <- function(directory, variable){
  
  files <- ListMassBinsFiles(directory)
  files <- files[grep(".rds",files)]
  files <- file.path(directory, files)
  files <- files[grep(variable,files)]
  tmp <- lapply(files,readRDS)
  
}

# Function to select elements of a list by their name or part of their name

filter_by_pattern <- function(pattern, your.list){
  
  names(your.list) %>% 
    str_detect(pattern) %>%
    keep(your.list, .)
  
}

# Function to convert monthly data to yearly (or whatever interval)

convert_timesteps <- function(dataframe, interval, func){
  
  monthly_matrix <- t(dataframe)
  
  n <- interval
  
  time_converted_matrix <- t(aggregate(monthly_matrix,list(rep(1:(nrow(monthly_matrix) %/% 
                             n + 1), each = n, len = nrow(monthly_matrix))), func))
  
  time_converted_matrix <- time_converted_matrix[-1,]
  
  return(time_converted_matrix)
  
}


remove_burn_in <- function(data, burnin) {
  
  data[,(burnin + 1):ncol(data)]
  
}

# Function to remove burn-in period and take the mean over a specified interval
# 'Return a list of functional group variables for each replicate

#' @param directory string - directory containing single replicate outputs from adaptor code
#' @param variable string - denoting variable (abundance or biomass) to retrieve
#' @param interval numeric - number of time steps you want to apply 'func' over
#' @param func function - to reduce multiple time steps - mean, sample, median etc
#' @param burnin numeric - number of time steps set as burnin (in years)


## TO DO: At the moment this spits out matrices with a weird first row
## that we don't want or need

clean_replicate_data <- function(directory, variable, interval, func, burnin) {
  
  log_biomass_replicates_all_timesteps <- get_rds_files(directory, variable)
  
   # Replace -9999 with NAs so you don't get strange values when you average the replicates
  
  log_biomass_replicates_NA <- list()

  for (i in seq_along(log_biomass_replicates_all_timesteps)) {

    log_biomass_replicates_NA[[i]] <-
      na_if(log_biomass_replicates_all_timesteps[[i]], -9999)

  }
  # 
  # Remove the burn-in period from each file
  
  remove_burn_in <- function(data, burnin) {
    
   data[,(burnin + 1):ncol(data)]
    
  }
  
  log_biomass_replicates_impact_timesteps_only <- lapply(log_biomass_replicates_NA, 
                                                         remove_burn_in, (burnin * 12) + 1 )

  
  # Function to convert monthly timesteps to yearly by taking the mean of a specified interval (12 to convert monthly to yearly)
  
  convert_timesteps <- function(dataframe, interval, func){
    
    monthly_matrix <- t(dataframe)
    
    n <- interval
    
    time_converted_matrix <- t(aggregate(monthly_matrix,list(rep(1:(nrow(monthly_matrix) %/% 
                              n + 1), each = n, len = nrow(monthly_matrix))), func))
    
    time_converted_matrix <- time_converted_matrix[-1,]
    
    return(time_converted_matrix)
    
  }
  
  # Loop through each replicate and convert biomass per month to mean annual biomass
  
  if (interval > 1) {
  
  time_converted_replicates <- list()
  
  for (i in seq_along(log_biomass_replicates_impact_timesteps_only)) {
    
    time_converted_replicates[[i]] <- 
      convert_timesteps(log_biomass_replicates_impact_timesteps_only[[i]], interval, func)
  }
  
  
  
  return(time_converted_replicates)
  
  } else if (interval == 1 ) {
    
  return(log_biomass_replicates_impact_timesteps_only)
    
  }
  
}

take_the_mean_of_replicates <- function(list, label) {
  
  if (names(list) == logical(0)) {
    
    category <- label
  
    } else {
  
  names <- names(list)
  category <- names
  
  }
  
  replicate_df <- do.call(cbind,list) 
  replicate_df <- as.data.frame(replicate_df[, -grep("year", colnames(replicate_df))])
  
  mean_value <-  replicate_df %>%
  mutate(label = rowMeans(replicate_df, na.rm = TRUE)) %>%
  mutate(year = c(1:nrow(replicate_df))) %>%
  select(label, year) %>%
  mutate(group = category)
  
  return(mean_value)
  
}

# Function to take the mean of a list of replicates

#' Return a matrix that is the mean of a list of proportion biomass replicates

#' @param list a list of proportion biomass replicates
#' @param identifier string, the name of the group (eg total, herbivore endotherm etc)


take_the_mean_of_replicates <- function(list, identifier) {
  
  list_df <- do.call(cbind,list) 
  
  list_df <- list_df[, -grep("year", colnames(list_df))]
  
  mean_df <-  list_df %>%
    mutate(relative_proportion_biomass = 
           rowMeans(list_df, na.rm = TRUE)) %>%
    mutate(year = c(1:nrow(list_df))) %>%
    dplyr::select(relative_proportion_biomass, year) %>%
    mutate(group = identifier)
  
  if (all(is.nan(as.matrix(mean_df$relative_proportion_biomass))) == TRUE) {
    
    message(message(paste("There is no biomass for", identifier, ", NaNs produced" )))
    
  } else {
    
  return(mean_df)
  
  }
}

# Take the mean biomass per group, per timestep, across all replicates

reduce_matrix_list <- function(list){
  
  mean_matrix <- Reduce("+", lapply(list, function(x) replace(x, is.na(x), 0)))/ 
                 Reduce("+", lapply(list, Negate(is.na)))
  
  return(mean_matrix)
  
}

# Visualise the bodymass over time of different virtual species (i.e. bodymass bins
# of a functional group)

# TODO: Work in progress.  Update it so it actually just takes output from
# adaptor code directly.  
# TODO: Species key also not working so doesn't return the third element
# TODO: Update it so it takes the functional group name as an argument and uses
# it to find files and name outputs
# TODO: Move this to process build model file when it is finished

#' @param data functional group yearly biomass dataframe returned by get_functional_groups.R 
#' @param groups
#' @param key
#' @return A dataframe of only groups with biomass with bodymass range included
#' @return A plot of each virtual species change in biomass over time
#' @return A list of actual species these virtual species might represent

visualise_groups <- function(data, groups, key, name, output_folder) {

# Remove functional groups with no biomass

groups_with_and_without_biomass <- data[[1]] # This will be different when takes adaptor output
groups_with_biomass <- groups_with_and_without_biomass[apply(groups_with_and_without_biomass, 1, function(y) !all(is.na(y))),]
groups_without_biomass <- groups_with_and_without_biomass[apply(groups_with_and_without_biomass, 1, function(y) all(is.na(y))),]

# Get id numbers of groups with biomass

groups_with_biomass_id <- rownames(groups_with_biomass)

# Get group attributes

groups_with_biomass_attributes <- groups_with_biomass_id[groups$group_id %in% 
                                  groups_with_biomass_id, ]

# Convert bodymass ranges from grams to kilograms

groups_with_biomass_attributes <- groups_with_biomass_attributes %>%
                                  mutate(mass_lower_kg = mass_lower/1000) %>%
                                  mutate(mass_upper_kg = mass_upper/1000) %>%
                                  mutate(mass_range_kg = paste(mass_lower_kg, 
                                                               "-", mass_upper_kg))

# Add bodymass range to model outputs
# TODO: Change so it merges or uses left join, so we don't accidentally end up 
# mismatching bodymass ranges with groups

groups_with_biomass_mass_range <- cbind(groups_with_biomass_attributes$mass_range_kg, 
                                        groups_with_biomass_attributes)

# Melt into long format ready to plot

groups_with_biomass_long <- melt(groups_with_biomass_mass_range, id.vars=c("mass_range_kg"))

names(groups_with_biomass_long) <- c("mass_range_kg", "year", "biomass")

plotName <- paste(name,".tiff",sep="")
tiff(file = (paste(output_folder,plotName, sep = "/")), units ="in", width=10, height=5, res=200)

groups_with_biomass_plot <- ggplot(data = groups_with_biomass_long, 
                            aes(x = year, y = biomass, group = mass_range_kg)) +
                            geom_line(aes(colour = mass_range_kg)) +
                            labs(x = "Time (years)", 
                                 y = "biomass") +
                            ggtitle("species biomass over time")

groups_with_biomass_plot

dev.off()

}

#' @param data matrix of rows (species) and columns (timesteps). Can use abundance
#' or biomass
#' @param number of burnin timesteps in months
#' @return saves a plot of total biomass in each bodymass bin, split by functional
#' groups
#' TODO: fix labels (of functional groups and plot names)
#' TODO: see if I can move this to the process output files.  Currently takes 
#' input from the adaptor code though.

plot_functional_groups <- function(data, burnin, output_folder, reference, scenario){
  
  require(dplyr)
  require(stringr)
  require(reshape)
  require(ggplot2)
  
  remove_burn_in <- function(data, burnin) {
    
  data[,(burnin + 1):ncol(data)]
    
  }
  
  data <- as.data.frame(data)
  is.na(data) <- !data
  
  new_data <- remove_burn_in(data, burnin) %>%
              dplyr::mutate(functional_group_index = row.names(.)) %>%
              dplyr::mutate(functional_group = substr(functional_group_index, 
                                            start = 1, stop = 2)) %>%
              dplyr::mutate(total_abundance = rowSums(.[,c(1:2400)])) %>%
              dplyr::select(functional_group, functional_group_index, total_abundance)
  
  plotName <- paste(reference, "_", scenario, "_name",".tiff",sep="")
  tiff(file = (paste(output_folder,plotName, sep = "/")), units ="in", width=10, height=5, res=200)
  
  
  plot <- ggplot(data = new_data, aes( x = functional_group_index, 
                                       y = total_abundance)) +
          geom_bar(stat = 'identity') +
          labs(title = 'test',
               x = "species", y = "abundance") +
          facet_wrap( ~ functional_group, nrow = 3)
  
  plot
  
  dev.off()
  
  
}

