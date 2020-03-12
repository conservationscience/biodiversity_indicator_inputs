
# https://github.com/conservationscience/model_outputs_to_indicator_inputs.git

# extract relevant data from each replicate
# NOTE - run process_species_list before this function
# NOTE - you have to create the [name of scenario] folder yourself
# buildmodel_folder should be:
# Indicators-Project/[name of region]/Inputs_to_adaptor_code/Madingley_simulation_outputs/[name of scenario]/[XX]_BuildModel
# output folder should be:
# Indicators-Project/[name of region]/Outputs_from_adaptor_code/[name of species list file]/[name of scenario]
# 
# 


process_buildmodel_folder <- function( buildmodel_folder, output_folder ) {
  
  # Check if the folder has already been processed
  
  folder_name <- basename(buildmodel_folder)
  output_folder_exists <- folder_name %in% list.dirs( output_folder, full.names = FALSE, recursive = FALSE )
  processing_complete <- "groups.csv" %in% list.files(file.path(output_folder, folder_name))

  if (!processing_complete) {
  
  # load comparable taxa
  comparable_taxa <- readRDS( file.path( dirname( output_folder ), "comparable_taxa.rds" ) )
  
  
  # add the [XX]_BuildModel folder to the output folder, so that batches of replicates are separated
  output_folder <- file.path( output_folder, basename( buildmodel_folder ) )
  
  # get name of the folder containing the replicates
  directories <- list.dirs( buildmodel_folder, full.names = FALSE, recursive = FALSE )
  replicate_folder = file.path( buildmodel_folder, directories[!(directories %in% c("input"))] )
  
  # get and save groups for this batch
  groups <- madingley_get_groups(
    read.csv(file.path(replicate_folder, "CohortFunctionalGroupDefinitions.csv")), # should not change between runs
    read.csv(file.path(buildmodel_folder, "input", "Model setup", "Ecological definition files", "MassBinDefinitions.csv") ), # could change if user changes mass bins
    read.csv(file.path(buildmodel_folder, "input", "Model setup", 'SimulationControlParameters.csv')) #Can change between runs
  )
  
  write.csv( groups, file.path( output_folder, "groups.csv") )
  saveRDS(   groups, file.path( output_folder, "groups.rds") )
  
  # match the groups to species for this batch
  species_and_groups_key <- madingley_get_species_and_groups_key( comparable_taxa, groups )
  
  write.csv( species_and_groups_key, file.path( output_folder, "species_and_groups_key.csv") )
  saveRDS  ( species_and_groups_key, file.path( output_folder, "species_and_groups_key.rds") )
  
  
  # Iterate over each MassBinFile, process both Abundance and Biomass data, and output results to output folder
  ListMassBinsFiles <- function(resultsDir){
    files<-dir(resultsDir)
    files<-files[grep("MassBinsOutputs",files)]
    files<-files[grep("Cell",files)]
    
    return(files)
  }
  
  massbin_files <- ListMassBinsFiles( replicate_folder )
  
  i <- 1
  for( file in massbin_files ) {
    
    input_file_path <- file.path( replicate_folder, file )
    biomass_output_file_path <- file.path( output_folder, sub( ".nc", "_biomass.rds", file ) )
    abundance_output_file_path <- file.path( output_folder, sub( ".nc", "_abundance.rds", file ) )
   
    cat( paste0( "processing file ", i, " of ", length( massbin_files ), "...\n") )
    
    log_biomass_through_time <- madingley_get_biomass_of_groups( input_file_path, groups )
    saveRDS( log_biomass_through_time, file = biomass_output_file_path )
    write.csv( log_biomass_through_time, file = sub( ".rds", ".csv", biomass_output_file_path ) )
    
    # remove the large variable before loading another one. Probably doesn't increase speed,
    # but shown here as an example of how you could keep your code using the minimum amount of 
    # memory necessary, in case you had variables that were say 500mb-1gb
    rm( log_biomass_through_time )
    
    log_abundance_through_time <- madingley_get_abundance_of_groups( input_file_path, groups )
    saveRDS( log_abundance_through_time, file = abundance_output_file_path )
    write.csv( log_abundance_through_time, file = sub( ".rds", ".csv", abundance_output_file_path ) )
    rm( log_abundance_through_time )
    
    i <- i + 1
  }
  
  rm(file)
  
  ListCellOutputFiles <- function(resultsDir){
    
    files <- dir(resultsDir)
    files <- files[grep("BasicOutputs",files)]
    files <- files[grep("Cell",files)]
    
    return(files)
    
  }
  
  basicoutput_files <- ListCellOutputFiles( replicate_folder )
  
  i <- 1
  for( file in basicoutput_files ) {
    
    input_file_path <- file.path( replicate_folder, file )
    autotroph_output_file_path <- file.path( output_folder, sub( ".nc", "_autotroph.rds", file ) )
    
    cat( paste0( "processing", folder_name, "file", i, " of ", length( basicoutput_files ), "...\n") )
    
    log_autotroph_biomass_through_time <- madingley_get_autotroph_biomass( input_file_path )
    saveRDS( log_autotroph_biomass_through_time, file = autotroph_output_file_path )
    write.csv( log_autotroph_biomass_through_time, file = sub( ".rds", ".csv", autotroph_output_file_path ) )
    
    # remove the large variable before loading another one. Probably doesn't increase speed,
    # but shown here as an example of how you could keep your code using the minimum amount of 
    # memory necessary, in case you had variables that were say 500mb-1gb
    rm( log_autotroph_biomass_through_time )
   
    i <- i + 1
  }
  cat( "done\n" )
  rm( i )
  
  } else {

    print(paste(folder_name, "folder has already been processed", sep = " "))

  }
}

