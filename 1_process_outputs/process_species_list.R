


# collect trait information from a list of species, and save it into output_folder
# output_folder should be:
# Indicators-Project/Serengeti/Outputs_from_adaptor_code/[name of species list file]
process_species_list <- function( species_list, databases, output_folder ) {
  
  #########################################
  # Get raw trait information from functionaltraits package
  #########################################
  
  #### WARNING: This function takes about three hours or more to run.
  #### TODO: Save the output of this function so that it doesn't have to be run twice
  #### TODO: Optimise the functionaltraits package so that if a large species list is provided,
  ####       it can break it down into smaller lists of 100 species and process them independently, 
  ####       in case network connection is lost. 
  ####       functionaltraits should save the raw data; the get_trait_data function should save the 
  ####       processed data. Should be implemented as a function in functionaltraits that wraps find_species_traits
  ####       also find_species_traits should be renamed or something. Should be removed, and instead
  ####       an option or something in Databases should be added to check synonyms as well.
  
  raw_trait_data_filename <- file.path( output_folder, "raw_trait_data.rds" )
  
  if( file.exists( raw_trait_data_filename ) ) {
    print( "reusing exising trait data: " )
    print( raw_trait_data_filename )
    raw_trait_data <- readRDS( raw_trait_data_filename )
  } else {
    
    if( !databases$ready() ) stop( "the databases are not ready. Try running databases$initialise() on the databases argument.")
    
    # check that all of the species names provided are unique
    # if there are duplicate species names in the list provided, then the user must
    # remove the duplicate names themselves. This is to ensure that this function returns the 
    # same number of species as it was provided with (which would likely be an assumption of users)
    
    # TODO: this check should actually be in functionaltraits::Databases::search(), not here.
    if( !all(!duplicated(species_list)) ) stop( "there were duplicated species names in the species list. Try removing them with unique()")
    
    search_results <- functionaltraits::find_species_traits( databases, species_list )
    
    raw_trait_data <- search_results$results
    
    saveRDS( raw_trait_data, raw_trait_data_filename )
    write.csv( raw_trait_data, file = file.path( output_folder, "raw_trait_data.csv" ) )
    
  }
    
    

  #########################################
  # Process it and save the outputs
  #########################################
  processed_trait_data <- madingley_process_trait_data( raw_trait_data )
  saveRDS  ( processed_trait_data, file = file.path( output_folder, "processed_trait_data.csv" ) )
  write.csv( processed_trait_data, file = file.path( output_folder, "processed_trait_data.csv" ) )
  
  rm( raw_trait_data_filename )
  
  
  
  #########################################
  # Work out which species we will be able to match with the model
  # and make a list of the species that we don't have enough information for
  #########################################
  
  
  ### TODO: What do we want to do with the species that weren't found by the taxonomic system?
  ### TODO: Should these unmatched species still be searched for in the databases?
  
  ### NOTE: more species are added to unmatched_taxa when some species are missing information
  
  indexes_of_comparable_taxa <- which(
    processed_trait_data$found == TRUE
    & !is.na( processed_trait_data$heterotroph_autotroph )
    & !is.na( processed_trait_data$nutrition_source )
    & !is.na( processed_trait_data$endo_ectotherm )
    & !is.na( processed_trait_data$bodymass )
  )
  
  comparable_taxa <- processed_trait_data[ indexes_of_comparable_taxa, ]
  incomparable_taxa <- processed_trait_data[ !(1:nrow( processed_trait_data ) %in% indexes_of_comparable_taxa), ]
  
  comparable_taxa$species_id <- 1:nrow(comparable_taxa)
  
  rm( indexes_of_comparable_taxa )
  
  #### Save outputs
  write.csv( comparable_taxa, file = file.path( output_folder, "comparable_taxa.csv" ) )
  saveRDS( comparable_taxa, file =  file.path( output_folder, "comparable_taxa.rds" ) )
  write.csv( incomparable_taxa, file = file.path( output_folder, "incomparable_taxa.csv" ) )
  saveRDS(  incomparable_taxa, file = file.path( output_folder, "incomparable_taxa.rds" ) )
}

