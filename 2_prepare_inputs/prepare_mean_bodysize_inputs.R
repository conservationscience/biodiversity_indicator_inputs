marine_data <- readRDS('N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_adaptor_code/map_of_life/Harvesting_carnivores/201_BuildModel/MassBinsOutputs_harvest_carnivores_201_0_Cell0_biomass.rds')

marine_data <- marine_data[,12000:15601]

saveRDS(marine_data, "N:\\Quantitative-Ecology\\Indicators-Project\\marine_data_201_0.rds")


inputs <-  "N:/Quantitative-Ecology/Simone/marine_data_201_0.rds"
outputs <- "N:/Quantitative-Ecology/Simone"


# Calculate total biomass over time

convert_timesteps <- function(dataframe, interval, func){
  
  monthly_matrix <- t(dataframe)
  
  n <- interval
  
  time_converted_matrix <- t(aggregate(monthly_matrix,list(rep(1:(nrow(monthly_matrix) %/%
                                                                    n + 1), each = n, len = nrow(monthly_matrix))), func, na.rm = TRUE))
  
  time_converted_matrix <- time_converted_matrix[-1,]
  
  time_converted_matrix[is.nan(time_converted_matrix)] = NA
  
  return(time_converted_matrix)
  
}

source("C:\\Users\\ssteven\\Desktop\\Serengeti-analysis\\model_outputs_to_indicator_inputs\\2_prepare_inputs\\prepare_proportion_total_biomass_inputs.R")

marine_data[marine_data == -9999] <- NA

mean_annual_biomass <- convert_timesteps(marine_data, 12, mean)


total_biomass_vector <- colSums(marine_data, na.rm = TRUE)


#total_biomass_vector <- total_biomass_monthly[seq(1, length(total_biomass_monthly), 12)]

timestep <- seq(1,length(total_biomass_vector), 1)

total_biomass <- as.data.frame(cbind(timestep, total_biomass_vector))

total_biomass <- total_biomass %>%
  dplyr::mutate(biomass_scaled = total_biomass_vector/
                  total_biomass_vector[1])

total_biomass <- total_biomass[1:200,] # remove recovery years

# Get abundance

abundance_data <- readRDS('N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_adaptor_code/map_of_life/Harvesting_carnivores/201_BuildModel/MassBinsOutputs_harvest_carnivores_201_0_Cell0_abundance.rds')

abundance_data <- abundance_data[,12000:15601]

abundance_data[abundance_data == -9999] <- NA

mean_annual_abundance <- convert_timesteps(abundance_data, 12, mean)

# Get mean body size

# Add body size

breaks <- read.csv('N:/Quantitative-Ecology/Simone/MassBinDefinitions.csv')

bodymass_bins_lower <- c(breaks[,1])
bodymass_bins_upper <- c(1000000000,
                         bodymass_bins_lower[bodymass_bins_lower !=
                                               min(bodymass_bins_lower)]) # Add a huge upper limit

bodymass <- rep(bodymass_bins_upper, 6)

x <- as.data.frame(cbind(bodymass_bins_upper, mean_annual_abundance))

y <- sapply(x, '*', x$bodymass_bins_upper)

z <- y[,2:ncol(y)]

mean_bodymass_vector <- colMeans(z, na.rm = TRUE)
mean_bodymass <- as.data.frame(cbind(c(1:301), mean_bodymass_vector))

mean_bodymass <- mean_bodymass %>%
  dplyr::mutate(bodymass_scaled = mean_bodymass_vector/
                  mean_bodymass_vector[1])