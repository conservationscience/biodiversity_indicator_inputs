
# Git repo "C:\\Users\\ssteven\\Dropbox\\Deakin\\Serengeti-analysis\\model_outputs_to_indicator_inputs\\1_process_outputs"





#' TODO: Automate the vertical lines representing impact start and end
#' 
#' This function plots the model outputs by functional group so you can have a 
#' quick look at the outputs in more detail than the log biomass plot that
#' summarises replicates
#' 
#' @param simulation_path A string that denotes the file path where the
#' processed abundance and biomass data lives.  It should be one directory
#' from a single simulation, that contains outputs from multiple replicates
#' 
#' @param simulation_number A string that denotes the number of the simulation
#' you are processing outputs for
#' 
#' @return Doesn't return anything to the environment, but will create 
#' a 'functional group plots' directory in the simulation directory, to
#' which it saves three different plots (heterotroph abundance and biomass and
#' autotroph biomass) per replicate (so if you have five replicates, you'll end
#' up with 15 plots in your new directory).  If the function has already been
#' run for that simulation, it won't process it again.
#' 
#' For the raw autotroph data that is read in, the groups are measured in  
#' biomass density [kg km-2]: The wet matter biomass density of the set of organisms in <group>. 

#' For all the heterotroph data read in, the group's abundance and biomass are 
#' measured in log values

# dev_mode <- TRUE
# 
# if(dev_mode == TRUE) {
# 
# burnin <- 1*12
# 
# simulation_path <- "N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Outputs_from_adaptor_code\\map_of_life\\Test_runs\\ae_BuildModel"
# output_path <- simulation_path
# simulation_number <- "ae"
#   
# burnin <- 0
# startimpact <- 12
# endimpact <- 24
# sample_year <- 1
# n <- 1
# 
# } else {
# 
# simulation_path <- "N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Outputs_from_adaptor_code\\map_of_life\\Harvesting_herbivores\\301_BuildModel"
# output_path <- simulation_path
# 
# burnin <- 1000 * 12
# startimpact <- 1100 * 12
# endimpact <- 1200 * 12
# sample_year <- 5*12
# n <- 120
# 
# }
# 
# plot_functional_groups(simulation_path, simulation_number)

plot_functional_groups <- function(simulation_path, simulation_number) {
  
  require(tidyverse)
  require(reshape2)
  require(ggplot2)

output_path <- file.path(simulation_path, paste(simulation_number,
                                                "_functional_group_plots",sep=""))

if( !dir.exists( file.path(output_path) ) ) {
  
  dir.create( file.path(output_path), recursive = TRUE )

# Find files

files <- list.files(simulation_path, recursive = FALSE) 

# Plot log abundance of functional groups ----

abundance_files <- paste(simulation_path, files[str_detect(files, "abundance.rds")], sep = "\\")
replicate_numbers <- 0:(length(abundance_files) - 1)

#abundance_files <- abundance_files[1:2] # TEMP CODE - subset so you don't read everything in

plot_names <- paste("Simulation", simulation_number, "replicate", replicate_numbers, 
                    "log abundance by heterotrophic functional group", sep = " ")

simulation_abundance_data <- lapply(abundance_files, readRDS)

#simulation_abundance_data <- lapply(simulation_abundance_data, exp)


# Remove burn in and sample every n years to reduce data for plotting

subset_simulation_data <- list()

for (i in seq_along(simulation_abundance_data)) {
  
  temp <- simulation_abundance_data[[i]][,burnin:ncol(simulation_abundance_data[[i]])]
  
  sample <- seq(1, ncol(temp), by = n)
  
  subset_simulation_data[[i]] <- temp[,sample]
  
}

rm(temp)

# Melt into long format

subset_simulation_data_long <- lapply(subset_simulation_data, melt)

# Replace no value with NA

for (i in seq_along(subset_simulation_data_long)) {
  
  subset_simulation_data_long[[i]][subset_simulation_data_long[[i]] == -9999] <- NA
  
}

# head(subset_simulation_data_long[[1]])

# Add grouping variables to prepare data for plotting

plot_data <- list()

for (i in seq_along(subset_simulation_data_long)) {
  
  names(subset_simulation_data_long[[i]]) <- c("functional_group", "timestep", "abundance")
  
  subset_simulation_data_long[[i]]$functional_group <- as.character(subset_simulation_data_long[[i]]$functional_group)
  
  subset_simulation_data_long[[i]] <- subset_simulation_data_long[[i]][!is.na(subset_simulation_data_long[[i]]$abundance),]
  
  plot_data[[i]] <- subset_simulation_data_long[[i]] %>%
    mutate(group = 
             str_sub(subset_simulation_data_long[[i]]$functional_group, 
                     start = 1, end= 2)) %>%
    # filter(timestep >= burnin) %>%
    mutate(functional_group_name = ifelse(group == 10, "herbivorous endotherms",
                                   ifelse(group == 11, "carnivorous  endotherms",
                                   ifelse(group == 12, "omnivorous  endotherms",
                                   ifelse(group == 13, "herbivorous ectotherms", # combine iteroparous and semelparous ectotherms
                                   ifelse(group == 14, "carnivorous ectotherms",
                                   ifelse(group == 15, "omnivorous ectotherms", 
                                   "NA"))))))) %>%
    arrange(functional_group, timestep) %>%
    mutate(bodymass_bin = str_sub(functional_group, start= -2)) %>%
    group_by(group, timestep) %>%
    mutate(mean_group_abundance = mean(abundance))
  
}

# dim(plot_data[[1]])
# 
# head(plot_data[[1]])

# Plot

abundance_plots <- list()

for (i  in seq_along(plot_data)) {
  
  abundance_plots[[i]] <- ggplot() +
    geom_path(aes(y = abundance, x = timestep, colour = bodymass_bin),
              data =  plot_data[[i]]) +
    geom_line(aes(y = mean_group_abundance, x = timestep), 
              data = plot_data[[i]], size = 1) +
    theme(legend.position = "none") +
    geom_vline(xintercept = 10, color = "red") +
    geom_vline(xintercept = 20 , color = "dark green") +
    facet_wrap(~ functional_group_name, ncol = 3) +
    ggtitle(plot_names[[i]])
  
  abundance_plots[[i]]
  
}

# abundance_plots[[1]]
# abundance_plots[[2]]
# abundance_plots[[3]]
# abundance_plots[[4]]
# abundance_plots[[5]]

for (i in seq_along(abundance_plots)) {
  
  ggsave(file.path(output_path, paste(plot_names[[i]],".pdf")), abundance_plots[[i]])

}

# Plot log biomass of functional groups ----

biomass_files <- paste(simulation_path, files[str_detect(files, "biomass.rds")], sep = "\\")


#biomass_files <- biomass_files[1:2] # TEMP CODE - subset so you don't read everything in

biomass_plot_names <- paste("Simulation", simulation_number, "replicate", replicate_numbers, 
                    "log biomass by heterotrophic functional group", sep = " ")

simulation_biomass_data <- lapply(biomass_files, readRDS)

#simulation_biomass_data <- lapply(simulation_biomass_data, exp)


# Remove burn in and sample every n years to reduce data for plotting

subset_biomass_simulation_data <- list()

for (i in seq_along(simulation_biomass_data)) {
  
  temp <- simulation_biomass_data[[i]][,burnin:ncol(simulation_biomass_data[[i]])]
  
  sample <- seq(1, ncol(temp), by = n)
  
  subset_biomass_simulation_data[[i]] <- temp[,sample]
  
}

rm(temp)

# Melt into long format

subset_biomass_simulation_data_long <- lapply(subset_biomass_simulation_data, melt)

# Replace no value with NA

for (i in seq_along(subset_biomass_simulation_data_long)) {
  
  subset_biomass_simulation_data_long[[i]][subset_biomass_simulation_data_long[[i]] == -9999] <- NA
  
}

# head(subset_biomass_simulation_data_long[[1]])

# Add grouping variables to prepare data for plotting

biomass_plot_data <- list()

for (i in seq_along(subset_biomass_simulation_data_long)) {
  
  names(subset_biomass_simulation_data_long[[i]]) <- c("functional_group", "timestep", "biomass")
  
  subset_biomass_simulation_data_long[[i]]$functional_group <- as.character(subset_biomass_simulation_data_long[[i]]$functional_group)
  
  subset_biomass_simulation_data_long[[i]] <- subset_biomass_simulation_data_long[[i]][!is.na(subset_biomass_simulation_data_long[[i]]$biomass),]
  
  biomass_plot_data[[i]] <- subset_biomass_simulation_data_long[[i]] %>%
    mutate(group = 
             str_sub(subset_biomass_simulation_data_long[[i]]$functional_group, 
                     start = 1, end= 2)) %>%
    # filter(timestep >= burnin) %>%
    mutate(functional_group_name = ifelse(group == 10, "herbivorous endotherms",
                                          ifelse(group == 11, "carnivorous  endotherms",
                                                 ifelse(group == 12, "omnivorous  endotherms",
                                                        ifelse(group == 13, "herbivorous ectotherms", # combine iteroparous and semelparous ectotherms
                                                               ifelse(group == 14, "carnivorous ectotherms",
                                                                      ifelse(group == 15, "omnivorous ectotherms", 
                                                                             "NA"))))))) %>%
    arrange(functional_group, timestep) %>%
    mutate(bodymass_bin = str_sub(functional_group, start= -2)) %>%
    group_by(group, timestep) %>%
    mutate(mean_group_biomass = mean(biomass))
  
}

# dim(biomass_plot_data[[1]])
# 
# head(biomass_plot_data[[1]])

# Plot

biomass_plots <- list()

for (i  in seq_along(biomass_plot_data)) {
  
  biomass_plots[[i]] <- ggplot() +
    geom_path(aes(y = biomass, x = timestep, colour = bodymass_bin),
              data =  biomass_plot_data[[i]]) +
    geom_line(aes(y = mean_group_biomass, x = timestep), 
              data = biomass_plot_data[[i]], size = 1) +
    theme(legend.position = "none") +
    geom_vline(xintercept = 10, color = "red") +
    geom_vline(xintercept = 20 , color = "dark green") +
    facet_wrap(~ functional_group_name, ncol = 3) +
    ggtitle(biomass_plot_names[[i]])
  
  biomass_plots[[i]]
  
}

# biomass_plots[[1]]
# biomass_plots[[2]]
# biomass_plots[[3]]
# biomass_plots[[4]]
# biomass_plots[[5]]

for (i in seq_along(biomass_plots)) {
  
  ggsave(file.path(output_path, paste(biomass_plot_names[[i]],".pdf")), biomass_plots[[i]])
  
}

# Plot autotroph biomass ----

autotroph_files <- paste(simulation_path, files[str_detect(files, "autotroph.rds")], 
                         sep = "\\")

auto_plot_names <- paste("Simulation", simulation_number, "replicate", replicate_numbers, 
                         "log biomass of autotrophs", sep = " ")

simulation_auto_data <- lapply(autotroph_files, readRDS)

simulation_auto_data <- lapply(simulation_auto_data, as.matrix)

# Remove burn in and sample every n years to reduce data for plotting

subset_simulation_auto_data <- list()

for (i in seq_along(simulation_auto_data)) {
  
  temp <- simulation_auto_data[[i]][,burnin:ncol(simulation_auto_data[[i]])]
  
  sample <- seq(1, ncol(temp), by = n)
  
  subset_simulation_auto_data[[i]] <- temp[,sample]
  
}

rm(temp)

# Melt into long format

subset_simulation_auto_data_long <- lapply(subset_simulation_auto_data, melt)

# Replace no value with NA

for (i in seq_along(subset_simulation_auto_data_long)) {
  
  subset_simulation_auto_data_long[[i]][subset_simulation_auto_data_long[[i]] == -9999] <- NA
  
}

# head(subset_simulation_auto_data_long[[1]])

# Name columns to prepare data for plotting

auto_plot_data <- list()

for (i in seq_along(subset_simulation_auto_data_long)) {
  
  names(subset_simulation_auto_data_long[[i]]) <- c("functional_group", 
                                                    "timestep", "biomass")
  
  auto_plot_data[[i]] <- subset_simulation_auto_data_long[[i]] %>%
                         mutate(new_timestep = as.numeric(substring(timestep,2))) 
  
  
}

# head(auto_plot_data[[1]])
# class(auto_plot_data[[1]]$new_timestep)

# Plot

autotroph_plots <- list()

for (i  in seq_along(auto_plot_data)) {
  
    autotroph_plots[[i]] <- ggplot() +
    geom_path(aes(y = biomass, x = new_timestep, colour = functional_group),
              data =  auto_plot_data[[i]]) +
    theme(legend.position = "right") +
    geom_vline(xintercept = 10, color = "red") +
    geom_vline(xintercept = 30 , color = "dark green") +
    ggtitle(auto_plot_names[[i]])
  
  autotroph_plots[[i]]
  
}

# 
# autotroph_plots[[1]]
# autotroph_plots[[2]]
# autotroph_plots[[3]]
# autotroph_plots[[4]]
# autotroph_plots[[5]]

for (i in seq_along(biomass_plots)) {
  
  ggsave(file.path(output_path, paste(auto_plot_names[[i]],".pdf")), autotroph_plots[[i]])
  
}

} else {
  
  print(paste("Functional group plots for Simulation", simulation_number, 
              "already exist"))
}

}
