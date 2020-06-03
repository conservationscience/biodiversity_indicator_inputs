
# Git repo "C:\\Users\\ssteven\\Dropbox\\Deakin\\Serengeti-analysis\\model_outputs_to_indicator_inputs\\1_process_outputs"


#' TODO: Automate the vertical lines representing impact start and end for autotrophs
#' TODO: IMPORTANT - add a plot of mean autototroph values
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
#' 
#' @param n Integer that denotes which timestep you want to sample. n = 1 will
#' include all timesteps, n = 12 will select the last timestep in each year, and so on.
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
# simulation_path <- "N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Outputs_from_adaptor_code\\map_of_life\\Test_runs\\ae_BuildModel"
# simulation_number <- "ae"
# burnin <- 0
# n <- 1
# 
# } else {
# 
# simulation_path <- "N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Outputs_from_adaptor_code\\map_of_life\\Harvesting_carnivores\\209_BuildModel"
# output_path <- simulation_path
# 
# burnin <- 1000 * 12
# startimpact <- 100
# endimpact <- 200
# sample_year <- 5*12
# n <- 96
# simulation_number <- "209"
# 
# }
# 
# plot_functional_groups(simulation_path, simulation_number, n)

plot_simulation_functional_groups <- function(simulation_path, simulation_number, n) {
  
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

abundance_files <- paste(simulation_path, files[str_detect(files, 
                                                           "abundance.rds")], 
                         sep = "\\")

replicate_numbers <- 0:(length(abundance_files) - 1)

plot_names <- paste("Simulation", simulation_number, "replicate", replicate_numbers, 
                    "log abundance by heterotrophic functional group", sep = " ")

simulation_abundance_data <- lapply(abundance_files, readRDS)


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

# Plot

abundance_plots <- list()

for (i  in seq_along(plot_data)) {
  
  abundance_plots[[i]] <- ggplot() +
    geom_path(aes(y = abundance, x = timestep, colour = bodymass_bin),
              data =  plot_data[[i]]) +
    geom_line(aes(y = mean_group_abundance, x = timestep), 
              data = plot_data[[i]], size = 1) +
    theme(legend.position = "none") +
    geom_vline(xintercept = max(plot_data[[i]]$timestep)/3, color = "red") +
    geom_vline(xintercept = ( max(plot_data[[i]]$timestep)/3) *2 ,
               color = "dark green") +
    facet_wrap(~ functional_group_name, ncol = 3) +
    ggtitle(plot_names[[i]])
  
  abundance_plots[[i]]
  
}

for (i in seq_along(abundance_plots)) {
  
  ggsave(file.path(output_path, paste(plot_names[[i]],".pdf")), abundance_plots[[i]])

}

# Plot log biomass of functional groups ----

biomass_files <- paste(simulation_path, files[str_detect(files, "biomass.rds")], 
                       sep = "\\")

biomass_plot_names <- paste("Simulation", simulation_number, "replicate", 
                            replicate_numbers, 
                    "log biomass by heterotrophic functional group", sep = " ")

simulation_biomass_data <- lapply(biomass_files, readRDS)

# Remove burn in and sample every n years to reduce data for plotting

subset_biomass_simulation_data <- list()

for (i in seq_along(simulation_biomass_data)) {
  
  temp <- simulation_biomass_data[[i]][,burnin:ncol(simulation_biomass_data[[i]])]
  
  temp[temp == -9999] <- NA
  
  sample <- seq(1, ncol(temp), by = n)
  
  subset_biomass_simulation_data[[i]] <- temp[,sample]
  
}

rm(temp)

# Melt into long format

subset_biomass_simulation_data_long <- lapply(subset_biomass_simulation_data, melt)

# Add grouping variables to prepare data for plotting

biomass_plot_data <- list()

for (i in seq_along(subset_biomass_simulation_data_long)) {
  
  names(subset_biomass_simulation_data_long[[i]]) <- c("functional_group", 
                                                       "timestep", "biomass")
  
  subset_biomass_simulation_data_long[[i]]$functional_group <- as.character(
    subset_biomass_simulation_data_long[[i]]$functional_group)
  
  subset_biomass_simulation_data_long[[i]] <- subset_biomass_simulation_data_long[[i]][!is.na(subset_biomass_simulation_data_long[[i]]$biomass),]
  
  biomass_plot_data[[i]] <- subset_biomass_simulation_data_long[[i]] %>%
    mutate(group = 
             str_sub(subset_biomass_simulation_data_long[[i]]$functional_group, 
                     start = 1, end= 2)) %>%
    mutate(functional_group_name = ifelse(group == 10, "herbivorous endotherms",
                                          ifelse(group == 11, "carnivorous  endotherms",
                                                 ifelse(group == 12, "omnivorous  endotherms",
                                                        ifelse(group == 13, "herbivorous ectotherms", # combine iteroparous and semelparous ectotherms
                                                               ifelse(group == 14, "carnivorous ectotherms",
                                                                      ifelse(group == 15, "omnivorous ectotherms", 
                                                                             "NA"))))))) %>%
    arrange(functional_group, timestep) %>%
    mutate(bodymass_bin = str_sub(functional_group, start = -2)) %>%
    group_by(group, timestep) %>%
    mutate(mean_group_biomass = mean(biomass))
  
}

# Plot

biomass_plots <- list()

for (i  in seq_along(biomass_plot_data)) {
  
  biomass_plots[[i]] <- ggplot() +
    geom_path(aes(y = biomass, x = timestep, colour = bodymass_bin),
              data =  biomass_plot_data[[i]]) +
    geom_line(aes(y = mean_group_biomass, x = timestep), 
              data = biomass_plot_data[[i]], size = 1) +
    theme(legend.position = "none") +
    geom_vline(xintercept = max(biomass_plot_data[[i]]$timestep)/3, 
               color = "red") +
    geom_vline(xintercept = (max(biomass_plot_data[[i]]$timestep)/3)*2 ,
               color = "dark green") +
    facet_wrap(~ functional_group_name, ncol = 3) +
    ggtitle(biomass_plot_names[[i]])
  
  biomass_plots[[i]]
  
}



for (i in seq_along(biomass_plots)) {
  
  ggsave(file.path(output_path, paste(biomass_plot_names[[i]],".pdf")), biomass_plots[[i]])
  
}

# Plot mean biomass ----

mean_biomass <- Reduce("+", lapply(subset_biomass_simulation_data, 
                        function(x) replace(x, is.na(x), 0))) / 
                Reduce("+", lapply(subset_biomass_simulation_data, Negate(is.na)))


mean_biomass_long <- melt(mean_biomass)
names(mean_biomass_long) <- c("functional_group", "timestep", "mean_biomass")


mean_biomass_plot_data <- mean_biomass_long %>%
        mutate(bodymass_bin = str_sub(functional_group, start = -2)) %>%
        mutate(group = str_sub(mean_biomass_long$functional_group, 
                               start = 1, end = 2)) %>%
        arrange(functional_group, timestep) %>%
        group_by(group, timestep) %>%
        dplyr::summarise(group_mean = mean(mean_biomass, na.rm = TRUE),
                         group_sd = sd(mean_biomass, na.rm = TRUE),
                         n = n(),
                         group_se = group_sd/sqrt(n),
                         group_lb = group_mean - (1.96 * group_se),
                         group_ub = group_mean + (1.96 * group_se)) %>%
        mutate(functional_group_name = ifelse(group == 10, "herbivorous endotherms",
                                       ifelse(group == 11, "carnivorous  endotherms",
                                       ifelse(group == 12, "omnivorous  endotherms",
                                       ifelse(group == 13, "herbivorous ectotherms", # combine iteroparous and semelparous ectotherms
                                       ifelse(group == 14, "carnivorous ectotherms",
                                       ifelse(group == 15, "omnivorous ectotherms", 
                                       "NA"))))))) 

mean_biomass_plot <- ggplot() +
                     geom_path(aes(y = group_mean, x = timestep, 
                                    colour = functional_group_name),
                                data =  mean_biomass_plot_data) +
                     theme(legend.position = "none") +
                     geom_vline(xintercept = max(mean_biomass_plot_data$timestep)/3,
                                color = "red") +
                     geom_vline(xintercept = (max(mean_biomass_plot_data$timestep)/3)*2,
                                color = "dark green") +
                     ggtitle(paste("Simulation", simulation_number, 
                                   "mean biomass by functional group")) +
                     theme(panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank()) + 
                     theme(legend.position="right") +
                     facet_wrap(~ functional_group_name)



mean_biomass_plot_95_CI <- mean_biomass_plot + 
                           geom_ribbon(aes(x = timestep,
                                           ymin = group_lb, 
                                           ymax = group_ub, 
                                           fill = functional_group_name), 
                                           alpha = 0.4, 
                                           data = mean_biomass_plot_data) 
mean_biomass_plot_95_CI

ggsave(file.path(output_path, paste("Simulation", simulation_number,
                                    "mean biomass by functional group.pdf", 
                                    sep = " ")), 
       mean_biomass_plot_95_CI)

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


# Name columns to prepare data for plotting

auto_plot_data <- list()

for (i in seq_along(subset_simulation_auto_data_long)) {
  
  names(subset_simulation_auto_data_long[[i]]) <- c("functional_group", 
                                                    "timestep", "biomass")
  
  auto_plot_data[[i]] <- subset_simulation_auto_data_long[[i]] %>%
                         mutate(new_timestep = as.numeric(substring(timestep,2))) 
                        
}

# Plot

autotroph_plots <- list()

for (i  in seq_along(auto_plot_data)) {
  
  range <- range(auto_plot_data[[i]]$new_timestep, na.rm = TRUE)
  range <- (c(range[1]: range[2]))
  start <- min(range) + 100*12 
  end <- max(range) - 100 *12
    
  autotroph_plots[[i]] <- ggplot() +
    geom_path(aes(y = biomass, x = new_timestep, colour = functional_group),
              data =  auto_plot_data[[i]]) +
    theme(legend.position = "right") +
    geom_vline(xintercept =  start, 
               color = "red") +
    geom_vline(xintercept = end , 
               color = "dark green") +
    ggtitle(auto_plot_names[[i]])
  
  autotroph_plots[[i]]
  
}

for (i in seq_along(autotroph_plots)) {
  
  ggsave(file.path(output_path, paste(auto_plot_names[[i]],".pdf")), 
         autotroph_plots[[i]])
  
}


print(paste("Plots for", simulation_number, 
            "saved to functional group plots folder"))
} else {
  
  print(paste("Functional group plots for Simulation", simulation_number, 
              "already exist"))
  }
}



#' TODO: Think about where all these outputs should be saved
#' 
#' This function summarises the biomass data of all replicates in a single 
#' scenario. It should produce a plot of mean biomass over time for the six functional
#' groups which it saves to the scenario folder, along with the biomass time
#' series plot data
#' 
#' @param processed_scenario_path A string that denotes the file path for
#' that scenario where all the individual simulation folders reside
#' 
#' @param scenario A string - the name of the scenario 
#' 
#' @param n Integer that denotes which timestep you want to sample. n = 1 will
#' include all timesteps, n = 12 will select the last timestep in each year, and so on.
#' 
#' @param burnin Integer that denotes how many monthly timesteps the burnin lasts
#' for, so it can be removed from the data (eg if the simulation has a burnin
#' period of 100 years, burnin should be 1000 * 12 = 12000)
#' 
#' @return should produce a plot of mean biomass over time for the six functional
#' groups which it saves to the scenario folder, along with the biomass time
#' series plot data

# processed_scenario_path <- "N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Outputs_from_adaptor_code\\map_of_life\\Test_runs"
# scenario <- "Test runs"
# n <- 1
# burnin <- 0
# 
# testruns <- plot_scenario_functional_groups(processed_scenario_path, scenario, burnin,n)
# testruns
# 
# processed_scenario_path <- "N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Outputs_from_adaptor_code\\map_of_life\\Harvesting_carnivores"
# scenario <- "Harvesting carnivores"
# n <- 96
# burnin <- 1000 * 12
# 
# carnivores <- plot_scenario_functional_groups(processed_scenario_path, scenario, burnin,n)
# carnivores
# 
# processed_scenario_path <- "N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Outputs_from_adaptor_code\\map_of_life\\Harvesting_herbivores"
# scenario <- "Harvesting herbivores"
# 
# herbivores <- plot_scenario_functional_groups(processed_scenario_path, scenario, burnin,n)
# herbivores
# 
# processed_scenario_path <- "N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Outputs_from_adaptor_code\\map_of_life\\Land_use"
# scenario <- "Land use"
# 
# landuse <- plot_scenario_functional_groups(processed_scenario_path, scenario, burnin,n)
# landuse
# 
# processed_scenario_path <- "N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Outputs_from_adaptor_code\\map_of_life\\Baseline"
# scenario <- "Baseline"
# 
# baseline <- plot_scenario_functional_groups(processed_scenario_path, scenario, burnin,n)
# baseline

plot_scenario_functional_groups <- function(processed_scenario_path, 
                                            scenario, 
                                            burnin, n) {
  
  require(tidyverse)
  require(reshape2)
  require(ggplot2)
  
output_path <- processed_scenario_path
  
scenario_simulation_paths <- list.dirs(processed_scenario_path, recursive = FALSE) 
  
  # Get the file paths for each simulation within the scenario
  
  replicate_paths <- list()
  replicate_numbers <- list()
  
  for (k in seq_along(scenario_simulation_paths)) { # For each simulation
    
    single_simulation_path <- scenario_simulation_paths[k]
    
    single_simulation_replicates <- list.files(single_simulation_path, 
                                               pattern = "biomass.rds") # list the biomass files for that simulation and store
    
    replicate_paths[[k]] <- file.path(single_simulation_path, 
                                      single_simulation_replicates)
    
    replicate_numbers[[k]] <- 0:(length( single_simulation_replicates) - 1)
    
  }
  
replicate_paths <- unlist(replicate_paths)

biomass_data <- lapply(replicate_paths, readRDS)

scenario_biomass_data <- lapply(replicate_paths, readRDS)

subset_scenario_biomass_data <- list()

for (i in seq_along(scenario_biomass_data)) {
  
  scenario_biomass_data[[i]][scenario_biomass_data[[i]] == -9999] <- NA
  
  temp <- scenario_biomass_data[[i]][,burnin:ncol(scenario_biomass_data[[i]])]
  
  sample <- seq(1, ncol(temp), by = n)
  
  subset_scenario_biomass_data[[i]] <- temp[,sample]
  
  }

rm(temp)

mean_scenario_biomass <- Reduce("+", lapply(subset_scenario_biomass_data, 
                         function(x) replace(x, is.na(x), 0))) / 
                         Reduce("+", lapply(subset_scenario_biomass_data, 
                                            Negate(is.na)))

mean_scenario_biomass_long <- melt(mean_scenario_biomass)

names(mean_scenario_biomass_long) <- c("functional_group", 
                                       "timestep", "mean_biomass")


mean_biomass_plot_data <- mean_scenario_biomass_long %>%
  mutate(bodymass_bin = str_sub(functional_group, start = -2)) %>%
  mutate(group = str_sub(mean_scenario_biomass_long$functional_group, 
                         start = 1, end = 2)) %>%
  arrange(functional_group, timestep) %>%
  group_by(group, timestep) %>%
  dplyr::summarise(group_mean = mean(mean_biomass, na.rm = TRUE),
                   group_sd = sd(mean_biomass, na.rm = TRUE),
                   n = n(),
                   group_se = group_sd/sqrt(n),
                   group_lb = group_mean - (1.96 * group_se),
                   group_ub = group_mean + (1.96 * group_se)) %>%
  mutate(functional_group_name = ifelse(group == 10, "herbivorous endotherms",
                                 ifelse(group == 11, "carnivorous  endotherms",
                                 ifelse(group == 12, "omnivorous  endotherms",
                                 ifelse(group == 13, "herbivorous ectotherms", # combine iteroparous and semelparous ectotherms
                                 ifelse(group == 14, "carnivorous ectotherms",
                                 ifelse(group == 15, "omnivorous ectotherms", 
                                        "NA"))))))) %>%
  mutate(scenario = scenario)

mean_biomass_plot <- ggplot() +
  geom_path(aes(y = group_mean, x = timestep, 
                colour = functional_group_name),
            data =  mean_biomass_plot_data,
            size = 1.5) +
  theme(legend.position = "none") +
  geom_vline(xintercept = max(mean_biomass_plot_data$timestep)/3,
             color = "red") +
  geom_vline(xintercept = (max(mean_biomass_plot_data$timestep)/3)*2,
             color = "dark green") +
  ggtitle(paste(scenario, 
                "mean biomass by functional group")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + 
  theme(legend.position="right") +
  facet_wrap(~ functional_group_name)



mean_biomass_plot_95_CI <- mean_biomass_plot + 
                           geom_ribbon(aes(x = timestep,
                                            ymin = group_lb, 
                                            ymax = group_ub, 
                                            fill = functional_group_name), 
                                        alpha = 0.4, 
                                        data = mean_biomass_plot_data) 



saveRDS(mean_biomass_plot_data, file.path(output_path, 
                                          paste(scenario, 
                                                "mean_biomass_plot_data.rds", 
                                                sep = "_")))                           
ggsave(file.path(output_path, paste(scenario, 
                                    "mean biomass by functional group.pdf",
                                    sep = " ")),
                 mean_biomass_plot_95_CI) 

mean_biomass_plot_95_CI  
   
}

#' This function summarises the biomass data of all replicates, for all 
#' scenarios, using the output data produced by the plot_scenario_functional_groups
#' function. It should produce a plot of mean biomass over time for six functional
#' groups, with a separate panel for each scenario, which it saves to the 
#' scenario folder, along with the biomass time series plot data
#' 
#' @param processed_scenario_paths A vector of strings, each element in the
#' vector is the file path for one scenario's processed outputs
#' 
#' @param output_path A string - the filepath for the directory you want to
#' save the plot
#' 
#' @return should produce a plot of mean biomass over time for the six functional
#' groups which it saves to the scenario folder, along with the biomass time
#' series plot data

# processed_outputs_path <- "N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Outputs_from_adaptor_code\\map_of_life"
# 
# processed_scenario_paths <- c("N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_adaptor_code/map_of_life/Baseline",             
#                               "N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_adaptor_code/map_of_life/Harvesting_carnivores",
#                               "N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_adaptor_code/map_of_life/Harvesting_herbivores",
#                               "N:/Quantitative-Ecology/Indicators-Project/Serengeti/Outputs_from_adaptor_code/map_of_life/Land_use") 


plot_all_functional_groups <- function(processed_scenario_paths, output_path) {
  
  require(tidyverse)
  require(reshape2)
  require(ggplot2)
  
  # Get the processed plot data for each scenario
  
  all_plot_data <- list()
  
  for (i in seq_along(processed_scenario_paths)) {

  files <- list.files(processed_scenario_paths[i], recursive = FALSE) 
    
  all_plot_data[[i]] <- readRDS(paste(processed_scenario_paths[i], 
                          files[str_detect(files, "mean_biomass_plot_data.rds")], 
                             sep = "\\"))
  
  }
  
all_plot_data_df <- do.call(rbind, all_plot_data)

all_scenarios_plot <- ggplot() +
                     geom_path(aes(y = group_mean, x = timestep, 
                                  colour = functional_group_name),
                              data =  all_plot_data_df,
                              size = 1) +
                     theme(legend.position = "none") +
                     geom_vline(xintercept = max(all_plot_data_df$timestep)/3,
                               color = "red") +
                     geom_vline(xintercept = (max(all_plot_data_df$timestep)/3)*2,
                               color = "dark green") +
                     ggtitle("Mean biomass over time by scenario") +
                     theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank()) + 
                     theme(legend.position="right") +
                     facet_wrap(~ scenario)



mean_biomass_plot_95_CI <- all_scenarios_plot + 
                           geom_ribbon(aes(x = timestep,
                                            ymin = group_lb, 
                                            ymax = group_ub, 
                                            fill = functional_group_name), 
                                        alpha = 0.3, 
                                        data = all_plot_data_df) 
                          
                          
                          
saveRDS(all_plot_data_df, file.path(output_path, 
                                          paste(scenario, 
                                                "all_scenarios_plot_data.rds", 
                                                sep = "_")))                           
                          
ggsave(file.path(output_path, "Mean biomass over time by scenario.pdf"),
                                 mean_biomass_plot_95_CI) 

mean_biomass_plot_95_CI 
  
  
}

