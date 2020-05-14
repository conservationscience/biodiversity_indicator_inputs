
# Git repo "C:\\Users\\ssteven\\Dropbox\\Deakin\\Serengeti-analysis\\model_outputs_to_indicator_inputs\\1_process_outputs"


library(reshape2)
library(ggplot2)
library(stringr)
library(dplyr)

#' TODO: Add in autotrophs
#' TODO: Add in bodymass labels to lines
#' TODO: Add save function
#' TODO: Generalise object names
#' TODO: Add a 'mean' line for each group (and test out just including mean with
#' geom_jitter behind?)
#' TODO: Add bodymass bins and colour by bin (eg so bodymass bin 75 is the same colour regardless of functional group)

burnin <- 1*12

simulation_path <- "N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Outputs_from_adaptor_code\\map_of_life\\Test_runs\\ae_BuildModel"

# burnin <- 1000 * 12
# 
# simulation_path <- "N:\\Quantitative-Ecology\\Indicators-Project\\Serengeti\\Outputs_from_adaptor_code\\map_of_life\\Harvesting_carnivores\\203_BuildModel"

files <- list.files(simulation_path, recursive = FALSE) 

abundance_files <- paste(simulation_path, files[str_detect(files, "biomass.rds")], sep = "\\")

abundance_files <- abundance_files[1:2] # TEMP CODE - subset so you don't read everything in

name_pt_1 <- basename(simulation_path)
names_pt_2 <- unlist(lapply(abundance_files, basename))
plot_names <- paste(name_pt_1, names_pt_2, sep = "_")

simulation_abundance_data <- lapply(abundance_files, readRDS)

# Melt into long format

simulation_abundance_data_long <- lapply(simulation_abundance_data, melt)

# Replace no value with NA

for (i in seq_along(simulation_abundance_data_long)) {
  
  simulation_abundance_data_long[[i]][simulation_abundance_data_long[[i]] == -9999] <- NA
  
}

plot_data <- list()

for (i in seq_along(simulation_abundance_data_long)) {
  
names(simulation_abundance_data_long[[i]]) <- c("functional_group", "timestep", "abundance")

simulation_abundance_data_long[[i]]$functional_group <- as.character(simulation_abundance_data_long[[i]]$functional_group)

simulation_abundance_data_long[[i]] <- simulation_abundance_data_long[[i]][!is.na(simulation_abundance_data_long[[i]]$abundance),]

plot_data[[i]] <- simulation_abundance_data_long[[i]] %>%
                  mutate(group = 
                         str_sub(simulation_abundance_data_long[[i]]$functional_group, 
                                 start = 1, end= 2)) %>%
                  filter(timestep >= burnin) %>%
                  mutate(functional_group_name = ifelse(group == 10, "herbivorous endotherms",
                         ifelse(group == 11, "carnivorous  endotherms",
                         ifelse(group == 12, "omnivorous  endotherms",
                         ifelse(group == 13, "herbivorous ectotherms", # combine iteroparous and semelparous ectotherms
                         ifelse(group == 14, "carnivorous ectotherms",
                         ifelse(group == 15, "omnivorous ectotherms", "NA"))))))) %>%
                  arrange(functional_group, timestep) %>%
                  filter(row_number() %% 360 == 0) %>%
                  mutate(years = seq(1:nrow(.)))
}

dim(plot_data[[1]])

head(plot_data[[1]])

# Plot

abundance_plots <- list()

for (i  in seq_along(plot_data)) {

abundance_plots[[i]] <- ggplot() +
             geom_path(aes(y = abundance, x = timestep, colour = functional_group),
                  data =  plot_data[[1]]) +
             theme(legend.position = "none") +
             geom_vline(xintercept = 1100 *12, color = "red") +
             geom_vline(xintercept = 1200 * 12, color = "dark green") +
             facet_wrap(~ functional_group_name, ncol = 3) +
             ggtitle(plot_names[[i]])

abundance_plots[[i]]

}

abundance_plots[[2]]

ggsave("C:\\Users\\ssteven\\Dropbox\\Deakin\\Serengeti-analysis\\herb303test1.pdf", abundance_plots[[1]],  device = "pdf")
ggsave("C:\\Users\\ssteven\\Dropbox\\Deakin\\Serengeti-analysis\\herb303test2.4.pdf", abundance_plots[[2]],  device = "pdf")
