# If you want to analyze the abundance of KEGG pathways instead of KO within the pathway, please set `ko_to_kegg` to TRUE.
# KEGG pathways typically have more descriptive explanations.

library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)

# Load necessary data: abundance data and metadata
abundance_file <- "../data/KO_metagenome_out/pred_metagenome_unstrat.tsv"
metadata <- read_delim(
  "../metadata/MetaLight_Sample_PatAge_V2.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

metadata <- metadata %>% mutate(age_group = ifelse(age_at_initial_pathologic_diagnosis<=50, "Young", "Old"))

# Run ggpicrust2 with input file path
results_file_input <- ggpicrust2(file = abundance_file,
                                 metadata = metadata,
                                 group = "age_group", ## For example dataset, group = "Environment"
                                 pathway = "KO",
                                 daa_method = "LinDA",
                                 ko_to_kegg = TRUE,
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")

# Run ggpicrust2 with imported data.frame
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)

# Run ggpicrust2 with input data
results_data_input <- ggpicrust2(data = abundance_data,
                                 metadata = metadata,
                                 group = "ShortLetterCode", # For example dataset, group = "Environment"
                                 pathway = "KO",
                                 daa_method = "LinDA",
                                 ko_to_kegg = TRUE,
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")

# Access the plot and results dataframe for the first DA method
example_plot <- results_file_input[[1]]$plot
example_results <- results_file_input[[1]]$results

# Use the example data in ggpicrust2 package
data(ko_abundance)
data(metadata)
results_file_input <- ggpicrust2(data = ko_abundance,
                                 metadata = metadata,
                                 group = "Environment",
                                 pathway = "KO",
                                 daa_method = "LinDA",
                                 ko_to_kegg = TRUE,
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")

# Analyze the EC or MetaCyc pathway
data(metacyc_abundance)
results_file_input <- ggpicrust2(data = metacyc_abundance,
                                 metadata = metadata,
                                 group = "Environment",
                                 pathway = "MetaCyc",
                                 daa_method = "LinDA",
                                 ko_to_kegg = FALSE,
                                 order = "group",
                                 p_values_bar = TRUE,
                                 x_lab = "description")
results_file_input[[1]]$plot
results_file_input[[1]]$results