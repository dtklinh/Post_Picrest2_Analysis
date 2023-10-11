library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)

# If you want to analyze KEGG pathway abundance instead of KO within the pathway, turn ko_to_kegg to TRUE.
# KEGG pathways typically have more explainable descriptions.

# Load metadata as a tibble
# data(metadata)
metadata <- read_delim("../metadata/TCMA/MetaLight_Sample_PatAge_V3.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
metadata <- metadata %>% mutate(age_group = ifelse(age_at_initial_pathologic_diagnosis<=50, "Young", "Old"))

# Load KEGG pathway abundance
# data(kegg_abundance)
kegg_abundance <- ko2kegg_abundance("../data/KO_metagenome_out/pred_metagenome_unstrat.tsv") 

## filter metadata
metadata <- metadata[metadata$SampleID %in% colnames(kegg_abundance),]

# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
# Please change group to "your_group_column" if you are not using example dataset
daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = "Sample", daa_method = "ALDEx2", select = NULL, p.adjust = "none", reference = NULL) 

# Filter results for ALDEx2_Welch's t test method
# Please check the unique(daa_results_df$method) and choose one
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Wilcoxon rank test", ]
daa_sub_method_results_df <- daa_sub_method_results_df[daa_sub_method_results_df$p_adjust<=0.05,]

# Annotate pathway results using KO to KEGG conversion

daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = TRUE)
daa_annotated_sub_method_results_df_2 <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df[81:nrow(daa_sub_method_results_df),], ko_to_kegg = TRUE)
daa_annotated_sub_method_results_df <- rbind(daa_annotated_sub_method_results_df_1, daa_annotated_sub_method_results_df_2)

# Generate pathway error bar plot
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = kegg_abundance, daa_results_df = daa_annotated_sub_method_results_df, Group = metadata$Sample, p_values_threshold = 0.05, order = "pathway_class", select = NULL, ko_to_kegg = TRUE, p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")

# If you want to analyze EC, MetaCyc, and KO without conversions, turn ko_to_kegg to FALSE.
######################################################################################
# Load metadata as a tibble
# data(metadata)
metadata <- read_delim("path/to/your/metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Load KO abundance as a data.frame
# data(ko_abundance)
ko_abundance <- read.delim("../data/KO_metagenome_out/pred_metagenome_unstrat.tsv", check.names = F)
ko_abundance <- ko_abundance %>% column_to_rownames("function")

metadata <- metadata[metadata$SampleID %in% colnames(ko_abundance),]

# Perform pathway DAA using ALDEx2 method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
daa_results_df <- pathway_daa(abundance = ko_abundance, metadata = metadata, group = "age_group", daa_method = "ALDEx2", select = NULL, reference = NULL, p.adjust = "none")

# Filter results for ALDEx2_Kruskal-Wallace test method
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Wilcoxon rank test", ]

# Annotate pathway results without KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = FALSE)

# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = ko_abundance, daa_results_df = daa_annotated_sub_method_results_df, Group = metadata$age_group, p_values_threshold = 0.05, order = "group",
                      select = daa_annotated_sub_method_results_df %>% arrange(p_adjust) %>% slice(1:20) %>% dplyr::select(feature) %>% pull(), 
                      ko_to_kegg = FALSE, 
                      p_value_bar = TRUE, 
                      colors = NULL, 
                      x_lab = "description")

plot(p)
# Workflow for MetaCyc Pathway and EC

# Load MetaCyc pathway abundance and metadata
data("metacyc_abundance")
data("metadata")

# Perform pathway DAA using LinDA method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment", daa_method = "LinDA")

# Annotate MetaCyc pathway results without KO to KEGG conversion
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_daa_results_df, ko_to_kegg = FALSE)

# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
pathway_errorbar(abundance = metacyc_abundance %>% column_to_rownames("pathway"), daa_results_df = metacyc_daa_annotated_results_df, Group = metadata$Environment, ko_to_kegg = FALSE, p_values_threshold = 0.05, order = "group", select = NULL, p_value_bar = TRUE, colors = NULL, x_lab = "description")

# Generate pathway heatmap
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
feature_with_p_0.05 <- metacyc_daa_results_df %>% filter(p_adjust < 0.05)
pathway_heatmap(abundance = metacyc_abundance %>% filter(pathway %in% feature_with_p_0.05$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment")

# Generate pathway PCA plot
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
pathway_pca(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment")

# Run pathway DAA for multiple methods
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
methods <- c("ALDEx2", "DESeq2", "edgeR")
daa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment", daa_method = method)
})

# Compare results across different methods
comparison_results <- compare_daa_results(daa_results_list = daa_results_list, method_names = c("ALDEx2_Welch's t test", "ALDEx2_Wilcoxon rank test", "DESeq2", "edgeR"))