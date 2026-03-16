library(tidyverse)
library(ggplot2)
library(here)
library(randomForest)

load_from_scratch <- TRUE

#model_path <- "/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/tmp/large_gene_stuff_ACTUAL_models/large_gene_model__A-06__Mycophenolate Mofetil__FALSE__1___FULL_MODEL.rds"
# load model_path from args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No model path provided. Please provide the model path as an argument.")
}
model_path <- args[1]
model_path_end <- str_split(model_path, "/")[[1]][length(str_split(model_path, "/")[[1]])]
held_out_patient <- str_split(model_path_end, "__")[[1]][2]
drug <- str_split(model_path_end, "__")[[1]][3]

# Load profiles for 27 PRISMA samples

if (!load_from_scratch) {
  gmgc_profile_raw <- read_tsv("/g/scb/zeller/SHARED/DATA/functional_profiles_gmgc_0_5_prevalence/27_PRISMA_samples/output/collated/collated.gene_counts.combined_scaled.txt.gz")
} else {
  gmgc_profile_raw <- read_tsv("/scratch/karcher/collated.gene_counts.combined_scaled.txt.gz")
}

gmgc_profile_raw[is.na(gmgc_profile_raw)] <- 0
gmgc_profile_all_genes <- gmgc_profile_raw %>%
    as.data.frame()
	
colSums <- c()
for (colIndex in 2:ncol(gmgc_profile_all_genes)) {
    colSums <- c(colSums, sum(gmgc_profile_all_genes[[colIndex]]))
}

gmgc_profile_all_genes <- gmgc_profile_all_genes %>%
    column_to_rownames('gene')
for (colIndex in 1:ncol(gmgc_profile_all_genes)) {
    gmgc_profile_all_genes[[colIndex]] <- gmgc_profile_all_genes[[colIndex]] / colSums[colIndex]
}

# load metadata
meta <- read_csv(here('data', 'MPA_Tox_df_for_Nic.csv'))
# load MMF models
model_object <- readRDS(model_path)[[1]]
model_features <- names(unlist(model_object$forest$xlevels))

# Keep only motus found in model
gmgc_profile_all_genes <- gmgc_profile_all_genes
gmgc_profile_all_genes <- gmgc_profile_all_genes[rownames(gmgc_profile_all_genes) %in% model_features, ]

# add motus missing in profiles and found in model
missing_features <- setdiff(model_features, rownames(gmgc_profile_all_genes))
if (length(missing_features) > 0) {
	for (feature in missing_features) {
		gmgc_profile_all_genes[feature, ] <- 0
	}
}
stopifnot(all(rownames(gmgc_profile_all_genes) %in% model_features)) # should be TRUE
stopifnot(all(model_features %in% rownames(gmgc_profile_all_genes))) # should be TRUE

pred <- predict(model_object, newdata = t(gmgc_profile_all_genes)) %>%
	enframe() %>%
	rename(
		sampleID = name, prediction = value
	) %>%
	mutate(
		raw_model_path = model_path_end,
		held_out_patient = held_out_patient,
		target = drug)

write_tsv(
	pred, 
	str_c("/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/tmp/predictions_for_27_samples/", 
		str_replace(model_path_end, "FULL_MODEL.rds", "27_samples_predictions"), 
		".tsv")
)
