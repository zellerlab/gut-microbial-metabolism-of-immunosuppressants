library(tidyverse)
library(here)
library(readxl)
# library(glmnet)
library(randomForest)
library(ggembl)
library(data.table)

source(here('scripts/Figure_6_utils.r'))
# Be careful here!
# this overwrites loadFromScratch, a variable defined in Figure_6_load_data_prep_env.r
# If you want to run this, please also define/overwrite scratchpath and make sure to put the corresponding gmgc file there.
loadFromScratch <- TRUE
scratchpath <- '/scratch/karcher/'
source(here('scripts/Figure_6_load_data_prep_env.r'))

#################################################################################
# Run this once to get .iv to send off with slurm (or your favourite HPC system :)
# BEFORE RUNNING, YOU HAVE TO THEN OUT-COMMENT
#################################################################################
# file_conn <- file(here("scripts/large_gene_model_lines.iv"), open = "wt")
# for (community in unique(tax_profiles_long$Stool_donor)) {
#     for (target in drugsToComputeMetricsFor) {
#         # Write the text to the file
#         writeLines(str_c("Rscript Figure_6_run_large_gene_models.r ", community, " ", str_replace_all(target, " ", "_"), " ", "FALSE", " ", 1), file_conn)
#     }
# }
# for (community in unique(tax_profiles_long$Stool_donor)) {
#     for (target in drugsToComputeMetricsFor) {
#         for (s in 1:100) {
#             writeLines(str_c("Rscript Figure_6_run_large_gene_models.r ", community, " ", str_replace_all(target, " ", "_"), " ", "TRUE", " ", s), file_conn)
#         }
#         # Write the text to the file
#     }
# }
# close(file_conn)
#################################################################################
# Run this once to get .iv to send off with slurm (or your favourite HPC system :)
# BEFORE RUNNING, YOU HAVE TO THEN OUT-COMMENT
#################################################################################

args <- commandArgs(trailingOnly = TRUE)
community = args[1]
target = str_replace_all(args[2], "_", " ")
permute = as.logical(args[3])
seed = as.numeric(args[4])

# ss <- gmgc_profile_all_genes_long %>%
#     group_by(gene) %>%
#     # a single gene has a dash in it, it looks like this can trip up randomForest
#     mutate(gene = str_replace_all(gene, "-", "_")) %>%
#     # This gene seems to also break the randomForest, for some reason
#     mutate(gene = ifelse(str_detect(gene, "GMGC10.295_880_869"), "fixed_gene_name", gene)) %>%
#     mutate(abundance = (abundance - mean(abundance)) / sd(abundance)) %>%
#     filter(!is.na(abundance)) %>%
#     ungroup()

ss <- gmgc_profile_all_genes_wide


if (!permute) {
    print("Calculating actual model performance")
    predictions_and_models_abundant_genes <- get_model_performance(
        gene_abundances = ss,
        targets = c(target),
        pivotWide = FALSE,
        stoolDonors = c(community),
        predictor_target_map = gene_target_map,
        numResamp = 1,
        mod = "abundant_genes_not_permuted",
        modelType = "RF",
        onlyReturnData = TRUE,
        how = 'all_predictors',
        preFilterCor = TRUE,
        permuteOutcome = FALSE)
    predictions_abundant_genes <- predictions_and_models_abundant_genes[[2]]
    write_tsv(predictions_abundant_genes, here("tmp/large_gene_stuff_real_models", str_c("large_gene_model_predictions__", community, "__", target, "__", permute, "__", as.character(seed), "___FULL_MODEL.tsv")))
} else {
    predictions_and_models_abundant_genes_random <- get_model_performance(
        gene_abundance = ss,
        seed = seed,
        targets = c(target),
        pivotWide = FALSE,
        stoolDonors = c(community),
        numResamp = 1,
        predictor_target_map = gene_target_map,
        mod = "abundant_genes_permuted",
        modelType = "RF",
        onlyReturnData = TRUE,
        how = 'all_predictors',
        preFilterCor = TRUE,
        permuteOutcome = TRUE)
    predictions_abundant_genes_permuted <- predictions_and_models_abundant_genes_random[[2]]
    write_tsv(predictions_abundant_genes_permuted, here("tmp/large_gene_stuff_permuted_models", str_c("large_gene_model_predictions__", community, "__", target, "__", permute, "__", seed, "___FULL_MODEL.tsv")))
}
