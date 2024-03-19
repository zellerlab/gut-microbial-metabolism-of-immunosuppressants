library(tidyverse)
library(here)
library(readxl)
# library(glmnet)
library(randomForest)
library(ggembl)
library(data.table)

## Disclaimer: This script is taking some time to run, as it is running a lot of models
## The most time consuming part - treaining and evaluating of the large gene models - is in a separate script. you'll see that part by inline comments
## The rest wil lrun for a day or so.

source(here('scripts/Figure_6_utils.r'))
# This runs for 10 minutes. The GMGC gene catalogue is large - sorry.
source(here('scripts/Figure_6_load_data_prep_env.r'))

predictions_and_models_gene_characterized <- get_model_performance(gene_abundances = gmgc_profile_gmgc_profile_homologues_to_characterized_enzymes_long %>%
    group_by(gene) %>%
    mutate(abundance = (abundance - mean(abundance)) / sd(abundance)) %>%
    filter(!is.na(abundance)) %>%
    ungroup(), numResamp = 1, targets = drugsToComputeMetricsFor, predictor_target_map = gene_target_map, mod = "characterized_enzymes", modelType = "RF", alsoReturnData = TRUE, how = 'only_testedactive')
# models_gene_characterized <- predictions_and_models_gene_characterized[[1]
predictions_gene_characterized <- predictions_and_models_gene_characterized[[2]]
# tdm_gene_characterized <- do.call('rbind', predictions_and_models_gene_characterized[[3]])
# directionality_gene_characterized <- do.call('rbind', predictions_and_models_gene_characterized[[4]])


# predictions_and_models_abundant_genes <- get_model_performance(gene_abundances =
#     gmgc_profile_all_genes_wide,
# numResamp = 1,
# pivotWide = FALSE,
# targets = drugsToComputeMetricsFor,
# predictor_target_map = gene_target_map,
# mod = "characterized_enzymes",
# modelType = "RF",
# alsoReturnData = TRUE,
# how = 'all_predictors',
# preFilterCor = TRUE
# )
# predictions_abundant_genes <- predictions_and_models_abundant_genes[[2]]
# To save time, I load here precomputed values (see Figure_6_get_large_gene_model_predictions.r, which also produces the baseline predictions, for this see below)
predictions_abundant_genes <- map(list.files(here('tmp/large_gene_stuff_real_models'), pattern = "__FALSE__[0-9]+.tsv", full.names = TRUE), \(x) fread(x, sep = '\t')) %>%
    do.call('rbind', .) %>%
    as_tibble()

predictions_and_models_strains_active <- get_model_performance(gene_abundances = tax_profiles_long %>%
    rename(gene = mOTUs_ID, stoolDonor = Stool_donor) %>%
    group_by(gene) %>%
    mutate(abundance = (abundance - mean(abundance)) / sd(abundance)) %>%
    filter(!is.na(abundance)) %>%
    ungroup(), numResamp = 1, targets = drugsToComputeMetricsFor, predictor_target_map = species_target_map, mod = "active_strains", modelType = "RF", alsoReturnData = TRUE, how = 'only_testedactive')
# models_strains_active <- predictions_and_models_strains_active[[1]]
predictions_strains_active <- predictions_and_models_strains_active[[2]]
# tdm_strains_active <- do.call('rbind', predictions_and_models_strains_active[[3]])
# directionality_strains_active <- do.call('rbind', predictions_and_models_strains_active[[4]])

predictions_and_models_motus_all <- get_model_performance(gene_abundances = tax_profiles_long %>%
    rename(gene = mOTUs_ID, stoolDonor = Stool_donor) %>%
    group_by(gene) %>%
    filter(mean(abundance > 0) > 0.2) %>%
    filter(max(abundance) > 1E-3) %>%
    mutate(abundance = (abundance - mean(abundance)) / sd(abundance)) %>%
    filter(!is.na(abundance)) %>%
    ungroup(),
numResamp = 1,
###############
###############
# In this function I'm looping over targets to get model performance
# Unfortunately I'm setting the seed only once, before the loop (my bad).
# Since results (performances/tpo features) are hence somewhat dependent on the order of the targets vector,
# I'm doing this here to ensure results are perfectly consistent with the previous runs.
# The model predictions_and_models_motus_all_supplement contains everything else
targets = DrugsToShowInMainFigure,
# targets = DrugsToShowInMainFigure,
predictor_target_map = species_target_map,
mod = "all_motus",
modelType = "RF",
alsoReturnData = TRUE,
how = 'all_predictors',
numFeaturesForTopDownSpearmanCorComparisons = 100) ############################
# Contains all large models with all features
models_motus_all <- predictions_and_models_motus_all[[1]]
# Contains predictions
predictions_motus_all <- predictions_and_models_motus_all[[2]]
# Contains predictions of top-down models
tdm_motus_all <- do.call('rbind', predictions_and_models_motus_all[[3]])
# Contains directionality of all features based on correlation with outcome
directionality_motus_all <- do.call(
    'rbind',
    map(predictions_and_models_motus_all[[4]], \(x) {
        x %>%
            as.data.frame() %>%
            rownames_to_column('gene')
    }))
# truncate_last_field <- function(string) {
#     sub("(_\\d{5})\\d*$", "\\1", string)
# }
directionality_motus_all <- directionality_motus_all %>%
    # mutate(gene = truncate_last_field(gene)) %>%
    group_by(gene, target, dataType, seed) %>%
    summarize(directionality = mean(directionality)) %>%
    ungroup() %>%
    mutate(directionality = case_when(
        # directionality > 0.9 ~ 'positive',
        # directionality < 0.1 ~ 'negative',
        directionality > 0.5 ~ 'positive',
        directionality <= 0.5 ~ 'negative',
        .default = NA
    ))
# Contains top down models themselfes
tdm_models_motus_all <- predictions_and_models_motus_all[[5]]
tdm_models_motus_all <- do.call('rbind', map2(tdm_models_motus_all, names(tdm_models_motus_all), \(x, nnn) {
    do.call('rbind', map2(x, names(x), \(x, nn) {
        importance(x) %>%
            as.data.frame() %>%
            rownames_to_column('feature') %>%
            as_tibble() %>%
            mutate(numTopFeatures = nn)
    })) %>%
        as.data.frame() %>%
        as_tibble() %>%
        mutate(heldOutCommunity = str_split_fixed(nnn, "__", n = 4)[, 1]) %>%
        mutate(target = str_split_fixed(nnn, "__", n = 4)[, 2]) %>%
        mutate(dataType = str_split_fixed(nnn, "__", n = 4)[, 3]) %>%
        mutate(seed = str_split_fixed(nnn, "__", n = 4)[, 4]) %>%
        mutate(seed = str_split_fixed(seed, "__", n = 4)[, 1])
}))
predictions_and_models_motus_all_supplement <- get_model_performance(gene_abundances = tax_profiles_long %>%
    rename(gene = mOTUs_ID, stoolDonor = Stool_donor) %>%
    group_by(gene) %>%
    filter(mean(abundance > 0) > 0.2) %>%
    filter(max(abundance) > 1E-3) %>%
    mutate(abundance = (abundance - mean(abundance)) / sd(abundance)) %>%
    filter(!is.na(abundance)) %>%
    ungroup(),
numResamp = 1,
targets = drugsToComputeMetricsFor[!drugsToComputeMetricsFor %in% DrugsToShowInMainFigure],
# targets = DrugsToShowInMainFigure,
predictor_target_map = species_target_map,
mod = "all_motus",
modelType = "RF",
alsoReturnData = TRUE,
how = 'all_predictors',
numFeaturesForTopDownSpearmanCorComparisons = 100) ############################
models_motus_all_supplement <- predictions_and_models_motus_all_supplement[[1]]
# Contains predictions
predictions_motus_all_supplement <- predictions_and_models_motus_all_supplement[[2]]
# Contains predictions of top-down models
tdm_motus_all_supplement <- do.call('rbind', predictions_and_models_motus_all_supplement[[3]])
# Contains directionality of all features based on correlation with outcome
directionality_motus_all_supplement <- do.call(
    'rbind',
    map(predictions_and_models_motus_all_supplement[[4]], \(x) {
        x %>%
            as.data.frame() %>%
            rownames_to_column('gene')
    }))
# truncate_last_field <- function(string) {
#     sub("(_\\d{5})\\d*$", "\\1", string)
# }
directionality_motus_all_supplement <- directionality_motus_all_supplement %>%
    # mutate(gene = truncate_last_field(gene)) %>%
    group_by(gene, target, dataType, seed) %>%
    summarize(directionality = mean(directionality)) %>%
    ungroup() %>%
    mutate(directionality = case_when(
        # directionality > 0.9 ~ 'positive',
        # directionality < 0.1 ~ 'negative',
        directionality > 0.5 ~ 'positive',
        directionality <= 0.5 ~ 'negative',
        .default = NA
    ))
# Contains top down models themselfes
tdm_models_motus_all_supplement <- predictions_and_models_motus_all_supplement[[5]]
tdm_models_motus_all_supplement <- do.call('rbind', map2(tdm_models_motus_all_supplement, names(tdm_models_motus_all_supplement), \(x, nnn) {
    do.call('rbind', map2(x, names(x), \(x, nn) {
        importance(x) %>%
            as.data.frame() %>%
            rownames_to_column('feature') %>%
            as_tibble() %>%
            mutate(numTopFeatures = nn)
    })) %>%
        as.data.frame() %>%
        as_tibble() %>%
        mutate(heldOutCommunity = str_split_fixed(nnn, "__", n = 4)[, 1]) %>%
        mutate(target = str_split_fixed(nnn, "__", n = 4)[, 2]) %>%
        mutate(dataType = str_split_fixed(nnn, "__", n = 4)[, 3]) %>%
        mutate(seed = str_split_fixed(nnn, "__", n = 4)[, 4]) %>%
        mutate(seed = str_split_fixed(seed, "__", n = 4)[, 1])
}))


all_predictions_and_models_strains_active_random <- list()
set.seed(1)
for (iteration in 1:100) {
    print(iteration)
    ss <- tax_profiles_long %>%
        rename(gene = mOTUs_ID, stoolDonor = Stool_donor) %>%
        group_by(gene) %>%
        mutate(abundance = (abundance - mean(abundance)) / sd(abundance)) %>%
        filter(!is.na(abundance)) %>%
        ungroup()
    predictions_and_models_strains_active_random <- get_model_performance(
        gene_abundances = ss,
        seed = iteration,
        targets = drugsToComputeMetricsFor,
        numResamp = 1,
        predictor_target_map = species_target_map,
        mod = "active_strains",
        modelType = "RF",
        alsoReturnData = TRUE,
        how = 'only_testedactive',
        permuteOutcome = TRUE)
    all_predictions_and_models_strains_active_random[[length(all_predictions_and_models_strains_active_random) + 1]] <- predictions_and_models_strains_active_random[[2]]
}
# For reloading later, if necessary, save here
# save(all_predictions_and_models_strains_active_random, file = "/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/data_tmp/all_predictions_and_models_strains_active_random.rdata")
all_predictions_and_models_motus_all_random <- list()
set.seed(1)
for (iteration in 1:100) {
    print(iteration)
    ss <- tax_profiles_long %>%
        rename(gene = mOTUs_ID, stoolDonor = Stool_donor) %>%
        group_by(gene) %>%
        filter(mean(abundance > 0) > 0.2) %>%
        filter(max(abundance) > 1E-3) %>%
        mutate(abundance = (abundance - mean(abundance)) / sd(abundance)) %>%
        filter(!is.na(abundance)) %>%
        ungroup()
    predictions_and_models_motus_all_random <- get_model_performance(
        gene_abundances = ss,
        seed = iteration,
        targets = drugsToComputeMetricsFor,
        # targets = c("Tacrolimus"),
        numResamp = 1,
        predictor_target_map = species_target_map,
        mod = "random_strains",
        modelType = "RF",
        alsoReturnData = TRUE,
        how = 'all_predictors',
        permuteOutcome = TRUE)
    all_predictions_and_models_motus_all_random[[length(all_predictions_and_models_motus_all_random) + 1]] <- predictions_and_models_motus_all_random[[2]]
}
# For reloading later, if necessary, save here
# save(all_predictions_and_models_motus_all_random, file = "/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/data_tmp/all_predictions_and_models_motus_all_random.rdata")

all_predictions_and_models_gene_characterized_random <- list()
set.seed(1)
for (iteration in 1:100) {
    print(iteration)
    ss <- gmgc_profile_gmgc_profile_homologues_to_characterized_enzymes_long %>%
        group_by(gene) %>%
        mutate(abundance = (abundance - mean(abundance)) / sd(abundance)) %>%
        filter(!is.na(abundance)) %>%
        ungroup()
    predictions_and_models_gene_characterized_random <- get_model_performance(gene_abundances = ss,
        numResamp = 1,
        seed = iteration,
        targets = drugsToComputeMetricsFor,
        # targets = DrugsToShowInMainFigure,
        predictor_target_map =
            gene_target_map,
        mod = "characterized_enzymes",
        modelType = "RF",
        alsoReturnData = TRUE,
        how = 'only_testedactive',
        permuteOutcome = TRUE)
    all_predictions_and_models_gene_characterized_random[[length(all_predictions_and_models_gene_characterized_random) + 1]] <- predictions_and_models_gene_characterized_random[[2]]
}
# For reloading later, if necessary, save here
# save(all_predictions_and_models_gene_characterized_random, file = "/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/tmp/all_predictions_and_models_gene_characterized_random.rdata")

##
## This is extremely slow, so I'm putting this script into a separate file if run sequentually
## See Figure_6_run_large_gene_models.r
#
# all_predictions_and_models_abundant_genes_random <- list()
# set.seed(1)
# for (iteration in 1:100) {
#     print(iteration)
#     ss <- gmgc_profile_all_genes_long %>%
#         group_by(gene) %>%
#         # a single gene has a dash in it, replace
#         mutate(gene = str_replace_all(gene, "-", "_")) %>%
#         # This gene seems to also break the randomForest, for some reason
#         mutate(gene = ifelse(str_detect(gene, "GMGC10.295_880_869"), "blaaaaaa", gene)) %>%
#         mutate(abundance = (abundance - mean(abundance)) / sd(abundance)) %>%
#         filter(!is.na(abundance)) %>%
#         ungroup()
#      ss <- gmgc_profile_all_genes_wide
#     predictions_and_models_abundant_genes_random <- get_model_performance(
#         gene_abundance = ss,
#         seed = iteration,
# pivotWide = FALSE,
#         targets = drugsToComputeMetricsFor,
#         numResamp = 1,
#         predictor_target_map = gene_target_map,
#         mod = "characterized_enzymes",
#         modelType = "RF",
#         alsoReturnData = TRUE,
#         how = 'all_predictors',
#         preFilterCor = TRUE,
#         permuteOutcome = TRUE)
#     all_predictions_and_models_abundant_genes_random[[length(all_predictions_and_models_abundant_genes_random) + 1]] <- predictions_and_models_abundant_genes_random[[2]]
# }
paths <- list.files(here('tmp/large_gene_stuff_permuted_models'), pattern = "__TRUE__[0-9]+.tsv", full.names = TRUE)
all_predictions_and_models_abundant_genes_random <- map(paths, \(x) {
    y <- fread(x, sep = '\t')
    seed <- str_split(x, "__")[[1]]
    seed <- seed[length(seed)]
    seed <- str_replace(seed, ".tsv", "")
    see <- seed
    y <- y %>%
        as.data.frame() %>%
        as_tibble()
    y <- y %>%
        mutate(seed = see)
    return(y)
})
all_predictions_and_models_abundant_genes_random <- do.call('rbind', all_predictions_and_models_abundant_genes_random) %>%
    as_tibble()


# motus
predictions_motus <- rbind(
    predictions_strains_active %>% mutate(features = "Identified active strains"),
    predictions_motus_all %>% mutate(features = "All species")
    predictions_motus_all_supplement %>% mutate(features = "All species")
) %>%
    mutate(condition = case_when(
        oxygen == "MA" ~ "Microaerobic",
        oxygen == "AA" ~ "Anaerobic"
    )) %>%
    select(-oxygen)

metrics_motus <- predictions_motus %>%
    group_by(target, features) %>%
    summarize(spearmanCore = cor(prediction, truth, method = "spearman"),
        pearsonCore = cor(prediction, truth),
        R2 = pearsonCore^2) %>%
    filter(!str_detect(target, 'IS'))

performance_by_feature_motus_all <- rbind(tdm_motus_all %>%
    group_by(numFeatures, target, dataType, seed) %>%
    summarize(spearmanCor = cor(prediction, truth, method = 'spearman')) %>%
    unnest() %>%
    arrange(target, numFeatures) %>%
    ungroup() %>%
    filter(!str_detect(target, "IS"))
    tdm_motus_all_supplement %>%
    group_by(numFeatures, target, dataType, seed) %>%
    summarize(spearmanCor = cor(prediction, truth, method = 'spearman')) %>%
    unnest() %>%
    arrange(target, numFeatures) %>%
    ungroup() %>%
    filter(!str_detect(target, "IS")))


numberFeaturesDataMotusAll <- performance_by_feature_motus_all %>%
    left_join(metrics_motus %>%
        filter(features == "All species") %>%
        select(target, features, spearmanCore),
    by = c('target')) %>%
    arrange(target, numFeatures) %>%
    mutate(close = spearmanCor > 0.9 * spearmanCore) %>%
    group_by(target) %>%
    nest() %>%
    mutate(data = map(data, \(x) {
        if (any(x$close)) {
            return(x %>% filter(close) %>% slice_head(n = 1))
        } else {
            return(x %>% ungroup() %>% slice_tail(n = 1))
        }
    })) %>%
    unnest() %>%
    rename(modelCloseToFullModelInPerformance = close)


plotObject <- ggplot(
    data = performance_by_feature_motus_all %>%
        left_join(numberFeaturesDataMotusAll %>%
            select(target, numFeatures) %>%
            rename(numFeaturesCloseEnough = numFeatures), by = c('target')) %>%
        filter(numFeatures <= numFeaturesCloseEnough + 3),
    aes(x = numFeatures, y = spearmanCor, group = 1)) +
    geom_rect(data =
        performance_by_feature_motus_all %>%
            ungroup() %>%
            select(target) %>%
            distinct() %>%
            mutate(maintext = ifelse(target %in% DrugsToShowInMainFigure, TRUE, FALSE)), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = maintext), inherit.aes = FALSE, show.legend = FALSE) +
    facet_wrap(. ~ target, nrow = 5, scales = "free_x") +
    # facet_wrap(. ~ target, nrow = 5) +
    geom_line() +
    geom_abline(data = metrics_motus %>% filter(features == "all species"), aes(intercept = spearmanCore, slope = 0, color = 'red'), show.legend = F) +
    theme_publication() +
    # scale_x_continuous(limits = c(0, 20), breaks = c(0, 5, 10, 15, 20)) +
    scale_fill_manual(values = c("TRUE" = "lightblue", "FALSE" = "white")) +
    ylab("Spearman correlation between measured\nand predicted metabolic activity") +
    xlab("Number of top features in model") +
    geom_rect(data = numberFeaturesDataMotusAll %>% filter(modelCloseToFullModelInPerformance), aes(xmin = numFeatures, xmax = +Inf, ymin = -Inf, ymax = Inf), color = 'darkgreen', alpha = 0.3) +
    NULL

ggsave(plot = plotObject, filename = here("results/plots/motus_all_performance_by_feature.pdf"), width = 7, height = 6.5)


predictions_motus$target <- as.factor(predictions_motus$target)
metrics_motus$target <- factor(metrics_motus$target, levels = levels(predictions_motus$target))
targetOrder <- levels(predictions_motus$target)

library(ggembl)
taxaScatters <- ggplot() +
    geom_point(data =
        predictions_motus %>%
            filter(features == 'All species') %>%
            filter(!str_detect(target, "IS")) %>%
            filter(target %in% DrugsToShowInMainFigure) %>%
            mutate(target = factor(target, levels = DrugsToShowInMainFigure)), aes(x = truth, y = prediction, color = condition), alpha = 1) +
    # geom_line(data = predictions_motus %>%
    #     filter(features == 'all species') %>%
    #     filter(!str_detect(target, "IS")) %>%
    #     filter(target %in% drugsToComputeMetricsFor), aes(x = truth, y = prediction, group = left_out_community), alpha = 0.3) +
    # theme_classic() +
    theme_presentation() +
    geom_text(data =
        metrics_motus %>%
            filter(features == 'All species') %>%
            filter(!str_detect(target, "IS")) %>%
            filter(target %in% DrugsToShowInMainFigure) %>%
            mutate(target = factor(target, levels = DrugsToShowInMainFigure)), aes(x = -Inf, y = Inf, label = str_c("Rho: ", round(spearmanCore, 3))), hjust = -0.05, vjust = 1.2, color = 'black', size = 8 * 0.3523) +
    scale_color_manual(values = aa_ma_colors) +
    facet_wrap(target ~ features, ncol = 1, scales = 'free', strip.position = "right") +
    scale_x_continuous(breaks = seq(floor(min(predictions_motus$truth)), ceiling(max(predictions_motus$truth)), by = 2)) +
    scale_y_continuous(breaks = seq(floor(min(predictions_motus$prediction)), ceiling(max(predictions_motus$prediction)), by = 1), position = "right") +
    # ggtitle("All species") +
    xlab("True drug metabolism [AUC]") +
    ylab("Predicted drug metabolism [AUC]") +
    theme(
        axis.ticks.length = unit(0.5, "mm"),
        strip.text = element_blank()
    ) +
    # xlim(c(0, 12)) +
    # scale_x_continuous(limits = c(0, 12), breaks = c(0, 3, 6, 9, 12)) +
    # scale_y_continuous(limits = c(0, 12), breaks = c(0, 3, 6, 9, 12)) +
    # geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "lightgrey") +
    NULL

ggsave(
    plot = taxaScatters,
    filename = here("results/plots/scatters_taxa.pdf"), width = 3.4, height = 8)

library(ggembl)
taxaScattersSupplement <- ggplot() +
    geom_point(data =
        predictions_motus %>%
            filter(features == 'All species') %>%
            filter(!str_detect(target, "IS")) %>%
            filter(target %in% drugsToComputeMetricsFor) %>%
            mutate(target = factor(target, levels = drugsToComputeMetricsFor)), aes(x = truth, y = prediction, color = condition), alpha = 1) +
    # geom_line(data = predictions_motus %>%
    #     filter(features == 'all species') %>%
    #     filter(!str_detect(target, "IS")) %>%
    #     filter(target %in% drugsToComputeMetricsFor), aes(x = truth, y = prediction, group = left_out_community), alpha = 0.3) +
    # theme_classic() +
    theme_presentation() +
    geom_text(data =
        metrics_motus %>%
            filter(features == 'All species') %>%
            filter(!str_detect(target, "IS")) %>%
            filter(target %in% drugsToComputeMetricsFor) %>%
            mutate(target = factor(target, levels = drugsToComputeMetricsFor)), aes(x = -Inf, y = Inf, label = str_c("Rho: ", round(spearmanCore, 3))), hjust = -0.05, vjust = 1.2, color = 'black', size = 8 * 0.3523) +
    scale_color_manual(values = aa_ma_colors) +
    facet_wrap(target ~ ., ncol = 3, scales = 'free', strip.position = "top") +
    scale_x_continuous(breaks = seq(floor(min(predictions_motus$truth)), ceiling(max(predictions_motus$truth)), by = 2)) +
    scale_y_continuous(breaks = seq(floor(min(predictions_motus$prediction)), ceiling(max(predictions_motus$prediction)), by = 1), position = "right") +
    # ggtitle("All species") +
    xlab("True drug metabolism [AUC]") +
    ylab("Predicted drug metabolism [AUC]") +
    theme(
        axis.ticks.length = unit(0.5, "mm"),
        # strip.text = element_blank()
    ) +
    # xlim(c(0, 12)) +
    # scale_x_continuous(limits = c(0, 12), breaks = c(0, 3, 6, 9, 12)) +
    # scale_y_continuous(limits = c(0, 12), breaks = c(0, 3, 6, 9, 12)) +
    # geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "lightgrey") +
    NULL

ggsave(
    plot = taxaScattersSupplement,
    filename = here("results/plots/scatters_taxa_supplement.pdf"), width = 5.35 * 1.075, height = 10 * 1.3)

plotObject <- evaluate_rf_feature_importances(
    # listOfModels = mOTUs_predictions_prev_mOTUs_of_interest[[1]],
    listOfModels = models_motus_all,
    modelsTopDown = tdm_models_motus_all,
    numberFeaturesData = numberFeaturesDataMotusAll,
    directionality = directionality_motus_all,
    returnData = FALSE,
    outputPath = '/g/scb/zeller/karcher/tmp/mipf.pdf',
    # width = 8,
    # height = 13,
    # testedStrainsAsReference = TRUE,
    importanceTypeToPlot = "PercIncMSE",
    yAxisLabel = "% increase MSE\n upon feature removal (OOB)",
    # plotHitsAsCategorical = TRUE,
    # plotHitsAsCategorical = FALSE,
    # prevalenceCutoff = prevCutoff,
    # yAxisLegend = "Strain activity (ranked)")
    yAxisLegend = "Strain activity",
    predictor_target_map = species_target_map,
    numFeaturesToShow = 100,
    # targetOrder = targetOrder[!str_detect(targetOrder, "IS")],
    targetOrder = c("Tacrolimus"),
    taxonBasedAnalysis = TRUE)

ggsave(plot = plotObject, filename = here('results/plots/feature_importances_taxa.pdf'), width = 4.25, height = 3)

plotObject <- evaluate_rf_feature_importances(
    # listOfModels = mOTUs_predictions_prev_mOTUs_of_interest[[1]],
    listOfModels = models_motus_all,
    modelsTopDown = tdm_models_motus_all,
    numberFeaturesData = numberFeaturesDataMotusAll,
    directionality = directionality_motus_all,
    returnData = FALSE,
    outputPath = '/g/scb/zeller/karcher/tmp/mipf.pdf',
    # width = 8,
    # height = 13,
    # testedStrainsAsReference = TRUE,
    importanceTypeToPlot = "PercIncMSE",
    yAxisLabel = "% increase MSE\n upon feature removal (OOB)",
    # plotHitsAsCategorical = TRUE,
    # plotHitsAsCategorical = FALSE,
    # prevalenceCutoff = prevCutoff,
    # yAxisLegend = "Strain activity (ranked)")
    yAxisLegend = "Strain activity",
    predictor_target_map = species_target_map,
    numFeaturesToShow = 100,
    # targetOrder = targetOrder[!str_detect(targetOrder, "IS")],
    targetOrder = DrugsToHighlightInMainFigure,
    taxonBasedAnalysis = TRUE)

ggsave(plot = plotObject, filename = here('results/plots/feature_importances_taxa_Tac_Eve_Sir_no_union.pdf'), width = 4.25, height = 8)

plotObject <- evaluate_rf_feature_importances(
    # listOfModels = mOTUs_predictions_prev_mOTUs_of_interest[[1]],
    listOfModels = models_motus_all,
    modelsTopDown = tdm_models_motus_all,
    numberFeaturesData = numberFeaturesDataMotusAll,
    directionality = directionality_motus_all,
    returnData = FALSE,
    outputPath = '/g/scb/zeller/karcher/tmp/mipf.pdf',
    # width = 8,
    # height = 13,
    # testedStrainsAsReference = TRUE,
    importanceTypeToPlot = "PercIncMSE",
    yAxisLabel = "% increase MSE\n upon feature removal (OOB)",
    # plotHitsAsCategorical = TRUE,
    # plotHitsAsCategorical = FALSE,
    # prevalenceCutoff = prevCutoff,
    # yAxisLegend = "Strain activity (ranked)")
    yAxisLegend = "Strain activity",
    predictor_target_map = species_target_map,
    numFeaturesToShow = 100,
    # targetOrder = targetOrder[!str_detect(targetOrder, "IS")],
    targetOrder = DrugsToHighlightInMainFigure,
    union = TRUE,
    taxonBasedAnalysis = TRUE)

ggsave(plot = plotObject, filename = here('results/plots/feature_importances_taxa_Tac_Eve_Sir_union.pdf'), width = 4.25, height = 8)

# Tac, Sir, Eve top feature comparison
# Get superset of top features
library(ggvenn)
tf <- evaluate_rf_feature_importances(
    # listOfModels = mOTUs_predictions_prev_mOTUs_of_interest[[1]],
    listOfModels = models_motus_all,
    modelsTopDown = tdm_models_motus_all,
    numberFeaturesData = numberFeaturesDataMotusAll,
    directionality = directionality_motus_all,
    returnData = TRUE,
    outputPath = '/g/scb/zeller/karcher/tmp/mipf.pdf',
    # width = 8,
    # height = 13,
    # testedStrainsAsReference = TRUE,
    importanceTypeToPlot = "PercIncMSE",
    yAxisLabel = "% increase MSE\n upon feature removal (OOB)",
    # plotHitsAsCategorical = TRUE,
    # plotHitsAsCategorical = FALSE,
    # prevalenceCutoff = prevCutoff,
    # yAxisLegend = "Strain activity (ranked)")
    yAxisLegend = "Strain activity",
    predictor_target_map = species_target_map,
    numFeaturesToShow = 100,
    # targetOrder = targetOrder[!str_detect(targetOrder, "IS")],
    targetOrder = DrugsToHighlightInMainFigure,
    # targetOrder = c("Tacrolimus")
    union = FALSE,
    taxonBasedAnalysis = TRUE)

topFeatureSuperset <- unique(tf$origM)

# Train and evaluate species models for this superset
predictions_and_models_motus_top <- get_model_performance(gene_abundances = tax_profiles_long %>%
    #filter(mOTUs_ID %in% topFeatureSuperset) %>%
    inner_join(data.frame(mOTUs_ID = topFeatureSuperset)) %>%
    rename(gene = mOTUs_ID, stoolDonor = Stool_donor) %>%
    group_by(gene) %>%
    filter(mean(abundance > 0) > 0.2) %>%
    filter(max(abundance) > 1E-3) %>%
    mutate(abundance = (abundance - mean(abundance)) / sd(abundance)) %>%
    filter(!is.na(abundance)) %>%
    ungroup(),
numResamp = 1,
targets = DrugsToShowInMainFigure,
predictor_target_map = species_target_map,
mod = "top_features",
modelType = "RF",
alsoReturnData = TRUE,
how = 'all_predictors',
numFeaturesForTopDownSpearmanCorComparisons = 100)
models_motus_top <- predictions_and_models_motus_top[[1]]
predictions_motus_top <- predictions_and_models_motus_top[[2]]
tdm_motus_top <- do.call('rbind', predictions_and_models_motus_top[[3]])
directionality_motus_top <- do.call(
    'rbind',
    map(predictions_and_models_motus_top[[4]], \(x) {
        x %>%
            as.data.frame() %>%
            rownames_to_column('gene')
    }))

directionality_motus_top <- directionality_motus_top %>%
    group_by(gene, target, dataType, seed) %>%
    summarize(directionality = mean(directionality)) %>%
    ungroup() %>%
    mutate(directionality = case_when(
        directionality > 0.5 ~ 'positive',
        directionality <= 0.5 ~ 'negative',
        .default = NA
    ))

# Contains top down models themselfes
tdm_models_motus_top <- predictions_and_models_motus_top[[5]]
tdm_models_motus_top <- do.call('rbind', map2(tdm_models_motus_top, names(tdm_models_motus_top), \(x, nnn) {
    do.call('rbind', map2(x, names(x), \(x, nn) {
        importance(x) %>%
            as.data.frame() %>%
            rownames_to_column('feature') %>%
            as_tibble() %>%
            mutate(numTopFeatures = nn)
    })) %>%
        as.data.frame() %>%
        as_tibble() %>%
        mutate(heldOutCommunity = str_split_fixed(nnn, "__", n = 4)[, 1]) %>%
        mutate(target = str_split_fixed(nnn, "__", n = 4)[, 2]) %>%
        mutate(dataType = str_split_fixed(nnn, "__", n = 4)[, 3]) %>%
        mutate(seed = str_split_fixed(nnn, "__", n = 4)[, 4]) %>%
        mutate(seed = str_split_fixed(seed, "__", n = 4)[, 1])
}))

# Evaluate feature importance values of models trained on superset
tf <- evaluate_rf_feature_importances(
    # listOfModels = mOTUs_predictions_prev_mOTUs_of_interest[[1]],
    listOfModels = models_motus_top,
    modelsTopDown = tdm_models_motus_top,
    numberFeaturesData = data.frame(target = c("Tacrolimus", "Sirolimus", "Everolimus"), numFeatures = length(topFeatureSuperset)),
    directionality = directionality_motus_top,
    returnData = TRUE,
    outputPath = '/g/scb/zeller/karcher/tmp/mipf.pdf',
    # width = 8,
    # height = 13,
    # testedStrainsAsReference = TRUE,
    importanceTypeToPlot = "PercIncMSE",
    yAxisLabel = "% increase MSE\n upon feature removal (OOB)",
    # plotHitsAsCategorical = TRUE,
    # plotHitsAsCategorical = FALSE,
    # prevalenceCutoff = prevCutoff,
    # yAxisLegend = "Strain activity (ranked)")
    yAxisLegend = "Strain activity",
    predictor_target_map = species_target_map,
    numFeaturesToShow = 100,
    # targetOrder = targetOrder[!str_detect(targetOrder, "IS")],
    targetOrder = DrugsToHighlightInMainFigure,
    # targetOrder = c("Tacrolimus")
    union = FALSE,
    taxonBasedAnalysis = TRUE)

## Scatter plot of feature importance values
library(GGally)

ggpairsplot <- tf %>%
        pivot_wider(id_cols = mOTUs_ID, names_from = target, values_from = importanceValue, values_fill = 0) %>%
        column_to_rownames(var = "mOTUs_ID") %>%
        as.data.frame()

plots <- list()
mi <- min(ggpairsplot)
ma <- max(ggpairsplot)
for (i1 in 1:dim(ggpairsplot)[2]) {
    for (i2 in 1:dim(ggpairsplot)[2]) {
        if (i1 > i2) {
            a <- colnames(ggpairsplot)[i1]
            b <- colnames(ggpairsplot)[i2]
            tmp <- ggpairsplot %>%
                select(all_of(c(a, b))) %>%
                as.data.frame() %>%
                as_tibble()
            co <- cor(tmp)[1,2] %>% round(3)
            plot <- ggplot(tmp, aes(x = !!sym(a), y = !!sym(b))) +
                geom_point(alpha = 0.5) +
                theme_presentation() +
                geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "lightgrey") +
                theme(axis.text = element_text(size = 6), axis.title = element_text(size = 6)) +
                xlab(a) +
                ylab(b) +
                annotate('text', x = 0, y = ma * 0.95, hjust = 0, label = str_c("Pearson's r: ", co), size = 2) +
                xlim(c(mi, ma)) +
                ylim(c(mi, ma)) +
                
                NULL
            plots[[length(plots) + 1]] <- plot
        }
    }
}
ggsave(plot = wrap_plots(plots, ncol = 3), filename = here('results/plots/ggpairsplot.pdf'), width = 8 * 0.6, height = 2.75 * 0.6)

# genes
predictions_genes <- rbind(
    predictions_gene_characterized %>% mutate(features = "Known relevant genes"),
    predictions_abundant_genes %>% mutate(features = "All genes")
) %>%
    mutate(condition = case_when(
        oxygen == "MA" ~ "Microaerobic",
        oxygen == "AA" ~ "Anaerobic"
    )) %>%
    select(-oxygen)

metrics_genes <- predictions_genes %>%
    group_by(target, features) %>%
    summarize(spearmanCore = cor(prediction, truth, method = "spearman"),
        pearsonCore = cor(prediction, truth),
        R2 = pearsonCore^2) %>%
    filter(!str_detect(target, "IS"))

# performance_by_feature_abundant_genes <- tdm_abundant_genes %>%
#     group_by(numFeatures, target, dataType, seed) %>%
#     summarize(spearmanCor = cor(prediction, truth, method = 'spearman')) %>%
#     unnest() %>%
#     arrange(target, numFeatures) %>%
#     filter(!str_detect(target, "IS"))


# plotObject <- ggplot(
#     data = performance_by_feature_abundant_genes,
#     aes(x = numFeatures, y = spearmanCor, group = 1)) +
#     facet_grid(. ~ target) +
#     geom_line() +
#     geom_abline(data = metrics_genes %>% filter(features == "abundant genes"), aes(intercept = spearmanCore, slope = 0, color = 'red')) +
#     theme_publication()

# ggsave(plot = plotObject, filename = "/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/plots_community_activity/abundant_genes_performance_by_feature.pdf", width = 10, height = 4)


predictions_genes$target <- as.factor(predictions_genes$target)
metrics_genes$target <- factor(metrics_genes$target, levels = levels(predictions_genes$target))
targetOrder <- levels(predictions_genes$target)

# library(ggembl)
# geneScatters <- ggplot() +
#     geom_point(data =
#         predictions_genes %>%
#             filter(features == 'All genes') %>%
#             inner_join(data.frame(target = unique(predictions_motus$target))) %>%
#             filter(target %in% drugsToComputeMetricsFor) %>%
#             mutate(target = factor(target, levels = drugsToComputeMetricsFor)), aes(x = truth, y = prediction, color = condition), alpha = 1) +
#     # geom_line(data = predictions_motus %>% filter(target == "Tacrolimus") %>% filter(features == 'all species'), aes(x = truth, y = prediction, group = left_out_community), alpha = 1) +
#     # theme_classic() +
#     theme_presentation() +
#     theme(axis.ticks.length = unit(0.5, "mm")) +
#     geom_text(data =
#         metrics_genes %>%
#             filter(features == 'abundant genes') %>%
#             inner_join(data.frame(target = unique(predictions_motus$target))) %>%
#             filter(target %in% drugsToComputeMetricsFor) %>%
#             mutate(target = factor(target, levels = drugsToComputeMetricsFor)), aes(x = -Inf, y = Inf, label = str_c("Spearman-cor: ", round(spearmanCore, 3))), hjust = -0.05, vjust = 1.2, color = 'black', size = 6 * 0.353) +
#     scale_color_manual(values = aa_ma_colors) +
#     facet_wrap(target ~ features, ncol = 1, scales = 'free', strip.position = "right") +
#     # ggtitle("Abundant genes") +
#     xlab("Truth") +
#     ylab("Prediction")

# ggsave(plot = geneScatters, filename = "/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/plots_community_activity/scatters_genes.pdf", width = 4.5, height = 6.5)
# ggsave(plot = geneScatters, filename = "/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/plots_community_activity/scatters_genes.pdf", width = 8, height = 12)

# plotObject <- evaluate_rf_feature_importances(
#     # listOfModels = mOTUs_predictions_prev_mOTUs_of_interest[[1]],
#     listOfModels = models_abundant_genes,
#     returnData = FALSE,
#     outputPath = '/g/scb/zeller/karcher/tmp/mipf.pdf',
#     # width = 8,
#     # height = 13,
#     # testedStrainsAsReference = TRUE,
#     importanceTypeToPlot = "PercIncMSE",
#     yAxisLabel = "% increase MSE\n upon feature removal (OOB)",
#     # plotHitsAsCategorical = TRUE,
#     # plotHitsAsCategorical = FALSE,
#     # prevalenceCutoff = prevCutoff,
#     # yAxisLegend = "Strain activity (ranked)")
#     yAxisLegend = "Gene activity",
#     predictor_target_map = gene_target_map,
#     numFeaturesToShow = 15,
#     targetOrder = targetOrder,
#     taxonBasedAnalysis = FALSE)

# ggsave(plot = plotObject, filename = "/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/plots_community_activity/feature_importances_genes.pdf", width = 7.5, height = 7)

# spearmanBars <- ggplot(data = metrics_genes, aes(x = target, y = spearmanCore, fill = features)) +
#     geom_bar(stat = "identity", position = "dodge") +
#     theme_presentation() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     ylab('Spearman correlation') +
#     xlab("Drug")
# ggsave(plot = spearmanBars, filename = "/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/plots_community_activity/spearman_bars_genes.pdf", width = 6.5, height = 3.5)

data <- rbind(
    do.call('rbind',
        map2(all_predictions_and_models_strains_active_random,
            1:length(all_predictions_and_models_strains_active_random), \(x, y) {
                x[[2]] %>%
                    # x %>%
                    mutate(seed = y)
            })) %>%
        group_by(target, data_type, seed) %>%
        summarize(spearmanCore = cor(prediction, truth, method = "spearman")) %>%
        mutate(`Feature set` = 'motus', features = "Identified active strains", `Label\npermuted` = TRUE),
    do.call('rbind',
        map2(all_predictions_and_models_motus_all_random,
            1:length(all_predictions_and_models_motus_all_random), \(x, y) {
                x[[2]] %>%
                    # x %>%
                    mutate(seed = y)
            })) %>%
        group_by(target, data_type, seed) %>%
        summarize(spearmanCore = cor(prediction, truth, method = "spearman")) %>%
        mutate(`Feature set` = 'motus', features = 'All species', `Label\npermuted` = TRUE),
    do.call('rbind',
        map2(all_predictions_and_models_gene_characterized_random,
            1:length(all_predictions_and_models_gene_characterized_random), \(x, y) {
                x %>%
                    mutate(seed = y)
            })) %>%
        group_by(target, data_type, seed) %>%
        summarize(spearmanCore = cor(prediction, truth, method = "spearman")) %>%
        mutate(`Feature set` = 'genes', features = 'Known relevant genes', `Label\npermuted` = TRUE),
    # do.call('rbind',
    #     map2(all_predictions_and_models_abundant_genes_random,
    #         1:length(all_predictions_and_models_abundant_genes_random), \(x, y) {
    #             x %>%
    #                 mutate(seed = y)
    #         })) %>%
    all_predictions_and_models_abundant_genes_random %>%
        mutate(seed = as.integer(seed)) %>%
        group_by(target, data_type, seed) %>%
        summarize(spearmanCore = cor(prediction, truth, method = "spearman")) %>%
        mutate(`Feature set` = 'genes', features = 'All genes', `Label\npermuted` = TRUE)
) %>%
    # mutate(features = factor(features, levels = (c("Identified active strains", "All species", "Known relevant genes", 'All genes')))) %>%
    select(-data_type) %>%
    # group_by(target, data_type, `Feature set`, features, `Label\npermuted`) %>%
    # summarize(tenPerc = quantile(spearmanCore, 0.05),
    #     ninetyPerc = quantile(spearmanCore, 0.95)) %>%
    # pivot_longer(c(tenPerc, ninetyPerc), names_to = "quantile", values_to = "value")
    filter(target %in% drugsToComputeMetricsFor) %>%
    mutate(target = factor(target, levels = drugsToComputeMetricsFor)) %>%
    identity() %>%
    mutate(features = factor(features, levels = (c("Identified active strains", "All species", "Known relevant genes", 'All genes'))))

data2 <- rbind(
    metrics_motus %>% mutate(`Feature set` = 'motus', `Label\npermuted` = FALSE),
    metrics_genes %>% mutate(`Feature set` = 'genes', `Label\npermuted` = FALSE)) %>%
    #    rbind(data) %>%
    select(-pearsonCore, -R2) %>%
    # inner_join(data, by = c("target", 'features'), suffix = c(".real", ".permuted")) %>%
    filter(target %in% drugsToComputeMetricsFor) %>%
    mutate(target = factor(target, levels = drugsToComputeMetricsFor)) %>%
    # mutate(spearmanCore.permuted = ifelse(spearmanCore.permuted < 0, 0, spearmanCore.permuted)) %>%
    # mutate(spearmanCorDiff = spearmanCore.real - spearmanCore.permuted)
    arrange(target) %>%
    identity() %>%
    mutate(features = factor(features, levels = (c("Identified active strains", "All species", "Known relevant genes", 'All genes'))))

n <- as.numeric(data$features)
nn <- data$features
n <- data.frame(n = n, nn = nn) %>%
    unique()
ve <- n$nn
names(ve) <- n$n

library(ggbeeswarm)

spearmanBarsAll <- ggplot(
    data =
        data2 %>%
            filter(target %in% DrugsToShowInMainFigure) %>%
            mutate(target = factor(target, levels = DrugsToShowInMainFigure)),
    aes(x = as.numeric(features), y = spearmanCore)) +
    theme(axis.text.x = element_blank()) +
    theme_presentation() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab('Spearman correlation') +
    facet_grid(target ~ ., switch = "y") +
    xlab("Feature sets") +
    # coord_flip() +
    ylab("Spearman correlation between prediction and truth") +
    geom_errorbar(data = data %>%
        group_by(target, data_type, `Feature set`, features, `Label\npermuted`) %>%
        summarize(fivePerc = quantile(spearmanCore, 0.05),
            ninetyFivePerc = quantile(spearmanCore, 0.95)), aes(x = as.numeric(features) + 0.3, ymin = fivePerc, ymax = ninetyFivePerc), width = 0.2, inherit.aes = F) +
    geom_beeswarm(data = data, aes(x = as.numeric(features), y = spearmanCore), inherit.aes = F, color = 'grey', alpha = 0.3, size = 1) +
    geom_point(color = 'red') +
    # scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "grey"))
    scale_y_continuous(breaks = seq(floor(min(data2$spearmanCore)), ceiling(max(data2$spearmanCore)), by = 0.2),
        labels = function(x) sprintf("%.1f", x), position = 'right') +
    scale_x_continuous(breaks = as.numeric(names(ve)), labels = unname(ve)) +
    NULL
ggsave(plot = spearmanBarsAll, filename = here("results/plots/spearman_bars_all.pdf"), width = 2.5, height = 8.2)
