loadFromScratch <- FALSE

getGeneTypeMapping <- function(gene) {
    return(case_when(
        str_detect(gene, "CLOSCI_00899") ~ "Steroids___CLOSCI_00899",
        str_detect(gene, "CLOSCI_00900") ~ "Steroids___CLOSCI_00900",
        str_detect(gene, "bu0920") ~ "MMF___bu0920",
        str_detect(gene, "bu0921") ~ "MMF___bu0921",
        str_detect(gene, "WP_005929331.1") ~ "Tacrolimus___WP_005929331.1",
        str_detect(gene, "WP_005934826.1") ~ "Tacrolimus___WP_005934826.1",
    ))
}

get_model_performance <- function(
    gene_abundances,
    numResamp = 1,
    seed = 1,
    # targets = unique(gene_target_map$target),
    targets = allTargets,
    predictor_target_map = NULL,
    mod = NULL,
    modelType = "RF",
    returnData = FALSE,
    pivotWide = TRUE,
    preFilterCor = FALSE,
    alsoReturnData = FALSE,
    onlyReturnData = FALSE,
    prevalenceCutoff = NULL,
    permuteOutcome = FALSE,
    stoolDonors = unique(tax_profiles_long$Stool_donor),
    numFeaturesForTopDownSpearmanCorComparisons = 1,
    how = 'only_testedactive',
    dataTypesToTest = c("drug_degradation")) {
    # if (numResamp == 1) {
    #    print("numResamp is 1, only for testing and quick iterations...")
    # }
    print(str_c('numResamp is ', numResamp))

    ps <- list()
    models <- list()
    topDownPredictionsAndTruthsAll <- list()
    directionalityAll <- list()
    modelsTopDownAll <- list()

    set.seed(seed)

    # for (dataType in c("drug_degradation", "metabolite_buildup")) {
    for (dataType in dataTypesToTest) {
        drug_data <- dataTypes[[dataType]]
        # for (predictor_target_index in 1:dim(predictor_target_map)[1]) {
        for (target in targets) {
            print(target)

            dd <- drug_data %>% filter(Drug == target) %>% select(-Drug, -Metabolized)
            if (permuteOutcome) {
                print("CAREFUL!! PERMUTING DEPENDNET VARIABLE!!")
                dd$effect_size <- sample(dd$effect_size)
            }
            if (pivotWide) {
                print("Pivoting wide...")
                if (how == 'only_testedactive') {
                    print("Careful! Taking only tested/active strains/genes!")
                    candidatePredictors <- predictor_target_map[predictor_target_map$target == target, ]
                    if (dim(candidatePredictors)[1] == 0) {
                        print(str_c("Skipping target ", target, " because no predictors available"))
                        next
                    }
                    data <- gene_abundances %>%
                        inner_join(candidatePredictors %>% select(predictor), by = c('gene' = 'predictor')) %>%
                        # pivot_wider(id_cols = c(stoolDonor, oxygen), names_from = gene, values_from = abundance) %>%
                        dcast(stoolDonor + oxygen ~ gene, value.var = "abundance") %>%
                        as_tibble() %>%
                        inner_join(dd, by = c("stoolDonor" = "Stool_donor", 'oxygen'))
                } else if (how == 'all_predictors') {
                    data <- gene_abundances %>%
                        # inner_join(candidatePredictors %>% select(predictor), by = c('gene' = 'predictor')) %>%
                        pivot_wider(id_cols = c(stoolDonor, oxygen), names_from = gene, values_from = abundance) %>%
                        inner_join(dd, by = c("stoolDonor" = "Stool_donor", 'oxygen'))
                } else {
                    # data <- gene_abundances %>%
                    #     # inner_join(candidatePredictors %>% select(predictor), by = c('gene' = 'predictor')) %>%
                    #     pivot_wider(id_cols = c(stoolDonor, oxygen), names_from = gene, values_from = abundance) %>%
                    #     inner_join(drug_data %>% filter(Drug == target) %>% select(-Drug, -Metabolized), by = c("stoolDonor" = "Stool_donor", 'oxygen'))
                    asdaasddas
                }
            } else {
                if (how == "only_testedactive") {
                    exit("This is not supported... exiting")
                }
                data <- gene_abundances %>%
                    # t() %>%
                    # as.data.frame() %>%
                    # rownames_to_column('stoolDonor') %>%
                    inner_join(dd, by = c("stoolDonor" = "Stool_donor", 'oxygen'))
            }
            data$label <- str_c(data$stoolDonor, "__", data$oxygen)
            data$stoolDonor <- NULL
            data$oxygen <- NULL
            print(str_c("Training RF for ", target, " and ", dataType, " and seed ", seed, " with ", dim(data)[1], " samples and ", dim(data)[2], " features"))

            if (dim(data)[1] == 0) {
                print(str_c("Skipping target ", target, " because no predictors available"))
                next
            }
            for (communityID in stoolDonors) {
                skip_to_next <<- FALSE
                trainData <- data %>%
                    ungroup() %>%
                    filter(!str_detect(label, communityID)) %>%
                    column_to_rownames('label')
                testData <- data %>%
                    ungroup() %>%
                    filter(str_detect(label, communityID)) %>%
                    column_to_rownames('label')


                es <- trainData$effect_size
                p <- map(trainData[1:(dim(trainData)[2] - 1)], \(x) cor(x, es, method = "spearman"))

                if (preFilterCor) {
                    # If this is et, we dont choose candidate predictors but instead use all genes but
                    # pre-selected using wilcox tests
                    print("Doing cor-based filtering of features in training fold... not cor.test anymore!")
                    # pInd <- which(map_dbl(p, \(x) x$p.value) < 0.001)
                    pInd <- head(order(map_dbl(p, \(x) x), decreasing = TRUE), 1000)
                    # pInd <- 1:10000
                    print(str_c("Keeping ", length(pInd), " features to train RF..."))
                    trainData <- trainData[, c(pInd, dim(trainData)[2])]
                    testData <- testData[, c(pInd, dim(testData)[2])]
                }

                # browser()

                if (modelType == "RF") {
                    if (any(is.na(trainData))) {
                        print("This actually happened... :(")
                        browser()
                    }
                    tryCatch(
                        modelFin <- randomForest(
                            as.formula(str_c("effect_size ~ ", str_c(colnames(trainData)[colnames(trainData) != "effect_size"], collapse = " + "))),
                            # as.formula(str_c("effect_size ~ ", str_c(colnames(trainData)[1:10000], collapse = " + "))),
                            data = trainData,
                            importance = TRUE)
                        ,
                        error = function(e) {
                            skip_to_next <<- TRUE
                        })
                    if (skip_to_next) {
                        print(str_c("Skipping ", communityID, " because of error in RF"))
                        next
                    }

                    # print("Doing cor.test-based testing of features in training fold for forward-selection...")

                    # corTestFeatureOrder <- order(map_dbl(p, \(x) x$p.value))
                    corTestFeatureOrder <- order(as.data.frame(importance(modelFin))[['%IncMSE']], decreasing = TRUE)
                    featuresOrdered <- colnames(trainData)[1:length(colnames(trainData))][corTestFeatureOrder]
                    directionality <- map_lgl(p, \(x) x > 0)
                    names(directionality) <- colnames(trainData)[1:(length(colnames(trainData)) - 1)]

                    modelsTopDown <- list()
                    # print("Top-down feature selection...")
                    # for (modelTopDownIndex in 1:min(10, length(featuresOrdered))) {
                    if (!permuteOutcome) {
                        for (modelTopDownIndex in 1:min(numFeaturesForTopDownSpearmanCorComparisons, length(featuresOrdered))) {
                            # print(str_c("effect_size ~ ", str_c(featuresOrdered[1:modelTopDownIndex], collapse = " + "), collapse = " + "))
                            tryCatch(model <- randomForest(
                                as.formula(str_c("effect_size ~ ", str_c(featuresOrdered[1:modelTopDownIndex], collapse = " + "), collapse = " + ")),
                                data = trainData,
                                importance = TRUE),
                            error = function(e) {
                                skip_to_next <<- TRUE
                            })
                            if (skip_to_next) {
                                next
                            }
                            modelsTopDown[[length(modelsTopDown) + 1]] <- model
                        }
                        names(modelsTopDown) <- 1:length(modelsTopDown)
                        modelsTopDownAll[[length(modelsTopDownAll) + 1]] <- modelsTopDown
                        names(modelsTopDownAll)[length(modelsTopDownAll)] <- str_c(communityID, "__", target, "__", dataType, "__", seed, "__")
                    }
                } else if (modelType == "LASSO") {
                    lambda.min <- cv.glmnet(
                        trainData %>% select(-effect_size) %>%
                            as.data.frame() %>%
                            as.matrix(),
                        y = trainData %>% pull(effect_size),
                        # lasso
                        alpha = 1,
                        # non-negative
                        lower.limits = 0,
                        nfolds = dim(trainData)[1])$lambda.min
                    model <- glmnet(
                        trainData %>% select(-effect_size) %>%
                            as.data.frame() %>%
                            as.matrix(),
                        y = trainData %>% pull(effect_size),
                        # lasso
                        alpha = 1,
                        #
                        lambda = lambda.min,
                        # non-negative
                        lower.limits = 0,
                    )
                } else {
                    sadassdsd
                }
                prediction <- predict(modelFin, testData %>% select(-effect_size) %>% as.matrix())
                if (!permuteOutcome) {
                    topDownPredictionsAndTruths <- do.call('rbind', map2(modelsTopDown, 1:length(modelsTopDown), function(x, l) {
                        predict(
                            x,
                            testData) %>%
                            as.data.frame() %>%
                            rename('prediction' = '.') %>%
                            mutate(numFeatures = l) %>%
                            mutate(left_out_community = communityID, target = target, dataType = dataType, seed = seed)
                    })) %>%
                        as.data.frame() %>%
                        rownames_to_column('meta') %>%
                        mutate(oxygen = case_when(
                            str_detect(meta, "__MA") ~ "MA",
                            str_detect(meta, "__AA") ~ "AA",
                            .default = NA
                        )) %>% left_join(
                            testData %>%
                                select(effect_size) %>%
                                mutate(oxygen = case_when(
                                    str_detect(rownames(.), "MA") ~ "MA",
                                    str_detect(rownames(.), "AA") ~ "AA",
                                    .default = NA
                                )), by = 'oxygen') %>%
                        rename(truth = effect_size)


                    topDownPredictionsAndTruthsAll[[length(topDownPredictionsAndTruthsAll) + 1]] <- topDownPredictionsAndTruths
                }

                directionalityAll[[length(directionalityAll) + 1]] <- as.data.frame(directionality) %>%
                    mutate(left_out_community = communityID, target = target, dataType = dataType, seed = seed)
                # truth <- testData$effect_size[1]
                prediction <- as.data.frame(prediction)
                truth <- testData %>% select(effect_size)
                stopifnot(rownames(prediction) == rownames(truth))

                ps[[length(ps) + 1]] <- cbind(prediction, truth)
                models[[length(models) + 1]] <- modelFin
                # if (modelType == "LASSO") {
                #     nonZeroCoef[[length(nonZeroCoef) + 1]] <- sum(coef(model) != 0)
                #     names(nonZeroCoef)[length(nonZeroCoef)] <- str_c(communityID, "__", target, "__", dataType, seed = seed)
                # }
                ##################### FIX THIS
                #####################
                #####################
                names(ps)[length(ps)] <- str_c(communityID, "__", target, "__", dataType, seed = seed)
                names(models)[length(models)] <- str_c(communityID, "__", target, "__", dataType, seed = seed)
            }
        }
    }
    pss <- map(ps, \(x) {
        rownames(x) <- str_replace(rownames(x), ".*__", "__")
        return(x)
    })
    psD <- do.call('rbind', pss)
    rownames(psD) <- str_replace(rownames(psD), "[0-9]+[.]__", "__")
    colnames(psD)[length(colnames(psD))] <- "truth"
    psD <- psD %>%
        as.data.frame() %>%
        rownames_to_column('meta') %>%
        pivot_longer(-c(truth, meta)) %>%
        # mutate(oxygen = str_split_fixed(meta, "__", n = 4)[, 4]) %>%
        # mutate(oxygen = ifelse(str_detect(name, "AA"), "AA", "MA")) %>%
        # mutate(meta = str_c(meta, "__", oxygen)) %>%
        select(-name) %>%
        rename(prediction = value) %>%
        as.data.frame() %>%
        column_to_rownames('meta')
    psD <- psD %>%
        as.data.frame() %>%
        rownames_to_column("meta") %>%
        mutate(left_out_community = str_split_fixed(meta, "__", n = 4)[, 1]) %>%
        mutate(target = str_split_fixed(meta, "__", n = 4)[, 2]) %>%
        mutate(data_type = str_split_fixed(meta, "__", n = 4)[, 3]) %>%
        mutate(data_type = str_replace(data_type, "[.]1", "")) %>%
        mutate(data_type = str_replace(data_type, "[0-9]{1}$", "")) %>%
        mutate(oxygen = str_split_fixed(meta, "__", n = 4)[, 4]) %>%
        mutate(seed = str_replace(meta, ".*buildup", "")) %>%
        mutate(seed = str_replace(seed, ".*degradation", "")) %>%
        mutate(seed = str_replace(seed, "[.]1", "")) %>%
        mutate(seed = str_replace(seed, "__.*", "")) %>%
        mutate(data_type = str_replace(data_type, "[0-9]{1}$", "")) %>%
        select(-meta) %>%
        as_tibble()

    for (dt in dataTypesToTest) {
        rf_eval_plot <- ggplot() +
            theme_classic() +
            facet_wrap(. ~ target, ncol = 3, scales = 'free') +
            geom_point(data = psD %>% filter(data_type == dt) %>% filter(seed == 1), aes(x = truth, y = prediction, color = oxygen), alpha = 0.3) +
            geom_line(data = psD %>% filter(data_type == dt) %>% filter(seed == 1), aes(x = truth, y = prediction, group = left_out_community)) +
            geom_text(data = psD %>%
                filter(data_type == dt) %>%
                group_by(target, data_type) %>%
                summarize(
                    core = cor(prediction, truth),
                    spearmanCore = cor(prediction, truth, method = "spearman"),
                    y = Inf,
                    x = -Inf) %>%
                mutate(core = round(core, 3)), aes(x = x, y = y, vjust = 1.0, hjust = 0, label = str_c("pearson: ", core, "\nSpearman: ", round(spearmanCore, 3), "\nR2: ", round(core^2, 3))), color = 'red', size = 3) +
            xlim(c(0, NA)) +
            ggtitle(str_c("Data type: ", dt, "\nModel: ", mod)) +
            xlab(str_c(effect_size_measure_map[[dt]], " (truth)")) +
            ylab(str_c(effect_size_measure_map[[dt]], " (prediction)"))
        # if (!onlyReturnData) {
        # ggsave(plot = rf_eval_plot, filename = str_c("/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/plots_community_activity/rf_eval_plot_example", dt, "__", mod, "__", prevalenceCutoff, ".pdf"), width = 10, height = 8)
        # ggsave(plot = rf_eval_plot, filename = str_c(here(str_c("results/plots/plots_community_activity/rf_eval_plot_example", dt, "__", mod, "__", prevalenceCutoff, ".pdf"))), width = 10, height = 8)#
        # }

    }

    library(patchwork)
    # if (!onlyReturnData) {
    #     ggsave(plot = (psD %>%
    #         group_by(target, data_type, seed) %>%
    #         summarize(core = cor(prediction, truth)) %>%
    #         mutate(R2 = core^2) %>%
    #         group_by(target, data_type) %>%
    #         summarize(fracPosCor = mean(core > 0)) %>%
    #         ggplot(aes(x = target, y = fracPosCor, color = data_type)) +
    #         theme_classic() +
    #         theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #         geom_point(position = position_dodge(width = 0.5)) +
    #         ylab("Fraction of models with\npositive slope")) + (psD %>%
    #         group_by(target, data_type, seed) %>%
    #         summarize(core = cor(prediction, truth)) %>%
    #         mutate(R2 = core^2) %>%
    #         ggplot(aes(x = target, y = R2, fill = data_type)) +
    #         theme_classic() +
    #         theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #         geom_boxplot()) + plot_layout(ncol = 1, heights = c(1, 2)), filename = str_c("/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/plots_community_activity/rf_eval_plot_boxplot__", mod, "__", prevalenceCutoff, ".pdf"), width = 6.5, height = 6.5)
    # }

    if (alsoReturnData || onlyReturnData) {
        return(list(models, psD, topDownPredictionsAndTruthsAll, directionalityAll, modelsTopDownAll))
    }
}

evaluate_rf_feature_importances <- function(
    listOfModels,
    modelsTopDown,
    numberFeaturesData,
    directionality = NULL,
    returnData = FALSE,
    outputPath = NULL,
    shadeByType = "drug_degradation",
    testedStrainsAsReference = FALSE,
    importanceTypeToPlot = "PercIncMSE",
    plotHitsAsCategorical = FALSE,
    yAxisLabel = NULL,
    prevalenceCutoff = NULL,
    yAxisLegend = NULL,
    numFeaturesToShow = 5,
    predictor_target_map = NULL,
    targetOrder = NULL,
    taxonBasedAnalysis = TRUE,
    union = FALSE) {

    prep <- function(union = FALSE, fs = NULL) {
        data <- modelsTopDown %>%
            left_join(numberFeaturesData %>% select(target, numFeatures), by = 'target') %>%
            filter(target %in% targetOrder)
        if (union) {
            data <- data %>%
                filter(numTopFeatures == 100) %>%
                # filter(mOTUs_ID %in% fs)
                inner_join(fs, by = c("feature" = 'mOTUs_ID'))
        } else {
            data <- data %>%
                filter(numTopFeatures == numFeatures)
        }
        data <- data %>%
            group_by(target, dataType) %>%
            nest() %>%
            mutate(data = map(data, \(x) {
                numModels <- length(unique(x$heldOutCommunity))
                # browser()
                y <- x %>%
                    pivot_wider(
                        id_cols = feature,
                        names_from = heldOutCommunity,
                        values_from = `%IncMSE`,
                        values_fill = 0) %>%
                    pivot_longer(-feature) %>%
                    group_by(feature) %>%
                    summarize(mean = mean(value)) %>%
                    rename(PercIncMSE = mean)
                x <- x %>%
                    group_by(feature) %>%
                    tally() %>%
                    mutate(fractionOfModels = n / numModels) %>%
                    select(-n)

                return(x %>% left_join(y, by = 'feature'))
            })) %>%
            unnest() %>%
            rename(importanceValue = PercIncMSE, mOTUs_ID = feature) %>%
            mutate(importanceType = 'PercIncMSE') %>%
            filter(fractionOfModels > 0.2) %>%
            filter(importanceType == importanceTypeToPlot) %>%
            left_join(species_motus_link_motus_extender %>% select(mOTUs_ID, Species), by = 'mOTUs_ID') %>%
            # I'm updating this with the species_motus link from motus-extender.
            # Make sure to understand which mOTUs is lost
            left_join(motus3.0_taxonomy %>% select(mOTUs_ID, Species), by = 'mOTUs_ID') %>%
            # equivalent to mutate(Species = ifelse(is.na(Species.x), Species.y, Species.x))
            # I find the beforementioned, commented-out line more readable as you don't need to know what coalesce does... but
            # but coalasce would work on more than 2 vectors also, which I guess would be useful to establish a sort of hierarchy
            # between the vectors
            mutate(Species = coalesce(Species.x, Species.y)) %>%
            select(-Species.x, -Species.y) %>%
            # Trying to get a hand on these fucking Species strings in motus
            group_by(mOTUs_ID, Species) %>%
            nest() %>%
            # mutate(SpeciesNameUnique = map_chr(Species, \(x) {
            #     tmp <- gsub(x = x, "\\[.*\\]", "")
            #     tmp <- gsub(x = tmp, "$", str_c(" ", round(runif(n = 1, 0, 10000000))))
            #     return(tmp)
            # })) %>%
            mutate(SpeciesNameUnique = Species) %>%
            unnest() %>%
            # This is already done above
            # group_by(target, mOTUs_ID, Species, dataType, importanceType) %>%
            # summarize(importanceValue = mean(importanceValue), fractionOfModels = mean(fractionOfModels)) %>%
            group_by(target, dataType, importanceType) %>%
            nest() %>%
            mutate(data = map(data, \(x) {
                return(x %>%
                    ungroup() %>%
                    arrange(desc(importanceValue)) %>%
                    head(numFeaturesToShow))

            })) %>%
            unnest()
        return(data)
    }

    if (union) {
        topFeatures <- data.frame(mOTUs_ID = prep() %>% pull(mOTUs_ID) %>% unique())
        data <- prep(union = TRUE, fs = topFeatures)
    } else {
        data <- prep()
    }
    data <- data %>%
        left_join(predictor_target_map %>% mutate(active = TRUE), by = c('target' = 'target', 'mOTUs_ID' = 'predictor')) %>%
        left_join(predictor_target_map %>% mutate(tested = TRUE) %>% select(-target) %>% distinct(), by = c('mOTUs_ID' = 'predictor')) %>%
        mutate(featureType = case_when(
            active ~ "active",
            tested ~ "inactive",
            .default = 'untested'
        )) %>%
        select(-active, -tested) %>%
        # mutate(m = mOTUs_ID) %>%
        mutate(origM = mOTUs_ID) %>%
        mutate(m = map_chr(mOTUs_ID, \(x) str_split(x, "_")[[1]][length(str_split(x, "_")[[1]])])) %>%
        mutate(Species = str_replace(Species, '^ ', "")) %>%
        mutate(Species = str_replace(Species, "\\[Eubacterium\\]", "Eubacterium")) %>%
        mutate(Species = str_replace(Species, "\\[Clostridium\\]", "Clostridium")) %>%
        mutate(Species = str_replace(Species, "\\[Ruminococcus\\]", "Ruminococcus")) %>%
        mutate(mOTUs_ID = map2_chr(Species, 1:length(Species), \(x, dummy) {
            x <- case_when(
                str_detect(x, " sp") & !str_detect(x, "sedis") ~ get_first_two_fields(x),
                str_detect(x, "sedis") ~ get_first_four_fields_but_truncate_first(x),
                TRUE ~ get_first_two_fields_but_truncate_first(x))
            # str_detect(x, " sp ") & !str_detect(x, "sedis") ~ 'a',
            # str_detect(x, "sedis") ~ 'b',
            # TRUE ~ 'c')
            return(x)
            #     return(str_c(x, str_c(rep(" ", dummy), collapse = ""), collapse = ""))
        })) %>%
        # mutate(mOTUs_ID = str_c(mOTUs_ID, " \n[", m, ']')) %>%
        mutate(mOTUs_ID = str_c(mOTUs_ID, " [", m, ']')) %>%
        identity()
    if (union) {
        data <- data %>%
            mutate(mOTUs_ID = factor(mOTUs_ID, levels = unique(data$mOTUs_ID)))
    }

    if (returnData) {
        return(data)
    }
    plots <- list()

    for (ta in targetOrder) {
        # for (ta in c("Tacrolimus")) {
        tmp <- data %>% filter(target == ta)
        tmp <- tmp %>%
            arrange(desc(importanceValue))
        if (!union) {
            tmp$mOTUs_ID <- factor(tmp$mOTUs_ID, levels = rev(tmp$mOTUs_ID))
        } else {
            tmp$mOTUs_ID <- factor(tmp$mOTUs_ID, levels = rev(levels(data$mOTUs_ID)))
        }

        tmp <- tmp %>%
            left_join(directionality %>% mutate(association = case_when(
                # Have to flip - since LOW AUC is HIGH metabolism.
                directionality == 'positive' ~ 'negative',
                directionality == 'negative' ~ 'positive',
                is.na(directionality) ~ 'untested'
            )), by = c("origM" = 'gene', 'target'))
        plo <- ggplot(tmp, aes(x = mOTUs_ID, y = importanceValue, fill = featureType, color = association)) +
            theme_presentation() +
            facet_wrap(~target, ncol = 1, strip.position = 'right') +
            geom_bar(stat = "identity") +
            coord_flip() +
            theme_presentation() +
            # theme(axis.text.x = element_text(angle = 90, hjust = -1, size = 6)) +
            # geom_text(data = tmp %>% filter(featureType == "inactive"), aes(x = mOTUs_ID, y = max(tmp$importanceValue) * 0.025, label = mOTUs_ID), color = 'lightgrey', angle = 90, hjust = 'left', size = 6 * 0.353, alpha = 0.6, lineheight = 0.75) +
            # geom_text(data = tmp %>% filter(featureType == "active"), aes(x = mOTUs_ID, y = max(tmp$importanceValue) * 0.025, label = mOTUs_ID), color = 'black', angle = 90, hjust = 'left', size = 6 * 0.353, alpha = 0.6, lineheight = 0.75) +
            # geom_text(data = tmp %>% filter(featureType == "untested"), aes(x = mOTUs_ID, y = max(tmp$importanceValue) * 0.025, label = mOTUs_ID), color = 'black', angle = 90, hjust = 'left', size = 6 * 0.353, alpha = 0.6, lineheight = 0.75) +
            scale_fill_manual(values = activity_colors, drop = FALSE) +
            scale_color_manual(values = c("positive" = 'orange', "negative" = 'blue')) +
            theme(
                # axis.title.x = element_blank(),
                # axis.text.x = element_blank(),
                # axis.ticks.x = element_blank(),
                axis.ticks.length = unit(0.5, "mm")) +
            ylab("") +
            # ggtitle(ta) +
            theme(axis.title.y = element_blank()) +
            theme(legend.position = 'none') +
            theme(
                axis.text = element_text(size = 6),
                axis.ticks = element_line(size = 0.5),
            ) +
            scale_x_discrete(drop = F)
        NULL
        plo2 <- ggplot() +
            theme_presentation() +
            geom_bar(data = tmp, aes(x = mOTUs_ID, y = 1), stat = 'identity', fill = NA, color = 'NA') +
            geom_bar(data = tmp, aes(x = mOTUs_ID, y = fractionOfModels), stat = 'identity', fill = 'lightgrey', color = NA) +
            theme(
                axis.title.y = element_blank(),
                # axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                # axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.length = unit(0.5, "mm")
            ) +
            ylab("") +
            scale_y_continuous(breaks = c(0, 1)) +
            coord_flip() +
            theme(axis.title.y = element_blank()) +
            theme(legend.position = 'none') +
            theme(
                axis.text = element_text(size = 6),
                axis.ticks = element_line(size = 0.5),
            ) +
            scale_x_discrete(drop = F)

        plots[[length(plots) + 1]] <- (plo | plo2) + plot_layout(nrow = 1, widths = c(5, 1))
    }
    return(wrap_plots(plots, ncol = 1) + plot_layout(guides = 'collect') & theme(legend.position = 'right'))
}
