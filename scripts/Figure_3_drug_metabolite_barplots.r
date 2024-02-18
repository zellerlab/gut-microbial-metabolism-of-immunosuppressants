library(tidyverse)
library(ggplot2)
library(patchwork)
library(readxl)
library(vegan)
library(RColorBrewer)
library(taxonomizr)
library(ape)
library(ggunileg)
library(ggrepel)
library(here)

source(here('data', 'utils.r'))

pseudoCount <- 1E-5

jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return(intersection / union)
}

buildup <- map(list.files(here('data'), full.names = TRUE)[str_detect(list.files(here('data'), full.names = TRUE), 'corrected_hits.tsv')], \(p) {
    pp <- str_split(p, "/")[[1]][length(str_split(p, "/")[[1]])]
    oxygen <- str_split_fixed(pp, "_", n = 20)[, 4]
    type <- str_split_fixed(pp, "_", n = 20)[, 3]
    return(read_tsv(p) %>%
        mutate(Oxygen = oxygen, Type = type) %>%
        mutate(DataType = 'buildup'))
})
buildup <- do.call('rbind', buildup) %>%
    rename(Drug = index, ParentCompound = parent_compound) %>%
    select(-Drug) %>%
    rename(
        # StrainCommunity = Strain,
        Drug = ParentCompound,
        metabolized = hit) %>%
    group_by(Drug, Strain, Oxygen, Type, DataType) %>%
    summarize(metabolized = any(metabolized))
degrade <- map2(c(
    7,
    8,
    11,
    12),
c(
    "C_AA",
    "C_MA",
    "SS_AA",
    "SS_MA"
), function(sheetnumber, tmp) {
    # tt <- read_xlsx('/g/scb/zeller/karcher/PRISMA/data/WGS_metadata/Tables_STM.xlsx', sheet = sheetnumber, skip = 1) %>% mutate(tmp = tmp)
    tt <- read_xlsx(here('data', "Supp_Tables_STM_240209.xlsx"), sheet = sheetnumber + 2, skip = 1) %>% mutate(tmp = tmp)
    if (any(colnames(tt) == "NCBI tax ID")) {
        tt <- tt %>%
            select(-`NCBI tax ID`)
    }
    if (any(colnames(tt) == "Stool_donor")) {
        tt <- tt %>%
            rename(Strain = Stool_donor)
    }
    return(tt)
})
# degrade <- map(1:12, function(sheetnumber) read_xlsx('/g/scb/zeller/karcher/PRISMA/data/WGS_metadata/Tables_STM.xlsx', sheet = sheetnumber, skip = 1))
degrade <- do.call('rbind', degrade) %>%
    mutate(Oxygen = str_split_fixed(tmp, "_", n = 2)[, 2]) %>%
    mutate(Type = str_split_fixed(tmp, "_", n = 2)[, 1]) %>%
    mutate(DataType = 'degradation')


allData <- rbind(buildup %>%
    # filter(!str_starts(Strain, "S9")) %>%
    # filter(!str_starts(Strain, "C9")) %>%
    mutate(isControl = str_detect(Strain, "Control") | str_detect(Strain, "^C[0-9]{3}")) %>%
    group_by(isControl) %>%
    nest() %>%
    mutate(data = map2(data, isControl, \(x, is) {
        if (is) {
            x <- x %>%
                group_by(Drug, Oxygen, Type) %>%
                summarize(metabolized = any(metabolized))
            x$Strain <- "Control"
            x$DataType <- "buildup"
            return(x)
        } else {
            return(x)
        }
    })) %>%
    unnest() %>%
    # filter(!str_starts(Strain, "Control")) %>%
    # Defining betabolite-side drug metabolism if at least one metabolite is a hit (meanHits > 0)
    # mutate(fc = meanRatioCompoundIntensityStartEnd, metabolized = meanHits > 0) %>%
    # # sometimes, meanRatioCompoundIntensityStartEnd can be < 1 (due to noise).
    # mutate(fc = ifelse(fc < 1, 1, fc)) %>%
    # rename(Drug = parentComp) %>%
    rename(StrainCommunity = Strain) %>%
    select(StrainCommunity, Drug, Oxygen, Type, DataType, metabolized) %>%
    filter(!str_detect(Drug, "IS_")) %>%
    group_by(Drug, Oxygen, Type, DataType),
# mutate(rankedEffectSize = rank(fc, ties = 'average')),
degrade %>%
    filter(!str_starts(Strain, "Control")) %>%
    rename(metabolized = Metabolized) %>%
    mutate(Percentage_of_degraded_drug = as.numeric(Percentage_of_degraded_drug)) %>%
    mutate(fc = Percentage_of_degraded_drug) %>%
    rename(StrainCommunity = Strain) %>%
    select(StrainCommunity, Drug, Oxygen, Type, DataType, fc, metabolized) %>%
    group_by(Drug, Oxygen, Type, DataType)) %>%
    #    mutate(rankedEffectSize = rank(-fc, ties = 'average'))) %>%
    group_by(Oxygen, Type, Drug) %>%
    nest()

# Write binary heatmap ordered by mean BUILDUP potential of C, AA
drugs_order <- read_tsv('/g/scb/zeller/karcher/maral_find_metabolite_buildup/results/drug_buildup_C_AA_TM_LMMs_FDR_corrected..hits') %>%
    group_by(ParentCompound) %>%
    summarize(hit = sum(hit)) %>%
    arrange(desc(hit)) %>%
    pull(ParentCompound)

for (type in c("C", "SS")) {
    plotObject <- allData %>%
        filter(Type == type) %>%
        unnest() %>%
        select(Drug, Oxygen, Type, StrainCommunity, DataType, metabolized) %>%
        rename(hit = metabolized) %>%
        pivot_wider(id_cols = c(Drug, Oxygen, Type, StrainCommunity), names_from = DataType, values_from = hit, values_fill = FALSE) %>%
        filter(!str_detect(Drug, "IS_")) %>%
        mutate(Drug = as.factor(Drug)) %>%
        # bar height should be equivalent to degradation
        filter(degradation) %>%
        filter(Drug %in% drugs_order) %>%
        mutate(type = case_when(
            buildup & degradation ~ "degradation + buildup",
            !buildup & degradation ~ "degradation + no buildup")) %>%
        # mutate(type = ifelse(type, 1, 0)) %>%
        mutate(type = factor(type, levels = c(
            "degradation + no buildup",
            "degradation + buildup"
        ))) %>%
        mutate(Drug = factor(Drug, levels = drugs_order))


    if (type == "SS") {
        plotObject <- plotObject %>%
            filter(!Drug %in% c("Methotrexate", "Azathioprine")) %>%
            mutate(Drug = factor(Drug, levels = levels(Drug)[!levels(Drug) %in% c("Methotrexate", "Azathioprine")]))
    }

    print(str_c(type, " ", dim(plotObject)[1], " drugs"))

    plotObject <- plotObject %>%
        ggplot(aes(x = Drug, fill = type, color = type)) +
        geom_bar(position = 'stack') +
        facet_grid(Oxygen ~ .) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, face = 'bold'), panel.border = element_rect(colour = 'black', fill = NA, size = 1), legend.position = 'bottom', axis.text = element_text(color = 'black')) +
        # ylab("Number of drugs") +
        xlab("Drug") +
        # ggtitle("Comparison of drug buildup and degradation") +
        scale_x_discrete(drop = F) +
        # scale_fill_manual(values = c("no degradation + buildup" = "white", "degradation + buildup" = "grey", "degradation" = "black")) +
        # scale_color_manual(values = c("no degradation + buildup" = "black", "degradation + buildup" = "black", "degradation" = 'black')),
        scale_fill_manual(values = c("degradation + no buildup" = "white", "degradation + buildup" = "#3182BD")) +
        scale_color_manual(values = c("degradation + no buildup" = "black", "degradation + buildup" = "black"))

    if (type == "C") {
        plotObject <- plotObject + ylab("Number of Communities")
    } else if (type == "SS") {
        plotObject <- plotObject + ylab("Number of Strains")
    }

    ggsave(plot = plotObject,
        filename = str_c(here('results', 'plots', str_c('/drug_degradation_buildup_barplot_', type, '.pdf'))), width = 11.3, height = 7, units = 'cm')
}
