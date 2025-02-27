library(here)
library(tidyverse)
# https://git.embl.de/grp-zeller/ggembl
library(ggembl)
library(readxl)
library(pheatmap)
library(scales)
library(GGally)

source(here('data', 'utils.r'))

# combs <- c("_C_AA_TM",
#     "_C_MA_TM",
#     "_SS_AA_TM",
#     "_SS_MA_TM")

params <- list(
    c("_C_AA_TM", "C", "AA", 9 + 2),
    c("_C_MA_TM", "C", "MA", 10 + 2),
    c("_SS_AA_TM", "SS", "AA", 13 + 3),
    c("_SS_MA_TM", "SS", "MA", 14 + 2)
)

allR <- list()

# for (cc in names(allR)) {
for (param in params) {

    # Write binary heatmap ordered by mean BUILDUP potential of C, AA
    cc <- param[1]
    typeE <- param[2]
    oxygen <- param[3]
    sheetNumber <- param[4]

    tested_strains <- read_xlsx(here('data/Supp_Tables_STM_240216.xlsx'), sheet = 5, skip = 2) %>%
        filter(!is.na(`NCBI tax ID`))

    if (oxygen == "MA") {
        tested_strains <- tested_strains %>%
            filter(str_detect(`Growth condition`, "Micro"))
    } else {
        tested_strains <- tested_strains %>%
            filter(str_detect(`Growth condition`, "Anaerobic") | str_detect(`Growth condition`, "anaerobic"))
    }

    tested_strains <- tested_strains %>%
        select(Name) %>%
        # filter(Name != "Actinomyces graevenitzii") %>%
        mutate(Name = ifelse(Name == "Bifidobacterium longum subsp. Infantis", "Bifidobacterium longum", Name)) %>%
        mutate(Name = ifelse(Name == "Escherichia coli BW25113", "Escherichia coli", Name)) %>%
        left_join(read_tsv(here('data/drug_degradation_potential_SS_ordered_strains.tsv'), col_names = F) %>%
            rename(Name = X1)) %>%
        pull(Name)


    if (cc == "_C_AA_TM") {
        sheetNumber <- 11
        m <- read_xlsx(here('data/Supp_Tables_STM_240216.xlsx'), sheet = sheetNumber, skip = 2) %>%
            rename(Strain = Stool_donor,
                ParentCompound = ParentDrug,
                pval = `p-value`, hit = Detected) %>%
            mutate(Oxygen = "AA", Type = "C") %>%
            select(Strain, ParentCompound, pval, hit, Oxygen, Type, Time, clean_median_intensity, Index)
        # distinct()
    } else if (cc == "_C_MA_TM") {
        sheetNumber <- 12
        m <- read_xlsx(here('data/Supp_Tables_STM_240216.xlsx'), sheet = sheetNumber, skip = 2) %>%
            rename(Strain = Stool_donor,
                ParentCompound = ParentDrug,
                pval = `p-value`, hit = Detected) %>%
            mutate(Oxygen = "MA", Type = "C") %>%
            select(Strain, ParentCompound, pval, hit, Oxygen, Type, Time, clean_median_intensity, Index)
        # distinct()
    } else if (cc == "_SS_AA_TM") {
        sheetNumber <- 15
        m <- read_xlsx(here('data/Supp_Tables_STM_240216.xlsx'), sheet = sheetNumber, skip = 2) %>%
            rename(ParentCompound = ParentDrug,
                Strain = Strain,
                pval = `p-value`,
                hit = Detected) %>%
            mutate(Oxygen = "AA", Type = "C") %>%
            select(Strain, ParentCompound, pval, hit, Oxygen, Type, Time, clean_median_intensity, Index)
        # distinct()
    } else if (cc == "_SS_MA_TM") {
        sheetNumber <- 16
        m <- read_xlsx(here('data/Supp_Tables_STM_240216.xlsx'), sheet = sheetNumber, skip = 2) %>%
            rename(ParentCompound = ParentDrug,
                Strain = Strain,
                pval = `p-value`,
                hit = Detected) %>%
            mutate(Oxygen = "MA", Type = "C") %>%
            select(Strain, ParentCompound, pval, hit, Oxygen, Type, Time, clean_median_intensity, Index)
        # distinct()
    }

    allR[[length(allR) + 1]] <- m
    names(allR)[length(allR)] <- cc

    m <- m %>% select(-Time, clean_median_intensity, Index) %>% distinct()

    drugs_order <- read_xlsx(here('data/Supp_Tables_STM_240216.xlsx'), sheet = 11, skip = 2) %>%
        rename(Strain = Stool_donor,
            ParentCompound = ParentDrug,
            pval = `p-value`, hit = Detected) %>%
        mutate(Oxygen = "AA", Type = "C") %>%
        select(Strain, ParentCompound, pval, hit, Oxygen, Type) %>%
        distinct() %>%
        filter(hit) %>%
        group_by(ParentCompound, Strain) %>%
        summarize(hit = any(hit)) %>%
        group_by(ParentCompound) %>%
        tally() %>%
        arrange(desc(n)) %>%
        pull(ParentCompound)
    strains_order <- read_tsv(str_c(here("data", str_c("drug_degradation_potential_", typeE, "_ordered_strains.tsv"))), col_names = F)$X1
    strains_order <- c(strains_order, "Control")
    hm <- m %>%
        mutate(Strain = ifelse(Strain == "Bifidobacterium longum subsp. infantis", "Bifidobacterium longum", Strain)) %>%
        ungroup() %>%
        group_by(Strain, ParentCompound) %>%
        summarize(hit = any(hit)) %>%
        select(Strain, ParentCompound, hit) %>%
        mutate(isControl = str_detect(Strain, "Control") | str_detect(Strain, "^C[0-9]{3}")) %>%
        filter(ParentCompound %in% drugs_order) %>%
        mutate(ParentCompound = factor(ParentCompound, levels = drugs_order)) %>%
        group_by(isControl) %>%
        nest() %>%
        mutate(data = map2(data, isControl, \(x, is) {
            if (is) {
                x <- x %>%
                    group_by(ParentCompound) %>%
                    summarize(hit = any(hit))
                x$Strain <- "Control"
                return(x)
            } else {
                return(x)
            }
        })) %>%
        select(-isControl) %>%
        unnest() %>%
        mutate(hit = ifelse(hit, 1, 0)) %>%
        ungroup()

    if (str_detect(cc, "_SS_")) {
        hm <- hm %>%
            complete(Strain = tested_strains, ParentCompound)
    }
    hm <- hm %>%
        mutate(Strain = factor(Strain, levels = strains_order)) %>%
        mutate(hit = ifelse(is.na(hit), 0, hit)) %>%
        arrange(Strain, ParentCompound) %>%
        pivot_wider(id_cols = Strain, names_from = ParentCompound, values_from = hit, values_fill = 0) %>%
        as.data.frame() %>%
        column_to_rownames("Strain")

    if (str_detect(cc, "_SS_") && str_detect(cc, "_AA_")) {
        h <- 6
    } else if (str_detect(cc, "_SS_") && str_detect(cc, "_MA_")) {
        h <- 2.75
    } else {
        h <- 4
    }


    pdf(here("results/plots", str_c("Heatmap", str_replace(cc, "_TM", ""), "_buildup.pdf")),
        width = 10.7,
        height = h)
    pheatmap::pheatmap(hm, cluster_rows = FALSE, cluster_cols = FALSE, scale = 'none', color = colorRampPalette(c("lightgrey", "black"))(2), fontsize = 6.4, legend = FALSE)
    dev.off()
}

# allR <- map(combs, \(x) {
#     read_tsv(here('data', str_c("metabolite_buildup", x, "_data_parsed.tsv")))
# })
# names(allR) <- combs

allRLong <- do.call('rbind', map2(allR, names(allR), function(x, y) x %>% mutate(tmp = y) %>%
    rename(
        AreaControlsubtracted = clean_median_intensity,
        Drug = Index) %>%
    select(
        Drug, Strain, Time, AreaControlsubtracted, tmp)
)) %>%
    mutate(Oxygen = str_split_fixed(tmp, "_", n = 4)[, 3]) %>%
    mutate(Type = str_split_fixed(tmp, "_", n = 5)[, 2]) %>%
    select(-tmp)

Color_code_community <- c("black",
    "#89CFF0",
    "#0000FF",
    "#7393B3",
    "#088F8F",
    "#0096FF",
    "#5F9EA0",
    "#0047AB",
    "#6495ED",
    "#6F8FAF",
    "#5D3FD3",
    "#C19A6B",
    "#954535",
    "#D27D2D",
    "#E97451",
    "#50C878",
    "forestgreen",
    "darkseagreen1",
    "#4CBB17",
    "#90EE90")

names(Color_code_community) <- c("Control",
    "A-01",
    "A-02",
    "A-03",
    "A-04",
    "A-05",
    "A-06",
    "A-07",
    "A-08",
    "A-09",
    "A-10",
    "P-01",
    "P-02",
    "P-03",
    "P-04",
    "T-01",
    "T-02",
    "T-03",
    "T-04",
    "T-05")

trapezoidal_rule <- function(x, y) {
    sum(diff(x) * (y[-1] + y[-length(y)]) / 2)
}


for (oxygen in c("AA", "MA")) {
    for (type in c("C", "SS")) {
        plotData <- (allRLong %>%
            ungroup() %>%
            filter(Drug %in% c("M_SIROLIMUS_1", "M_EVEROLIMUS_1", "M_TACROLIMUS_8")) %>%
            # select(Drug, Strain, data, Oxygen, Type) %>%
            filter(Oxygen == oxygen, Type == type) %>%
            unnest() %>%
            # select(Drug, Strain, Oxygen, Type, Time, AreaControlsubtracted) %>%
            # filter(Time == 12) %>%
            # group_by(Drug, Strain, Oxygen, Type) %>%
            # summarize(AreaControlsubtracted = mean(AreaControlsubtracted)) %>%
            group_by(Drug, Strain, Oxygen, Time, Type) %>%
            # 1. take mean ion intensity over pools
            summarize(AreaControlsubtracted = mean(AreaControlsubtracted)) %>%
            group_by(Drug, Strain, Oxygen, Type) %>%
            # nest() %>%
            # 2. Compute area under the curve (trapezoidal rule)
            # mutate(Time = as.numeric(factor(Time, levels = sort(unique(Time))))) %>%
            summarize(AreaControlsubtracted = trapezoidal_rule(Time, AreaControlsubtracted)) %>%
            pivot_wider(id_cols = c(Strain, Oxygen, Type), names_from = Drug, values_from = AreaControlsubtracted) %>%
            ungroup() %>%
            select(-Oxygen, -Type))
        plots <- list()
        for (metabolitePair in list(
            c("M_EVEROLIMUS_1", "M_SIROLIMUS_1"),
            c("M_EVEROLIMUS_1", "M_TACROLIMUS_8"),
            c("M_SIROLIMUS_1", "M_TACROLIMUS_8")
        )) {
            p1 <- metabolitePair[1]
            p2 <- metabolitePair[2]
            p1AxisLabel <- str_replace(p1, pattern = "M_", replacement = "")
            # Make p1AxisLabel lowercase
            p1AxisLabel <- str_to_sentence(str_to_lower(p1AxisLabel))
            p1AxisLabel <- str_replace(p1AxisLabel, pattern = "_[0-9]", replacement = " M1")
            # Make p2AxisLabel lowercase
            # p2AxisLabel <- str_to_lower(p2AxisLabel)
            p2AxisLabel <- str_replace(p2, pattern = "M_", replacement = "")
            p2AxisLabel <- str_to_sentence(str_to_lower(p2AxisLabel))
            p2AxisLabel <- str_replace(p2AxisLabel, pattern = "_[0-9]", replacement = " M1")
            plotObject <- plotData %>%
                ggplot(data = .) +
                theme_classic() +
                ggtitle(str_c("Oxygen: ", oxygen, "\n", "Type: ", type)) +
                xlab(p1AxisLabel) +
                ylab(p2AxisLabel)
            # corValue <- round(cor(plotData[[p1]], plotData[[p2]], use = "complete.obs"), 2)
            # corString <- str_c("Pearson's r = ", corValue)
            plotObject <- plotObject +
                # annotate('text', x = 10, y = max(plotData[[p2]]) * 0.95, label = corString, size = 3, hjust = 0) +
                theme_publication()
            if (type == "C") {
                plotObject <- plotObject +
                    geom_point(aes_string(x = p1, y = p2, color = "Strain"), show.legend = FALSE) +
                    scale_color_manual(values = Color_code_community) +
                    NULL
            } else {
                plotObject <- plotObject +
                    geom_point(aes_string(x = p1, y = p2)) +
                    # scale_color_manual(values = Color_code_community) +
                    NULL
            }
            f <- function(x) trans_format("log10", math_format(10^.x))(x)

            # plots[[length(plots) + 1]] <- plotObject + theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), axis.text = element_text(size = 6)) + scale_x_continuous(labels = scientific, breaks = breaks_func) + scale_y_continuous(labels = scientific, breaks = breaks_func)
            plots[[length(plots) + 1]] <- plotObject + theme(panel.border = element_rect(colour = 'black', fill = NA, size = 1), axis.text = element_text(size = 6)) + scale_x_continuous(labels = scientific, breaks = pretty_breaks(n = 3), limits = c(1, max(plotData[[p1]]) * 1.05)) + scale_y_continuous(labels = scientific, breaks = pretty_breaks(n = 3), limits = c(1, max(plotData[[p2]]) * 1.05))
        }
        library(patchwork)
        # ggsave(plot = wrap_plots(plots, nrow = 1) + plot_layout(guides = "collect"), filename = str_c("/g/scb/zeller/karcher/maral_find_metabolite_buildup/other_plots/M1_SIR_EVE_TAC__", oxygen, "__", type, ".pdf"), width = 5.75, height = 2)
        ggsave(plot = wrap_plots(plots, nrow = 1) + plot_layout(guides = "collect"), filename = str_c(here('results', 'plots', str_c("M1_SIR_EVE_TAC__", oxygen, "__", type, ".pdf"))), width = 5.75, height = 2)
    }
}
