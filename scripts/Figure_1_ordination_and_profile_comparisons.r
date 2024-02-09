library(here)
library(tidyverse)
library(vegan)
library(ggsignif)
library(lmerTest)
library(lme4)
library(RColorBrewer)

source(here('data/utils.r'))

ref_profiles_wide <- read_tsv(here('data/profiles/reference_profiles_WGS_wide.tsv'))
overnight_culture_profiles_long <- read_tsv(here('data/profiles/overnight_cultures_WGS_profiles_long.tsv'))

profilesTogetherLong <- rbind(
    ref_profiles_wide %>%
        select(-c(species, kingdom, taxon)) %>%
        pivot_longer(-c(mOTU_ID, phylum, class, order, family, genus)) %>%
        rename(sampleID = name, relAb = value) %>%
        mutate(oxygen = "NA", profileType = 'reference'),
    overnight_culture_profiles_long %>%
        ungroup() %>%
        rename(mOTU_ID = mOTUs_ID) %>%
        select(stoolDonor, oxygen, phylum, class, order, family, genus, mOTU_ID, relAb) %>%
        mutate(profileType = 'overnight_culture') %>%
        rename(sampleID = stoolDonor)
)

profilesTogether <- profilesTogetherLong %>%
    mutate(log10relAb = log10(relAb + pseudoCount)) %>%
    select(sampleID, oxygen, mOTU_ID, log10relAb, profileType) %>%
    pivot_wider(id_cols = c(sampleID, oxygen, profileType), names_from = mOTU_ID, values_from = log10relAb, values_fill = log10(pseudoCount))


####################
# Ordination
####################


# This runs for a few minutes
pwdistancesBig <- vegan::vegdist(profilesTogether %>%
    mutate(profilesTogether = str_c(sampleID, as.character(oxygen), profileType, sep = "___")) %>%
    select(-c(sampleID, oxygen, profileType)) %>%
    as.data.frame() %>%
    column_to_rownames('profilesTogether'), method = 'euclidean')

pcoaBig <- cmdscale(pwdistancesBig, k = 2)
pcoaBig <- pcoaBig %>%
    as.data.frame()
colnames(pcoaBig) <- c("PCo 1", "PCo 2")
pcoaBig <- pcoaBig %>%
    rownames_to_column('tmpSampleID') %>%
    left_join(
        profilesTogether %>%
            ungroup() %>%
            mutate(tmpSampleID = str_c(sampleID, as.character(oxygen), profileType, sep = "___")) %>%
            select(sampleID, oxygen, profileType, tmpSampleID) %>%
            distinct(),
        by = "tmpSampleID"
    )

pcoaBig <- pcoaBig %>%
    mutate(sampleType = case_when(
        str_detect(sampleID, "A-") ~ "healthy adults",
        str_detect(sampleID, "T-") ~ "transplant patients",
        str_detect(sampleID, "P-") ~ "healthy children"
    ))


meta <- data.frame(raw = rownames(as.matrix(pwdistancesBig))) %>%
    mutate(profileType = str_split_fixed(raw, "___", n = 3)[, 3]) %>%
    column_to_rownames('raw')

stopifnot(all(rownames(meta) == rownames(as.matrix(pwdistancesBig))))

# Create a named vector with the colors
# ordination_color_vector <- c("healthy adults" = "#ACB78E", "healthy children" = "#F88379", "transplant patients" = "#B0C4DE")
# Update the color vector with the new blue
ordination_color_vector <- c("healthy adults" = "#19984A", "healthy children" = "#317EC2", "transplant patients" = "#E7872B")



p1Ordination <- ggplot() +
    geom_point(data = pcoaBig %>%
        filter(profileType == "reference"), aes(x = `PCo 1`, y = `PCo 2`), color = 'grey', alpha = 0.25) +
    geom_line(data = pcoaBig %>%
        filter(profileType != "reference"), aes(x = `PCo 1`, y = `PCo 2`, group = sampleID), color = 'black', alpha = 0.5) +
    geom_point(data = pcoaBig %>%
        filter(profileType != "reference"), aes(x = `PCo 1`, y = `PCo 2`, color = sampleType, shape = oxygen), alpha = 1) +
    theme_classic() +
    # scale_alpha(range  = c(0.25, 1)) +
    guides(alpha = 'none') +
    scale_color_manual(values = ordination_color_vector)

ggsave(plot = p1Ordination, filename = here('results/plots/Figure_1_Ordination.pdf'), width = 6.75, height = 4.75)


####################
# Shannon div comparisons
####################

tmp <- profilesTogetherLong %>%
    select(sampleID, oxygen, profileType, genus, relAb) %>%
    group_by(sampleID, oxygen, profileType, genus) %>%
    # mutate(relAb = 10^relAb - pseudoCount) %>%
    summarize(relAb = sum(relAb)) %>%
    group_by(sampleID, oxygen, profileType) %>%
    summarize(`Shannon` = diversity(relAb, 'shannon')) %>%
    pivot_longer(-c(sampleID, oxygen, profileType))

tmp <- tmp %>%
    rename(`Diversity index` = name) %>%
    mutate(donorType = case_when(
        str_detect(sampleID, "A-") ~ "healthy adults",
        str_detect(sampleID, "T-") ~ "transplant patients",
        str_detect(sampleID, "P-") ~ "healthy children",
        # )) %>%
        .default = "reference")) %>%
    mutate(donorType = factor(donorType, levels = c(
        "healthy adults",
        "healthy children",
        "transplant patients",
        "reference"
    )))



comparisons = list(
    c("healthy adults", "healthy children"),
    c("healthy adults", "transplant patients"),
    c("healthy children", "transplant patients"),
    c("healthy adults", "reference"),
    c("healthy children", "reference"),
    c("transplant patients", "reference")
)

diffResults <- c()
for (compPair in comparisons) {
    a <- compPair[[1]]
    b <- compPair[[2]]

    tmp2 <- tmp %>% filter(donorType %in% c(a, b))
    lmmResultsShannon <- summary(lmer(value ~ donorType + (1 | oxygen), data = tmp2 %>% filter(`Diversity index` == "Shannon")))
    # lmmResultsRichness <- summary(lmer(value ~ donorType + (1 | oxygen), data = tmp2 %>% filter(`Diversity index` == "richness")))

    diffResults[[length(diffResults) + 1]] <- lmmResultsShannon$coefficients[rownames(lmmResultsShannon$coefficients) != '(Intercept)', 'Pr(>|t|)']
    names(diffResults)[length(diffResults)] <- str_c(str_c(compPair[[1]], compPair[[2]], sep = "___"), "Shannon", sep = "___")

    # diffResults[[length(diffResults) + 1]] <- lmmResultsRichness$coefficients[rownames(lmmResultsShannon$coefficients) != '(Intercept)', 'Pr(>|t|)']
    # names(diffResults)[length(diffResults)] <- str_c(str_c(compPair[[1]], compPair[[2]], sep = "___"), "richness", sep = "___")
}

bla <- tmp %>%
    group_by(`Diversity index`) %>%
    summarize(maxVal = max(value))

maxVals <- bla$maxVal
names(maxVals) <- bla$`Diversity index`

annotations <- data.frame(diffResults, check.names = FALSE) %>%
    t() %>%
    data.frame(check.names = FALSE) %>%
    rename(values = '.') %>%
    rownames_to_column('raw') %>%
    as_tibble() %>%
    mutate(comparisons = map(raw, \(x){
        tmp <- str_split_fixed(x, "___", n = 3)[, 1:2]
        return(c(tmp[[1]], tmp[[2]]))
    })) %>%
    mutate(`Diversity index` = map_chr(raw, \(x) str_split_fixed(x, "___", n = 3)[, 3])) %>%
    mutate(comparisonsLong = map_chr(comparisons, \(x) str_c(x[[1]], x[[2]], sep = "___"))) %>%
    mutate(label = round(values, 4)) %>%
    mutate(start = map_chr(comparisons, \(x) x[[1]])) %>%
    mutate(end = map_chr(comparisons, \(x) x[[2]])) %>%
    group_by(`Diversity index`) %>%
    nest() %>%
    mutate(data = map2(data, `Diversity index`, \(x, di) {
        desiredOrder <- map_chr(comparisons, \(x) str_c(x[[1]], x[[2]], sep = "___"))
        x <- x[match(desiredOrder, x$comparisonsLong), ]
        # maxVal  <- max(x$values)
        maxVal <- maxVals[[di]]
        x$y <- map_dbl((1:dim(x)[1]) / dim(x)[1], \(x) {
            return(maxVal + maxVal * x * 0.75)
        })
        return(x)
    })) %>%
    unnest()


(
    # mutate(SampleType = ifelse(is.na(donorType), "reference", "Stool donor")) %>%
    ggplot(data = tmp) +
        facet_wrap(~`Diversity index`, scales = "free") +
        #    theme_publication() +
        # theme(
        #     axis.text.x = element_blank(),
        #     axis.ticks.x = element_blank(),
        #     axis.title.x = element_blank()) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_boxplot(aes(x = donorType, y = value)) +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
        # geom_signif(
        #     comparisons = list(
        #         c("healthy adults", "healthy children"),
        #         c("healthy adults", "transplant patients"),
        #         c("healthy children", "transplant patients"),
        #         c("healthy adults", "reference"),
        #         c("healthy children", "reference"),
        #         c("transplant patients", "reference")
        #     ),
        #     step_increase = 0.15,
        #     map_signif_level = FALSE,
        #     size = 0.2,
        #     textsize = 2.25
        # ) +
        geom_signif(data = annotations %>% select(start, end, label, y, `Diversity index`), aes(
            xmin = start,
            xmax = end,
            annotations = label,
            y_position = y), manual = TRUE, textsize = 2.25, size = 0.2) +
        ylab("Shannon diversity") +
        xlab("Sample type")) %>%
    # ggsave(plot = ., filename = "/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/diversity_measures_naive_wilcox_test.pdf", width = 4, height = 3.5)
    # ggsave(plot = ., filename = "/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/diversity_measures_lmm.pdf", width = 4, height = 3.5)
    ggsave(plot = ., filename = here('results', 'plots', 'shannon_comparison_lmm.pdf'), width = 2, height = 3.75)


#########################
# Family level barplot
#########################

get_family_level_barplot <- function(pObj, dataB, taxLevel = 'Family', levelsToShow = NULL, pc = pseudoCount) {
    # browser()
    dataB <- dataB %>%
        # mutate(relAb = (10^relAb) - pc) %>%
        mutate(taxa = .data[[taxLevel]]) %>%
        mutate(taxa = as.character(taxa)) %>%
        mutate(taxa = ifelse(taxa %in% levelsToShow, taxa, "other")) %>%
        # mutate(taxa = factor(taxa, levels = c(levelsToShow, "other", 'unclassified'))) %>%
        group_by(taxa, sampleID, oxygen) %>%
        summarize(relAb = sum(relAb))
    dataC <- dataB %>%
        ungroup() %>%
        group_by(sampleID, oxygen) %>%
        summarize(relAb = 1 - sum(relAb)) %>%
        mutate(taxa = "unclassified")
    dataB <- rbind(dataB, dataC) %>%
        mutate(taxa = factor(as.character(taxa), levels = c(levelsToShow, "other", 'unclassified')))

    dataB <- dataB %>%
        mutate(sampleID = str_c(sampleID, oxygen, sep = "__")) %>%
        mutate(oxygen = str_split_fixed(sampleID, "__", n = 2)[, 2]) %>%
        mutate(sampleID = str_split_fixed(sampleID, "__", n = 2)[, 1]) %>%
        mutate(oxygen = case_when(
            oxygen == "MA" ~ "Microaerobic",
            oxygen == "AA" ~ "Anaerobic",
        ))

    pObj <- pObj +
        geom_bar(data = dataB,
            aes(x = sampleID, y = relAb, fill = taxa), position = 'stack', stat = 'identity') +
        theme_classic() +
        facet_grid(. ~ oxygen)
    return(pObj)
}
# for (taxL in c("Phylum", "Class", "Order", "Family", "Genus")) {
for (taxL in c("phylum", "class", "order", "family", "genus")) {
    numTaxa <- 10
    levelsToShow <- profilesTogetherLong %>%
        # mutate(relAb = 10^relAb - pseudoCount) %>%
        group_by(sampleID, .data[[taxL]]) %>%
        summarize(relAb = sum(relAb)) %>%
        group_by(.data[[taxL]]) %>%
        summarize(m = mean(relAb)) %>%
        filter(!str_detect(.data[[taxL]], "incertae")) %>%
        arrange(desc(m)) %>%
        head(numTaxa) %>%
        pull(.data[[taxL]])
    getPalette <- colorRampPalette(brewer.pal(9, "Paired"))
    set.seed(1231321)
    colors <- sample(getPalette(numTaxa))
    colors <- c(colors, "#808080", "#D3D3D3")

    p2 <- get_family_level_barplot(
        ggplot(),
        profilesTogetherLong %>%
            filter(profileType == "overnight_culture") %>%
            group_by(sampleID, oxygen, profileType) %>%
            nest() %>%
            ungroup() %>%
            # if you want a bit more resolution
            # slice_sample(n = 10) %>%
            unnest(),
        taxL,
        levelsToShow = levelsToShow) +
        scale_fill_manual(values = colors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5)) +
        ggtitle("Overnight Cultures") +
        guides(fill = guide_legend(ncol = 2)) +
        theme_publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(plot = p2,
        # filename = str_c('/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/AAAKVT7HR_', taxL, '_level_barplot.pdf'), width = 6.5, height = 2)
        filename = str_c(here('results', 'plots', str_c(taxL, '_level_barplot.pdf'))), width = 6.5, height = 2)

}
