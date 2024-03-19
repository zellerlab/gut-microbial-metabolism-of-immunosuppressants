library(here)

pseudoCount <- 1E-5
veryLargeGeneModel <- TRUE

get_first_two_fields <- function(string) {
    fields <- strsplit(string, " ")[[1]]
    return(paste(fields[1:2], collapse = " "))
}
get_first_two_fields_but_truncate_first <- function(string) {
    fields <- strsplit(string, " ")[[1]]
    fields[1] <- str_c(str_split(fields[1], "")[[1]][1], ".", collapse = "")
    return(paste(fields[1:2], collapse = " "))
}

get_first_four_fields_but_truncate_first <- function(string) {
    fields <- strsplit(string, " ")[[1]]
    fields[2] <- "sp."
    return(paste(fields[1:2], collapse = " "))
}


aa_ma_colors <- c("Anaerobic" = "#66B3FF", "Microaerobic" = "#99CC00")
activity_colors <- c("inactive" = "#171f96", "active" = "#e0c941", 'untested' = "grey")

# load(url("https://github.com/AlessioMilanese/motus_taxonomy/blob/master/data/motus_taxonomy_3.0.1.Rdata?raw=true"))
# motus3.0_taxonomy <- read.table('/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/motus3.taxonomy.tsv', sep = "\t")
motus3.0_taxonomy <- read.table(here("data/motus3.taxonomy.tsv"), sep = "\t", )
motus3.0_taxonomy$mOTUs_ID <- str_replace(motus3.0_taxonomy$mOTUs_ID, "_v3_", "_v31_")
motus3.0_taxonomy <- motus3.0_taxonomy %>%
    mutate(Species = sub("^\\S+\\s*", "", Species))

# diamond_with_characterized_enzymes <- read_tsv('/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/all_genes_against_GMGC10.human-gut.95nr.protein.DIAMOND.out', col_names = F)
diamond_with_characterized_enzymes <- read_tsv(here('data/all_genes_against_GMGC10.human-gut.95nr.protein.DIAMOND.out'), col_names = F)
colnames(diamond_with_characterized_enzymes) <- vector <- c("qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "sstart", "send", "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", "gapopen", "gaps", "btop", "cigar", "stitle", "qcovhsp", "scovhsp", "qtitle", "qqual")

# These are putative homologues found using characterized enzymes
diamond_with_characterized_enzymes <- diamond_with_characterized_enzymes %>%
    select(qseqid, qlen, length, sseqid, slen, evalue, bitscore, pident) %>%
    # I think length == length of alignment
    mutate(fractionQueryAligned = length / qlen) %>%
    filter(fractionQueryAligned > 0.7 & pident > 40) %>%
    mutate(func = map_chr(qseqid, getGeneTypeMapping)) %>%
    # This (for now) removes the 2/4 F. prausnitzii genes that are not on Maral's slides for some reason.
    filter(!is.na(func))


offst <- 2
community_drug_degradation_AA <- read_xlsx(here('data/Supp_Tables_STM_240216.xlsx'), sheet = 7 + offst, skip = 2) %>%
    # select(Stool_donor, Drug, Percentage_of_degraded_drug) %>%
    select(Stool_donor, Drug, AUC, Metabolized) %>%
    distinct() %>%
    rename(effect_size = AUC) %>%
    mutate(effect_size = ifelse(effect_size > 12, 12, effect_size)) %>%
    mutate(oxygen = "AA")
community_drug_degradation_MA <- read_xlsx(here('data/Supp_Tables_STM_240216.xlsx'), sheet = 8 + offst, skip = 2) %>%
    # select(Stool_donor, Drug, Percentage_of_degraded_drug) %>%
    select(Stool_donor, Drug, AUC, Metabolized) %>%
    distinct() %>%
    rename(effect_size = AUC) %>%
    mutate(effect_size = ifelse(effect_size > 12, 12, effect_size)) %>%
    mutate(oxygen = "MA")
strain_drug_degradation_AA <- read_xlsx(here('data/Supp_Tables_STM_240216.xlsx'), sheet = 11 + offst, skip = 2) %>%
    # select(Stool_donor, Drug, Percentage_of_degraded_drug) %>%
    select(Strain, Drug, AUC, Metabolized) %>%
    distinct() %>%
    rename(effect_size = AUC) %>%
    mutate(effect_size = ifelse(effect_size > 12, 12, effect_size)) %>%
    mutate(oxygen = "AA")
strain_drug_degradation_MA <- read_xlsx(here('data/Supp_Tables_STM_240216.xlsx'), sheet = 12 + offst, skip = 2) %>%
    # select(Stool_donor, Drug, Percentage_of_degraded_drug) %>%
    select(Strain, Drug, AUC, Metabolized) %>%
    distinct() %>%
    rename(effect_size = AUC) %>%
    mutate(effect_size = ifelse(effect_size > 12, 12, effect_size)) %>%
    mutate(oxygen = "MA")

dataTypes <- list()
dataTypes[['drug_degradation']] <- rbind(community_drug_degradation_AA, community_drug_degradation_MA)
# dataTypes[['metabolite_buildup']] <- rbind(community_metabolite_buildup_AA, community_metabolite_buildup_MA)
dataTypesStrains <- list()
dataTypesStrains[['drug_degradation']] <- rbind(
    strain_drug_degradation_AA %>%
        mutate(Strain = ifelse(Strain == "Bifidobacterium longum subsp. infantis", "Bifidobacterium longum", Strain)) %>%
        mutate(Strain = ifelse(Strain == "Actinomyces graevenitzi", "Actinomyces graevenitzii", Strain)) %>%
        mutate(Strain = ifelse(Strain == "Bacteroides xylanisolvens", "Bacteroides ovatus", Strain)) %>%
        mutate(Strain = ifelse(Strain == "Bacteroides ovatus", "Bacteroides ovatus/Bacteroides xylanisolvens", Strain)),
    strain_drug_degradation_MA %>%
        mutate(Strain = ifelse(Strain == "Bifidobacterium longum subsp. infantis", "Bifidobacterium longum", Strain)) %>%
        mutate(Strain = ifelse(Strain == "Actinomyces graevenitzi", "Actinomyces graevenitzii", Strain)) %>%
        mutate(Strain = ifelse(Strain == "Bacteroides xylanisolvens", "Bacteroides ovatus", Strain)) %>%
        mutate(Strain = ifelse(Strain == "Bacteroides ovatus", "Bacteroides ovatus/Bacteroides xylanisolvens", Strain))
) %>%
    filter(!str_detect(Strain, "Control")) %>%
    # Merge B. ovatus and B. xylanisolvens and get the mean of their effect sizes
    group_by(Strain, Drug, oxygen) %>%
    summarize(effect_size = mean(effect_size), Metabolized = any(Metabolized))

if (!loadFromScratch) {
    # gmgc_profile_raw <- read_tsv('/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/collated.gene_counts.combined_scaled.txt')
    gmgc_profile_raw <- read_tsv(here('data/collated.gene_counts.combined_scaled.txt'))
} else {
    print("Loading data from scratch...")
    gmgc_profile_raw <- read_tsv(str_c(scratchpath, 'collated.gene_counts.combined_scaled.txt'))
}

gmgc_profile_raw[is.na(gmgc_profile_raw)] <- 0

gmgc_profile_all_genes <- gmgc_profile_raw %>%
    as.data.frame()

## Both apply and prop.table are very slow, presumably because we have such insane amounts of rows
## Solution: Calculate colSums, subset, calcualte proportions (of remaining rows)
# gmgc_profile <- as.data.frame(apply(gmgc_profile, 2, \(x) x / sum(x)))
# gmgc_profile <- prop.table(gmgc_profile)
## FOr some reason this is just much faster than colSums
colSums <- c()
for (colIndex in 2:ncol(gmgc_profile_all_genes)) {
    colSums <- c(colSums, sum(gmgc_profile_all_genes[[colIndex]]))
}

# In this case, apply seems 'fast' enough (fast is relative, still takes a few minutes)
meanAbundance <- apply(gmgc_profile_all_genes[, -1], 1, \(x) mean(x))
prevalence <- apply(gmgc_profile_all_genes[, -1], 1, \(x) sum(x > 0))

keepIndex <- prevalence > 35 & meanAbundance > 10

gmgc_profile_homologues_to_characterized_enzymes <- gmgc_profile_all_genes %>%
    inner_join(diamond_with_characterized_enzymes %>% rename(gene = sseqid) %>% select(gene) %>% distinct()) %>%
    column_to_rownames('gene')
for (colIndex in 1:ncol(gmgc_profile_homologues_to_characterized_enzymes)) {
    gmgc_profile_homologues_to_characterized_enzymes[[colIndex]] <- gmgc_profile_homologues_to_characterized_enzymes[[colIndex]] / colSums[colIndex]
}


gmgc_profile_all_genes <- gmgc_profile_all_genes %>%
    column_to_rownames('gene')
for (colIndex in 1:ncol(gmgc_profile_all_genes)) {
    gmgc_profile_all_genes[[colIndex]] <- gmgc_profile_all_genes[[colIndex]] / colSums[colIndex]
}

# This controls whether to keep all genes for the full model or not
if (!veryLargeGeneModel) {
    gmgc_profile_all_genes <- gmgc_profile_all_genes[keepIndex, ]
}


# Make sure all gmgc genes are part of the profile
# Turns out this is not true since we searched against the non-rare GMGC.
# TODO: Repeat search against the non-rare GMGC
# stopifnot(all((diamond_with_characterized_enzymes %>% rename(gene = sseqid) %>% select(gene) %>% distinct()  %>% pull(gene)) %in% rownames(gmgc_profile)))
print(str_c("Number of hits in gmgc: ", diamond_with_characterized_enzymes %>% rename(gene = sseqid) %>% select(gene) %>% distinct() %>% pull(gene) %>% length()))
print(str_c("Number of those I can find in profiles: ", dim(gmgc_profile_homologues_to_characterized_enzymes)[1]))

sampleID_donorID_map <- read.delim(text = "sampleID	ExpID	stoolDonor
MB001	1	A-01
MB002	2	A-02
MB003	3	A-03
MB005	4	A-04
MB010	5	A-05
MB006	6	A-06
MB007	7	A-07
MB008	8	A-08
MB009	9	A-09
MB021	21	A-10
MB015	15	P-01
MB016	16	P-02
MB017	17	P-03
MB018	18	P-04
MB013	13	T-01
MB014	14	T-02
MB019	19	T-03
MB020	20	T-04
MB011	22	T-05", sep = "\t") %>%
    as_tibble()


gmgc_profile_gmgc_profile_homologues_to_characterized_enzymes_long <- gmgc_profile_homologues_to_characterized_enzymes %>%
    rownames_to_column('gene') %>%
    pivot_longer(-gene, names_to = 'sampleIDGenecore', values_to = 'abundance') %>%
    mutate(sampleIDGenecore = str_replace(sampleIDGenecore, ".*lane1", "")) %>%
    # left_join(read_tsv('/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/MB002A_Samples_for_WGS.tsv',
    left_join(read_tsv(here('data/MB002A_Samples_for_WGS.tsv'),
        col_types = readr::cols(`Sample ID (Genecore)` = col_character())) %>%
        select(Samples, `Sample ID (Genecore)`, `Anaerob (AA)/Microaerob (MA)`) %>%
        rename(
            sampleIDGenecore = `Sample ID (Genecore)`,
            sampleID = Samples,
            oxygen = `Anaerob (AA)/Microaerob (MA)`
        ), by = 'sampleIDGenecore') %>%
    anti_join(data.frame(sampleID = c("MB006S", "MB010 160321"))) %>%
    mutate(sampleID = ifelse(sampleID == "MB006A", "MB006", sampleID)) %>%
    mutate(sampleID = ifelse(sampleID == "MB010 290720", "MB010", sampleID)) %>%
    # filter(oxygen == "AA") %>%
    select(-sampleIDGenecore) %>%
    inner_join(sampleID_donorID_map, by = 'sampleID')

gmgc_profile_all_genes_wide <- gmgc_profile_all_genes
colnames(gmgc_profile_all_genes_wide) <- str_replace(colnames(gmgc_profile_all_genes_wide), ".*lane1", "")
#tmp <- read_tsv('/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/MB002A_Samples_for_WGS.tsv',
tmp <- read_tsv(here('data/MB002A_Samples_for_WGS.tsv'),
    col_types = readr::cols(`Sample ID (Genecore)` = col_character())) %>%
    select(Samples, `Sample ID (Genecore)`, `Anaerob (AA)/Microaerob (MA)`) %>%
    rename(
        sampleIDGenecore = `Sample ID (Genecore)`,
        sampleID = Samples,
        oxygen = `Anaerob (AA)/Microaerob (MA)`
    ) %>%
    .[match(colnames(gmgc_profile_all_genes_wide), .$sampleIDGenecore), ]
stopifnot(all(colnames(gmgc_profile_all_genes_wide) == tmp$sampleIDGenecore))
colnames(gmgc_profile_all_genes_wide) <- str_c(tmp$sampleID, tmp$oxygen, sep = "___")
gmgc_profile_all_genes_wide <- gmgc_profile_all_genes_wide[, !colnames(gmgc_profile_all_genes_wide) %in% c("MB006S___AA", "MB010 160321___MA")]
colnames(gmgc_profile_all_genes_wide)[colnames(gmgc_profile_all_genes_wide) == "MB006A___AA"] <- "MB006___AA"
colnames(gmgc_profile_all_genes_wide)[colnames(gmgc_profile_all_genes_wide) == "MB010 290720___MA"] <- "MB010___MA"
colnames(gmgc_profile_all_genes_wide)[colnames(gmgc_profile_all_genes_wide) == "MB010 290720___AA"] <- "MB010___AA"
for (i in 1:dim(sampleID_donorID_map)[1]) {
    sampleID <- sampleID_donorID_map$sampleID[i]
    stoolDonor <- sampleID_donorID_map$stoolDonor[i]
    w <- map_chr(colnames(gmgc_profile_all_genes_wide), \(x) str_split(x, "___")[[1]][1]) == sampleID
    stopifnot(sum(w) == 2)
    o <- map_chr(colnames(gmgc_profile_all_genes_wide), \(x) str_split(x, "___")[[1]][2])[w]
    stopifnot(length(unique(o)) == 2)
    colnames(gmgc_profile_all_genes_wide)[w] <- str_c(stoolDonor, o, sep = "___")
}
gmgc_profile_all_genes_wide <- as.data.frame(t(as.data.frame(gmgc_profile_all_genes_wide)))
gmgc_profile_all_genes_wide$stoolDonor <- str_split_fixed(rownames(gmgc_profile_all_genes_wide), "___", 2)[, 1]
gmgc_profile_all_genes_wide$oxygen <- str_split_fixed(rownames(gmgc_profile_all_genes_wide), "___", 2)[, 2]

colnames(gmgc_profile_all_genes_wide) <- str_replace_all(colnames(gmgc_profile_all_genes_wide), "-", "_")
colnames(gmgc_profile_all_genes_wide) <- str_replace_all(colnames(gmgc_profile_all_genes_wide), "\\(", "_")
colnames(gmgc_profile_all_genes_wide) <- str_replace_all(colnames(gmgc_profile_all_genes_wide), "\\)", "_")


# WGS taxonomic profiles
tax_profiles_long <- rbind(
    # read_xlsx('/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/Supp_Tables_STM_240216.xlsx', sheet = 18 + offst, skip = 2) %>%
    read_xlsx(here('data/Supp_Tables_STM_240216.xlsx'), sheet = 18 + offst, skip = 2) %>%
        pivot_longer(-c(phylum, class, order, family, genus, species, mOTUs_ID), names_to = "Stool_donor", values_to = "abundance") %>%
        mutate(oxygen = "AA"),
    # read_xlsx('/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/Supp_Tables_STM_240216.xlsx', sheet = 19 + offst, skip = 2) %>%
    read_xlsx(here('data/Supp_Tables_STM_240216.xlsx'), sheet = 19 + offst, skip = 2) %>%
        pivot_longer(-c(phylum, class, order, family, genus, species, mOTUs_ID), names_to = "Stool_donor", values_to = "abundance") %>%
        mutate(oxygen = "MA")
)

# species_motus_link_motus_extender <- read_tsv('/g/scb/zeller/karcher/maral_dig_into_metabolizing_enzymes/45_genomes_mOTUs_Species_name.tsv') %>%
species_motus_link_motus_extender <- read_tsv(here('data/45_genomes_mOTUs_Species_name.tsv')) %>%
    group_by(mOTUs_ID) %>%
    summarize(Species = str_c(Species, collapse = "/")) %>%
    filter(Species != "Actinomyces graevenitzii")


tax_profiles_long_only_motus_of_interest <- tax_profiles_long %>%
    # inner_join(species_motus_link %>% select(mOTUs_ID), by = "mOTUs_ID")
    inner_join(species_motus_link_motus_extender %>% select(mOTUs_ID), by = "mOTUs_ID")


gene_target_map <- c(
    "WP_005929331.1" = "Tacrolimus",
    "WP_005934826.1" = "Tacrolimus",
    "WP_005929331.1" = "Everolimus",
    "WP_005934826.1" = "Everolimus",
    "WP_005929331.1" = "Sirolimus",
    "WP_005934826.1" = "Sirolimus",
    "bu0920" = "Mycophenolate Mofetil",
    "bu0921" = "Mycophenolate Mofetil",
    "CLOSCI_00899" = "Betamethasone",
    "CLOSCI_00899" = "Budesonid",
    "CLOSCI_00899" = "Cortisone",
    "CLOSCI_00899" = "Dexamethasone",
    "CLOSCI_00899" = "Hydrocortisone",
    "CLOSCI_00899" = "Methylprednisolone",
    "CLOSCI_00899" = "Prednisolone",
    "CLOSCI_00899" = "Prednisone",
    "CLOSCI_00900" = "Betamethasone",
    "CLOSCI_00900" = "Budesonid",
    "CLOSCI_00900" = "Cortisone",
    "CLOSCI_00900" = "Dexamethasone",
    "CLOSCI_00900" = "Hydrocortisone",
    "CLOSCI_00900" = "Methylprednisolone",
    "CLOSCI_00900" = "Prednisolone",
    "CLOSCI_00900" = "Prednisone"
)

gene_target_map <- data.frame(gene = names(gene_target_map), target = gene_target_map) %>%
    as_tibble() %>%
    rename(predictor = gene, target = target) %>%
    # For genes, we need to extend the gene_target map to the GMGC gene IDs
    left_join(diamond_with_characterized_enzymes %>%
        select(qseqid, sseqid) %>%
        rename(baitGene = qseqid, gene = sseqid), by = c('predictor' = 'baitGene')) %>%
    select(gene, target) %>%
    distinct() %>%
    rename(predictor = gene)

# As species_target_map, use this
species_target_map <- dataTypesStrains[['drug_degradation']] %>%
    group_by(Strain, Drug) %>%
    summarize(Metabolized = any(Metabolized)) %>%
    filter(Metabolized) %>%
    select(-Metabolized) %>%
    rename(predictor = Strain, target = Drug) %>%
    inner_join(species_motus_link_motus_extender, by = c('predictor' = 'Species')) %>%
    mutate(predictor = mOTUs_ID) %>%
    select(-mOTUs_ID)

effect_size_measure_map <- list('drug_degradation' = "Area under Curve", 'metabolite_buildup' = "Median intensity of metabolite")

# allTargets <- unique(gene_target_map$target)
allTargets <- unique(dataTypesStrains$drug_degradation$Drug)
allTargets <- allTargets[!str_detect(allTargets, "IS_")]

# drugsToComputeMetricsFor <- c("Tacrolimus", "Everolimus", "Sirolimus", "Methylprednisolone", "Mycophenolate Mofetil")
drugsToComputeMetricsFor <- dataTypes[['drug_degradation']] %>% pull(Drug) %>% unique() %>% .[!str_detect(., "IS_")]
DrugsToHighlightInMainFigure <- c("Tacrolimus", "Everolimus", "Sirolimus")
DrugsToShowInMainFigure <- c("Tacrolimus", "Everolimus", "Sirolimus", "Methylprednisolone", "Mycophenolate Mofetil")
