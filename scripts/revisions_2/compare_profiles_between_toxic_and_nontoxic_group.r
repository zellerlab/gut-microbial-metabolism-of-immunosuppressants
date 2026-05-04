# Load libraries
library(tidyverse)
library(here)

# Load profiles
source(here('data/utils.r'))
prof <- read_tsv(here('data/MPA_Tox_PRISMA_motu_profiles.tsv'))

# Load metadata
meta <- read_csv(here('data/MPA_Tox_df_for_Nic.csv')) %>%
    select(v65_pat_id, v62_visit_number, GIT_tox, MPA_MMF_Ratio) %>%
    rename(
        PSN = v65_pat_id,
        visit = v62_visit_number,
        tox = GIT_tox,
        mpa_mmf_ratio = MPA_MMF_Ratio
    ) %>%
    mutate(tox = ifelse(tox, "toxic", 'non-toxic'))

# Subset profiles
ss <- read_tsv(here(
    "data/MPA_Tox_PRISMA_sample_list.tsv"
)) %>%
    filter(!is.na(sampleID)) %>%
    select(PSN, sampleID, visit) %>%
    distinct()

prof <- prof %>% 
    group_by(sampleID) %>% 
    inner_join(
        ss
    ) %>%
    inner_join(meta, by = c("PSN", "visit"))  %>%
    ungroup() %>%
    group_by(sampleID) %>%
    filter(sum(count) > 1000) %>%
    ungroup()

# Filter taxa and other preprocessing
prof <- prof %>%
    mutate(
        phylum = str_split_fixed(motu_raw, "[|]", n = 10)[, 2],
        order = str_split_fixed(motu_raw, "[|]", n = 10)[, 3],
        class = str_split_fixed(motu_raw, "[|]", n = 10)[, 4],
        family = str_split_fixed(motu_raw, "[|]", n = 10)[, 5],
        genus = str_split_fixed(motu_raw, "[|]", n = 10)[, 6],
        species = str_split_fixed(motu_raw, "[|]", n = 10)[, 7],
        mOTU = str_split_fixed(motu_raw, "[|]", n = 10)[, 8]
    ) %>%
    group_by(
        species, mOTU
    ) %>%
    filter(
        sum(count > 0) >= 5
    ) %>%
    group_by(
        sampleID
    ) %>%
    mutate(
        rel_ab = count/sum(count)
    )  %>%
    ungroup()

min_nonzero_relab <- min(prof %>% filter(rel_ab > 0) %>% pull(rel_ab))
prof <- prof %>%
    #select(-count) %>%
    mutate(
        rel_ab_log10 = log10(rel_ab + min_nonzero_relab)
    )   

# Ordination
pairwise_distances <- prof %>%
    select(sampleID, species, mOTU, rel_ab_log10) %>%
    pivot_wider(names_from = c(species, mOTU), values_from = rel_ab_log10) %>%
    column_to_rownames("sampleID") %>%
    as.matrix() %>%
    dist(method = "euclidean")

pcoa <- cmdscale(pairwise_distances, k = 2, eig = TRUE)$points %>%
    as.data.frame() %>%
    rownames_to_column("sampleID") %>%
    inner_join(
        prof %>% select(sampleID, tox) %>% distinct(),
        by = "sampleID"
    )

pcoa_plot <- ggplot(pcoa, aes(x = V1, y = V2, color = tox)) +
    geom_point(alpha = 0.5) +
    labs(x = "PCoA 1", y = "PCoA 2", color = "Toxicity") +
    theme_classic() +
    theme(
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
    ) 

ggsave(
    here('figures/revisions_2/pcoa_plot.pdf'),
    pcoa_plot,
    width = 4,
    height = 3
)

get_n_distinct_colors <- function(n) {
  # sample evenly around HCL color wheel
  hues  <- seq(0, 360, length.out = n + 1)[- (n + 1)]
  # fix chroma & luminance for consistency
  hcl(h = hues, c = 60, l = 70)
}

get_color_scale_per_phylum <- function(unique_phyla) {
    n <- length(unique_phyla)
    colors <- RColorBrewer::brewer.pal(n, "Set3")
    names(colors) <- unique_phyla
    return(colors)
}

#barplot_level <- 'phylum'
barplot_level <- 'genus'



# phylum barplots
phylum_data <- prof %>%
    group_by(sampleID, .data[[barplot_level]], tox, mpa_mmf_ratio) %>%
    summarise(rel_ab = sum(rel_ab)) %>%
    ungroup()

phylum_data[[barplot_level]] <- ifelse(phylum_data[[barplot_level]] == "", "unclassified", phylum_data[[barplot_level]])

# differential abundance using wilcox.test
da_results <- prof %>%
    #group_by(.data[[barplot_level]], tox) %>%
    #summarise(rel_ab = sum(rel_ab)) %>%
    #pivot_wider(names_from = tox, values_from = rel_ab) %>%
    #ungroup() %>%
    group_by(.data[[barplot_level]]) %>%
    mutate(mpa_mmf_ratio_elevated = ifelse(
        mpa_mmf_ratio > 4, "elevated", "not_elevated"
    )) %>%
    summarise(
        #p_value = wilcox.test(rel_ab[tox =='toxic'], rel_ab[tox == 'non-toxic'])$p.value,
        p_value = wilcox.test(rel_ab[mpa_mmf_ratio_elevated == 'elevated'], rel_ab[!mpa_mmf_ratio_elevated != "not_elevated"])$p.value,
    ) %>%
    arrange(p_value)

da_results$p.adj <- p.adjust(da_results$p_value, method = "BH")
da_results$p_value <- NULL

p <- ggplot(
    phylum_data %>%
        inner_join(
            da_results %>%
                filter(p.adj < 0.1) %>%
                select(.data[[barplot_level]])
        ) %>%
        mutate(mpa_mmf_ratio_elevated = ifelse(
            mpa_mmf_ratio > 4, "elevated", "not_elevated"
        )) %>%                        
        mutate(
            rel_ab_log10 = log10(rel_ab + min_nonzero_relab)
        ),
    aes(x = mpa_mmf_ratio_elevated, y = rel_ab_log10, color = tox)
) + 
    geom_point(position = position_jitter(width = 0.2, height = 0), alpha = 0.5) +
    facet_wrap(.~ .data[[barplot_level]], scales = "free_y") +
    labs(x = "mpa/mmf ratio", y = "Log10 Relative Abundance", fill = "Toxicity") +
    theme_classic() +
    theme(
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
    )

ggsave(
    here(str_c('figures/revisions_2/', barplot_level, "_differential_abundance.pdf")),
    p,
    width = 8,
    height = 4
)

common_taxa <- phylum_data %>%
    group_by(.data[[barplot_level]]) %>%
    filter(.data[[barplot_level]] != "unclassified") %>%
    #summarise(common_taxon_bool = sum(rel_ab > 0) > 10) %>%
    summarise(
        twelve_most_common_taxa = mean(rel_ab),
    ) %>%
    arrange(desc(
        twelve_most_common_taxa
    )) %>%
    head(10) %>%
    pull(.data[[barplot_level]]) 


phylum_data[[barplot_level]] <- ifelse (phylum_data[[barplot_level]] %in% common_taxa, phylum_data[[barplot_level]], "other")

color_scale <- get_color_scale_per_phylum(unique(phylum_data[[barplot_level]]))
if (
    "other" %in% names(color_scale)
) {
    color_scale["other"] <- "lightgrey"
}
color_scale["unclassified"] <- "grey"

total_depth_per_sample <- prof %>%
    group_by(sampleID, tox) %>%
    summarise(total_depth = sum(count))
#barplot_level <- 'class'

p <- ggplot(phylum_data, aes(x = sampleID, y = rel_ab, fill = .data[[barplot_level]])) +
    geom_bar(stat = "identity") +
    #facet_grid(~ tox, scales = "free_x") +
    facet_grid(.~ tox, scale="free_x", space="free") +
    labs(x = "Sample ID", y = "Relative Abundance", fill = "Phylum") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
    ) +
    scale_fill_manual(values = color_scale) +
    # color the depth labels by tox
    scale_color_manual(values = c("outlier" = "firebrick", "non-outlier" = "steelblue4")) +
    geom_text(
        data = total_depth_per_sample %>%
            inner_join(
                pcoa %>% 
                mutate(outlier_in_pcoa = case_when(
                    pcoa$V1 < 0 | pcoa$V2 > 0.2 ~ TRUE,
                    TRUE ~ FALSE
                )) %>%
                select(
                    sampleID, outlier_in_pcoa
                ) %>%
                mutate(outlier_in_pcoa = ifelse(outlier_in_pcoa, "outlier", "non-outlier"))
            ),
        aes(x = sampleID, y = 1.2, label = total_depth, color = outlier_in_pcoa),
        inherit.aes = FALSE,
        size = 3,
        angle = 90
    ) + 
    ylim(0, 1.3) +
    NULL

ggsave(
    here(str_c('figures/revisions_2/', barplot_level, "_barplpot.pdf")),
    p,
    width = 8,
    height = 4
)
