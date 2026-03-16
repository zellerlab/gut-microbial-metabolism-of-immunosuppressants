library(tidyverse)
library(ggembl)


blast_result_paths <- list.files(
    '/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/blast_results',
    full.names = TRUE
)

motus_species_map <- read_tsv(
            '/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/45_genomes_mOTUs_Species_name.tsv'
        ) %>%
        rename(
            strain = strainID
        )

map(
    blast_result_paths,
    ~ read_tsv(.x, col_names = F, show_col_types = FALSE) %>%
        mutate(
            type = ifelse(
                str_detect(.x, "all"), "all", 'reviewed'
            ),
        )
) -> blast_results

bind_rows(
    blast_results
) %>%
    rename(
        strain_query = X1,
        uniprot_entry = X2,
        percent_identity = X3,
        alignment_length = X4,
        e_value = X11
    ) %>%
    select(
        !contains("X")
    ) %>%
    mutate(
        strain = str_split_fixed(strain_query, "___", 2)[,1],
        query_id = str_split_fixed(strain_query, "___", 2)[,2]
    ) %>%
    left_join(
        read_tsv(
            '/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/45_genome_protein_lengths.tsv',
            col_names = F
        ) %>%
            rename(
                query_id = X1,
                protein_length = X2
            ) %>%
            mutate(
                query_id = str_split_fixed(query_id, " ", 2)[,1]
            ),
        by = "query_id"
    ) %>%
    filter(e_value < 1E-5) -> data_filtered


for (pid_cutoff in c(50, 60, 70, 80, 90, 95)) {
data_filtered %>%
    arrange(desc(percent_identity)) %>%
    left_join(
        read_tsv(
            '/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/uniprot_data/uniprotkb_beta_glucuronidase_AND_taxono_2026_03_13_all.tsv'
        ) %>%
        select(
            entry, organism
        ) %>%
        rename(
            uniprot_entry = entry,
            uniprot_organism = organism
        )
    ) %>%
    mutate(
        alignment_length_rel = alignment_length / protein_length
    ) %>%
    filter(percent_identity > pid_cutoff, alignment_length_rel > 0.8) %>%
    select(
        strain, 
        query_id,
        percent_identity,
        alignment_length_rel,
        uniprot_entry, 
        uniprot_organism) %>%
    mutate(
        strain = str_replace_all(strain, "_", "")
    ) %>%
    left_join(
        motus_species_map
    ) %>%
    distinct(strain, query_id, .keep_all = TRUE) -> data_final

# Script will initially fail because it can't find strain_order - run this once before so things are consistent.
# strain_order <- data_final %>%
#     group_by(Species) %>%
#     summarise(
#         max_identity = max(percent_identity),
#         median_identity = median(percent_identity),
#         ) %>%
#     arrange(desc(median_identity)) %>%
#     pull(Species) %>%
#     c(
#         ., motus_species_map$Species[!motus_species_map$Species %in% .]
#     )


p <- ggplot() +
    geom_point(
        data = data_final,
        aes(x = Species, y = percent_identity),
        position = position_jitter(width = 0.2, height = 0),
        alpha = 0.2
    ) +
    geom_text(
        data = data_final %>% 
            group_by(Species) %>%
            summarize(total_hits = n()) %>%
            mutate(
                Species = factor(Species, levels = strain_order)
            ),
        aes(x = Species, y = 110, label = total_hits)
    ) +
    # show unseen factor levels as empty space:
    scale_x_discrete(drop = FALSE, limits = strain_order) +
    theme_presentation() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, 110) +
    NULL

ggsave(
    p,
    filename = str_c('/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/figures/revisions_2/blast_results_summary_', pid_cutoff, '.pdf'),
    width = 10,
    height = 6
)
}

