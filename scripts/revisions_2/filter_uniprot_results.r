library(tidyverse)
data <- read_tsv('/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/uniprot_data/uniprotkb_beta_glucuronidase_AND_taxono_2026_03_13.tsv')

data %>%
    filter(
        str_detect(
            `Protein names`,
            'Beta-glucuronidase'
        )
    ) %>%
    select(Entry, Reviewed, `Protein names`, `Gene Names`, Organism, Sequence) %>%
    rename(
        entry = Entry,
        reviewed = Reviewed,
        protein_names = `Protein names`,
        gene_names = `Gene Names`,
        organism = Organism,
        sequence = Sequence
    ) -> data

data_reviewed <- data %>%
    inner_join(
        data.frame(
            reviewed = "reviewed"
        )
    )

data_all <- data

data_reviewed %>%
    write_tsv(
        '/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/uniprot_data/uniprotkb_beta_glucuronidase_AND_taxono_2026_03_13_reviewed.tsv'
    )
data_all %>%
    write_tsv(
        '/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/uniprot_data/uniprotkb_beta_glucuronidase_AND_taxono_2026_03_13_all.tsv'
    )
# also write them out as fastas
write_fasta <- function(df, path) {
    fasta_lines <- c()
    for (i in seq_len(nrow(df))) {
        entry <- df$entry[i]
        protein_names <- df$protein_names[i]
        gene_names <- df$gene_names[i]
        organism <- df$organism[i]
        sequence <- df$sequence[i]
        
        header <- paste0(">", entry)
        fasta_lines <- c(fasta_lines, header, sequence)
    }
    writeLines(fasta_lines, con = path)
}

write_fasta(
    data_reviewed,
    '/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/uniprot_data/uniprotkb_beta_glucuronidase_AND_taxono_2026_03_13_reviewed.fasta'
)
write_fasta(
    data_all,
    '/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/uniprot_data/uniprotkb_beta_glucuronidase_AND_taxono_2026_03_13_all.fasta'
)   