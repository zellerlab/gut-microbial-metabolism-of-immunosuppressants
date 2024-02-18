library(tidyverse)
library(ggembl)
library(patchwork)
library(ggrepel)
library(ggplot2)
library(gggenes)
library(here)

# source('../data/utils.r')
source(here("data", 'utils.r'))

if (rstudioapi::isAvailable()) {
    print("We are within Rstudio, so mounted.")
    wdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
    print("We are on the server side, probably.")
    wdir <- getwd()
}
setwd(wdir)

# I used to blast against CL03T12C37_GCF_018292165.1_ASM1829216v1_genomic.fna, which is the wrong strain
# blast <- read_tsv("/g/scb/zeller/karcher/maral_genomic_screen_hit_ID/blast_OLD.out", col_names = F)
# In order to get locus tags that are compatible with what we show in fiogures, I'll now blast against the ATCC 8492 genome
# blastn -outfmt 6 -subject /g/scb/zeller/karcher/maral_genomic_screen_hit_ID/GCF_000154205.1_ASM15420v1_genomic.fna  -query /g/scb/zeller/karcher/maral_genomic_screen_hit_ID/all.fas > blast_NEW.out
blast <- read_tsv(here('data', 'blast.out'), col_names = F)
blast <- blast %>%
    select(X1, X2, X3, X11, X9, X10, X12) %>%
    rename(insert = X1,
        subjectContig = X2,
        meanSeqID = X3,
        start = X9,
        end = X10,
        eval = X11,
        # -outfmt "6 sstrand"...
        strandedness = X12) %>%
    filter(eval < 1E-30)
# Keep only forward reads for insertion position info
# blast <- blast %>% filter(str_detect(string = insert, "_F"))
# blast$insert <- map_chr(blast$insert, function(x) str_c(str_split(x, "_")[[1]][1:(length(str_split(x, "_")[[1]])-2)], collapse = "_"))
blast$whichRead <- map_chr(blast$insert, function(x) {
    if (str_detect(x, "_Fwd")) {
        return("Forward")
    } else {
        return('Reverse')
    }
})
blast <- blast %>%
    group_by(insert, whichRead) %>%
    do(., {
        tmp <- .
        tmp <- tmp[which(tmp$eval == min(tmp$eval)), ]
        if (dim(tmp)[1] > 1) {
            tmp <- tmp[1, ]
        }
        print(dim(tmp))
        tmp
    })
blast$replicate <- map_chr(blast$insert, function(x) {
    a <- str_count(x, pattern = "_")
    if (a == 3) {
        return(str_split(x, "_")[[1]][2])
    } else {
        return(str_c(str_split(x, "_")[[1]][3], collapse = "_"))
    }
})

blast$ID <- map_chr(blast$insert, function(x) {
    a <- str_count(x, pattern = "_")
    if (a == 3) {
        return(str_split(x, "_")[[1]][1])
    } else {
        return(str_c(str_split(x, "_")[[1]][1:2], collapse = "_"))
    }
})

blast$insertID <- map2_chr(blast$insert, blast$replicate, function(x, y) {
    a <- str_count(x, pattern = "_")
    if (a == 3) {
        return(str_c(str_split(x, "_")[[1]][1], "_", y, collapse = ""))
    } else {
        return(str_c(str_c(str_split(x, "_")[[1]][1:2], collapse = "_"), "_", y, collapse = ""))
    }
})

# control is somwhere completely different...
blast <- blast %>% filter(!str_detect(insertID, "ctrl"))

mi <- min(c(blast$start, blast$end))
ma <- max(c(blast$start, blast$end))


# gff <- read_tsv("CL03T12C37_GCF_018292165.1_ASM1829216v1_genomic_OLD.gff", col_names = F)
# gff <- read_tsv('/g/scb/zeller/karcher/maral_genomic_screen_hit_ID/GCF_000154205.1_ASM15420v1_genomic.gff', col_names = F, comment = '#')
gff <- read_tsv(here('data', 'GCF_000154205.1_ASM15420v1_genomic.gff'), col_names = F, comment = '#')
gff <- gff %>%
    select(X1, X3, X4, X5, X7, X9) %>%
    rename(subjectContig = X1,
        type = X3,
        start = X4,
        end = X5,
        strand = X7,
        info = X9) %>%
    inner_join(blast %>% ungroup() %>% select(subjectContig) %>% distinct()) %>%
    filter(start > (mi - 5000) & end < (ma + 5000)) %>%
    filter(type == "CDS") %>%
    mutate(infoParsed = map_chr(info, function(x) {
        tmp <- str_split(x, 'product=')[[1]][2]
        tmp <- str_split(tmp, ";")[[1]][1]
        return(tmp)
    })) %>%
    mutate(infoParsed2 = map_chr(info, function(x) {
        tmp <- str_split(x, 'RefSeq:')[[1]][2]
        tmp <- str_split(tmp, ";")[[1]][1]
        return(tmp)
    })) %>%
    mutate(infoParsed3 = map_chr(info, function(x) {
        tmp <- str_split(x, 'locus_tag=')[[1]][2]
        tmp <- str_split(tmp, ";")[[1]][1]
        return(tmp)
    }))
gff$xReal <- apply(gff %>% select(start, end), 1, mean)
gff$yOffset <- 1:dim(gff)[1]

library(patchwork)

p1 <- ggplot() +
    geom_gene_arrow(data = blast %>%
        filter(insertID != "ctrl") %>%
        # filter(!ID %in% c("P57")) %>%
        # filter(insertID == "P19_A" | insertID == "P14_L19_A" | insertID == "P14_L18_A"),
        filter(insertID == "P19_A" | insertID == "P14_L19_A" | insertID == "P14_L18_A"),
    aes(xmin = start, xmax = end, y = insertID, fill = ID, orientation = T, strand = 'forward')) +
    # xlim(c(4505000, 4526500)) +
    xlim(c(0, 20000)) +
    theme_embl() +
    geom_vline(aes(xintercept = gff$start), color = 'red', alpha = 0.3) +
    geom_vline(aes(xintercept = gff$end), color = 'red', alpha = 0.3) +
    geom_rect(data = gff,
        aes(xmin = start, xmax = end, ymin = 0, ymax = 3), fill = 'red', alpha = 0.15)


p2 <- ggplot() +
    geom_rect(data = gff,
        aes(xmin = start, xmax = end, ymin = yOffset, ymax = yOffset + 1), fill = "red", color = 'red', alpha = 0.3) +
    geom_text(data = gff,
        aes(x = xReal,
            label = infoParsed,
            y = yOffset + 0.5), force_pull = 50) +
    # xlim(c(4505000, 4526500)) +
    xlim(c(0, 20000)) +
    theme_embl()

p3 <- ggplot() +
    geom_rect(data = gff,
        aes(xmin = start, xmax = end, ymin = yOffset, ymax = yOffset + 1), fill = "red", color = 'red', alpha = 0.3) +
    geom_text(data = gff,
        aes(x = xReal,
            label = infoParsed2,
            y = yOffset + 0.5), force_pull = 50) +
    xlim(c(0, 20000)) +
    theme_embl()

p4 <- ggplot() +
    geom_rect(data = gff,
        aes(xmin = start, xmax = end, ymin = yOffset, ymax = yOffset + 1), fill = "red", color = 'red', alpha = 0.3) +
    geom_text(data = gff,
        aes(x = xReal,
            label = infoParsed3,
            y = yOffset + 0.5), force_pull = 50) +
    xlim(c(0, 20000)) +
    theme_embl()

# ggsave(plot = p2 / p3 / p4 / p1 + plot_layout(nrow = 4, heights = c(0.5, 0.5, 0.5, 0.1)), filename = 'gof_plot_raw.pdf', width = 8, height = 9)
ggsave(plot = p2 / p3 / p4 / p1 + plot_layout(nrow = 4, heights = c(0.5, 0.5, 0.5, 0.1)), filename = here('results', 'plots', 'gof_plot_raw.pdf'), width = 8, height = 9)
# ggplot(example_genes, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
#   geom_gene_arrow() +
#   facet_wrap(~ molecule, scales = "free", ncol = 1) +
#   scale_fill_brewer(palette = "Set3")
