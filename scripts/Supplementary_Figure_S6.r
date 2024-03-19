library(here)
library(tidyverse)
# https://git.embl.de/grp-zeller/ggembl
library(ggembl)
library(readxl)
library(pheatmap)
library(scales)
library(GGally)

source(here('data', 'utils.r'))

# _C_AA_TM
sheetNumber <- 7 + 2
community_drug_abundances_AA <- read_xlsx(here('data/Supp_Tables_STM_240216.xlsx'), sheet = sheetNumber, skip = 2) %>%
    select(Stool_donor, Drug, AUC) %>%
    mutate(Oxygen = "AA", Type = "C") %>%
    distinct()

# _C_MA_TM
sheetNumber <- 8 + 2
community_drug_abundances_MA <- read_xlsx(here('data/Supp_Tables_STM_240216.xlsx'), sheet = sheetNumber, skip = 2) %>%
    select(Stool_donor, Drug, AUC) %>%
    mutate(Oxygen = "MA", Type = "C") %>%
    distinct()

community_drug_abundances <- do.call('rbind', list(community_drug_abundances_AA, community_drug_abundances_MA))

community_drug_abundances <- community_drug_abundances %>%
    # inner_join(data.frame(Drug = c("Tacrolimus", "Mycophenolate Mofetil", "Methylprednisolone"))) %>%
    # cap AUC at 12...
    mutate(AUC = ifelse(AUC > 12, 12, AUC)) %>%
    pivot_wider(id_cols = c(Stool_donor, Oxygen), names_from = Drug, values_from = AUC) %>%
    filter(!str_detect(Stool_donor, "Control"))

plots <- list()
community_drug_abundances <- community_drug_abundances %>%
    mutate(rownames = str_c(Stool_donor, Oxygen, sep = '___')) %>%
    select(-Stool_donor, -Oxygen) %>%
    select(all_of(c("Tacrolimus", "Mycophenolate Mofetil", "Methylprednisolone", "rownames"))) %>%
    as.data.frame() %>%
    column_to_rownames("rownames")
mi <- min(community_drug_abundances)
ma <- max(community_drug_abundances)
for (i1 in 1:dim(community_drug_abundances)[2]) {
    for (i2 in 1:dim(community_drug_abundances)[2]) {
        if (i1 > i2) {
            a <- colnames(community_drug_abundances)[i1]
            b <- colnames(community_drug_abundances)[i2]
            tmp <- community_drug_abundances %>%
                select(all_of(c(a, b))) %>%
                as.data.frame()
            tmp$Oxygen <- str_split_fixed(rownames(tmp), "___", 2)[, 2]
            cors <- tmp %>%
                group_by(Oxygen) %>%
                summarize(co = cor(!!sym(a), !!sym(b), method = "spearman"))

            plot <- ggplot(tmp, aes(x = !!sym(a), y = !!sym(b))) +
                geom_point(alpha = 0.5) +
                theme_presentation() +
                geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "lightgrey") +
                theme(axis.text = element_text(size = 6), axis.title = element_text(size = 8)) +
                xlab(a) +
                ylab(b) +
                facet_grid(. ~ Oxygen) +
                # annotate('text', x = 0, y = ma * 0.95, hjust = 0, label = str_c("Pearson's r: ", co), size = 2) +
                geom_text(data = cors, aes(label = str_c("Rho: ", round(co, 3))), x = 1, y = ma * 0.95, hjust = 0, inherit.aes = FALSE, size = 2.5) +
                xlim(c(mi, ma)) +
                ylim(c(mi, ma)) +
                NULL
            plots[[length(plots) + 1]] <- plot
        }
    }
}
ggsave(plot = wrap_plots(plots, ncol = 3), filename = here('results/plots/Tac_MMF_MP_AUC_scatters.pdf'), width = 8 * 1, height = 2.75 * 0.675)
