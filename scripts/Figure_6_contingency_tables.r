library(tidyverse)
library(here)
library(ggembl)

cont <- read_csv(here("data", "MMF_strain_activity_predictions_contingency_tables.tsv")) %>%
    mutate(truth = case_when(
        kind %in% c("TP", "FN") ~ "positive",
        kind %in% c("FP", "TN") ~ "negative"
    )) %>%
    mutate(predicted = case_when(
        kind %in% c("TP", "FP") ~ "positive",
        kind %in% c("FN", "TN") ~ "negative"
    )) %>%
    rename(count = outcome) %>%
    mutate(truth = factor(truth, levels = c("negative", "positive")),
        predicted = factor(predicted, levels = c("positive", "negative"))) %>%
    mutate(features = ifelse(features == "sequence-based homology", "sequence-based", "structurally refined")) %>%
    mutate(features = factor(features, levels = c('sequence-based', 'structurally refined'))) %>%
    rename(Predicted = predicted, Measured = truth)

plotObject <- ggplot(cont, aes(x = Predicted, y = Measured)) +
    geom_tile(aes(fill = kind), show.legend = F) +
    scale_fill_manual(values = c(
        "TP" = "#56B203",
        "FP" = 'white',
        "FN" = 'white',
        "TN" = '#56B203'
    )) +
    geom_text(aes(label = count), size = 9 * 0.35) +
    theme_presentation() +
    theme(
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 7)) +
    facet_wrap(. ~ features) +
    geom_hline(yintercept = 1.5, color = "black", size = 0.5) +
    geom_vline(xintercept = 1.5, color = "black", size = 0.5) +
    NULL

ggsave(here('results/plots/contingency_tables.pdf'), plotObject, width = 2.5, height = 1.25)
