library(tidyverse)
library(readxl)
library(pheatmap)
library(data.table)
library(here)

for (id in c("SS_AA", "SS_MA")) {
    if (id == "SS_MA") {
        uniq_min_df_to_plot <- read_xlsx(here('data/Supp_Tables_STM_240216.xlsx'), sheet = 14, skip = 2) %>%
            select(-any_of(c("Time", "Pool", "clean_median_intensity"))) %>%
            distinct()
    } else {
        uniq_min_df_to_plot <- read_xlsx(here('data/Supp_Tables_STM_240216.xlsx'), sheet = 13, skip = 2) %>%
            select(-any_of(c("Time", "Pool", "clean_median_intensity"))) %>%
            distinct()
    }

    uniq_min_df_to_plot <- uniq_min_df_to_plot %>%
        select(Strain, Drug, FC) %>%
        rename(Sample = Strain) %>%
        mutate(Sample = ifelse(Sample == 'Actinomyces naeslundi', "Actinomyces naeslundii", Sample)) %>%
        identity()

    #### ADJUSTING DF FOT PLOTTING HEATMAP
    # uniq_min_df_to_plot <- uniq_min_df_to_plot %>% rename(Sample = abb_strain_name)
    # uniq_min_df_to_plot = subset(unique(uniq_min_df_to_plot))
    perc_degr_tmp <- uniq_min_df_to_plot
    uniq_min_df_to_plot$Sample[grepl('Control_', uniq_min_df_to_plot$Sample)] = 'Sterile control'
    uniq_min_df_to_plot <- uniq_min_df_to_plot %>% group_by(Drug, Sample) %>% summarize(FC = median(FC), .groups = "drop")

    uniq_min_df_to_plot$Drug <- factor(uniq_min_df_to_plot$Drug, levels = c("Azathioprine", "Mycophenolate Mofetil", "Sirolimus", "Tacrolimus",
        "Methotrexate", "Cyclosporine", "Everolimus", "Hydrocortisone", "Prednisolone", "Prednisone", "Cortisone",
        "Methylprednisolone", "Betamethasone", "Dexamethasone", "Budesonid", "Mycophenolic Acid",
        "Simvastatin", "Diltiazem", "Amlodipin", "Candesartan", "Ramipril", "Ranitidine", "Atenolol", "Cinacalcet", "Furosemide", "IS_WARFARIN", "IS_CAFFEINE", "IS_IPRIFLAVONE", "IS_LISINOPRIL", "IS_SULFAMETHOXAZOLE"))

    uniq_min_df_to_plot$Sample <- factor(uniq_min_df_to_plot$Sample, levels = c(unique(uniq_min_df_to_plot$Sample)[unique(uniq_min_df_to_plot$Sample) != "Sterile control"], "Sterile control"))

    uniq_min_df_to_plot <- as.data.table(uniq_min_df_to_plot)
    uniq_min_df_to_plot[FC > 1, FC := 1]
    uniq_min_df_to_plot[, Perc_degr := (100 - (FC * 100))]
    uniq_min_df_to_plot[, FC := NULL]

    perc_degr_tmp <- as.data.table(perc_degr_tmp)
    perc_degr_tmp[FC > 1, FC := 1]
    perc_degr_tmp[, Perc_degr := (100 - (FC * 100))]
    perc_degr_tmp[, FC := NULL]

    cast_heat_data = dcast(uniq_min_df_to_plot[, .(Sample, Drug, Perc_degr)], Sample ~ Drug, value.var = "Perc_degr")

    setDF(cast_heat_data)
    rownames(cast_heat_data) = cast_heat_data$Sample
    cast_heat_data = cast_heat_data[, -c(1)]

    # remove IS drugs
    interesting_drugs <- grep("IS_", colnames(cast_heat_data), fixed = T, invert = T, value = T)
    cast_heat_data_red <- cast_heat_data[, colnames(cast_heat_data) %in% interesting_drugs]
    rownames(cast_heat_data_red) <- str_replace_all(rownames(cast_heat_data_red), "  ", " ")
    # cast_heat_data_red <- left_join(
    if (str_detect(id, "MA")) {
        bla <- read_tsv(here('data/itol_tree_leaf_order_MA.txt'), col_names = F)
    } else if (str_detect(id, "AA")) {
        bla <- read_tsv(here('data/itol_tree_leaf_order_AA.txt'), col_names = F)
    }
    cast_heat_data_red <- inner_join(
        bla %>%
            # mutate(X1 = map_chr(X1, \(x){
            #     a <- str_split(x, " ")[[1]][1]
            #     a <- str_split(a, "")[[1]][1]
            #     a <- str_to_upper(a)
            #     b <- str_split(x, " ")[[1]][2]
            #     return(str_c(a, ". ", b, sep = ''))
            # })) %>%
            identity() %>%
            rbind(., data.frame(X1 = "Sterile control")) %>%
            rename(species = X1) %>%
            mutate(species = ifelse(species == 'B. longum', "B. longum subsp. infantis", species)),
        cast_heat_data_red %>%
            rownames_to_column('species') %>%
            # mutate(species = ifelse(species == 'A. naeslundi', "A. naeslundii", species)) %>%
            as_tibble()
    ) %>%
        column_to_rownames('species') %>%
        as.data.frame() %>%
        as.matrix()


    library(viridis)
    # pdf("/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/Heatmap_SS_MA.pdf", width = 19, height = 14)
    if (id == "SS_MA") {
        # pdf(str_c("/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/Heatmap_", id, ".pdf"), width = 19, height = 7)
        pdf(here('results/plots/', str_c('Heatmap_', id, '_degradation.pdf')), width = 19, height = 7)
    } else if (id == "SS_AA") {
        # pdf(str_c("/g/scb/zeller/karcher/PRISMA/plots/AAAKVT7HR/Heatmap_", id, ".pdf"), width = 19, height = 14)
        pdf(here('results/plots/', str_c('Heatmap_', id, '_degradation.pdf')), width = 19, height = 14)
    } else {
        asdadad
    }

    pheatmap(cast_heat_data_red, cluster_rows = F, cluster_cols = F, angle_col = 45, fontsize = 16, color = viridis(12))
    dev.off()

}
