# Load library --------------------------------------------------------------------
library(tidyverse)
library(ggembl)
library(readxl)
library(data.table)
library(pheatmap)
library(viridis)
library(here)

# Load data --------------------------------------------------------------------
C_AA_df <- read_xlsx(here("data/Supp_Tables_STM_240329.xlsx"), sheet = 8, skip = 2)
C_MA_df <- read_xlsx(here("data/Supp_Tables_STM_240329.xlsx"), sheet = 9, skip = 2)


# Define output path --------------------------------------------------------------------
output_path <- here("results/plots/")

# Drug order
drug_order <- factor(c("Mycophenolate Mofetil", "Cortisone", "Prednisone", "Hydrocortisone",
                       "Prednisolone", "Methylprednisolone", "Everolimus", "Methotrexate",
                       "Azathioprine", "Sirolimus", "Budesonid", "Dexamethasone",
                       "Betamethasone", "Tacrolimus", "Cyclosporine", "Mycophenolic Acid",
                       "Ranitidine", "Simvastatin", "Amlodipin", "Atenolol",
                       "Candesartan", "Cinacalcet", "Diltiazem", "Furosemide",
                       "Ramipril", "IS_WARFARIN", "IS_CAFFEINE", "IS_IPRIFLAVONE", "IS_LISINOPRIL", "IS_SULFAMETHOXAZOLE"))

# Plot heatmap communities-parents drugs anaerobic--------------------------------------------------------------------
# Initial data selection and cleaning
df_heatmap_C_AA <- C_AA_df %>%
  select(Drug, FC, Stool_donor) %>%
  distinct() %>%
  mutate(Stool_donor = if_else(grepl('Control_', Stool_donor), 'Control', Stool_donor)) %>%
  group_by(Drug, Stool_donor) %>%
  summarize(FC = median(FC), .groups = "drop") %>%
  mutate(Drug = factor(Drug, levels = drug_order),
         Stool_donor = fct_relevel(Stool_donor, "Control", after = Inf)) %>%
  mutate(FC = if_else(FC > 1, 1, FC),
         Perc_degr = 100 - (FC * 100)) %>%
  select(-FC)

# Prepare data for heatmap
df_heatmap_C_AA <- df_heatmap_C_AA %>%
  spread(key = Drug, value = Perc_degr) %>%
  remove_rownames() %>%
  column_to_rownames(var = "Stool_donor") %>%
  select(-matches("^IS_"))

# Define the color palette and breaks for a cohesive 0 to 100 scale
color_palette <- viridis(12)
num_colors <- length(color_palette)
scale_breaks <- seq(from = 0, to = 100, length.out = num_colors + 1)

# Generating the heatmap
pdf(str_c(output_path, "/", 'FC_Heatmap_C_AA.pdf'), width = 8, height = 3.5)
pheatmap(df_heatmap_C_AA, cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 45, fontsize = 9, color = color_palette, breaks = scale_breaks)
dev.off()

# Plot heatmap communities-parents drugs microaerobic--------------------------------------------------------------------

# Initial data selection and cleaning
df_heatmap_C_MA <- C_MA_df %>%
  select(Drug, FC, Stool_donor) %>%
  distinct() %>%
  mutate(Stool_donor = if_else(grepl('Control_', Stool_donor), 'Control', Stool_donor)) %>%
  group_by(Drug, Stool_donor) %>%
  summarize(FC = median(FC), .groups = "drop") %>%
  mutate(Drug = factor(Drug, levels = drug_order),
         Stool_donor = fct_relevel(Stool_donor, "Control", after = Inf)) %>%
  mutate(FC = if_else(FC > 1, 1, FC),
         Perc_degr = 100 - (FC * 100)) %>%
  select(-FC)

# Prepare data for heatmap
df_heatmap_C_MA <- df_heatmap_C_MA %>%
  spread(key = Drug, value = Perc_degr) %>%
  remove_rownames() %>%
  column_to_rownames(var = "Stool_donor") %>%
  select(-matches("^IS_"))

# Define the color palette and breaks for a cohesive 0 to 100 scale
color_palette <- viridis(12)
num_colors <- length(color_palette)
scale_breaks <- seq(from = 0, to = 100, length.out = num_colors + 1)

# Generating the heatmap
pdf(str_c(output_path, "/", 'FC_Heatmap_C_MA.pdf'), width = 8, height = 3.5)
pheatmap(df_heatmap_C_MA, cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 45, fontsize = 9, color = color_palette, breaks = scale_breaks)
dev.off()
