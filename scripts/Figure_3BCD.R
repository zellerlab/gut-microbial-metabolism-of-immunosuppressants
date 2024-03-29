# Load library --------------------------------------------------------------------
library(tidyverse)
library(ggembl)
library(readxl)
library(here)

# Load data --------------------------------------------------------------------
C_AA_df <- read_xlsx(here("data/Supp_Tables_STM_240329.xlsx"), sheet = 8, skip = 2)
C_AA_TM_df <- read_xlsx(here("data/Supp_Tables_STM_240329.xlsx"), sheet = 14, skip = 2)

# Define output path --------------------------------------------------------------------
output_path <- here("results/plots/")

# Define function for clean y-axis labels --------------------------------------------------------------------
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# Define custom theme -------------------------------------------------------------------
my_theme <- theme_presentation() + theme(
  axis.title.x = element_text(size = 6, face = "bold"),
  axis.title.y = element_text(size = 6, face = "bold"),
  axis.text = element_text(size = 6),  
  plot.title = element_text(size = 6, face = "bold"),
  strip.text = element_text(size = 6, face = "bold"),  
  legend.title = element_text(size = 6, face = "bold", hjust = 0.5),
  legend.text = element_text(size = 6),
  legend.key.size = unit(0.1, "cm"),
  axis.ticks.length = unit(0.05, "cm"),
  legend.margin = margin(0, 0, 0, 0),
  axis.ticks = element_line(size = 0.5),  # Set the axis ticks line size
  panel.border = element_rect(size = 0.5)  # Set the panel (frame) line size
)


# Define color scheme for communities
color_scheme <- c("Control" = "black", "A-01" = "#89CFF0", "A-02" = "#0000FF", "A-03" = "#7393B3", "A-04" = "#088F8F", "A-05" = "#0096FF", "A-06" = "#5F9EA0", "A-07"= "#0047AB", "A-08" = "#6495ED", "A-09" = "#6F8FAF","A-10"= "#5D3FD3",
                  "P-01" = "#C19A6B", "P-02"= "#954535", "P-03"="#D27D2D", "P-04"="#E97451","T-01"="#50C878", "T-02"="forestgreen", "T-03"="darkseagreen1", "T-04"="#4CBB17", "T-05"="#90EE90")

# Plot individual communities - parents drugs --------------------------------------------------------------------

# Filter data for drugs of interest
#Tacrolimus, Sirolimus, Everolimus (Figure 3B-D) 

my_df_selected <- C_AA_df %>%
  filter(Drug %in% c("Tacrolimus", "Sirolimus", "Everolimus")) %>%
  mutate(
    Drug = factor(Drug, levels = c("Tacrolimus", "Sirolimus", "Everolimus")))

my_df_selected$Stool_donor = as.character(my_df_selected$Stool_donor)
my_df_selected$Stool_donor[grep('Control',my_df_selected$Stool_donor)] = 'Control'
my_df_selected$Stool_donor = factor(my_df_selected$Stool_donor, levels = c(grep('Control',unique(my_df_selected$Stool_donor), value = T, invert = T),'Control'))


# Create the plot
p <- ggplot(my_df_selected, aes(x=Time, y=(Scaled_Area)*100, color=Stool_donor, group=Stool_donor, fill=Stool_donor)) +
  stat_summary(fun.data='mean_se', geom="smooth", se=TRUE, alpha=0.08, size=0.7) +
  facet_wrap(~Drug, nrow=1, scales="free") +
  scale_x_continuous(name="Time [h]", breaks=seq(0, 12, 2)) +
  scale_y_continuous(name="Remaining Drug [%]", limits=c(0, 150), breaks=seq(0, 150, 25)) +
  scale_colour_manual(values=color_scheme) +
  scale_fill_manual(values=color_scheme) + 
  expand_limits(y=0) +
  labs(color="Fecal \nCommunity") +
  my_theme +
  guides(color=guide_legend(ncol=2), fill="none")  

p

# Save the plot
ggsave(filename = "Tac_Eve_Sir_C_AA.pdf", plot = p, device = "pdf", width = 4.87, height = 1.39, path = output_path)


# Plot individual communities-metabolites drugs --------------------------------------------------------------------

# Filter data for metabolites of interest
#Tacrolimus M1, Sirolimus M1, Everolimus M1 (Figure 3B-D) 

my_df_selected <- C_AA_TM_df %>%
  filter(Index %in% c("M_TACROLIMUS_8", "M_SIROLIMUS_1", "M_EVEROLIMUS_1"))

my_df_selected$Stool_donor = as.character(my_df_selected$Stool_donor)
my_df_selected$Stool_donor[grep('Control',my_df_selected$Stool_donor)] = 'Control'
my_df_selected$Stool_donor = factor(my_df_selected$Stool_donor, levels = c(grep('Control',unique(my_df_selected$Stool_donor), value = T, invert = T),'Control'))

my_df_selected$Index <- as.character(my_df_selected$Index)
my_df_selected$Index[my_df_selected$Index=="M_TACROLIMUS_8"] <- "Tacrolimus M1"
my_df_selected$Index[my_df_selected$Index=="M_MMF_9"] <- "Mycophenolic Acid"
my_df_selected$Index[my_df_selected$Index=="M_SIROLIMUS_1"] <- "Sirolimus M1"
my_df_selected$Index[my_df_selected$Index=="M_EVEROLIMUS_1"] <- "Everolimus M1"
my_df_selected$Index[my_df_selected$Index=="M_METHYLPREDNISOLONE_1"] <- "Hydroxymethylandrostadienedione"

my_df_selected$Index <- factor(my_df_selected$Index, levels = c("Tacrolimus M1", "Sirolimus M1", "Everolimus M1"))

p <-ggplot(my_df_selected, aes(x= Time, y= clean_intensity, group = Stool_donor, color = Stool_donor, fill = Stool_donor), size=0.7)+
  stat_summary(fun.data='mean_se', geom="smooth", se=TRUE, alpha=0.08, size=0.7) +
  facet_wrap(~Index, nrow=1, scales="free") +
  scale_x_continuous(name =  "Time [h]", breaks = seq(0,12,2))+
  scale_y_continuous(name = "Ion intensity", label=scientific_10, n.breaks = 4) +
  scale_colour_manual(values=color_scheme) +
  scale_fill_manual(values=color_scheme) + 
  expand_limits(y=0)+
  labs(colour="Fecal \ncommunities")+
  my_theme+
  guides(color=guide_legend(ncol=2), fill="none")  

p
ggsave(filename = "M1_Tac_Eve_Sir_C_AA.pdf", plot = p, device = "pdf", width = 5.4, height = 1.39, path = output_path)


