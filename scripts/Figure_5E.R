# Load library --------------------------------------------------------------------
library(tidyverse)
library(ggembl)
library(readxl)
library(here)

# Load data --------------------------------------------------------------------
df <- read_xlsx(here("data/Supp_Tables_STM_240329.xlsx"), sheet = 19, skip = 2)

# Define output path --------------------------------------------------------------------
output_path <- here("results/plots/")

# Visualizing enzyme expression activity --------------------------------------------------------------------
df_EE <- df %>% filter(c(Strain != "B. uniformis RS05310"))
  
# Define color scheme
colors <- c("B. uniformis WT" = "#3690c0", 
            "B. uniformis RS05305" = "black", 
            "Control" = "grey40",
            "P2E3" = "#fcbba1", 
            "P1E4" = "#fc9272", 
            "P5E4" = "#fb6a4a",
            "P2E5" = "#ef3b2c", 
            "P4E5" = "#cb181d", 
            "P1E6" = "#99000d")

# Define custom theme
my_theme <- theme_presentation() +
  theme(
    axis.title = element_text(size = 6, face = "bold"),
    axis.text = element_text(size = 6, face = "plain"),
    strip.text = element_text(size = 6, face = "bold"),
    strip.text.y = element_text(size = 6, angle = 360, hjust = 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    plot.title = element_blank(),
    legend.key.size = unit(0.1, "cm"),
    axis.ticks.length = unit(0.05, "cm"),
    legend.margin = margin(0)
  )

# Reorder strains
df_EE$Strain <- factor(df_EE$Strain, levels = c("P2E3", "P1E4","P5E4","P2E5","P4E5","P1E6", "B. uniformis RS05305","B. uniformis WT", "Control"))

# Define function for clean y-axis labels
scientific_10 <- function(x) {
  parse(text=gsub("e", "%*% 10^", scales::scientific_format()(x)))
}

#Create plot 
p <- ggplot(df_EE, aes(x = Time, y = correctedArea, color = Strain, group = Strain,fill = Strain)) +
  #geom_line(linewidth = 0.7)+
  stat_summary(fun.data = 'mean_se', geom = "smooth", se = TRUE, alpha = 0.2, linewidth = 0.7) +
  facet_wrap(~Drug, nrow = 1, scales = "free") +
  scale_x_continuous(name = "Time [h]", breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(name = "Ion intensity", label=scientific_10, n.breaks = 6) +
  scale_colour_manual(values = c(colors))+
  scale_fill_manual(values = c(colors))+
  my_theme

# Print plot
print(p)

# Save plot
ggsave(filename = paste("Enzyme_expression.pdf", sep = ""),
       plot = p,
       device = "pdf",
       path = file.path(output_path), width = 4.2, height = 1.47)
    
# Visualizing esterase activity --------------------------------------------------------------------

df_EA <- df %>% filter(c(Strain == "B. uniformis RS05310" | Strain == "B. uniformis RS05305" |Strain == "Control" |Strain == "B. uniformis WT"))
  
# Define color scheme
colors <- c("B. uniformis WT" = "#3690c0", 
            "B. uniformis RS05305" = "#fdbb84",
            "B. uniformis RS05310" = "#e34a33",
            "Control" = "grey40")

# Reorder strains
df_EA$Strain <- factor(df_EA$Strain, levels = c("B. uniformis WT","B. uniformis RS05310","B. uniformis RS05305","Control"))

# Define function for clean y-axis labels
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

#Create plot 
p <- ggplot(df_EA, aes(x = Time, y = correctedArea, color = Strain, group = Strain, fill = Strain)) +
  stat_summary(fun.data = 'mean_se', geom = "smooth", se = TRUE, alpha = 0.2, linewidth = 0.7) +
  facet_wrap(~Drug, nrow = 1, scales = "free") +
  scale_x_continuous(name = "Time [h]", breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(name = "Ion intensity", label=scientific_10, n.breaks = 6) +
  scale_colour_manual(values = c(colors))+
  scale_fill_manual(values = c(colors))+
  my_theme

# Print plot
print(p)

# Save plot
ggsave(filename = paste("Esterase_comparison.pdf", sep = ""),
       plot = p,
       device = "pdf",
       path = file.path(output_path), width = 4.2, height = 1.47)
