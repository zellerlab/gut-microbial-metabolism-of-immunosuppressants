# Load library --------------------------------------------------------------------
library(tidyverse)
library(ggembl)
library(readxl)
library(here)

# Load data --------------------------------------------------------------------
df <- read_xlsx(here("data/Supp_Tables_STM_240329.xlsx"), sheet = 18, skip = 2)

# Adjust classes
str(df)
df$Side <- as.factor(df$Side)
df$Replicate <- as.factor(df$Replicate)
df$Drug <- as.factor(df$Drug)

# Define output path --------------------------------------------------------------------
output_path <- here("results/plots/")

# Define colors
colors <- c("#A6CEE3","#4682B4")

# Define custom theme
my_theme <- theme_publication()+
  theme(axis.title.x = element_text(size = 6, face = "bold"),
        axis.title.y = element_text(size = 6, face = "bold"),
        axis.text.x= element_text(size = 6, face = "plain"),
        axis.text.y= element_text(size =6, face = "plain"),
        plot.title= element_text(size =6, face = "bold"),
        legend.title = element_text(size = 6, face = "bold", hjust = 0.5),
        axis.ticks.length = unit(0.05, "cm"),    
        legend.margin = margin(0),
        legend.text = element_text(size = 6))


# Plot parent drug sterile --------------------------------------------------------------------
df_ps <- df %>% filter(Condition  == "parent_sterile")

p <- ggplot(df_ps, aes(x=Time, y=Concentration_µM, color=Side)) + 
  geom_point(aes(color=Side), size = 0.7, alpha=0.5, pch=16) +
  stat_summary(aes(y = Concentration_µM,group = Side), fun=mean, geom="line", linewidth=0.7) + 
  labs(title= "Methylprednisolone sterile", 
       x = "Time [h]", 
       y = "Concentration [µM]", 
       color = "Epithelial side")+
  scale_color_manual(values = colors, labels = c("Apical", "Basolateral"), 
                     name = "Epithelial side") + 
  scale_x_continuous(breaks = seq(0, 5, 1)) + 
  scale_y_continuous(limits = c(0, 2.5)) + 
  my_theme

# Print plot
print(p)

# Save plot
ggsave(filename = paste("Parent_sterile.pdf", sep = ""),
       plot = p,
       device = "pdf",
       path = file.path(output_path), width = 2.4, height = 1.45)

# Plot parent drug sterile --------------------------------------------------------------------
df_ms <- df %>% filter(Condition  == "metab_sterile")

p <- ggplot(df_ms, aes(x=Time, y=Concentration_µM, color=Side)) + 
  geom_point(aes(color=Side), size = 0.7, alpha=0.5, pch=16) +
  stat_summary(aes(y = Concentration_µM,group = Side), fun=mean, geom="line", linewidth=0.7) + 
  labs(title= "Hydroxymethylandrostadienedione sterile", 
       x = "Time [h]", 
       y = "Concentration [µM]", 
       color = "Epithelial side")+
  scale_color_manual(values = colors, labels = c("Apical", "Basolateral"), 
                     name = "Epithelial side") + 
  scale_x_continuous(breaks = seq(0, 5, 1)) + 
  scale_y_continuous(limits = c(0, 2.5)) + 
  theme_publication()+
  my_theme

# Print plot
print(p)

# Save plot
ggsave(filename = paste("Metab_sterile.pdf", sep = ""),
       plot = p,
       device = "pdf",
       path = file.path(output_path), width = 2.4, height = 1.45)

# Plot parent drug bacteria exposed --------------------------------------------------------------------
df_pb <- df %>% filter(Condition  == "parent_bacteria")

p <- ggplot(df_pb, aes(x=Time, y=Concentration_µM, color=Side)) + 
  geom_point(aes(color=Side), size = 0.7, alpha=0.5, pch=16) +
  stat_summary(aes(y = Concentration_µM,group = Side), fun=mean, geom="line", linewidth=0.7) + 
  labs(title= "Methylprednisolone bacteria exposed", 
       x = "Time [h]", 
       y = "Concentration [µM]", 
       color = "Epithelial side")+
  scale_color_manual(values = colors, labels = c("Apical", "Basolateral"), 
                     name = "Epithelial side") + 
  scale_x_continuous(breaks = seq(0, 5, 1)) + 
  scale_y_continuous(limits = c(0, 2.5)) + 
  my_theme

# Print plot
print(p)

# Save plot
ggsave(filename = paste("Parent_bact.pdf", sep = ""),
       plot = p,
       device = "pdf",
       path = file.path(output_path), width = 2.4, height = 1.45)

# Plot parent drug bacteria exposed --------------------------------------------------------------------
df_mb <- df %>% filter(Condition  == "metab_bacteria")

p <- ggplot(df_mb, aes(x=Time, y=Concentration_µM, color=Side)) + 
  geom_point(aes(color=Side), size = 0.7, alpha=0.5, pch=16) +
  stat_summary(aes(y = Concentration_µM,group = Side), fun=mean, geom="line", linewidth=0.7) +  
  labs(title= "Hydroxymethylandrostadienedione bacteria produced", 
       x = "Time [h]", 
       y = "Concentration [µM]", 
       color = "Epithelial side")+
  scale_color_manual(values = colors, labels = c("Apical", "Basolateral"), 
                     name = "Epithelial side") + 
  scale_x_continuous(breaks = seq(0, 5, 1)) + 
  scale_y_continuous(limits = c(0,0.6), breaks = seq(0, 0.6, 0.1) ) + 
  my_theme

# Print plot
print(p)

# Save plot
ggsave(filename = paste("Metab_bact.pdf", sep = ""),
       plot = p,
       device = "pdf",
       path = file.path(output_path), width = 2.4, height = 1.45)
