#plotting data and analysis for species-level data
library(tidyverse)
library(ggnewscale)

setwd("/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/")
df<-readRDS("Data/Trait_Data/species_level_data_final.rds")


source(file ="R/Figures/ggplot_colors.R")

df<-df %>%
  filter(!duplicated(Genus.Species)) %>%
  rowwise() %>%
  mutate(Q_index = mean(c(QWD, QWD_HL2,
                          QWD_HW2,
                          QWD_EL2,
                          QWD_ML2, na.rm =T)),
         Q_index2 = mean(c(QWD, QWD_HL,
                          QWD_HW,
                          QWD_EL,
                          QWD_ML, na.rm =T))
         )



Figure2_df<-df %>%
  dplyr::select(Genus.Species,diet, nesting_strata, morphism, status,
                QWD, QWD_EL, QWD_HL, QWD_HW, QWD_ML) %>%
  pivot_longer(cols = c(QWD, QWD_EL, QWD_HL, QWD_HW, QWD_ML), names_to = "type",
               values_to = "Average") %>%
  mutate(type = as.factor(type))

#let's change the levels to the factor "type"
levels(Figure2_df$type)<-c("Weber's length",
                           "Eye length",
                           "Head length",
                           "Head width",
                           "Mandible length")


##Figure 1
ggplot(data = df %>% drop_na(diet)) +
  geom_point(aes(x=Q_index, y = median_size), alpha = 0.8,
             pch =21,
             color = "black", fill = "darkgrey",
             size = 3) +
  scale_y_log10() +
  labs(x = "Average Queen-Worker Dimorphism (%)",
       y = "Colony size (log-scale)") +
  theme_bw() +
  theme(axis.title = element_text(size = 15))


#Figure 2
status<-ggplot(Figure2_df %>% filter(!is.na(status))) +
  geom_boxplot(aes(x = status, y = Average, fill = type,
  )) +
  scale_fill_cvi_d("my_favourite_colours") +
  scale_colour_cvi_d("my_favourite_colours") +
  labs(title = "Native vs. Nonnative",
       y = "QWD", fill = "QWD type", color = "QWD type") +
  theme_bw()+
  geom_hline(yintercept = 0, lty = 2, color = "black") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14))
status

#Figure 3 and 4 were made in the `environmental_models.R` file

setwd("/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/")

