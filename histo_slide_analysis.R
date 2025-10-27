library(tidyverse)
library(boot)
library(MASS)
library(pscl)
library(car)
library(cowplot)


### Loading data from wildcaught fish

cilia_AGZC_CALL <- read.csv("histology_data_cilia_CALL_AGZC.csv")
blinding_key <- read.csv("Histology_blinding_key_CALL_AGZC_labCilia.csv") %>%
  mutate(X.Slide = Blinded.sample.ID)

# Lengths of each sensory bed
cilia_lengths <- read.csv("AGZC_CALL_cilia_bed_lengths.csv") %>%
  dplyr::select(-Notes)

# Lengths of entire olfactory rosette, including sensory and non-sensory tissue
rosette_lengths <- read.csv("AGZC_CALL_olfactory_rosette_lengths.csv") %>%
  dplyr::select(-Notes) %>%
  group_by(Real.sample.ID) %>%
  mutate(total_rosette_length = sum(rosette_length_um, na.rm = T))

# Attaching true IDs to blinded measurements
cilia_unblinded <- left_join(cilia_AGZC_CALL, blinding_key) %>%
  mutate(sampleset = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][1]),
         population = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][2]),
         indiv = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][3]),
         damagedBool = X..Present < 100) %>%
  left_join(cilia_lengths)

summ_cilia <- cilia_unblinded %>% group_by(sampleset, population, indiv) %>% filter(!is.na(X..Present)) %>%
  dplyr::summarize(mean_present = mean(X..Present/100, na.rm = T),
                   total_bed_length = sum(Length_um, na.rm = T),
                   total_present = sum(X..Present/100 * Length_um, na.rm = T),
                   weighted_mean = weighted.mean(X..Present/100, w = Length_um, na.rm = T),
                   total_detached = sum(X..Detached/100 * Length_um, na.rm = T),
                   num_beds = n()) %>%
  arrange(population, indiv)

# Test for differences by sample set before lumping them together
wilcox.test(filter(summ_cilia, sampleset == "M3")$weighted_mean, filter(summ_cilia, sampleset == "M4")$weighted_mean)
# Non-parametric test for difference by population
wilcox.test(filter(summ_cilia, population == "CALL")$weighted_mean, filter(summ_cilia, population == "AGZC")$weighted_mean)


# Plots

cilia_barplot <- ggplot(summ_cilia, aes(x = population, y = weighted_mean)) +
  labs(x = "Population", y = "% Olfactory Surface Ciliated") +
  scale_x_discrete(labels = c("Upstream", "Downstream")) +
  theme_bw() +
  geom_boxplot(fill = "cadetblue")
cilia_barplot
ggsave("cilia_present_wildfish_boxplot.pdf",
       cilia_barplot, device = "pdf", width = 2.8, height = 3)


# Diagnostic plots and tests ruling out confounders
cilia_bed_histogram <- ggplot(cilia_unblinded, aes(x = X..Present, fill = population)) +
  theme_bw() +
  geom_histogram(position = "dodge")
cilia_bed_histogram

ggplot(summ_cilia, aes(x = total_bed_length, fill = population)) +
  theme_bw() +
  geom_histogram(position = "dodge")
ks.test(filter(summ_cilia, population == "CALL")$total_bed_length,filter(summ_cilia, population == "AGZC")$total_bed_length)

ggplot(summ_cilia, aes(x = num_beds, fill = population)) +
  theme_bw() +
  geom_histogram(position = "dodge")
ks.test(filter(summ_cilia, population == "CALL")$num_beds,filter(summ_cilia, population == "AGZC")$num_beds)

ggplot(cilia_unblinded, aes(x = Length_um, fill = population)) +
  theme_bw() +
  geom_histogram(position = "dodge")
ks.test(filter(cilia_unblinded, population == "CALL")$Length_um,filter(cilia_unblinded, population == "AGZC")$Length_um)

cilia_length_plot <- ggplot(cilia_unblinded, aes(x = Length_um, y = X..Present, color = population)) +
  facet_grid(~sampleset) +
  theme_bw() +
  geom_jitter()
cilia_length_plot

ks.test(filter(summ_cilia, sampleset == "M3")$weighted_mean,filter(summ_cilia, sampleset == "M4")$weighted_mean)


#### AB/PAS Counts of Goblet Cells ####

abpas_wild <- read.csv("histology_data_ABPAS_AGZC_CALL.csv") %>%
  dplyr::select(Slide:Notes)
blinding_key <- read.csv("Histology_blinding_key_CALL_AGZC_labCilia.csv") %>%
  mutate(Slide = Blinded.sample.ID)
abpas_wild_unblinded <- left_join(abpas_wild, blinding_key) %>%
  mutate(population = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][2]),
         indiv = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][3])) %>%
  left_join(filter(summ_cilia, sampleset == "M3")) %>%
  left_join(rosette_lengths)

summ_abpas_wild <- abpas_wild_unblinded %>%
  group_by(Slide, Real.sample.ID, population, indiv) %>%
  dplyr::summarize(sum_goblet_cells = sum(Goblet.Cells, na.rm = T),
                   total_observed_length_um = sum(rosette_length_um, na.rm = T),
                   goblets_per_um = sum_goblet_cells/total_observed_length_um) 

# Not normal, running non-parametric test for difference between populations
wilcox.test(filter(summ_abpas_wild, population == "AGZC")$sum_goblet_cells, filter(summ_abpas_wild, population == "CALL")$sum_goblet_cells)
wilcox.test(filter(summ_abpas_wild, population == "AGZC")$goblets_per_um, filter(summ_abpas_wild, population == "CALL")$goblets_per_um)

wild_abpas_barplot <- ggplot(summ_abpas_wild, aes(x = population, y = goblets_per_um)) +
  labs(x = "Population", y = "Goblet Cells per µm") +
  scale_x_discrete(labels = c("Upstream", "Downstream")) +
  theme_bw() +
  geom_boxplot(fill = "cadetblue")
wild_abpas_barplot
ggsave("abpas_count_wildfish_boxplot.pdf",
       wild_abpas_barplot, device = "pdf", width = 2.8, height = 3)




#########
# Lab fish from chemical exposure trial
#########


# Ciliation in lab treated fish

cilia_lab <- read.csv("histology_data_cilia_lab.csv") %>%
  dplyr::select(Slide:Notes)
lab_blinding_key <- read.csv("Histology_blinding_key_CALL_AGZC_labCilia.csv") %>%
  mutate(Slide = Blinded.sample.ID)

# Lengths of each individual sensory bed in the olfactory rosette
lab_cilia_lengths <- read.csv("lab_cilia_bed_lengths.csv") %>%
  dplyr::select(-Notes)

lab_cilia_unblinded <- left_join(cilia_lab, lab_blinding_key) %>%
  left_join(lab_cilia_lengths, by =) %>%
  mutate(month = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][1]),
         treatment = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][2]),
         sex = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][3]),
         indiv = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][4]))

lab_summ_cilia <- lab_cilia_unblinded %>% 
  group_by(Slide, Real.sample.ID, month, treatment, sex, indiv) %>%
  dplyr::summarize(mean_present = mean(X..Present/100, na.rm = T),
                   total_bed_length = sum(Length_um, na.rm = T),
                   total_present = sum(X..Present/100 * Length_um, na.rm = T),
                   weighted_mean = weighted.mean(X..Present/100, w = Length_um, na.rm = T),
                   total_detached = sum(X..Detached/100 * Length_um, na.rm = T),
                   num_beds = n()) %>%
  mutate(treatment = factor(treatment, levels = c("NE","HA","AM","CU","EV")))

lab_cilia_barplot <- ggplot(lab_summ_cilia, aes(x = treatment, y = weighted_mean)) +
  labs(x = "Treatment", y = "% Olfactory Surface Ciliated") +
  scale_x_discrete(labels = c("Control", "Humic Acid", expression(NH[4]~Cl), expression(Cu~Cl[2]), expression(HA+NH[4]~Cl+CuCl[2]))) +
  #scale_y_continuous(breaks = c(.25, .5, .75, 1), limits = c(.1, 1)) + # To make comparable to scale in main text figure
  theme_bw() +
  geom_boxplot(fill = "darkseagreen", outlier.shape = NA) +
  geom_jitter(width = 0.25, color = "gray50")
lab_cilia_barplot

ggsave("cilia_present_labfish_boxplot.png",
       lab_cilia_barplot, device = "png", width = 2.8, height = 3)


summary(aov(weighted_mean ~ treatment, data = lab_summ_cilia))
plot(aov(weighted_mean ~ treatment, data = lab_summ_cilia))
# Not normal, running non-parametric (shows no significant difference)
kruskal.test(x = lab_summ_cilia$weighted_mean, g = lab_summ_cilia$treatment)



### ABPAS in lab treatments


abpas_lab <- read.csv("histology_data_ABPAS_lab.csv") %>%
  dplyr::select(Slide:Notes)
blinding_key <- read.csv("Histology_blinding_key_lab_ABPAS.csv") %>%
  mutate(Slide = Blinded.sample.ID)
# Lengths of entire olfactory rosette, including sensory and non-sensory tissue
lab_rosette_lengths <- read.csv("lab_olfactory_rosette_lengths.csv") %>%
  dplyr::select(-Notes)


abpas_lab_unblinded <- left_join(abpas_lab, blinding_key) %>%
  left_join(lab_rosette_lengths) %>%
  mutate(month = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][1]),
         treatment = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][2]),
         sex = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][3]),
         indiv = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][4]))


summ_abpas_lab <- abpas_lab_unblinded %>%
  group_by(Slide, Real.sample.ID, month, treatment, sex, indiv) %>%
  dplyr::summarize(sum_goblet_cells = sum(Goblet_cells, na.rm = T),
                   total_observed_length_um = sum(rosette_length_um, na.rm = T),
                   goblets_per_um = sum_goblet_cells/total_observed_length_um) %>%
  mutate(treatment = factor(treatment, levels = c("NE","HA","AM","CU","EV")))

# Running non-parametric test for difference in goblet density between treatments
kruskal.test(x = summ_abpas_lab$goblets_per_um, g = summ_abpas_lab$treatment)


lab_abpas_barplot <- ggplot(summ_abpas_lab, aes(x = treatment, y = goblets_per_um)) +
  labs(x = "Treatment", y = "Goblet Cells per µm") +
  #lims(y = c(0,.045)) + # To make comparable to scale in main text figure
  scale_x_discrete(labels = c("Control", "Humic Acid", expression(NH[4]~Cl), expression(Cu~Cl[2]), expression(HA+NH[4]~Cl+CuCl[2]))) +
  theme_bw() +
  geom_boxplot(fill = "darkseagreen", outlier.shape = NA) +
  geom_jitter(width = 0.25, color = "gray50")
lab_abpas_barplot
ggsave("abpas_count_labfish_boxplot.png",
       lab_abpas_barplot, device = "png", width = 2.8, height = 3)

# Assembling Fig. S15
all_lab <- plot_grid(lab_cilia_barplot, lab_abpas_barplot, nrow = 2, labels = "AUTO")
ggsave("lab_cilia_abpas_boxplot.png", 
       all_lab, width = 6.5, height = 6)








