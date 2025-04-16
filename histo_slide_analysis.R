library(tidyverse)
library(lme4)
library(broom)
library(boot)
library(MASS)
library(pscl)
library(car)

malcol=rgb(0/255,0/255,139/255)
hetcol=rgb(65/255,105/255,225/255)
bircol=rgb(255/255,0/255,0/255)


cilia_AGZC_CALL <- read.csv("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/histology_data_cilia_CALL_AGZC.csv")
blinding_key <- read.csv("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/Histology_blinding_key.csv") %>%
  mutate(X.Slide = Blinded.sample.ID)

cilia_lengths <- read.csv("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/AGZC_CALL_cilia_bed_lengths.csv") %>%
  dplyr::select(-Notes)

rosette_lengths <- read.csv("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/AGZC_CALL_olfactory_rosette_lengths.csv") %>%
  dplyr::select(-Notes) %>%
  group_by(Real.sample.ID) %>%
  mutate(total_rosette_length = sum(Length_um, na.rm = T))


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
# Straightforward non-parametric approach
wilcox.test(filter(summ_cilia, population == "CALL")$weighted_mean, filter(summ_cilia, population == "AGZC")$weighted_mean)

# Proportional odds logistic regression (lets us try to account for sample set)
odds <- polr(as.factor(weighted_mean) ~ population + sampleset, data = summ_cilia)
Anova(odds)

# Permutation approach
set.seed(25)
scrambled_pop = permute::shuffleSet(summ_cilia$population, nset = 100000, control = permute::how(block = summ_cilia$sampleset))
testp <- apply(scrambled_pop, 1, function(p) {
  permd = summ_cilia$population[p]
  scrambled_cilia = cbind(summ_cilia, "scrambled_pop" = permd)
  median(filter(scrambled_cilia, scrambled_pop == "AGZC")$weighted_mean, na.rm = T) - median(filter(scrambled_cilia, scrambled_pop == "CALL")$weighted_mean)
}) 
median(filter(summ_cilia, population == "AGZC")$weighted_mean) - median(filter(summ_cilia, population == "CALL")$weighted_mean)
testd <- ecdf(testp)
plot(testd)
testd(median(filter(summ_cilia, population == "AGZC")$weighted_mean) - median(filter(summ_cilia, population == "CALL")$weighted_mean))


# What's the difference? Would be that the rank-based methods above don't take into account
# the distance between variables, just whether the expected rank of one is higher than of the other.
# The permutation test takes into account the distance between the points and their summary stats, 
# and asks if we were likely to see similar differences in median by chance

# Plots

cilia_barplot <- ggplot(summ_cilia, aes(x = population, y = weighted_mean)) +
  labs(x = "Population", y = "% olfactory surface ciliated") +
  scale_x_discrete(labels = c("Upstream", "Downstream")) +
  theme_bw() +
  geom_boxplot(fill = "cadetblue")
cilia_barplot
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/cilia_present_wildfish_boxplot.png",
       cilia_barplot, device = "png", width = 2.8, height = 3)

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

### Old garbage models
#cilia_generalized_mixed <- glmer(X..Present/100 ~ population + (1 | Real.sample.ID) + Length_um + sampleset, family = "binomial", data = cilia_unblinded)
#summary(cilia_generalized_mixed)
#plot(cilia_generalized_mixed)
#confint(cilia_generalized_mixed)

#cilia_damagedBool_mixed <- glmer(damagedBool ~ population + (1 | Real.sample.ID) + Length_um + sampleset, family = "binomial", data = cilia_unblinded)
#summary(cilia_damagedBool_mixed)
#plot(cilia_damagedBool_mixed)
#confint(cilia_damagedBool_mixed)

#cilia_zeroinflated <- zeroinfl(100 - X..Present ~ Length_um + sampleset + population | population, dist = "negbin", link = "logit", data = drop_na(cilia_unblinded))
#summary(cilia_zeroinflated)
#plot(cilia_zeroinflated)



#### AB/PAS Counts of Goblet Cells ####

abpas_wild <- read.csv("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/histology_data_ABPAS_AGZC_CALL.csv") %>%
  dplyr::select(Slide:Notes)
blinding_key <- read.csv("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/Histology_blinding_key.csv") %>%
  mutate(Slide = Blinded.sample.ID)
abpas_wild_unblinded <- left_join(abpas_wild, blinding_key) %>%
  mutate(population = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][2]),
         indiv = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][3])) %>%
  left_join(filter(summ_cilia, sampleset == "M3")) %>%
  left_join(rosette_lengths)

summ_abpas_wild <- abpas_wild_unblinded %>%
  group_by(Slide, Real.sample.ID, population, indiv, total_bed_length, total_rosette_length) %>%
  dplyr::summarize(goblet_cells = mean(Goblet.Cells, na.rm = T)) %>%
  dplyr::mutate(goblets_per_um = goblet_cells/total_rosette_length,
                prop_ciliarybed = total_bed_length/total_rosette_length) %>%
  drop_na() #%>%
  filter(!Slide %in% c("M3-19", "M3-6", "M3-16"))

#Clearly not normal, running non-parametric
kruskal.test(x = summ_abpas_wild$goblets_per_um, g = summ_abpas_wild$population)
wilcox.test(filter(summ_abpas_wild, population == "AGZC")$goblet_cells, filter(summ_abpas_wild, population == "CALL")$goblet_cells)
wilcox.test(filter(summ_abpas_wild, population == "AGZC")$goblets_per_um, filter(summ_abpas_wild, population == "CALL")$goblets_per_um)

#Permutation Based P-value
set.seed(36)
scrambled_pop <- permute::shuffleSet(summ_abpas_wild$population, nset = 100000)
testp_goblets <- apply(scrambled_pop, 1, function(p) {
  permd = summ_abpas_wild$population[p]
  scrambled_goblets = cbind(summ_abpas_wild, "scrambled_pop" = permd)
  median(filter(scrambled_goblets, scrambled_pop == "AGZC")$goblets_per_um, na.rm = T) - median(filter(scrambled_goblets, scrambled_pop == "CALL")$goblets_per_um)
}) 
median(filter(summ_abpas_wild, population == "AGZC")$goblets_per_um) - median(filter(summ_abpas_wild, population == "CALL")$goblets_per_um)
testd_goblets <- ecdf(testp_goblets)
plot(testd_goblets)
testd_goblets(median(filter(summ_cilia, population == "AGZC")$weighted_mean) - median(filter(summ_cilia, population == "CALL")$weighted_mean))


ggplot(summ_abpas_wild, aes(x = population, y = goblets_per_um)) +
  theme_bw() +
  geom_jitter(width = 0.2, height = 0, color = "darkgrey") +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)))

ggplot(summ_abpas_wild, aes(x = population, y = prop_ciliarybed)) +
  theme_bw() +
  geom_jitter(width = 0.2, height = 0, color = "darkgrey") +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)))

wild_abpas_barplot <- ggplot(summ_abpas_wild, aes(x = population, y = goblets_per_um)) +
  labs(x = "Population", y = "Goblet cells per µm") +
  scale_x_discrete(labels = c("Upstream", "Downstream")) +
  theme_bw() +
  geom_boxplot(fill = "cadetblue")
wild_abpas_barplot
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/abpas_count_wildfish_boxplot.png",
       wild_abpas_barplot, device = "png", width = 2.8, height = 3)





### CASIII for apoptosis ####


cas3_wild <- read.csv("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/histology_data_CASIII_AGZC_CALL.csv") %>%
  dplyr::select(Slide:Notes)
blinding_key <- read.csv("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/Histology_blinding_key.csv") %>%
  mutate(Slide = Blinded.sample.ID)
cas3_wild_unblinded <- left_join(cas3_wild, blinding_key) %>%
  mutate(population = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][2]),
         indiv = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][3])) %>%
  left_join(cilia_lengths)


summ_cas3_wild <- cas3_wild_unblinded %>%
  group_by(Slide, Real.sample.ID, population, indiv) %>%
  dplyr::summarize(average_positive = mean(CASIII_positive, na.rm = T),
                   total_positive = sum(CASIII_positive, na.rm = T),
                   total_bed_length = sum(Length_um, na.rm = T),
                   total_bed_area = sum(Area_CASIII_um2, na.rm = T),
                   weighted_mean = weighted.mean(CASIII_positive, w = Area_CASIII_um2, na.rm = T),
                   num_beds = n()) %>%
  dplyr::mutate(positive_per_um2 = total_positive/total_bed_area) %>%
  drop_na()#%>%



t.test(filter(summ_cas3_wild, population == "AGZC")$weighted_mean, filter(summ_cas3_wild, population == "CALL")$weighted_mean)
t.test(filter(summ_cas3_wild, population == "AGZC")$positive_per_um2, filter(summ_cas3_wild, population == "CALL")$positive_per_um2)
#Clearly not normal, running non-parametric
kruskal.test(x = summ_cas3_wild$positive_per_um2, g = summ_cas3_wild$population)
wilcox.test(filter(summ_cas3_wild, population == "AGZC")$positive_per_um2, filter(summ_cas3_wild, population == "CALL")$positive_per_um2)


ggplot(summ_cas3_wild, aes(x = population, y = positive_per_um2)) +
  theme_bw() +
  geom_jitter(width = 0.2, height = 0, color = "darkgrey") +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)))

wild_cas3_barplot <- ggplot(summ_cas3_wild, aes(x = population, y = positive_per_um2)) +
  labs(x = "Population", y = "Caspase 3+ cells per µm^2") +
  scale_x_discrete(labels = c("Upstream", "Downstream")) +
  theme_bw() +
  geom_boxplot(fill = "cadetblue")
wild_cas3_barplot
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/cas3_count_wildfish_boxplot.png",
       wild_cas3_barplot, device = "png", width = 2.8, height = 3)





#########
# Lab Fish
#########


cilia_lab <- read.csv("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/histology_data_cilia_lab.csv") %>%
  dplyr::select(Slide:Notes)
blinding_key <- read.csv("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/Histology_blinding_key.csv") %>%
  mutate(Slide = Blinded.sample.ID)

cilia_unblinded <- left_join(cilia_lab, blinding_key) %>%
  mutate(month = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][1]),
         treatment = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][2]),
         sex = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][3]),
         indiv = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][4]))

summ_cilia <- cilia_unblinded %>%
  group_by(Slide, Real.sample.ID, month, treatment, sex, indiv) %>%
  dplyr::summarize(mean_present = mean(X..Present, na.rm = T),
                   mean_detached = mean(X..Detached, na.rm = T),
                   num_beds = n()) %>%
  mutate(treatment = factor(treatment, levels = c("NE","HA","AM","CU","EV")))

lab_cilia_barplot <- ggplot(summ_cilia, aes(x = treatment, y = mean_present)) +
  labs(x = "Treatment", y = "% olfactory surface ciliated") +
  scale_x_discrete(labels = c("Control", "Humic Acid", expression(NH[4]~Cl), expression(Cu~Cl[2]), expression(HA+NH[4]~Cl+CuCl[2]))) +
  theme_bw() +
  geom_boxplot(fill = "darkseagreen", outlier.shape = NA) +
  geom_jitter(width = 0.25, color = "gray50")
lab_cilia_barplot


ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/cilia_present_labfish_boxplot.png",
       lab_cilia_barplot, device = "png", width = 2.8, height = 3)

summary(glm(mean_present/100 ~ treatment, family = "quasibinomial", data = summ_cilia))

summary(aov(mean_present ~ treatment, data = summ_cilia))
plot(aov(mean_present ~ treatment, data = summ_cilia))
#Clearly not normal, running non-parametric
kruskal.test(x = summ_cilia$mean_present, g = summ_cilia$treatment)



### ABPAS


abpas_lab <- read.csv("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/histology_data_ABPAS_lab.csv") %>%
  dplyr::select(Slide:Notes)
blinding_key <- read.csv("~/Documents/schumer_lab/water_chemistry/Histology_blinding_key_casIII.csv") %>%
  mutate(Slide = Blinded.sample.ID)
abpas_lab_unblinded <- left_join(abpas_lab, blinding_key) %>%
  mutate(month = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][1]),
         treatment = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][2]),
         sex = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][3]),
         indiv = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][4]))

summ_abpas_lab <- abpas_lab_unblinded %>%
  group_by(Slide, Real.sample.ID, month, treatment, sex, indiv) %>%
  dplyr::summarize(goblet_cells = mean(Goblet_cells, na.rm = T)) %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("NE","HA","AM","CU","EV")))

summ_abpas_lab <- abpas_lab_unblinded %>%
  group_by(Slide, Real.sample.ID, month, treatment, sex, indiv) %>%
  dplyr::summarize(goblet_cells = mean(Goblet_cells, na.rm = T)) %>%
  #dplyr::mutate(goblets_per_um = goblet_cells/total_rosette_length,
  #              prop_ciliarybed = total_bed_length/total_rosette_length) %>%
  drop_na()


summary(glm(Goblet_cells ~ treatment, family = "quasipoisson", data = summ_cas3))

#Clearly not normal, running non-parametric
kruskal.test(x = summ_abpas_lab$goblet_cells, g = summ_abpas_lab$treatment)

ggplot(summ_abpas_lab, aes(x = treatment, y = goblet_cells)) +
  theme_bw() +
  geom_jitter(width = 0.2, height = 0, color = "darkgrey") +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)))

lab_abpas_barplot <- ggplot(summ_abpas_lab, aes(x = treatment, y = goblet_cells)) +
  labs(x = "Treatment", y = "Goblet cells per rosette") +
  scale_x_discrete(labels = c("Control", "Humic Acid", expression(NH[4]~Cl), expression(Cu~Cl[2]), expression(HA+NH[4]~Cl+CuCl[2]))) +
  theme_bw() +
  geom_boxplot(fill = "darkseagreen", outlier.shape = NA) +
  geom_jitter(width = 0.25, color = "gray50")
lab_abpas_barplot
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/abpas_count_labfish_boxplot.png",
       lab_abpas_barplot, device = "png", width = 2.8, height = 3)

all_lab <- plot_grid(lab_cilia_barplot, lab_abpas_barplot, nrow = 2, labels = "AUTO")
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/lab_cilia_abpas_boxplot.png", 
       all_lab, width = 6.5, height = 6)






cas3_lab <- read.csv("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/histology_data_CASIII_lab.csv") %>%
  dplyr::select(Slide:Notes)
blinding_key <- read.csv("~/Documents/schumer_lab/water_chemistry/Histology_blinding_key_casIII.csv") %>%
  mutate(Slide = Blinded.sample.ID)
cas3_unblinded <- left_join(cas3_lab, blinding_key) %>%
  mutate(month = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][1]),
         treatment = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][2]),
         sex = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][3]),
         indiv = sapply(Real.sample.ID, function(x) str_split(x, "-")[[1]][4]))

summ_cas3 <- cas3_unblinded %>%
  group_by(Slide, Real.sample.ID, month, treatment, sex, indiv) %>%
  dplyr::summarize(cas3_positive = mean(CASIII_positive, na.rm = T)) %>%
  mutate(treatment = factor(treatment, levels = c("NE","HA","AM","CU","EV")))

summ_cas3_filtered <- filter(summ_cas3, Real.sample.ID != "Jul-NE-M-1")


summary(glm(cas3_positive ~ treatment, family = "quasipoisson", data = summ_cas3))

## LAB-2023-19 is an outlier in the controls, unfortunately...
summary(aov(cas3_positive ~ treatment, data = summ_cas3))
plot(aov(cas3_positive ~ treatment, data = summ_cas3))
t.test(filter(summ_cas3, treatment == "CU")$cas3_positive, filter(summ_cas3, treatment == "NE")$cas3_positive)
t.test(filter(summ_cas3, treatment == "EV")$cas3_positive, filter(summ_cas3, treatment == "NE")$cas3_positive)
#Clearly not normal, running non-parametric
kruskal.test(x = summ_cas3$cas3_positive, g = summ_cas3$treatment)
wilcox.test(filter(summ_cas3, treatment == "CU")$cas3_positive, filter(summ_cas3, treatment == "NE")$cas3_positive)
wilcox.test(filter(summ_cas3, treatment == "EV")$cas3_positive, filter(summ_cas3, treatment == "NE")$cas3_positive)

ggplot(summ_cas3, aes(x = treatment, y = cas3_positive)) +
  theme_bw() +
  geom_jitter(width = 0.2, height = 0, color = "darkgrey") +
  stat_summary(fun = mean, geom = "pointrange", fun.max = function(x) mean(x) + 2 * sd(x) / sqrt(length(x)), fun.min = function(x) mean(x) - 2 * sd(x) / sqrt(length(x)))

cas3_barplot <- ggplot(summ_cas3, aes(x = treatment, y = cas3_positive)) +
  labs(x = "Treatment", y = "Caspase 3+ cells") +
  scale_x_discrete(labels = c("None","HA","NH3","Cu","A+C+H")) +
  theme_bw() +
  geom_boxplot(fill = "darkseagreen")
cas3_barplot
ggsave("~/Documents/Presentations/Evolution_2024_plots/cas3_count_labfish_boxplot.png",
       cas3_barplot, device = "png", width = 2.8, height = 3)





