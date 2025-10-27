library(tidyverse)
library(cowplot)

### Script to parse & visualize data on water chemistry across treatments during chemical exposure histology trials

### Starting with ammonia (measured as ammonium) levels

nh3_data <- read.csv("05082024_schumer_chemtrial_NH4.csv") %>%
  mutate(treatment = sapply(SampleID, function(x) strsplit(x, " ")[[1]][1]),
         month = sapply(SampleID, function(x) strsplit(x, " ")[[1]][2]),
         tank = sapply(SampleID, function(x) strsplit(x, " ")[[1]][3])) %>%
  filter(SampleID != "NH4L")
nh3_data$treatment <- factor(nh3_data$treatment, levels = c("Co","AC","EV")) # Co = Control, "AC" = ammonium chloride, "EV" = everything (AC + HA + CuCl2)

nh3_plot <- ggplot(nh3_data, aes(x = treatment, y = Concentration)) +
  theme_bw() +
  theme(axis.title.y = element_text(hjust = -1)) +
  geom_boxplot(fill = "gold", outlier.shape = NA) +
  labs(x = "Treatment", y = "Ammonia (mg/L N)") +
  geom_jitter(width = 0.25, color = "gray50") +
  scale_x_discrete(labels = c("Control", expression(NH[4]~Cl), expression(HA+NH[4]~Cl+CuCl[2])))
nh3_plot

### Then levels of dissolved organic carbon (DOC)

DOC_data <- read_tsv("20240417_olfact_chemtrials_DOCnewcalibration_cleaned.txt") %>%
  distinct() %>%
  rename(`NPOC` = `Mean Conc.`) %>%
  filter(Timepoint != 0)
DOC_data$Treatment <- factor(DOC_data$Treatment, levels = c("Control", "HA", "All")) # HA = Humic acid, "All" = everything

DOC_plot <- ggplot(DOC_data, aes(x = Treatment, y = NPOC)) +
  theme_bw() +
  geom_boxplot(position = "dodge", fill = "burlywood") +
  scale_x_discrete(labels = c("Control", "Humic Acid", expression(HA+NH[4]~Cl+CuCl[2]))) +
  labs(x = "Day", y = "DOC (mg/L C)") +
  geom_jitter(width = 0.25, color = "gray40")
DOC_plot

### Last, copper (Cu) concentrations measured by ICP/MS
# Unlike other 2 parameters, Cu was measured across 4 treatments
# because commercial humic acid has been shown 
# to contain substantial heavy metal content

cu_data <- read.csv("EMF_ICPMS_Results_Ben_Michael_Moran20240507.csv") %>%
  drop_na() %>%
  mutate(treatment = sapply(Label, function(x) strsplit(x, " ")[[1]][1]),
         month = sapply(Label, function(x) strsplit(x, " ")[[1]][2]),
         tank = sapply(Label, function(x) strsplit(x, " ")[[1]][3])) %>%
  filter(!grepl("Std", Label),
         !grepl("QC", Label),
         !grepl("Ref", Label),
         !grepl("Blank", Label))
cu_data$treatment = factor(cu_data$treatment, levels = c("Co","HA","Cu","Ev")) 

cu_plot <- ggplot(cu_data, aes(x = treatment, y = Value)) +
  theme_bw() +
  geom_boxplot(fill = "goldenrod", outlier.shape = NA) +
  #scale_y_log10(breaks = c(0.125, 0.5, 2.0, 10.0)) +
  labs(x = "Treatment", y = "Copper (μg/L N)") +
  geom_jitter(width = 0.25, color = "gray50") +
  scale_x_discrete(labels = c("Control", "Humic Acid", expression(CuCl[2]), expression(HA+NH[4]~Cl+CuCl[2])))
cu_plot

### Combined plot

waterQC <- plot_grid(nh3_plot, DOC_plot, cu_plot, nrow = 3, rel_heights = c(1,1,2), labels = "AUTO")
waterQC
ggsave("chemtrials_waterQC.png",
       waterQC, device = "png", width = 6.5, height = 7.5)
