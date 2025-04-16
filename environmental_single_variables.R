library(tidyverse)
library(scales)
library(cowplot)
library(ggstance)
library(car)
library(basemaps)
library(ggmap)
library(sf)
library(sp)
library(gstat)
library(vegan)
set_null_device("cairo")

##### Plots ######
malcol=rgb(0/255,0/255,139/255)
hetcol=rgb(65/255,105/255,225/255)
bircol=rgb(255/255,0/255,0/255)
viridis_scale = "rocket"
hex <- hue_pal()(4)
names(hex) = c("Calnali", "Conzintla", "Huazalingo", "Pochula")

river_colors <- c('#8B1E3F','#89BD9E','#EBB65C','#DB4C40')
names(river_colors) <- c("Calnali","Conzintla", "Huazalingo","Pochula")

is_all_numeric <- function(x) {
  !any(is.na(suppressWarnings(as.numeric(na.omit(x))))) & is.character(x)
}

site_distances <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/site_riverkm_distances.csv") %>%
  dplyr::select(-Drainage)

water_chemistry <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/waterchem_with_averagedsonde_labmetals_data.csv") %>%
  filter(Site.Code != "N/A", !Date %in% c("15-Jun-22","16-Jun-22")) %>%
  mutate(day = sapply(Date, function(x) strsplit(x, split = "-")[[1]][1]),
         month = sapply(Date, function(x) strsplit(x, split = "-")[[1]][2]),
         year = sapply(Date, function(x) strsplit(x, split = "-")[[1]][3]),
         dist_from_hybridstart = ifelse(Drainage == "Calnali", Lon - -98.585, ifelse(Drainage == "Conzintla", Lon - -98.62412, ifelse(Drainage == "Pochula", Lon - -98.57946, ifelse(Drainage == "Huazalingo", Lon - -98.63713, NA)))),
         Turbidity = ifelse(Turbidity_instrument == "orion", Turbidity, NA),
         Drainage = ifelse(Site.Code %in% c("PEZM", "CAPS", "PLAZ","CALL", "TLCD"), "Calnali - Downstream", ifelse(Drainage == "Calnali", "Calnali - Upstream", Drainage))) %>%
  filter(year %in% c(22,23,24)) %>%
  filter(!(month == "Sep")) %>%
  filter(Drainage %in% c("Calnali - Upstream", "Calnali - Downstream", "Pochula", "Huazalingo", "Conzintla")) %>%
  left_join(site_distances)
water_chemistry <- as.data.frame(apply(water_chemistry, 2,
                                       FUN = function(x) ifelse(grepl("<",x), as.numeric(str_remove(x, "<"))/2, ifelse(x == "", NA, x)))) %>%
  mutate_if(is_all_numeric,as.numeric)



###### Individual Variable Plots ######

kruskal.test(Cu_t ~ Drainage, water_chemistry)

kruskal.test(Pb_t ~ Drainage, water_chemistry)

kruskal.test(Cd_t ~ Drainage, water_chemistry)

kruskal.test(Al_t ~ Drainage, water_chemistry)

kruskal.test(Zn_t ~ Drainage, water_chemistry)

kruskal.test(Ni_t ~ Drainage, water_chemistry)


ammonia_model <- aov(log(NH3 + .001) ~ Drainage * month, data = water_chemistry)
plot(ammonia_model)
shapiro.test(ammonia_model$residuals)
Anova(ammonia_model)
TukeyHSD(ammonia_model, which = "Drainage")

no2_model <- aov(log(NO2 + .001) ~ Drainage * month, data = water_chemistry)
plot(no2_model)
shapiro.test(no2_model$residuals)
Anova(no2_model)
TukeyHSD(no2_model, which = "Drainage")

fDOM_model <- aov(log(fDOM_QSU) ~ Drainage * month, data = water_chemistry)
plot(fDOM_model)
shapiro.test(fDOM_model$residuals)
Anova(fDOM_model)
TukeyHSD(fDOM_model, which = "Drainage")

turb_model <- aov(log(Turbidity) ~ Drainage * month, data = water_chemistry)
plot(turb_model)
shapiro.test(turb_model$residuals)
Anova(turb_model)
TukeyHSD(turb_model, which = "Drainage")

kruskal.test(Mn_t ~ Drainage, water_chemistry)


kruskal.test(As_t ~ Drainage, water_chemistry)


kruskal.test(Ca_t ~ Drainage, water_chemistry)

cond_model <- aov(log(Conductivity) ~ Drainage * month, data = water_chemistry)
plot(cond_model)
shapiro.test(cond_model$residuals)
Anova(cond_model)
TukeyHSD(cond_model)

water_chemistry_nosplit <- water_chemistry %>%
  mutate(Drainage = ifelse(Drainage == "Calnali - Downstream", "Calnali", ifelse(Drainage == "Calnali - Upstream", "Calnali", Drainage)))

DOMplot <- ggplot(water_chemistry_nosplit,
                  aes(x = dist_from_first_site_km, y = fDOM_QSU)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.position = "inside",
        legend.position.inside = c(.2, .8)) +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Dissolved Organic Matter (QSU)") +
  geom_smooth(aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
DOMplot
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/fDOM_riverkm_multiriver.pdf",
       DOMplot, device = "pdf", width = 7, height = 5.5)

nh3plot <- ggplot(filter(water_chemistry_nosplit, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla"), !(month == "Jun")),
                  aes(x = dist_from_first_site_km, y = NH3 + .001)) +
  coord_cartesian(x = c(0, 20), y = c(0.001, 30)) +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 24),
        legend.position = "inside",
        legend.text = element_text(size = 20),
        legend.position.inside = c(.2, .8)) +
  labs(x = "River Distance (km)", y = "Ammonia (mg/L N)") +
  #lims(y = c(0,5)) +
  scale_color_manual(values = river_colors) +
  scale_y_log10(labels = label_number(accuracy = .01)) +
  #annotate("text", x = 12, y = 5, label = "Outlier at 31 mg/L N", size = 6) +
  geom_smooth(aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
nh3plot
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/nh3_riverkm_multiriver.pdf",
       nh3plot, device = "pdf", width = 7, height = 5.5)

water_chem_censored_raw <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/waterchem_with_averagedsonde_labmetals_data.csv") %>%
  filter(Site.Code != "N/A", !Date %in% c("15-Jun-22","16-Jun-22")) %>%
  mutate(day = sapply(Date, function(x) strsplit(x, split = "-")[[1]][1]),
         month = sapply(Date, function(x) strsplit(x, split = "-")[[1]][2]),
         year = sapply(Date, function(x) strsplit(x, split = "-")[[1]][3]),
         dist_from_hybridstart = ifelse(Drainage == "Calnali", Lon - -98.585, ifelse(Drainage == "Conzintla", Lon - -98.62412, ifelse(Drainage == "Pochula", Lon - -98.57946, ifelse(Drainage == "Huazalingo", Lon - -98.63713, NA)))),
         Turbidity = ifelse(Turbidity_instrument == "orion", Turbidity, NA)) %>%
  filter(year %in% c(22,23,24)) %>%
  filter(!(month == "Sep")) %>%
  filter(Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")) %>%
  left_join(site_distances)
water_chem_censored_values <- as.data.frame(apply(water_chem_censored_raw, 2,
                                       FUN = function(x) ifelse(grepl("<",x), as.numeric(str_remove(x, "<"))/2, ifelse(x == "", NA, x)))) %>%
  mutate_if(is_all_numeric,as.numeric)
water_chem_censored_values <- as.data.frame(apply(water_chem_censored_values, 2,
                                       FUN = function(x) {
                                         num_cen <- sum(grepl("<",x))
                                         DL <- as.numeric(str_remove(x, "<"))
                                         ifelse(grepl("<",x), num_cen, ifelse(x == "", NA, x))
                                         })) %>%
  mutate_if(is_all_numeric,as.numeric)

Cu_censored <- water_chemistry_nosplit %>%
  filter(Cu_t == 0.00025) %>%
  group_by(round(dist_from_first_site_km, 1)) %>%
  mutate(num_cen = n(),
         id = row_number(round(dist_from_first_site_km, 1)),
         dummy_val = .0005/2 - (num_cen/2 - row_number(dist_from_first_site_km) + 0.5) * .0005/8)

Cuplot <- ggplot(filter(water_chemistry_nosplit, Cu_t != .00025),
                 aes(x = dist_from_first_site_km, y = Cu_t * 1000)) +
  coord_cartesian(x = c(0, 20), y = c(0, 2.5)) +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.position = "inside",
        legend.position.inside = c(.2, .8)) +
  scale_color_manual(values = river_colors, drop = F) +
  labs(x = "River Distance (km)", y = "Total Copper (ug/L)") +
  geom_smooth(data = water_chemistry_nosplit, aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River", override.aes = list(shape = 16))) +
  geom_hline(yintercept = 0.5, lty = "dashed", color = "darkgray") +
  geom_point(data = Cu_censored, aes(x = dist_from_first_site_km, y = dummy_val * 1000, color = Drainage), shape = 16, size = 2) +
  annotate("rect", xmin = -5, xmax = 25, ymin = 0, ymax = .0005 * 1000, alpha = .5,fill = "white") +
  geom_jitter(aes(color = Drainage), size = 2)
  Cuplot
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/Cu_riverkm_multiriver.pdf",
       Cuplot, device = pdf, width = 7, height = 5.5)
#### Supplemental Figures #####

turbplot <- ggplot(water_chemistry_nosplit,
                  aes(x = dist_from_first_site_km, y = Turbidity)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(
        legend.position = "inside",
        legend.position.inside = c(.2, .8)) +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Turbidity (NTU)") +
  geom_smooth(aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
turbplot

no2plot <- ggplot(water_chemistry_nosplit,
                   aes(x = dist_from_first_site_km, y = NO2)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(.2, .8)) +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Nitrite (mg/L N)") +
  geom_smooth(aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
no2plot

calnali_diff_nonmetals <- plot_grid(turbplot, no2plot, nrow = 2, labels = "AUTO")
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/Turb_NO2_riverkm_multiriver.pdf",
       calnali_diff_nonmetals, device = "pdf", width =6.5, height = 8)


condplot <- ggplot(water_chemistry_nosplit,
                   aes(x = dist_from_first_site_km, y = Conductivity)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
    legend.position = "none",
    legend.position.inside = c(.9, .2)) +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Conductivity\n(μS/cm)") +
  geom_smooth(aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
condplot

alkplot <- ggplot(water_chemistry_nosplit,
                  aes(x = dist_from_first_site_km, y = Alkalinity)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
    legend.position = "none",
    legend.position.inside = c(.9, .2)) +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = expression(atop("Alkalinity","(mg"~CaCO[3]~"Eq.)"))) +
  geom_smooth(aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
alkplot

hardplot <- ggplot(water_chemistry_nosplit,
                  aes(x = dist_from_first_site_km, y = Total_Hardness)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
    legend.position = "none",
    legend.position.inside = c(.9, .2)) +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = expression(atop("Hardness","(mg"~CaCO[3]~"Eq.)"))) +
  geom_smooth(aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
hardplot

pHplot <- ggplot(water_chemistry_nosplit,
                 aes(x = dist_from_first_site_km, y = pH)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.position.inside = c(.9, .2)) +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "pH") +
  geom_smooth(aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage))
pHplot

nonmetals_other <- plot_grid(condplot, alkplot, hardplot, pHplot, nrow = 4, labels = "AUTO", rel_heights = c(1,1,1,1.5), align = "v")
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/cond_alk_hard_pH_riverkm_multiriver.png",
       nonmetals_other, device = "png", width =6.5, height = 8)



Fe_censored <- water_chemistry_nosplit %>%
  filter(Fe_t == 0.005) %>%
  group_by(round(dist_from_first_site_km, 1)) %>%
  mutate(num_cen = n(),
         id = row_number(round(dist_from_first_site_km, 1)),
         dummy_val = 0.01/2 - (num_cen/2 - row_number(dist_from_first_site_km) + 0.5) * 0.01/8)



Feplot <- ggplot(filter(water_chemistry_nosplit, Fe_t != 0.005),
                 aes(x = dist_from_first_site_km, y = Fe_t * 1000)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(.25, .72)) +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Iron (μg/L)") +
  geom_smooth(data = water_chemistry_nosplit, aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_hline(yintercept = 0.01 * 1000, lty = "dashed", color = "darkgray") +
  geom_point(data = Fe_censored, aes(x = dist_from_first_site_km, y = dummy_val * 1000, color = Drainage), shape = 16, size = 2) +
  annotate("rect", xmin = -5, xmax = 25, ymin = 0, ymax = .01 * 1000, alpha = .5,fill = "white") +
  geom_jitter(aes(color = Drainage)) #+
Feplot


# No censored values for Aluminum
Alplot <- ggplot(filter(water_chemistry_nosplit),
                 aes(x = dist_from_first_site_km, y = Al_t * 1000)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.position.inside = c(.2, .8)) +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Aluminum (μg/L)") +
  geom_smooth(aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
Alplot

Pb_censored <- water_chemistry_nosplit %>%
  filter(Pb_t == 0.00005/2) %>%
  group_by(round(dist_from_first_site_km, 1)) %>%
  mutate(num_cen = n(),
         id = row_number(round(dist_from_first_site_km, 1)),
         dummy_val = 0.00005/2 - (num_cen/2 - row_number(dist_from_first_site_km) + 0.5) * 0.00005/8)

Pbplot <- ggplot(filter(water_chemistry_nosplit, Pb_t != 0.00005/2),
                 aes(x = dist_from_first_site_km, y = Pb_t * 1000)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.position.inside = c(.2, .8)) +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Lead (μg/L)") +
  geom_smooth(data = water_chemistry_nosplit, aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_point(data = Pb_censored, aes(x = dist_from_first_site_km, y = dummy_val * 1000, color = Drainage), shape = 16, size = 2) +
  geom_hline(yintercept = 0.00005 * 1000, lty = "dashed", color = "darkgray") +
  annotate("rect", xmin = -5, xmax = 25, ymin = 0, ymax = 0.00005 * 1000, alpha = .5,fill = "white") +
  geom_jitter(aes(color = Drainage)) #+
Pbplot

Ti_censored <- water_chem_censored_raw %>%
  filter(grepl('<', Ti_t)) %>%
  group_by(round(dist_from_first_site_km, 1)) %>%
  mutate(num_cen = n(),
         id = row_number(round(dist_from_first_site_km, 1)),
         dummy_val = as.numeric(str_remove(Ti_t, "<"))/2 - (num_cen/2 - row_number(dist_from_first_site_km) + 0.5) * as.numeric(str_remove(Ti_t, "<"))/8)

Tiplot <- ggplot(filter(water_chem_censored_raw, !grepl("<", Ti_t)),
                 aes(x = dist_from_first_site_km, y = as.numeric(Ti_t) * 1000)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.position.inside = c(.2, .8)) +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Titanium (μg/L)") +
  geom_smooth(data = water_chemistry_nosplit, aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_point(data = Ti_censored, aes(x = dist_from_first_site_km, y = dummy_val * 1000, color = Drainage), shape = 16, size = 2) +
  geom_hline(yintercept = 0.0003 * 1000, lty = "dashed", color = "darkgray") +
  annotate("rect", xmin = -5, xmax = 25, ymin = 0, ymax = 0.0003 * 1000, alpha = .5,fill = "white") +
  geom_jitter(aes(color = Drainage)) #+
Tiplot

Zn_censored <- water_chem_censored_raw %>%
  filter(grepl('<', Zn_t)) %>%
  group_by(round(dist_from_first_site_km, 1)) %>%
  mutate(num_cen = n(),
         id = row_number(round(dist_from_first_site_km, 1)),
         dummy_val = as.numeric(str_remove(Zn_t, "<"))/2 - (num_cen/2 - row_number(dist_from_first_site_km) + 0.5) * as.numeric(str_remove(Zn_t, "<"))/8)


Znplot <- ggplot(filter(water_chem_censored_raw, !grepl("<",Zn_t)),
                 aes(x = dist_from_first_site_km, y = as.numeric(Zn_t) * 1000)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        #axis.title.x = element_blank(),
        legend.position = "none",
        legend.position.inside = c(.2, .8)) +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Zinc (μg/L)") +
  geom_smooth(data = water_chemistry_nosplit, aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_point(data = Zn_censored, aes(x = dist_from_first_site_km, y = dummy_val * 1000, color = Drainage), shape = 16, size = 2) +
  geom_hline(yintercept = max(as.numeric(str_remove(Zn_censored$Zn_t, "<"))) * 1000, lty = "dashed", color = "darkgray") +
  annotate("rect", xmin = -5, xmax = 25, ymin = 0, ymax = max(as.numeric(str_remove(Zn_censored$Zn_t, "<"))) * 1000, alpha = .5,fill = "white") +
  geom_jitter(aes(color = Drainage)) #+
Znplot


Co_censored <- water_chem_censored_raw %>%
  filter(grepl('<', Co_t)) %>%
  group_by(round(dist_from_first_site_km, 1)) %>%
  mutate(num_cen = n(),
         id = row_number(round(dist_from_first_site_km, 1)),
         dummy_val = as.numeric(str_remove(Co_t, "<"))/2 - (num_cen/2 - row_number(dist_from_first_site_km) + 0.5) * as.numeric(str_remove(Co_t, "<"))/8)

Coplot <- ggplot(filter(water_chem_censored_raw, !grepl("<",Co_t)),
                 aes(x = dist_from_first_site_km, y = as.numeric(Co_t) * 1000)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.title.x = element_blank(),
    legend.position = "none",
    legend.position.inside = c(.2, .8)) +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Cobalt (μg/L)") +
  geom_smooth(data = water_chemistry_nosplit, aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_point(data = Co_censored, aes(x = dist_from_first_site_km, y = dummy_val * 1000, color = Drainage), shape = 16, size = 2) +
  geom_hline(yintercept = max(as.numeric(str_remove(Co_censored$Co_t, "<"))) * 1000, lty = "dashed", color = "darkgray") +
  annotate("rect", xmin = -5, xmax = 25, ymin = 0, ymax = max(as.numeric(str_remove(Co_censored$Co_t, "<"))) * 1000, alpha = .5,fill = "white") +
  geom_jitter(aes(color = Drainage)) #+
Coplot

metals_highcalnali <- plot_grid(Feplot, Alplot, Pbplot, Tiplot, Znplot, Coplot, nrow = 3, labels = "AUTO", rel_heights = c(1,1,1,1.25), align = "v")
metals_highcalnali
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/Fe_Al_Ti_Zn_Co_riverkm_multiriver.png",
       metals_highcalnali, device = "png", width =6.5, height = 8)

Cd_censored <- water_chem_censored_raw %>%
  filter(grepl('<', Cd_t)) %>%
  group_by(round(dist_from_first_site_km, 1)) %>%
  mutate(num_cen = n(),
         id = row_number(round(dist_from_first_site_km, 1)),
         dummy_val = as.numeric(str_remove(Cd_t, "<"))/2 - (num_cen/2 - row_number(dist_from_first_site_km) + 0.5) * as.numeric(str_remove(Cd_t, "<"))/8)


Cdplot <- ggplot(filter(water_chem_censored_raw, !grepl("<",Cd_t)),
                 aes(x = dist_from_first_site_km, y = as.numeric(Cd_t) * 1000)) +
  coord_cartesian(x = c(0, 20), y = c(0.001, 0.025)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
    legend.position = "none",
    legend.direction = "horizontal") +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Cadmium (μg/L)") +
  geom_smooth(data = water_chemistry_nosplit, aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_point(data = Cd_censored, aes(x = dist_from_first_site_km, y = dummy_val * 1000, color = Drainage), shape = 16, size = 2) +
  geom_hline(yintercept = max(as.numeric(str_remove(Cd_censored$Cd_t, "<"))) * 1000, lty = "dashed", color = "darkgray") +
  annotate("rect", xmin = -5, xmax = 25, ymin = 0, ymax = max(as.numeric(str_remove(Cd_censored$Cd_t, "<"))) * 1000, alpha = .5,fill = "white") +
  geom_jitter(aes(color = Drainage)) #+
Cdplot

Cd_censored <- water_chem_censored_raw %>%
  filter(grepl('<', Cd_t)) %>%
  group_by(round(dist_from_first_site_km, 1)) %>%
  mutate(num_cen = n(),
         id = row_number(round(dist_from_first_site_km, 1)),
         dummy_val = as.numeric(str_remove(Cd_t, "<"))/2 - (num_cen/2 - row_number(dist_from_first_site_km) + 0.5) * as.numeric(str_remove(Cd_t, "<"))/8)

# No censored points in Mn plot

Mnplot <- ggplot(filter(water_chemistry_nosplit),
                 aes(x = dist_from_first_site_km, y = Mn_t * 1000)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.position.inside = c(.2, .8)) +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Manganese (μg/L)") +
  geom_smooth(aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
Mnplot

As_censored <- water_chem_censored_raw %>%
  filter(grepl('<', As_t)) %>%
  group_by(round(dist_from_first_site_km, 1)) %>%
  mutate(num_cen = n(),
         id = row_number(round(dist_from_first_site_km, 1)),
         dummy_val = as.numeric(str_remove(As_t, "<"))/2 - (num_cen/2 - row_number(dist_from_first_site_km) + 0.5) * as.numeric(str_remove(As_t, "<"))/8)


Asplot <- ggplot(filter(water_chem_censored_raw, !grepl("<",As_t)),
                 aes(x = dist_from_first_site_km, y = as.numeric(As_t) * 1000)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.position.inside = c(.2, .8)) +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Arsenic (μg/L)") +
  geom_smooth(data = water_chemistry_nosplit, aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_point(data = As_censored, aes(x = dist_from_first_site_km, y = dummy_val * 1000, color = Drainage), shape = 16, size = 2) +
  geom_hline(yintercept = max(as.numeric(str_remove(As_censored$As_t, "<"))) * 1000, lty = "dashed", color = "darkgray") +
  annotate("rect", xmin = -5, xmax = 25, ymin = 0, ymax = max(as.numeric(str_remove(As_censored$As_t, "<"))) * 1000, alpha = .5,fill = "white") +
  geom_jitter(aes(color = Drainage)) #+
Asplot

# No censored values for Ca

Caplot <- ggplot(filter(water_chemistry_nosplit),
                 aes(x = dist_from_first_site_km, y = Ca_t)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.position.inside = c(.2, .8)) +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Calcium (mg/L)") +
  geom_smooth(aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
Caplot

# No censored values for Mg

Mgplot <- ggplot(filter(water_chemistry_nosplit),
                 aes(x = dist_from_first_site_km, y = Mg_t)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(
        legend.position = "none",
        legend.direction = "horizontal") +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Magnesium (mg/L)") +
  geom_smooth(aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
Mgplot

# No censored values for Si

Siplot <- ggplot(filter(water_chemistry_nosplit),
                 aes(x = dist_from_first_site_km, y = Si_t)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.direction = "horizontal") +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Silicon (mg/L)") +
  geom_smooth(aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
Siplot

dummy_Siplot <- ggplot(filter(water_chemistry_nosplit),
                 aes(x = dist_from_first_site_km, y = Si_t)) +
  coord_cartesian(x = c(0, 20)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal") +
  scale_color_manual(values = river_colors) +
  labs(x = "River Distance (km)", y = "Silicon (mg/L)") +
  geom_smooth(aes(color = Drainage, group = Drainage), se = F) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
dummy_Siplot

get_legend_35 <- function(plot, legend_number = 1) {
  # find all legend candidates
  legends <- get_plot_component(plot, "guide-box", return_all = TRUE)
  # find non-zero legends
  idx <- which(vapply(legends, \(x) !inherits(x, "zeroGrob"), TRUE))
  # return either the chosen or the first non-zero legend if it exists,
  # and otherwise the first element (which will be a zeroGrob) 
  if (length(idx) >= legend_number) {
    return(legends[[idx[legend_number]]])
  } else if (length(idx) >= 0) {
    return(legends[[idx[1]]])
  } else {
    return(legends[[1]])
  }
}

leg <- get_legend_35(dummy_Siplot)
leg_plot <- ggdraw(leg)
leg_plot
metals_other <- plot_grid(Cdplot, Mnplot, Asplot, Caplot, Mgplot, Siplot, nrow = 3, labels = "AUTO", rel_heights = c(1,1,1.25), align = "v")
metals_other_wleg <- plot_grid(metals_other, leg, nrow = 2, rel_heights = c(15, 1))
metals_other_wleg
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/Cd_Mn_As_Ca_Mg_riverkm_multiriver.png",
       metals_other_wleg, device = "png", width =6.5, height = 8)


test <- map(names(water_chemistry_nosplit)[grepl("_t", names(water_chemistry_nosplit))], function(name) {
  ggplot(filter(water_chemistry_nosplit),
         aes_string(x = "dist_from_first_site_km", y = name)) +
    coord_cartesian(x = c(0, 20)) +
    theme_bw() +
    theme(
      legend.position = "none",
      legend.direction = "horizontal") +
    scale_color_manual(values = river_colors) +
    labs(x = "River Distance (km)", y = name) +
    geom_smooth(aes(color = Drainage, group = Drainage), se = F) +
    guides(color=guide_legend(title="River")) +
    geom_jitter(aes(color = Drainage)) #+
})
test

sample_sizes <- water_chemistry_nosplit %>%
  group_by(Site.Code) %>%
  summarize(across(DO_sat:Zr_t, function(x) sum(!is.na(x))))
write.csv(sample_sizes, "Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/waterchem_samplesizes.csv", sep = ",")
