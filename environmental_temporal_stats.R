library(tidyverse)
library(ggbiplot)
library(scales)
library(cowplot)
library(factoextra)
library(car)
library(basemaps)
library(vegan)


malcol=rgb(0/255,0/255,139/255)
hetcol=rgb(65/255,105/255,225/255)
bircol=rgb(255/255,0/255,0/255)
viridis_scale = "rocket"
#hex <- hue_pal()(4)
#hex <- hue_pal()(3)
hex <- c('#01161E','#4D91A8','#D295BF')
names(hex) = c("Feb", "May", "Sep")

section_colors <- c('#3C153B','#8B1E3F','#89BD9E','#EBB65C','#DB4C40')
names(section_colors) <- c("Calnali - Downstream","Calnali - Upstream", "Conzintla", "Huazalingo","Pochula")

river_colors <- c('#8B1E3F','#89BD9E','#EBB65C','#DB4C40')
names(river_colors) <- c("Calnali","Conzintla", "Huazalingo","Pochula")

section_shapes <- c(17, 16, 16, 16, 16)
names(section_colors) <- c("Calnali - Downstream","Calnali - Upstream", "Conzintla", "Huazalingo","Pochula")

is_all_numeric <- function(x) {
  !any(is.na(suppressWarnings(as.numeric(na.omit(x))))) & is.character(x)
}

site_distances <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/site_riverkm_distances.csv") %>%
  dplyr::select(-Drainage)



######## Variation within the Calnali #######


water_chemistry_Calnali_censored_raw <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/waterchem_with_averagedsonde_labmetals_data.csv") %>%
  filter(Site.Code != "N/A", Drainage == "Calnali", !Date %in% c("15-Jun-22","16-Jun-22")) %>%
#  filter(Site.Code != "N/A", !) %>%
  mutate(day = sapply(Date, function(x) strsplit(x, split = "-")[[1]][1]),
         month = sapply(Date, function(x) strsplit(x, split = "-")[[1]][2]),
         year = sapply(Date, function(x) strsplit(x, split = "-")[[1]][3]),
         dist_from_hybridstart = ifelse(Drainage == "Calnali", Lon - -98.585, ifelse(Drainage == "Conzintla", Lon - -98.62412, ifelse(Drainage == "Pochula", Lon - -98.57946, ifelse(Drainage == "Huazalingo", Lon - -98.63713, NA)))),
         Turbidity = ifelse(Turbidity_instrument == "orion", Turbidity, NA),
         Drainage = ifelse(Site.Code %in% c("PEZM", "CAPS", "PLAZ","CALL", "TLCD"), "Calnali - Downstream", ifelse(Drainage == "Calnali", "Calnali - Upstream", Drainage)),
         month = ifelse(month == "Jun", "May", month)) %>%
  filter(year %in% c(22,23,24)) %>%
  left_join(site_distances)
water_chemistry_Calnali <- as.data.frame(apply(water_chemistry_Calnali_censored_raw, 2,
                                       FUN = function(x) ifelse(grepl("<",x), as.numeric(str_remove(x, "<"))/2, ifelse(x == "", NA, x)))) %>%
  mutate_if(is_all_numeric,as.numeric)

pca_nonmetals_Calnali <- water_chemistry_Calnali %>%
  #filter(Drainage %in% c("Calnali", "Calnali - Upstream", "Calnali - Downstream", "Conzintla", "Huazalingo", "Pochula")) %>%
  dplyr::select(-matches("_d$")) %>%
  dplyr::select(-matches("_t$"))%>%
  dplyr::select(-c(DO_sat, Temp, fDOM_RFU, DOC_lab, N_lab, NO3, Calcium_Hardness, Color_Apparent, Color_True, Sulfite, Sulfide,
                   Turbidity_instrument, Fe_colorimeter, Cu_free_colorimeter, SO4,
                   Cu_total_colorimeter, Ni_colorimeter, Mn_colorimeter, Weather, Observations, day, year, dist_from_hybridstart, dist_from_first_site_km, dist_from_prev_site_km, elevation, graphing_elevation)) %>%
  rename("fDOM" = fDOM_QSU, "Oxygen" = DO_conc, "Hardness" = Total_Hardness) %>%
  drop_na() %>%
  select_if(function(x) (any(x != 0)))

wc.pca.nonmetals.Calnali <- prcomp(pca_nonmetals_Calnali[,8:16], center = TRUE,scale. = TRUE)
summary(wc.pca.nonmetals.Calnali)
nonmetal_pca_Calnali <- fviz_pca_biplot(wc.pca.nonmetals.Calnali, axes = c(1,2), repel = T, label = "var", alpha.var = 0.20, mean.point = FALSE, col.var = "black", habillage = pca_nonmetals_Calnali$month, title = NULL)+
  theme(legend.position = "inside",
        legend.position.inside = c(.87, .23),
        legend.background = element_rect()) +
  scale_color_manual(values = hex, labels = c("Feb","May-Jun","Sep")) +
  scale_shape_discrete(labels = c("Feb","May-Jun","Sep")) +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title="Month")) +
  labs(x = bquote("PC 1 ("*.(round(summary(wc.pca.nonmetals.Calnali)$importance[2,1],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(wc.pca.nonmetals.Calnali)$importance[2,2],3)*100)*"%)"))
nonmetal_pca_Calnali
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/nonmetal_pca_Calnali_monthly.pdf",
       nonmetal_pca_Calnali, device = "pdf", width = 6.5 , height = 3)


aggregate(PC1 ~ month, data = data.frame(PC1 = wc.pca.nonmetals.Calnali$x[,1], month = pca_nonmetals_Calnali$month), mean)

pca_metals_Calnali <- water_chemistry_Calnali %>%
  filter(!(Site.Code %in% c("TLMC", "COAC")),
         !(Drainage %in% c("Tanchachin","Claro"))) %>%
  dplyr::select(c(Date, Time, Site.Name, Drainage, Site.Code, Lon, Lat, matches("_t$")), month) %>%
  dplyr::select(-c(Te_t, Th_t, W_t, Zr_t, Bi_t, Sn_t, Tl_t, Ag_t, Be_t)) %>% 
  rename_with(function(x) gsub("_t", "", x)) %>% # Removing metals never observed or only observed once
  drop_na() %>%
  dplyr::select_if(function(x) (any(x != 0)))

wc.pca.metals.Calnali <- prcomp(pca_metals_Calnali[,8:37], center = TRUE,scale. = TRUE)
summary(wc.pca.metals.Calnali)
metal_pca_Calnali <- fviz_pca_biplot(wc.pca.metals.Calnali, axes = c(1,2), repel = T, label = "var", alpha.var = 0.20, mean.point = FALSE, col.var = "black", habillage = pca_metals_Calnali$month, title = NULL)+
  theme(legend.position = "inside",
        legend.position.inside = c(.13, .83),
        legend.background = element_rect()) +
  scale_color_manual(values = hex, labels = c("Feb","May-Jun","Sep")) +
  scale_shape_discrete(labels = c("Feb","May-Jun","Sep")) +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title="Month")) +
  labs(x = bquote("PC 1 ("*.(round(summary(wc.pca.metals.Calnali)$importance[2,1],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(wc.pca.metals.Calnali)$importance[2,2],3)*100)*"%)"))
metal_pca_Calnali
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/metal_pca_Calnali_monthly.pdf",
       metal_pca_Calnali, device = "pdf", width = 6.5 , height = 6.5 * (22.3/38.1))

aggregate(PC1 ~ month, data = data.frame(PC1 = wc.pca.metals.Calnali$x[,1], month = pca_metals_Calnali$month), mean)

Calnali_monthly_PCAs <- plot_grid(nonmetal_pca_Calnali, metal_pca_Calnali, nrow = 2, ncol = 1, rel_height = c(1,4), labels = c("A","B"))
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/Calnali_monthly_PCAs.pdf",
       Calnali_monthly_PCAs, width = 6.5, height = 8)
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/Calnali_monthly_PCAs.png",
       Calnali_monthly_PCAs, width = 6.5, height = 8)

water_chemistry_alltime <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/waterchem_with_averagedsonde_labmetals_data.csv") %>%
  filter(Site.Code != "N/A", !Date %in% c("15-Jun-22","16-Jun-22")) %>%
  mutate(day = sapply(Date, function(x) strsplit(x, split = "-")[[1]][1]),
         month = sapply(Date, function(x) strsplit(x, split = "-")[[1]][2]),
         year = sapply(Date, function(x) strsplit(x, split = "-")[[1]][3]),
         dist_from_hybridstart = ifelse(Drainage == "Calnali", Lon - -98.585, ifelse(Drainage == "Conzintla", Lon - -98.62412, ifelse(Drainage == "Pochula", Lon - -98.57946, ifelse(Drainage == "Huazalingo", Lon - -98.63713, NA)))),
         Turbidity = ifelse(Turbidity_instrument == "orion", Turbidity, NA),
         Drainage = ifelse(Site.Code %in% c("PEZM", "CAPS", "PLAZ","CALL", "TLCD"), "Calnali - Downstream", ifelse(Drainage == "Calnali", "Calnali - Upstream", Drainage))) %>%
  filter(year %in% c(22,23,24)) %>%
  left_join(site_distances)
water_chemistry_alltime <- as.data.frame(apply(water_chemistry_alltime, 2,
                                               FUN = function(x) ifelse(grepl("<",x), as.numeric(str_remove(x, "<"))/2, ifelse(x == "", NA, x)))) %>%
  mutate_if(is_all_numeric,as.numeric)

pca_metals_alltime <- water_chemistry_alltime %>%
  filter(!(Site.Code %in% c("TLMC", "COAC")),
         !(Drainage %in% c("Tanchachin","Claro"))) %>%
  dplyr::select(c(Date, Time, Site.Name, Drainage, Site.Code, Lon, Lat, matches("_t$")), day, month, year) %>%
  dplyr::select(-c(Te_t, Th_t, W_t, Zr_t, Bi_t, Sn_t, Tl_t, Ag_t, Be_t)) %>% # Removing metals never observed or only observed once
  rename_with(function(x) gsub("_t", "", x)) %>% 
  mutate(sampdate = paste(month, year),
         samp_datetime = parse_date_time(paste(Time, Date), orders = "H:M d-b-y")) %>%
  drop_na() %>%
  dplyr::select_if(function(x) (any(x != 0)))

wc.pca.metals.alltime <- prcomp(pca_metals_alltime[,8:37], center = TRUE,scale. = TRUE)
summary(wc.pca.metals.alltime)
#plot(wc.pca.metals.alltime$sdev^2 / sum(wc.pca.metals.alltime$sdev^2), xlab="Principal component", ylab="Proportion of variance explained", type='b')
metal_pca_alltime_monthwise <- fviz_pca_biplot(wc.pca.metals.alltime, axes = c(1,2), repel = T, label = "var", alpha.var = 0.20, mean.point = FALSE, col.var = "black", habillage = pca_metals_alltime$month, title = NULL)+
  theme(legend.position = "inside",
        legend.position.inside = c(.13, .23),
        legend.background = element_rect()) +
  scale_color_manual(values = hex, labels = c("Feb","May-Jun","Sep")) +
  scale_shape_discrete(labels = c("Feb","May-Jun","Sep")) +
  guides(color=guide_legend(title="Month"), shape = guide_legend(title="Month"))  +
  labs(x = bquote("PC 1 ("*.(round(summary(wc.pca.metals.alltime)$importance[2,1],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(wc.pca.metals.alltime)$importance[2,2],3)*100)*"%)"))
metal_pca_alltime_monthwise

metal_pca_alltime_drainagewise <- fviz_pca_biplot(wc.pca.metals.alltime, axes = c(1,2), repel = T, label = "var", alpha.var = 0.20, mean.point = FALSE, col.var = "black", habillage = pca_metals_alltime$Drainage, title = NULL)+
  theme(legend.position = "inside",
        legend.position.inside = c(.15, .23),
        legend.background = element_rect()) +
  scale_color_manual(values = section_colors) +
  scale_shape_manual(values = section_shapes) +
  guides(color=guide_legend(title="Drainage"), shape = guide_legend(title="Drainage"))+
  labs(x = bquote("PC 1 ("*.(round(summary(wc.pca.metals.alltime)$importance[2,1],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(wc.pca.metals.alltime)$importance[2,2],3)*100)*"%)"))
metal_pca_alltime_drainagewise

alltime_monthly_PCAs <- plot_grid(metal_pca_alltime_monthwise, metal_pca_alltime_drainagewise, nrow = 2, ncol = 1, labels = c("A","B"))
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/alltime_monthly_PCAs.pdf",
       alltime_monthly_PCAs, width = 6.5, height = 8)
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/alltime_monthly_PCAs.png",
       alltime_monthly_PCAs, width = 6.5, height = 8)


time_sampling <- data.frame(pca_metals_alltime, PC1 = wc.pca.metals.alltime$x[,1], PC2 = wc.pca.metals.alltime$x[,2])
ggplot(filter(time_sampling, month == "Sep", year == "22"), aes(x = samp_datetime, y = PC2, color = Drainage)) + geom_point()
ggplot(filter(time_sampling, month == "Sep", year == "24"), aes(x = samp_datetime, y = PC2, color = Drainage)) + geom_point()


aggregate(PC1 ~ month, data = data.frame(PC1 = wc.pca.metals.alltime$x[,1], month = pca_metals_alltime$month), mean)
aggregate(PC2 ~ month, data = data.frame(PC2 = wc.pca.metals.alltime$x[,2], month = pca_metals_alltime$month), mean)

####### Individual variables in the Calnali


ammonia_model <- aov(log(NH3 + .001) ~ dist_from_first_site_km * month, data = water_chemistry_Calnali)
plot(ammonia_model)
shapiro.test(ammonia_model$residuals)
Anova(ammonia_model)

no2_model <- aov(log(NO2 + .001) ~ dist_from_first_site_km * month, data = water_chemistry_Calnali)
plot(no2_model)
shapiro.test(no2_model$residuals)
Anova(no2_model)
TukeyHSD(no2_model)

fDOM_model <- aov(log(fDOM_QSU) ~ dist_from_first_site_km * month, data = water_chemistry_Calnali)
plot(fDOM_model)
shapiro.test(fDOM_model$residuals)
Anova(fDOM_model)

turb_model <- aov(log(Turbidity) ~ dist_from_first_site_km * month, data = water_chemistry_Calnali)
plot(turb_model)
shapiro.test(turb_model$residuals)
Anova(turb_model)

cond_model <- aov(log(Conductivity) ~ dist_from_first_site_km * month, data = water_chemistry_Calnali)
plot(cond_model)
shapiro.test(cond_model$residuals)
Anova(cond_model)

Cu_model <- aov(log(Cu_t) ~ dist_from_first_site_km * month, data = water_chemistry_Calnali)
plot(Cu_model)
shapiro.test(Cu_model$residuals)
Anova(Cu_model)

Pb_model <- aov(log(Pb_t) ~ dist_from_first_site_km * month, data = water_chemistry_Calnali)
plot(Pb_model)
shapiro.test(Pb_model$residuals)
Anova(Pb_model)
TukeyHSD(Pb_model)

Cd_model <- aov(log(Cd_t) ~ dist_from_first_site_km * month, data = water_chemistry_Calnali)
plot(Cd_model)
shapiro.test(Cd_model$residuals)
Anova(Cd_model)

Mn_model <- aov(log(Mn_t) ~ dist_from_first_site_km * month, data = water_chemistry_Calnali)
plot(Mn_model)
shapiro.test(Mn_model$residuals)
Anova(Mn_model)

As_model <- aov(log(As_t) ~ dist_from_first_site_km * month, data = water_chemistry_Calnali)
plot(As_model)
shapiro.test(As_model$residuals)
Anova(As_model)

Ca_model <- aov(Ca_t ~ dist_from_first_site_km * month, data = water_chemistry_Calnali)
plot(Ca_model)
shapiro.test(Ca_model$residuals)
Anova(Ca_model)

Calnh3plot <- ggplot(water_chemistry_Calnali,
                     aes(x = dist_from_first_site_km, y = NH3)) +
  coord_cartesian(x = c(0, 16)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(hjust = 0),
        legend.position = "none") +
  labs(x = "River Distance (km)", y = "Ammonia (mg/L N)") +
  #lims(y = c(0,5)) +
  scale_color_manual(values = hex) +
  scale_y_log10() +
  #annotate("text", x = 12, y = 5, label = "Outlier at 31 mg/L N", size = 6) +
  geom_smooth(aes(color = month, group = month), se = F) +
  guides(color=guide_legend(title="Month")) +
  geom_jitter(aes(color = month)) #+
Calnh3plot



Calno2plot <- ggplot(water_chemistry_Calnali,
                     aes(x = dist_from_first_site_km, y = NO2)) +
  coord_cartesian(x = c(0, 16)) +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 24),
        legend.position = "inside",
        legend.text = element_text(size = 20),
        legend.position.inside = c(.2, .8)) +
  labs(x = "River Distance (km)", y = "Nitrite (mg/L N)") +
  #lims(y = c(0,5)) +
  scale_color_manual(values = hex) +
  scale_y_log10() +
  #annotate("text", x = 12, y = 5, label = "Outlier at 31 mg/L N", size = 6) +
  geom_smooth(aes(color = month, group = month), se = F) +
  guides(color=guide_legend(title="Month")) +
  geom_jitter(aes(color = month)) #+
Calno2plot


CalDOMplot <- ggplot(water_chemistry_Calnali,
                     aes(x = dist_from_first_site_km, y = fDOM_QSU)) +
  coord_cartesian(x = c(0, 16)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(.2, .65),
        legend.background = element_blank()) +
  scale_color_manual(values = hex) +
  labs(x = "River Distance (km)", y = "DOM (QSU)") +
  geom_smooth(aes(color = month, group = month), se = F) +
  guides(color=guide_legend(title="Month")) +
  geom_jitter(aes(color = month)) #+
CalDOMplot


CalTurbplot <- ggplot(water_chemistry_Calnali,
                     aes(x = dist_from_first_site_km, y = Turbidity)) +
  coord_cartesian(x = c(0, 16)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = hex) +
  labs(x = "River Distance (km)", y = "Turbidity (NTU)") +
  geom_smooth(aes(color = month, group = month), se = F) +
  guides(color=guide_legend(title="Month")) +
  geom_jitter(aes(color = month)) #+
CalTurbplot


CalCondplot <- ggplot(water_chemistry_Calnali,
                      aes(x = dist_from_first_site_km, y = Conductivity)) +
  coord_cartesian(x = c(0, 16)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(hjust = -2),
        legend.position = "none") +
  scale_color_manual(values = hex) +
  labs(x = "River Distance (km)", y = "Conductivity (μS/cm)") +
  geom_smooth(aes(color = month, group = month), se = F) +
  guides(color=guide_legend(title="Month")) +
  geom_jitter(aes(color = month)) #+
CalCondplot

CalCu_censored <- water_chemistry_Calnali_censored_raw %>%
  filter(grepl('<', Cu_t)) %>%
  group_by(round(dist_from_first_site_km, 1)) %>%
  mutate(num_cen = n(),
         id = row_number(round(dist_from_first_site_km, 1)),
         dummy_val = as.numeric(str_remove(Cu_t, "<"))/2 - (num_cen/2 - row_number(dist_from_first_site_km) + 0.5) * as.numeric(str_remove(Cu_t, "<"))/8)

CalCuplot <- ggplot(filter(water_chemistry_Calnali_censored_raw, !grepl("<",Cu_t)),
                    aes(x = dist_from_first_site_km, y = as.numeric(Cu_t) * 1000)) +
  coord_cartesian(x = c(0, 16)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = hex) +
  labs(x = "River Distance (km)", y = "Copper (μg/L)") +
  geom_smooth(data = water_chemistry_Calnali, aes(color = month, group = month), se = F) +
  guides(color=guide_legend(title="Month")) +
  geom_point(data = CalCu_censored, aes(x = dist_from_first_site_km, y = dummy_val * 1000, color = month), shape = 16, size = 2) +
  geom_hline(yintercept = max(as.numeric(str_remove(CalCu_censored$Cu_t, "<"))) * 1000, lty = "dashed", color = "darkgray") +
  annotate("rect", xmin = -5, xmax = 25, ymin = 0, ymax = max(as.numeric(str_remove(CalCu_censored$Cu_t, "<"))) * 1000, alpha = .5,fill = "white") +
  geom_jitter(aes(color = month)) #+
CalCuplot

CalPb_censored <- water_chemistry_Calnali_censored_raw %>%
  filter(grepl('<', Pb_t)) %>%
  group_by(round(dist_from_first_site_km, 1)) %>%
  mutate(num_cen = n(),
         id = row_number(round(dist_from_first_site_km, 1)),
         dummy_val = as.numeric(str_remove(Pb_t, "<"))/2 - (num_cen/2 - row_number(dist_from_first_site_km) + 0.5) * as.numeric(str_remove(Pb_t, "<"))/8)

CalPbplot <- ggplot(filter(water_chemistry_Calnali_censored_raw, !grepl("<",Pb_t)),
                    aes(x = dist_from_first_site_km, y = as.numeric(Pb_t) * 1000)) +
  coord_cartesian(x = c(0, 16)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = hex) +
  labs(x = "River Distance (km)", y = "Lead (μg/L)") +
  geom_smooth(data = water_chemistry_Calnali, aes(color = month, group = month), se = F) +
  guides(color=guide_legend(title="Month")) +
  geom_point(data = CalPb_censored, aes(x = dist_from_first_site_km, y = dummy_val * 1000, color = month), shape = 16, size = 2) +
  geom_hline(yintercept = max(as.numeric(str_remove(CalPb_censored$Pb_t, "<"))) * 1000, lty = "dashed", color = "darkgray") +
  annotate("rect", xmin = -5, xmax = 25, ymin = 0, ymax = max(as.numeric(str_remove(CalPb_censored$Pb_t, "<"))) * 1000, alpha = .5,fill = "white") +
  geom_jitter(aes(color = month)) #+
CalPbplot

CalCd_censored <- water_chemistry_Calnali_censored_raw %>%
  filter(grepl('<', Cd_t)) %>%
  group_by(round(dist_from_first_site_km, 1)) %>%
  mutate(num_cen = n(),
         id = row_number(round(dist_from_first_site_km, 1)),
         dummy_val = as.numeric(str_remove(Cd_t, "<"))/2 - (num_cen/2 - row_number(dist_from_first_site_km) + 0.5) * as.numeric(str_remove(Cd_t, "<"))/8)

CalCdplot <- ggplot(filter(water_chemistry_Calnali_censored_raw, !grepl("<",Cd_t)),
                    aes(x = dist_from_first_site_km, y = as.numeric(Cd_t) * 1000)) +
  coord_cartesian(x = c(0, 16)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = hex) +
  labs(x = "River Distance (km)", y = "Cadmium (μg/L)") +
  geom_smooth(data = water_chemistry_Calnali, aes(color = month, group = month), se = F) +
  guides(color=guide_legend(title="Month")) +
  geom_point(data = CalCd_censored, aes(x = dist_from_first_site_km, y = dummy_val * 1000, color = month), shape = 16, size = 2) +
  geom_hline(yintercept = max(as.numeric(str_remove(CalCd_censored$Cd_t, "<"))) * 1000, lty = "dashed", color = "darkgray") +
  annotate("rect", xmin = -5, xmax = 25, ymin = 0, ymax = max(as.numeric(str_remove(CalCd_censored$Cd_t, "<"))) * 1000, alpha = .5,fill = "white") +
  geom_jitter(aes(color = month)) #+
CalCdplot

# None below the detection limit
CalMnplot <- ggplot(water_chemistry_Calnali,
                      aes(x = dist_from_first_site_km, y = Mn_t * 1000)) +
  coord_cartesian(x = c(0, 16)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(hjust = 0),
        legend.position = "none") +
  scale_color_manual(values = hex) +
  labs(x = "River Distance (km)", y = "Manganese (μg/L)") +
  geom_smooth(aes(color = month, group = month), se = F) +
  guides(color=guide_legend(title="Month")) +
  geom_jitter(aes(color = month)) #+
CalMnplot

CalAs_censored <- water_chemistry_Calnali_censored_raw %>%
  filter(grepl('<', As_t)) %>%
  group_by(round(dist_from_first_site_km, 1)) %>%
  mutate(num_cen = n(),
         id = row_number(round(dist_from_first_site_km, 1)),
         dummy_val = as.numeric(str_remove(As_t, "<"))/2 - (num_cen/2 - row_number(dist_from_first_site_km) + 0.5) * as.numeric(str_remove(As_t, "<"))/8)

CalAsplot <- ggplot(filter(water_chemistry_Calnali_censored_raw, !grepl("<", As_t)),
                    aes(x = dist_from_first_site_km, y = as.numeric(As_t) * 1000)) +
  coord_cartesian(x = c(0, 16)) +
  theme_bw() +
  theme(axis.title.y = element_text(hjust = 0),
        legend.position = "none") +
  scale_color_manual(values = hex) +
  labs(x = "River Distance (km)", y = "Arsenic (μg/L)") +
  geom_smooth(data = water_chemistry_Calnali, aes(color = month, group = month), se = F) +
  guides(color=guide_legend(title="Month")) +
  geom_point(data = CalAs_censored, aes(x = dist_from_first_site_km, y = dummy_val * 1000, color = Drainage), shape = 16, size = 2) +
  geom_hline(yintercept = max(as.numeric(str_remove(CalAs_censored$As_t, "<"))) * 1000, lty = "dashed", color = "darkgray") +
  annotate("rect", xmin = -5, xmax = 25, ymin = 0, ymax = max(as.numeric(str_remove(CalAs_censored$As_t, "<"))) * 1000, alpha = .5,fill = "white") +
  geom_jitter(aes(color = month)) #+
CalAsplot


CalCaplot <- ggplot(water_chemistry_Calnali,
                    aes(x = dist_from_first_site_km, y = Ca_t)) +
  coord_cartesian(x = c(0, 16)) +
  theme_bw() +
  theme(axis.title.y = element_text(hjust = 0),
        legend.position = "none") +
  scale_color_manual(values = hex) +
  labs(x = "River Distance (km)", y = "Calcium (mg/L)") +
  geom_smooth(aes(color = month, group = month), se = F) +
  guides(color=guide_legend(title="Month")) +
  geom_jitter(aes(color = month)) #+
CalCaplot

Calnali_monthly_line_plots <- plot_grid(Calnh3plot, CalDOMplot, CalTurbplot, CalCondplot, CalCuplot, CalPbplot, CalCdplot, CalMnplot, CalAsplot, CalCaplot, 
          ncol = 2, labels = "AUTO")
Calnali_monthly_line_plots
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/Calnali_monthly_line_plots.png", 
       Calnali_monthly_line_plots, device = "png", width = 6.5, height = 9)

sample_sizes_alltime <- water_chemistry_alltime %>%
  group_by(Site.Code) %>%
  summarize(across(DO_sat:Zr_t, function(x) sum(!is.na(x))))
write.csv(sample_sizes_alltime, "Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/waterchem_samplesizes_alltime.csv", sep = ",")

