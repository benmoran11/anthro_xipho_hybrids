library(tidyverse)
library(cowplot)
library(vegan)

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
  left_join(site_distances)
#set.seed(64)
water_chemistry <- as.data.frame(apply(water_chemistry, 2,
                                       FUN = function(x) ifelse(grepl("<",x), as.numeric(str_remove(x, "<"))/2, ifelse(x == "", NA, x)))) %>%
  #FUN = function(x) ifelse(grepl("<",x), runif(1, min = 0, max = as.numeric(str_remove(x, "<"))), ifelse(x == "", NA, x)))) %>%
  #FUN = function(x) ifelse(grepl("<",x), 0, ifelse(x == "", NA, x)))) %>%
  mutate_if(is_all_numeric,as.numeric)

sites_drainages <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/all_sites_coords_PAPA_HZNP_Conz.csv")

## And now land use

#landuse_fullwatershed <- readxl::read_xlsx("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/LandUse_Sentinel_APXC_full_2020s_Spe18Spa19min20_202412141134378024739_PixelCounts_Watersheds.xlsx",
#                                           na = "NA") %>%
landuse_fullwatershed <- read.table("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/Subcatchment_LandUse_full.csv", sep = ",", header = T) %>%
  left_join(site_distances) %>%
  mutate(prop_misclassified_total = Misclassified_developed/Total_pixels,
         prop_misclassified_developed = Misclassified_developed/developed) %>%
  left_join(sites_drainages)
#landuse_500m <- readxl::read_xlsx("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/LandUse_Sentinel_APXC_full_2020s_Spe18Spa19min20_202412141134378024739_PixelCounts_Watersheds_half_km_buffer.xlsx",
#                                  na = "NA") %>%
landuse_500m <- read.table("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/Subcatchment_LandUse_500m.csv", sep = ",", header = T) %>%
  left_join( site_distances) %>%
  mutate(prop_misclassified_total = Misclassified_developed/Total_pixels,
         prop_misclassified_developed = Misclassified_developed/developed) %>%
  left_join(sites_drainages)
#landuse_100m <- readxl::read_xlsx("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/LandUse_Sentinel_APXC_full_2020s_Spe18Spa19min20_202412141134378024739_PixelCounts_Watersheds_100mBuffer.xlsx",
#                                  na = "NA") %>%
landuse_100m <- read.table("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/Subcatchment_LandUse_100m.csv", sep = ",", header = T) %>%
  left_join( site_distances) %>%
  mutate(prop_misclassified_total = Misclassified_developed/Total_pixels,
         prop_misclassified_developed = Misclassified_developed/developed) %>%
  left_join(sites_drainages)
#landuse_50m <- readxl::read_xlsx("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/LandUse_Sentinel_APXC_full_2020s_Spe18Spa19min20_202412141134378024739_PixelCounts_Watersheds_50mBuffer.xlsx",
#                                 na = "NA") %>%
landuse_50m <- read.table("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/Subcatchment_LandUse_50m.csv", sep = ",", header = T) %>%
  left_join( site_distances) %>%
  mutate(prop_misclassified_total = Misclassified_developed/Total_pixels,
         prop_misclassified_developed = Misclassified_developed/developed) %>%
  left_join(sites_drainages)


# Loading in water chem data

pca_ready_siteavgs <- water_chemistry %>%
  filter(!(Site.Code %in% c("TLMC", "COAC")),
         !(Drainage %in% c("Tanchachin","Claro"))) %>%
  dplyr::select(-matches("_d$")) %>%
  dplyr::select(-c(Date, Time, Lon, Lat, Site.Upstream, DO_sat, Temp, fDOM_RFU, DOC_lab, N_lab, NO3, Calcium_Hardness, Color_Apparent, Color_True, Sulfite, Sulfide,
                   Turbidity_instrument, Fe_colorimeter, Cu_free_colorimeter,
                   Cu_total_colorimeter, Ni_colorimeter, Mn_colorimeter, Weather, Observations, day, month, year, dist_from_hybridstart, dist_from_first_site_km, dist_from_prev_site_km, elevation, graphing_elevation)) %>%
  rename("fDOM" = fDOM_QSU, "Oxygen" = DO_conc, "Hardness" = Total_Hardness) %>%
  #dplyr::select(c(Date, Time, Site.Name, Drainage, Site.Code, Lon, Lat, matches("_t$")), month) %>%
  dplyr::select(-c(Te_t, Th_t, W_t, Zr_t, Bi_t, Sn_t, Tl_t, Ag_t, Be_t)) %>% 
  rename_with(function(x) gsub("_t", "", x)) %>% # Removing metals never observed or only observed once
  drop_na() %>%
  dplyr::select_if(function(x) (any(x != 0))) %>%
  group_by(Site.Name, Drainage, Site.Code) %>%
  summarize(across(Oxygen:Zn, mean))

all.pca <- prcomp(pca_ready_siteavgs[,4:44], center = TRUE,scale. = TRUE)
major_axes <- cbind(pca_ready_siteavgs,all.pca$x[, 1:4])


pca_ready_siteavgs_nosplit <- pca_ready_siteavgs %>%
  mutate(Drainage = ifelse(Drainage == "Calnali - Downstream", "Calnali", ifelse(Drainage == "Calnali - Upstream", "Calnali", Drainage)))
major_axes_nosplit <- major_axes %>% 
  mutate(Drainage = ifelse(Drainage == "Calnali - Downstream", "Calnali", ifelse(Drainage == "Calnali - Upstream", "Calnali", Drainage)))
landuse_chem_full <- left_join(landuse_fullwatershed, pca_ready_siteavgs_nosplit) %>%
  left_join(major_axes_nosplit) %>%
  mutate(cum_total_km2 = Cum_total / 10000)
landuse_chem_500m <- left_join(landuse_500m, pca_ready_siteavgs_nosplit) %>%
  left_join(major_axes_nosplit)%>%
  mutate(cum_total_km2 = Cum_total / 10000)
landuse_chem_100m <- left_join(landuse_100m, pca_ready_siteavgs_nosplit) %>%
  left_join(major_axes_nosplit)%>%
  mutate(cum_total_km2 = Cum_total / 10000)
landuse_chem_50m <- left_join(landuse_50m, pca_ready_siteavgs_nosplit) %>%
  left_join(major_axes_nosplit)%>%
  mutate(cum_total_km2 = Cum_total / 10000)



developed_plot <- ggplot(filter(landuse_fullwatershed, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                         aes(x = elevation, y = corrected_prop_developed * 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside = c(.8, .8)) +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("Full Watershed","% Developed"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) +
  ggtitle("Individual Site")
developed_plot

developed_plot_500m <- ggplot(filter(landuse_500m, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                              aes(x = elevation, y = corrected_prop_developed * 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside = c(.8, .8)) +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("500 m Buffer","% Developed"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
developed_plot_500m

developed_plot_100m <- ggplot(filter(landuse_100m, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                              aes(x = elevation, y = corrected_prop_developed * 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside = c(.8, .8)) +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("100 m Buffer","% Developed"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
developed_plot_100m

developed_plot_50m <- ggplot(filter(landuse_50m, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                              aes(x = elevation, y = corrected_prop_developed * 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside = c(.8, .8)) +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("% Developed", "50 m Buffer"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
developed_plot_50m

#### Cumulative Input now

developed_plot_cum_full <- ggplot(filter(landuse_fullwatershed, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                             aes(x = elevation, y = corrected_cumprop_developed * 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside = c(.8, .8)) +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("Full Watershed","% Developed"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) +
  ggtitle("Cumulative")
developed_plot_cum_full

developed_plot_cum_500m <- ggplot(filter(landuse_500m, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                                  aes(x = elevation, y = corrected_cumprop_developed * 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside = c(.8, .8)) +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("500 m Buffer","% Developed"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
developed_plot_cum_500m


developed_plot_cum_100m <- ggplot(filter(landuse_100m, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                                  aes(x = elevation, y = corrected_cumprop_developed * 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside = c(.8, .8)) +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("100 m Buffer","% Developed"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
developed_plot_cum_100m

developed_plot_cum_50m <- ggplot(filter(landuse_50m, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                                  aes(x = elevation, y = corrected_cumprop_developed * 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside = c(.8, .8)) +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("50 m Buffer","% Developed"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
developed_plot_cum_50m

dummy_developed_plot_cum_50m <- ggplot(filter(landuse_50m, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                                 aes(x = elevation, y = corrected_cumprop_developed * 100)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("50 m Buffer","% Developed"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
dummy_developed_plot_cum_50m


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

leg <- get_legend_35(dummy_developed_plot_cum_50m)
all_developed_plots <- plot_grid(developed_plot, developed_plot_cum_full,
          developed_plot_500m, developed_plot_cum_500m,
          developed_plot_100m, developed_plot_cum_100m,
          developed_plot_50m, developed_plot_cum_50m,
          nrow = 4, labels = "AUTO")
all_developed_plots_wleg <- plot_grid(all_developed_plots, leg, nrow = 2, rel_heights = c(15, 1))
all_developed_plots_wleg
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/developed_land_megaplot.pdf",
       all_developed_plots_wleg, device = "pdf", width =6.5, height = 8)





##### Now cleared land


herb_plot <- ggplot(filter(landuse_fullwatershed, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                         aes(x = elevation, y = corrected_prop_herb * 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside = c(.8, .8)) +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("Full Watershed","% Cleared"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) +
  ggtitle("Individual Site")
herb_plot

herb_plot_500m <- ggplot(filter(landuse_500m, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                              aes(x = elevation, y = corrected_prop_herb * 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside = c(.8, .8)) +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("500 m Buffer","% Cleared"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
herb_plot_500m

herb_plot_100m <- ggplot(filter(landuse_100m, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                              aes(x = elevation, y = corrected_prop_herb * 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside = c(.8, .8)) +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("100 m Buffer","% Cleared"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
herb_plot_100m

herb_plot_50m <- ggplot(filter(landuse_50m, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                             aes(x = elevation, y = corrected_prop_herb * 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside = c(.8, .8)) +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("50 m Buffer","% Cleared"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
herb_plot_50m

#### Cumulative Input now

herb_plot_cum_full <- ggplot(filter(landuse_fullwatershed, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                                  aes(x = elevation, y = corrected_cumprop_herb * 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside = c(.8, .8)) +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("Full Watershed","% Cleared"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) +
  ggtitle("Cumulative")
herb_plot_cum_full

herb_plot_cum_500m <- ggplot(filter(landuse_500m, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                                  aes(x = elevation, y = corrected_cumprop_herb * 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside = c(.8, .8)) +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("500 m Buffer","% Cleared"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
herb_plot_cum_500m


herb_plot_cum_100m <- ggplot(filter(landuse_100m, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                                  aes(x = elevation, y = corrected_cumprop_herb * 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside = c(.8, .8)) +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("100 m Buffer","% Cleared"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
herb_plot_cum_100m

herb_plot_cum_50m <- ggplot(filter(landuse_50m, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                                 aes(x = elevation, y = corrected_cumprop_herb * 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside = c(.8, .8)) +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("50 m Buffer","% Cleared"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
herb_plot_cum_50m

dummy_herb_plot_cum_50m <- ggplot(filter(landuse_chem_50m, Drainage %in% c("Calnali", "Pochula", "Huazalingo", "Conzintla")),
                                       aes(x = elevation, y = corrected_cumprop_herb * 100)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_reverse()+
  scale_color_manual(values = river_colors) +
  labs(x = "Elevation (m)", y = expression(atop("50 m Buffer","% Cleared"))) +
  geom_line(aes(color = Drainage, group = Drainage)) +
  guides(color=guide_legend(title="River")) +
  geom_jitter(aes(color = Drainage)) #+
dummy_herb_plot_cum_50m


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

leg <- get_legend_35(dummy_herb_plot_cum_50m)
all_herb_plots <- plot_grid(herb_plot, herb_plot_cum_full,
                                 herb_plot_500m, herb_plot_cum_500m,
                                 herb_plot_100m, herb_plot_cum_100m,
                                 herb_plot_50m, herb_plot_cum_50m,
                                 nrow = 4, labels = "AUTO")
all_herb_plots_wleg <- plot_grid(all_herb_plots, leg, nrow = 2, rel_heights = c(15, 1))
all_herb_plots_wleg
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/herb_land_megaplot.pdf",
       all_herb_plots_wleg, device = "pdf", width =6.5, height = 8)















### Land Use Stats by drainage

landuse_chem_full_split <- landuse_chem_full %>%
  mutate(Drainage = ifelse(Site.Code %in% c("PLAZ", "PEAT", "CAPS", "CALL", "TLCD", "PEZM", "CHAF"), "Calnali - Downstream", ifelse(Drainage == "Calnali", "Calnali - Upstream", Drainage)))
landuse_chem_500m_split <- landuse_chem_500m %>%
  mutate(Drainage = ifelse(Site.Code %in% c("PLAZ", "PEAT", "CAPS", "CALL", "TLCD", "PEZM", "CHAF"), "Calnali - Downstream", ifelse(Drainage == "Calnali", "Calnali - Upstream", Drainage)))
landuse_chem_100m_split <- landuse_chem_100m %>%
  mutate(Drainage = ifelse(Site.Code %in% c("PLAZ", "PEAT", "CAPS", "CALL", "TLCD", "PEZM", "CHAF"), "Calnali - Downstream", ifelse(Drainage == "Calnali", "Calnali - Upstream", Drainage)))
landuse_chem_50m_split <- landuse_chem_50m %>%
  mutate(Drainage = ifelse(Site.Code %in% c("PLAZ", "PEAT", "CAPS", "CALL", "TLCD", "PEZM", "CHAF"), "Calnali - Downstream", ifelse(Drainage == "Calnali", "Calnali - Upstream", Drainage)))

Anova(aov(corrected_cumprop_developed ~ Drainage, landuse_chem_full_split))
ggplot(landuse_chem_full_split, aes(x = Drainage, y = corrected_prop_developed)) +
  theme_bw() +
  geom_boxplot()
ggplot(landuse_chem_full_split, aes(x = Drainage, y = corrected_cumprop_developed)) +
  theme_bw() +
  geom_boxplot()
ggplot(landuse_chem_full_split, aes(x = Drainage, y = corrected_prop_developed)) +
  theme_bw() +
  geom_jitter()+
  stat_summary(fun.data = "mean_se")
ggplot(landuse_chem_full_split, aes(x = Drainage, y = corrected_cumprop_developed)) +
  theme_bw() +
  geom_jitter()+
  stat_summary(fun.data = "mean_se")







### Water chemistry statistics



Cu_land_model <- lm(Cu ~ corrected_cumprop_developed + 
                      corrected_prop_developed + 
                      corrected_cumprop_herb +
                      corrected_prop_herb, data = landuse_chem_full)
shapiro.test(Cu_land_model$residuals)
summary(Cu_land_model)
vif(Cu_land_model)
plot(Cu_land_model)

fDOM_land_model <- lm(log(fDOM) ~ corrected_cumprop_developed + 
                        corrected_prop_developed + 
                        corrected_cumprop_herb +
                        corrected_prop_herb, data = landuse_chem_full)
shapiro.test(fDOM_land_model$residuals)
summary(fDOM_land_model)
vif(fDOM_land_model)
plot(fDOM_land_model)

NH3_land_model <- lm(log(NH3 + .001) ~ corrected_cumprop_developed + 
                       corrected_prop_developed + 
                       corrected_cumprop_herb +
                       corrected_prop_herb, data = landuse_chem_full)
shapiro.test(NH3_land_model$residuals)
summary(NH3_land_model)
vif(NH3_land_model)
plot(NH3_land_model)

PC1_land_model <- lm(PC1 ~ corrected_cumprop_developed + 
                       corrected_prop_developed + 
                       corrected_cumprop_herb +
                       corrected_prop_herb + cum_total_km2, data = landuse_chem_full)
shapiro.test(PC1_land_model$residuals)
summary(PC1_land_model)
vif(PC1_land_model)
plot(PC1_land_model)

PC1_land_model_500m_buff <- lm(PC1 ~ corrected_cumprop_developed + 
                                corrected_prop_developed + 
                                corrected_cumprop_herb +
                                corrected_prop_herb + cum_total_km2, data = landuse_chem_500m)
shapiro.test(PC1_land_model_500m_buff$residuals)
summary(PC1_land_model_500m_buff)
vif(PC1_land_model_500m_buff)
plot(PC1_land_model_500m_buff)

PC1_land_model_100m_buff <- lm(PC1 ~ corrected_cumprop_developed + 
                                 corrected_prop_developed + 
                                 corrected_cumprop_herb +
                                 corrected_prop_herb + cum_total_km2, data = landuse_chem_100m)
shapiro.test(PC1_land_model_100m_buff$residuals)
summary(PC1_land_model_100m_buff)
vif(PC1_land_model_100m_buff)
plot(PC1_land_model_100m_buff)

PC1_land_model_50m_buff <- lm(PC1 ~ corrected_cumprop_developed + 
                       corrected_prop_developed + 
                       corrected_cumprop_herb +
                       corrected_prop_herb + cum_total_km2, data = landuse_chem_50m)
shapiro.test(PC1_land_model_50m_buff$residuals)
summary(PC1_land_model_50m_buff)
vif(PC1_land_model_50m_buff)
plot(PC1_land_model_50m_buff)

PC2_land_model <- lm(PC2 ~ corrected_cumprop_developed + 
                       corrected_prop_developed + 
                       corrected_cumprop_herb +
                       corrected_prop_herb + cum_total_km2, data = landuse_chem_full)
shapiro.test(PC2_land_model$residuals)
summary(PC2_land_model)
vif(PC2_land_model)
plot(PC2_land_model)

PC2_land_model_500m_buff <- lm(PC2 ~ corrected_cumprop_developed + 
                                 corrected_prop_developed + 
                                 corrected_cumprop_herb +
                                 corrected_prop_herb + cum_total_km2, data = landuse_chem_500m)
shapiro.test(PC2_land_model_500m_buff$residuals)
summary(PC2_land_model_500m_buff)
vif(PC2_land_model_500m_buff)
plot(PC2_land_model_500m_buff)

PC2_land_model_100m_buff <- lm(PC2 ~ corrected_cumprop_developed + 
                                 corrected_prop_developed + 
                                 corrected_cumprop_herb +
                                 corrected_prop_herb + cum_total_km2, data = landuse_chem_100m)
shapiro.test(PC2_land_model_100m_buff$residuals)
summary(PC2_land_model_100m_buff)
vif(PC2_land_model_100m_buff)
plot(PC2_land_model_100m_buff)

PC2_land_model_50m_buff <- lm(PC2 ~ corrected_cumprop_developed + 
                                corrected_prop_developed + 
                                corrected_cumprop_herb +
                                corrected_prop_herb + cum_total_km2, data = landuse_chem_50m)
shapiro.test(PC2_land_model_50m_buff$residuals)
summary(PC2_land_model_50m_buff)
vif(PC2_land_model_50m_buff)
plot(PC2_land_model_50m_buff)

PC3_land_model <- lm(PC3 ~ corrected_cumprop_developed + 
                       corrected_prop_developed + 
                       corrected_cumprop_herb +
                       corrected_prop_herb, data = landuse_chem_full)
shapiro.test(PC3_land_model$residuals)
summary(PC3_land_model)
vif(PC3_land_model)
plot(PC3_land_model)

PC3_land_model_50m_buff <- lm(PC3 ~ corrected_cumprop_developed + 
                                corrected_prop_developed + 
                                corrected_cumprop_herb +
                                corrected_prop_herb, data = landuse_chem_50m)
shapiro.test(PC3_land_model_50m_buff$residuals)
summary(PC3_land_model_50m_buff)
vif(PC3_land_model_50m_buff)
plot(PC3_land_model_50m_buff)

PC4_land_model <- lm(PC4 ~ corrected_cumprop_developed + 
                       corrected_prop_developed + 
                       corrected_cumprop_herb +
                       corrected_prop_herb, data = landuse_chem_full)
shapiro.test(PC4_land_model$residuals)
summary(PC4_land_model)
vif(PC4_land_model)
plot(PC4_land_model)

PC4_land_model_50m_buff <- lm(PC4 ~ corrected_cumprop_developed + 
                                corrected_prop_developed + 
                                corrected_cumprop_herb +
                                corrected_prop_herb, data = landuse_chem_50m)
shapiro.test(PC4_land_model_50m_buff$residuals)
summary(PC4_land_model_50m_buff)
vif(PC3_land_model_50m_buff)
plot(PC4_land_model_50m_buff)

### Statistical Tests for association between chemical PCs and land use

PC_land_model<-lm(corrected_cumprop_developed ~ PC1 + PC2 + PC3 + PC4, 
                  data = landuse_chem_50m)
shapiro.test(PC_land_model$residuals)
summary(PC_land_model)
vif(PC_land_model)
plot(PC_land_model)



ggplot(landuse_chem_full, aes(x = corrected_cumprop_developed, y = PC1, color = Drainage)) +
  theme_bw() +
  geom_point()

ggplot(landuse_chem_500m, aes(x = corrected_cumprop_developed, y = PC2, color = Drainage)) +
  theme_bw() +
  geom_point()

ggplot(landuse_chem_100m, aes(x = corrected_cumprop_developed, y = PC3, color = Drainage)) +
  theme_bw() +
  geom_point()

ggplot(landuse_chem_50m, aes(x = corrected_cumprop_developed, y = PC4, color = Drainage)) +
  theme_bw() +
  geom_point()

cor.test(x = landuse_chem_full$corrected_cumprop_developed, y = landuse_chem_full$Cu)
cor.test(x = landuse_chem_500m$corrected_cumprop_developed, y = landuse_chem_full$Cu)
cor.test(x = landuse_chem_100m$corrected_cumprop_developed, y = landuse_chem_full$Cu)
cor.test(x = landuse_chem_50m$corrected_cumprop_developed, y = landuse_chem_full$Cu)

plotty <- landuse_chem_full %>% dplyr::select(corrected_prop_developed, 
                                              corrected_cumprop_developed,
                                              corrected_prop_herb,
                                              corrected_cumprop_herb,
                                              Cum_total,
                                              PC1,
                                              PC2)

plot(plotty)

plotty_50m <- landuse_chem_50m %>% dplyr::select(corrected_prop_developed, 
                                              corrected_cumprop_developed,
                                              corrected_prop_herb,
                                              corrected_cumprop_herb,
                                              Cum_total,
                                              PC1,
                                              PC2)

plot(plotty_50m)


### What if we do a site-wise PCA with land use variables?

chem_land_pca_ready <- landuse_chem_full_split %>%
  select(Drainage, Site.Code, corrected_prop_developed, corrected_prop_herb, corrected_prop_forest, corrected_cumprop_developed, corrected_cumprop_herb, corrected_cumprop_forest, Oxygen:Zn) %>%
  mutate(prop_developed_50m = landuse_chem_50m_split$corrected_prop_developed, prop_herb_50m = landuse_chem_50m_split$corrected_prop_herb, prop_forest_50m = landuse_chem_50m_split$corrected_prop_forest, cumprop_developed_50m = landuse_chem_50m_split$corrected_cumprop_developed, cumprop_herb_50m = landuse_chem_50m_split$corrected_cumprop_herb, cumprop_forest_50m = landuse_chem_50m_split$corrected_cumprop_forest) %>%
  drop_na()
landchem.pca <- prcomp(chem_land_pca_ready[,3:55], center = TRUE,scale. = TRUE)

summary(landchem.pca)
landchem_pca_vars <- get_pca_var(landchem.pca)
plot(landchem.pca$sdev^2 / sum(landchem.pca$sdev^2), xlab="Principal component", ylab="Proportion of variance explained", type='b')

landchem_pca <- fviz_pca_biplot(landchem.pca, axes = c(1,2), repel = T, label = "var", alpha.var = 0.20, col.var = "black", mean.point = FALSE, habillage = chem_land_pca_ready$Drainage, title = NULL)  +
  theme(legend.position = "inside",
        legend.position.inside = c(.8, .8),
        legend.background = element_rect()) +
  scale_color_manual(values = section_colors) +
  scale_fill_manual(values = section_colors) +
  scale_shape_manual(values = section_shapes) +
  labs(x = bquote("PC 1 ("*.(round(summary(landchem.pca)$importance[2,1],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(landchem.pca)$importance[2,2],3)*100)*"%)"))
landchem_pca



### Trying out RDA
landuse_vars <- cbind(filter(landuse_fullwatershed, Site.Code %in% pca_ready_siteavgs$Site.Code)[c(25,27,29,31,33,42)], 
                           filter(landuse_50m, Site.Code %in% pca_ready_siteavgs$Site.Code)[c(26,28,32,34)])
rownames(landuse_vars) <- filter(landuse_fullwatershed, Site.Code %in% pca_ready_siteavgs$Site.Code)$Site.Code
colnames(landuse_vars) <- c("corrected_prop_developed", "corrected_prop_herb", "corrected_cum_total",  "corrected_cumprop_developed", "corrected_cumprop_herb","Drainage", "corrected_prop_developed_50m", "corrected_prop_herb_50m", "corrected_cumprop_developed_50m", "corrected_cumprop_herb_50m")
landuse_vars <- landuse_vars %>%
  mutate(log_cum_area = log(corrected_cum_total))
constraining_vars <-  landuse_vars %>%
  select(Drainage, everything()) %>%
  mutate(across(starts_with("corrected"), scale))
#constraining_vars <- filter(landuse_fullwatershed, Site.Code %in% pca_ready_siteavgs$Site.Code)[c(25,26,27,31,32,33)]
#rownames(constraining_vars) <- filter(landuse_fullwatershed, Site.Code %in% pca_ready_siteavgs$Site.Code)$Site.Code
#colnames(constraining_vars) <- c("corrected_prop_developed", "corrected_prop_forest", "corrected_prop_herb", "corrected_cumprop_developed", "corrected_cumprop_forest", "corrected_cumprop_herb")

allchem <- pca_ready_siteavgs[4:44]
rownames(allchem) <- pca_ready_siteavgs$Site.Code

enviro_vars<- pca_ready_siteavgs[4:44] %>%
  scale()
rownames(enviro_vars) <- pca_ready_siteavgs$Site.Code

variable_variables <- allchem %>%
  summarize(across(1:41, function(x) length(unique(x)) >= 20))
enviro_vars <- enviro_vars[,as.logical(variable_variables[1,])] %>%
#  scale()

enviro_vars <- pca_ready_siteavgs %>%
  ungroup() %>%
  select(fDOM, NH3, PO4, Turbidity, Al, Cd, Cu, Fe, Pb, P, K, Mn, Conductivity, Sr, S, Na, Rb, Mg, Co, Ca, B, As) %>%
  scale()
rownames(enviro_vars) <- pca_ready_siteavgs$Site.Code

chem_land.rda <- rda(formula = enviro_vars ~ corrected_cum_total +
                       Drainage + 
                       corrected_prop_developed +
                       corrected_prop_herb + 
                       corrected_cumprop_developed + 
                       corrected_cumprop_herb + 
                       corrected_prop_developed_50m + 
                       corrected_prop_herb_50m + 
                       corrected_cumprop_developed_50m + 
                       corrected_cumprop_herb_50m, 
                     data = landuse_vars)
fwd.sel <- ordiR2step(rda(enviro_vars ~ 1, data = landuse_vars),
                      scope = formula(chem_land.rda), direction = "forward", R2scope = TRUE, Pin = 0.99,
                      pstep = 10000, trace = T)
fwd.sel$call

plot(chem_land.rda, choices = c(1,2), scaling = 3)
summary(chem_land.rda)
RsquareAdj(chem_land.rda)
screeplot(chem_land.rda)
tryit <- scores(chem_land.rda, display = c("sp", "lc"), scaling = 3, tidy = T)
anova.cca(chem_land.rda)
