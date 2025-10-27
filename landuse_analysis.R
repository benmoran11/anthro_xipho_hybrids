library(tidyverse)
library(cowplot)
library(vegan)

# Colors to be used for different drainages - keeping Calnali together
river_colors <- c('#8B1E3F','#89BD9E','#EBB65C','#DB4C40')
names(river_colors) <- c("Calnali","Conzintla", "Huazalingo","Pochula")

# Utility to check for numeric values - useful for parsing cases where numerics  
# not initially parsed successfully because of "<" signs present in ICP/MS data
is_all_numeric <- function(x) {
  !any(is.na(suppressWarnings(as.numeric(na.omit(x))))) & is.character(x)
}

# Loading site distances along river
site_distances <- read.csv("site_riverkm_distances.csv") %>%
  dplyr::select(-Drainage)

# And water chemistry data
water_chemistry <- read.csv("waterchem_with_averagedsonde_labmetals_data.csv") %>%
  filter(Site.Code != "N/A", !Date %in% c("15-Jun-22","16-Jun-22")) %>%
  mutate(day = sapply(Date, function(x) strsplit(x, split = "-")[[1]][1]),
         month = sapply(Date, function(x) strsplit(x, split = "-")[[1]][2]),
         year = sapply(Date, function(x) strsplit(x, split = "-")[[1]][3]),
         dist_from_hybridstart = ifelse(Drainage == "Calnali", Lon - -98.585, ifelse(Drainage == "Conzintla", Lon - -98.62412, ifelse(Drainage == "Pochula", Lon - -98.57946, ifelse(Drainage == "Huazalingo", Lon - -98.63713, NA)))),
         Turbidity = ifelse(Turbidity_instrument == "orion", Turbidity, NA), # removing measurements from a different turbidity meter - not comparable to other measurements
         Drainage = ifelse(Site.Code %in% c("PEZM", "CAPS", "PLAZ","CALL", "TLCD"), "Calnali - Downstream", ifelse(Drainage == "Calnali", "Calnali - Upstream", Drainage))) %>%
  filter(year %in% c(22,23,24)) %>% # removing early years with more limited measurements - large number of NAs limit PCA variables
  filter(!(month == "Sep")) %>%  # removing data from the rainy season in this analysis
  left_join(site_distances)

water_chemistry <- as.data.frame(apply(water_chemistry, 2,
                                       FUN = function(x) ifelse(grepl("<",x), as.numeric(str_remove(x, "<"))/2, ifelse(x == "", NA, x)))) %>%
  #FUN = function(x) ifelse(grepl("<",x), runif(1, min = 0, max = as.numeric(str_remove(x, "<"))), ifelse(x == "", NA, x)))) %>%
  #FUN = function(x) ifelse(grepl("<",x), 0, ifelse(x == "", NA, x)))) %>%
  #FUN = function(x) ifelse(grepl("<",x), as.numeric(str_remove(x, "<")), ifelse(x == "", NA, x)))) %>%
  # ^ each of the above lines was used to test whether the choice of arbitrary value to substitute for ICP/MS measurements
  # below the detection limit (default 1/2 of detection limit) affected statistical outcomes; change the line that is not
  # commented out to recreate the results of Tables S10-12
  mutate_if(is_all_numeric,as.numeric)

# Loading site coordinate data
sites_drainages <- read.csv("all_sites_coords_PAPA_HZNP_Conz.csv")

## And now land use

# Land use stats within the entire upstream subcatchment
landuse_fullwatershed <- read.table("Subcatchment_LandUse_full.csv", sep = ",", header = T) %>%
  left_join(site_distances) %>%
  mutate(prop_misclassified_total = Misclassified_developed/Total_pixels,
         prop_misclassified_developed = Misclassified_developed/developed) %>%
  left_join(sites_drainages)

# Land use stats within a 500 m buffer around streams
landuse_500m <- read.table("Subcatchment_LandUse_500m.csv", sep = ",", header = T) %>%
  left_join( site_distances) %>%
  mutate(prop_misclassified_total = Misclassified_developed/Total_pixels,
         prop_misclassified_developed = Misclassified_developed/developed) %>%
  left_join(sites_drainages)

# Land use stats within a 100 m buffer around streams
landuse_100m <- read.table("Subcatchment_LandUse_100m.csv", sep = ",", header = T) %>%
  left_join( site_distances) %>%
  mutate(prop_misclassified_total = Misclassified_developed/Total_pixels,
         prop_misclassified_developed = Misclassified_developed/developed) %>%
  left_join(sites_drainages)

# Land use stats within a 50 m buffer around streams
landuse_50m <- read.table("Subcatchment_LandUse_50m.csv", sep = ",", header = T) %>%
  left_join( site_distances) %>%
  mutate(prop_misclassified_total = Misclassified_developed/Total_pixels,
         prop_misclassified_developed = Misclassified_developed/developed) %>%
  left_join(sites_drainages)


# Summarizeing water chem data as PCs and attaching to land use data

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


# Plotting developed land within each subcatchment for different buffer sizes
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

#### Now cumulative land use in all upstream subcatchments

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

# Combining into one mega-plot
leg <- get_legend_35(dummy_developed_plot_cum_50m)
all_developed_plots <- plot_grid(developed_plot, developed_plot_cum_full,
          developed_plot_500m, developed_plot_cum_500m,
          developed_plot_100m, developed_plot_cum_100m,
          developed_plot_50m, developed_plot_cum_50m,
          nrow = 4, labels = "AUTO")
all_developed_plots_wleg <- plot_grid(all_developed_plots, leg, nrow = 2, rel_heights = c(15, 1))
all_developed_plots_wleg
ggsave("developed_land_megaplot.pdf",
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
ggsave("herb_land_megaplot.pdf",
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


PC1_land_model <- lm(PC1 ~ corrected_cumprop_developed + 
                       corrected_prop_developed + 
                       corrected_cumprop_herb +
                       corrected_prop_herb + cum_total_km2, data = landuse_chem_full)
shapiro.test(PC1_land_model$residuals)
summary(PC1_land_model)
vif(PC1_land_model)
plot(PC1_land_model)

PC2_land_model <- lm(PC2 ~ corrected_cumprop_developed + 
                       corrected_prop_developed + 
                       corrected_cumprop_herb +
                       corrected_prop_herb + cum_total_km2, data = landuse_chem_full)
shapiro.test(PC2_land_model$residuals)
summary(PC2_land_model)
vif(PC2_land_model)
plot(PC2_land_model)

PC1_land_model_50m_buff <- lm(PC1 ~ corrected_cumprop_developed + 
                       corrected_prop_developed + 
                       corrected_cumprop_herb +
                       corrected_prop_herb + cum_total_km2, data = landuse_chem_50m)
shapiro.test(PC1_land_model_50m_buff$residuals)
summary(PC1_land_model_50m_buff)
vif(PC1_land_model_50m_buff)
plot(PC1_land_model_50m_buff)

PC2_land_model_50m_buff <- lm(PC2 ~ corrected_cumprop_developed + 
                                corrected_prop_developed + 
                                corrected_cumprop_herb +
                                corrected_prop_herb + cum_total_km2, data = landuse_chem_50m)
shapiro.test(PC2_land_model_50m_buff$residuals)
summary(PC2_land_model_50m_buff)
vif(PC2_land_model_50m_buff)
plot(PC2_land_model_50m_buff)

### Plots of association between chemical PCs and land use


ggplot(landuse_chem_full, aes(x = corrected_cumprop_developed, y = PC1, color = Drainage)) +
  theme_bw() +
  geom_point()

ggplot(landuse_chem_full, aes(x = corrected_cumprop_developed, y = PC2, color = Drainage)) +
  theme_bw() +
  geom_point()

ggplot(landuse_chem_50m, aes(x = corrected_cumprop_developed, y = PC1, color = Drainage)) +
  theme_bw() +
  geom_point()

ggplot(landuse_chem_50m, aes(x = corrected_cumprop_developed, y = PC2, color = Drainage)) +
  theme_bw() +
  geom_point()
