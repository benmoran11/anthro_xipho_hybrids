library(tidyverse)
library(ggbiplot)
library(ggforce)
library(scales)
library(cowplot)
library(factoextra)
library(car)
library(vegan)

section_colors <- c('#3C153B','#8B1E3F','#89BD9E','#EBB65C','#DB4C40')
names(section_colors) <- c("Calnali - Downstream","Calnali - Upstream", "Conzintla", "Huazalingo","Pochula")

river_colors <- c('#8B1E3F','#89BD9E','#EBB65C','#DB4C40')
names(river_colors) <- c("Calnali","Conzintla", "Huazalingo","Pochula")

section_shapes <- c(17, 16, 16, 16, 16)
names(section_colors) <- c("Calnali - Downstream","Calnali - Upstream", "Conzintla", "Huazalingo","Pochula")

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

##### Plots ######
malcol=rgb(0/255,0/255,139/255)
hetcol=rgb(65/255,105/255,225/255)
bircol=rgb(255/255,0/255,0/255)
viridis_scale = "rocket"
hex <- hue_pal()(4)
names(hex) = c("Calnali", "Conzintla", "Huazalingo", "Pochula")

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
water_chemistry <- as.data.frame(apply(water_chemistry, 2,
                                       FUN = function(x) ifelse(grepl("<",x), as.numeric(str_remove(x, "<"))/2, ifelse(x == "", NA, x)))) %>%
                                       #FUN = function(x) ifelse(grepl("<",x), 0, ifelse(x == "", NA, x)))) %>%
  mutate_if(is_all_numeric,as.numeric)

pca_ready_nonmetals <- water_chemistry %>%
  filter(!(Site.Code %in% c("TLMC", "COAC")),
         !(Drainage %in% c("Tanchachin","Claro"))) %>%
  dplyr::select(-matches("_d$")) %>%
  dplyr::select(-matches("_t$"))%>%
  dplyr::select(-c(DO_sat, Temp, fDOM_RFU, DOC_lab, N_lab, NO3, Calcium_Hardness, Color_Apparent, Color_True, Sulfite, Sulfide,
                   Turbidity_instrument, Fe_colorimeter, Cu_free_colorimeter, 
                   Cu_total_colorimeter, Ni_colorimeter, Mn_colorimeter, Weather, Observations, day, month, year, dist_from_hybridstart, dist_from_first_site_km, dist_from_prev_site_km, elevation, graphing_elevation, Site.Upstream)) %>%
  rename("fDOM" = fDOM_QSU, "Oxygen" = DO_conc, "Hardness" = Total_Hardness) %>%
  drop_na() %>%
  select_if(function(x) (any(x != 0)))

#set.seed(15)
#wc.nmds <- metaMDS(pca_ready_nonmetals[,8:18])
#stressplot(wc.nmds)
#wc.nmds.scores <- as.data.frame(scores(wc.nmds)$sites)
#nonmetals_nmds_plotting_dataset <- cbind(pca_ready_nonmetals, wc.nmds.scores)
#ggplot(nonmetals_nmds_plotting_dataset, aes(x = NMDS1, y = NMDS2, color = Drainage)) + geom_point()

wc.pca <- prcomp(pca_ready_nonmetals[,8:18], center = TRUE,scale. = TRUE)
summary(wc.pca)
nonmetals_pca_vars <- get_pca_var(wc.pca)
plot(wc.pca$sdev^2 / sum(wc.pca$sdev^2), xlab="Principal component", ylab="Proportion of variance explained", type='b')

nonmetal_pca <- fviz_pca_biplot(wc.pca, axes = c(1,2), repel = T, label = "var", alpha.var = 0.20, col.var = "black", mean.point = FALSE, habillage = pca_ready_nonmetals$Drainage, title = NULL)  +
  theme(legend.position = "none",
        legend.position.inside = c(.35, .2),
        legend.background = element_rect()) +
  scale_color_manual(values = section_colors) +
  scale_fill_manual(values = section_colors) +
  scale_shape_manual(values = section_shapes) +
  labs(x = bquote("PC 1 ("*.(round(summary(wc.pca)$importance[2,1],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(wc.pca)$importance[2,2],3)*100)*"%)"))
nonmetal_pca
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/pca_nonmetals.pdf",
       nonmetal_pca, device = "pdf", width = 4 * (42.5/28.7), height = 4)

dummy_nonmetal_pca <- fviz_pca_biplot(wc.pca, axes = c(1,2), repel = T, label = "var", alpha.var = 0.20, col.var = "black", mean.point = FALSE, habillage = pca_ready_nonmetals$Drainage, title = NULL)  +
  theme(legend.position = "inside",
        legend.position.inside = c(.5, .5),
        legend.background = element_rect()) +
  scale_color_manual(values = section_colors) +
  scale_shape_manual(values = section_shapes) +
  labs(x = bquote("PC 1 ("*.(round(summary(wc.pca)$importance[2,1],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(wc.pca)$importance[2,2],3)*100)*"%)"))
leg <- get_legend_35(dummy_nonmetal_pca)
leg_plot <- ggdraw(leg)
leg_plot
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/pca_legend.pdf",
       leg_plot, device = "pdf", width = 2, height = 2)


pca_ready_metals <- water_chemistry %>%
  filter(!(Site.Code %in% c("TLMC", "COAC")),
         !(Drainage %in% c("Tanchachin","Claro"))) %>%
  dplyr::select(c(Date, Time, Site.Name, Drainage, Site.Code, Lon, Lat, matches("_t$")), month) %>%
  dplyr::select(-c(Te_t, Th_t, W_t, Zr_t, Bi_t, Sn_t, Tl_t, Ag_t, Be_t)) %>% 
  rename_with(function(x) gsub("_t", "", x)) %>% # Removing metals never observed or only observed once
  drop_na() %>%
  dplyr::select_if(function(x) (any(x != 0)))


wc.pca.metals <- prcomp(pca_ready_metals[,8:37], center = TRUE,scale. = TRUE)
summary(wc.pca.metals)
metals_contrib <- get_pca_var(wc.pca.metals)$contrib
plot(wc.pca.metals$sdev^2 / sum(wc.pca.metals$sdev^2), xlab="Principal component", ylab="Proportion of variance explained", type='b')

metal_pca <- fviz_pca_biplot(wc.pca.metals, axes = c(1,2), repel = T, label = "var", alpha.var = 0.20, mean.point = FALSE, col.var = "black", habillage = pca_ready_metals$Drainage, title = NULL)+
  theme(legend.position = "none",
        legend.position.inside = c(.13, .83),
        legend.background = element_rect()) + 
  scale_color_manual(values = section_colors) +
  scale_shape_manual(values = section_shapes) +
  scale_x_reverse() +
  labs(x = bquote("PC 1 ("*.(round(summary(wc.pca.metals)$importance[2,1],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(wc.pca.metals)$importance[2,2],3)*100)*"%)"))
metal_pca
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/pca_metals.pdf",
       metal_pca, device = "pdf", width = 4 * (30.6/22.6), height = 4)

metal_pca_supp <- fviz_pca_ind(wc.pca.metals, axes = c(1,2), repel = T, alpha.var = 0.20, geom = "point", mean.point = FALSE, col.var = "black", habillage = pca_ready_metals$Drainage, title = "")+
  theme(legend.position = "inside",
        legend.position.inside = c(.20, .87),
        legend.background = element_rect()) +
  scale_color_manual(values = section_colors) +
  scale_shape_manual(values = section_shapes) +
  scale_x_reverse() +
  geom_mark_ellipse(data = wc.pca.metals$x[c(47,62),],aes(x = PC1, y = PC2), expand = .015) +
  geom_mark_ellipse(data = wc.pca.metals$x[c(40,52),],aes(x = PC1, y = PC2), expand = .015) +
  geom_mark_ellipse(data = rbind(wc.pca.metals$x[35,],wc.pca.metals$x[35,]), aes(x = PC1, y = PC2), expand = .011) +
  geom_text(label = c("VNNO"), x = 1.25, y = 2) +
  geom_text(label = c("ATMP"), x = -0.75, y = -2.25) +
  geom_text(label = c("ZACA"), x = 2, y = -2) +
  guides(color=guide_legend(title="Drainage"), shape = guide_legend(title="Drainage")) +
  labs(x = bquote("PC 1 ("*.(round(summary(wc.pca.metals)$importance[2,1],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(wc.pca.metals)$importance[2,2],3)*100)*"%)"))
metal_pca_supp
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/pca_metals_highlight_outliers.pdf",
       metal_pca_supp, device = "pdf", width = 6.5 , height = 6.5* (24.8/30.1))
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/pca_metals_highlight_outliers.png",
       metal_pca_supp, device = "png", width = 6.5 , height = 6.5* (24.8/30.1))

##### Attempt at improved model #######
# Analyzing all variables together


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
summary(all.pca)
all_contrib <- get_pca_var(all.pca)$contrib
plot(all.pca$sdev^2 / sum(all.pca$sdev^2), xlab="Principal component", ylab="Proportion of variance explained", type='b')


#nonmetals_contrib <- get_pca_var(wc.pca)$contrib

siteavgs_pca <- fviz_pca_biplot(axes = c(1,2), all.pca, repel = T, label = "var", alpha.var = 0.20, col.var = "black", mean.point = FALSE, habillage = pca_ready_siteavgs$Drainage, title = NULL)  +
  theme(legend.position = "none",
        legend.position.inside = c(.35, .2),
        legend.background = element_rect())
siteavgs_pca

#scatter3d(x = all.pca$x[,1], y = all.pca$x[,2], z = all.pca$x[,3])

major_axes <- cbind(pca_ready_siteavgs,all.pca$x[, 1:4])

#sites_upstream <- c("CALM", "CTMX", "XOCH", "PLAZ", "AGZC", "VNNO", "CAPAC", "PILO", "CALL", "CAPS")
#sites_upstream <- c("CALM", "AGZC", "CAPAC", "PILO", "PLAZ", "CALL", "CAPS")
sites_upstream <- c("CALM", "CTMX","PLAZ","VNNO","CALL")
genetics_siteavgchemPCs <- structure_stats %>%
  mutate(Drainage = ifelse(Site.Code %in% c("PEZM", "CAPS", "CALL", "PLAZ", "TLCD"), "Calnali - Downstream", ifelse(Drainage == "Calnali", "Calnali - Upstream", as.character(Drainage)))) %>%
  #filter(Site.Code %in% c("AGZC","XOCH","MXLA","PILO","CRZN","CALM","PLAZ","CALL","CAPS","TLCD","PEZM")) %>%
  #filter(Site.Code %in% c("AGZC","PILO","CALM","PLAZ","CALL","CAPS","TLCD","PEZM")) %>%
  filter(Site.Code %in% c("CRZN", "XOCH", "AGZC", "CALL","CAPS")) %>%
  left_join(site_distances) %>%
  left_join(major_axes) %>%
  left_join(downsampled_diptests) %>%
  mutate(discrepancy = statistic.D - mean_dipstat,
         Site.Upstream = sites_upstream) %>%
  right_join(structure_stats %>% filter(Site.Code %in% sites_upstream) %>% rename(`D.Upstream` = statistic.D, `Site.Upstream` = Site.Code) %>% select(Site.Upstream, D.Upstream)) %>%
  dplyr::select(Site.Code,mean_dipstat,mean_dipstat_replace,everything())
  
# Linear model of Dip statistic as a function of the main PCs and the 1-lag (from site upstream) D statistic
# This is garbage, cutting after documenting it
PC_dip_model<-lm(mean_dipstat ~ PC1 + PC2 + D.Upstream, 
                 data = genetics_siteavgchemPCs)
shapiro.test(PC_dip_model$residuals)
summary(PC_dip_model)
#summary(step(PC_dip_model))
vif(PC_dip_model)
plot(PC_dip_model)


# Supplemental Plot of PCs

siteavgs_pca12 <- fviz_pca_ind(axes = c(1,2), all.pca, repel = T, label = "none", alpha.var = 0.20, col.var = "black", mean.point = FALSE, habillage = pca_ready_siteavgs$Drainage, title = NULL)  +
  theme(legend.position = "none",
        legend.position.inside = c(.35, .2),
        legend.background = element_rect()) +
  scale_color_manual(values = section_colors) +
  scale_shape_manual(values = section_shapes) +
  ggtitle(NULL) +
  labs(x = bquote("PC 1 ("*.(round(summary(all.pca)$importance[2,1],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(all.pca)$importance[2,2],3)*100)*"%)"))
siteavgs_pca12

siteavgs_vars12 <- fviz_pca_var(axes = c(1,2), all.pca, repel = T, label = "var", alpha.var = 0.20, col.var = "black", mean.point = FALSE, title = NULL)  +
  theme(legend.position = "none",
        legend.position.inside = c(.35, .2),
        legend.background = element_rect())  +
  ggtitle(NULL) +
  labs(x = bquote("PC 1 ("*.(round(summary(all.pca)$importance[2,1],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(all.pca)$importance[2,2],3)*100)*"%)"))
siteavgs_vars12

siteavgs_vars13 <- fviz_pca_var(axes = c(1,3), all.pca, repel = T, label = "var", alpha.var = 0.20, col.var = "black", mean.point = FALSE, title = NULL)  +
  theme(legend.position = "none",
        legend.position.inside = c(.35, .2),
        legend.background = element_rect())  +
  ggtitle(NULL) +
  labs(x = bquote("PC 1 ("*.(round(summary(all.pca)$importance[2,1],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(all.pca)$importance[2,3],3)*100)*"%)"))
siteavgs_vars13

siteavgs_pca13 <- fviz_pca_ind(axes = c(1,3), all.pca, repel = T, label = "none", alpha.var = 0.20, col.var = "black", mean.point = FALSE, habillage = pca_ready_siteavgs$Drainage, title = NULL)  +
  theme(legend.position = "inside",
        legend.position.inside = c(.27, .27),
        legend.background = element_rect()) +
  scale_color_manual(values = section_colors) +
  scale_shape_manual(values = section_shapes) +
  ggtitle(NULL) +
  labs(x = bquote("PC 1 ("*.(round(summary(all.pca)$importance[2,1],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(all.pca)$importance[2,3],3)*100)*"%)"))
siteavgs_pca13

siteavgs_pca14 <- fviz_pca_biplot(axes = c(1,4), all.pca, repel = T, label = "var", alpha.var = 0.20, col.var = "black", mean.point = FALSE, habillage = pca_ready_siteavgs$Drainage, title = NULL)  +
  theme(legend.position = "none",
        legend.position.inside = c(.35, .2),
        legend.background = element_rect()) +
  scale_color_manual(values = section_colors) +
  scale_shape_manual(values = section_shapes) +
  ggtitle(NULL) +
  labs(x = bquote("PC 1 ("*.(round(summary(all.pca)$importance[2,1],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(all.pca)$importance[2,4],3)*100)*"%)"))
siteavgs_pca14

siteavgs_pca23 <- fviz_pca_ind(axes = c(2,3), all.pca, repel = T, label = "none", alpha.var = 0.20, col.var = "black", mean.point = FALSE, habillage = pca_ready_siteavgs$Drainage, title = NULL)  +
  theme(legend.position = "none",
        legend.position.inside = c(.85, .2),
        legend.background = element_rect()) +
  scale_color_manual(values = section_colors) +
  scale_shape_manual(values = section_shapes) +
  ggtitle(NULL) +
  labs(x = bquote("PC 1 ("*.(round(summary(all.pca)$importance[2,2],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(all.pca)$importance[2,3],3)*100)*"%)"))
siteavgs_pca23

siteavgs_pca24 <- fviz_pca_biplot(axes = c(2,4), all.pca, repel = T, label = "var", alpha.var = 0.20, col.var = "black", mean.point = FALSE, habillage = pca_ready_siteavgs$Drainage, title = NULL)  +
  theme(legend.position = "inside",
        legend.position.inside = c(.5, .5),
        legend.background = element_rect()) +
  scale_color_manual(values = section_colors) +
  scale_shape_manual(values = section_shapes) +
  ggtitle(NULL) +
  labs(x = bquote("PC 1 ("*.(round(summary(all.pca)$importance[2,2],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(all.pca)$importance[2,4],3)*100)*"%)"))
siteavgs_pca24

leg <- get_legend(siteavgs_pca24)

legbox <- ggdraw(leg)

siteavgs_pca34 <- fviz_pca_biplot(axes = c(3,4), all.pca, repel = T, label = "var", alpha.var = 0.20, col.var = "black", mean.point = FALSE, habillage = pca_ready_siteavgs$Drainage, title = NULL)  +
  theme(legend.position = "none",
        legend.position.inside = c(.35, .2),
        legend.background = element_rect()) +
  scale_color_manual(values = section_colors) +
  scale_shape_manual(values = section_shapes) +
  ggtitle(NULL) +
  labs(x = bquote("PC 1 ("*.(round(summary(all.pca)$importance[2,3],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(all.pca)$importance[2,4],3)*100)*"%)"))
siteavgs_pca34

all_siteavged_PCs <- plot_grid(siteavgs_pca12,   siteavgs_pca13, siteavgs_pca23, siteavgs_vars12, nrow = 2, labels = "AUTO")
all_siteavged_PCs
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/all_siteavged_PCs.pdf",
       all_siteavged_PCs, width = 6.5, height = 6.5)
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/all_siteavged_PCs.png",
       all_siteavged_PCs, width = 6.5, height = 6.5)

# Various plots to visualize relationships in data
ggplot(genetics_siteavgchemPCs, aes(x = statistic.D, y = mean_dipstat)) + geom_point() + geom_text(aes(label = Site.Code)) + geom_abline(slope = 1, intercept = 0)
ggplot(genetics_siteavgchemPCs, aes(x = log(nobs), y = discrepancy)) + geom_point() + geom_text(aes(label = Site.Code)) + geom_abline(slope = 1, intercept = 0)
ggplot(genetics_siteavgchemPCs, aes(x = statistic.D, y = D.Upstream)) + geom_point() + geom_text(aes(label = Site.Code)) + geom_abline(slope = 1, intercept = 0)
ggplot(genetics_siteavgchemPCs, aes(x = PC1, y = PC2, shape = Drainage, color = mean_dipstat)) +
  geom_point()+
  geom_text(aes(label = Site.Code))
ggplot(genetics_siteavgchemPCs, aes(x = PC1, y = PC3, shape = Drainage, color = mean_dipstat)) +
  geom_point()+
  geom_text(aes(label = Site.Code))
ggplot(genetics_siteavgchemPCs, aes(x = PC2, y = PC3, shape = Drainage, color = mean_dipstat)) +
  geom_point()+
  geom_text(aes(label = Site.Code))
ggplot(genetics_siteavgchemPCs, aes(x = PC3, y = mean_dipstat, shape = Drainage)) +
  geom_point()#+
  geom_text(aes(label = Site.Code))
  


####### Trying to correct for spatial autocorrelation #######

site_distance_matrix <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/site_distance_matrix.csv", row.names = 1) %>%
  as.matrix()

for(i in 1:nrow(site_distance_matrix)) {
  for(j in 1:i) {
    site_distance_matrix[j,i] = site_distance_matrix[i,j]
  }
}

mds_site_distance <- cmdscale(site_distance_matrix, eig = TRUE, k = 2)
df_site_distances <- data.frame(Site.Code = row.names(site_distance_matrix),
                                x = mds_site_distance$points[,1],
                                y = mds_site_distance$points[,2])

ggplot(df_site_distances, aes(x = x, y = y)) + geom_point() + geom_text(aes(label = Site.Code))




###### Redo with different substitutions for detection limit

# Substituting 0
water_chemistry_sub0 <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/waterchem_with_averagedsonde_labmetals_data.csv") %>%
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
water_chemistry_sub0 <- as.data.frame(apply(water_chemistry_sub0, 2,
                                       FUN = function(x) ifelse(grepl("<",x), 0, ifelse(x == "", NA, x)))) %>%
  mutate_if(is_all_numeric,as.numeric)


pca_ready_metals_sub0 <- water_chemistry_sub0 %>%
  filter(!(Site.Code %in% c("TLMC", "COAC")),
         !(Drainage %in% c("Tanchachin","Claro"))) %>%
  dplyr::select(c(Date, Time, Site.Name, Drainage, Site.Code, Lon, Lat, matches("_t$")), month) %>%
  dplyr::select(-c(Te_t, Th_t, W_t, Zr_t, Bi_t, Sn_t, Tl_t, Ag_t, Be_t)) %>% 
  rename_with(function(x) gsub("_t", "", x)) %>% # Removing metals never observed or only observed once
  drop_na() %>%
  dplyr::select_if(function(x) (any(x != 0)))

wc.pca.metals_sub0 <- prcomp(pca_ready_metals_sub0[,8:37], center = TRUE,scale. = TRUE)
summary(wc.pca.metals_sub0)
metals_contrib_sub0 <- get_pca_var(wc.pca.metals_sub0)$contrib
plot(wc.pca.metals_sub0$sdev^2 / sum(wc.pca.metals_sub0$sdev^2), xlab="Principal component", ylab="Proportion of variance explained", type='b')

metal_pca_sub0 <- fviz_pca_biplot(wc.pca.metals_sub0, axes = c(1,2), repel = T, label = "var", alpha.var = 0.20, mean.point = FALSE, col.var = "black", habillage = pca_ready_metals$Drainage, title = NULL)+
  theme(legend.position = "none",
        legend.position.inside = c(.13, .83),
        legend.background = element_rect()) + scale_x_reverse() +
  scale_color_manual(values = section_colors) +
  scale_shape_manual(values = section_shapes) +
  labs(x = bquote("PC 1 ("*.(round(summary(wc.pca.metals_sub0)$importance[2,1],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(wc.pca.metals_sub0)$importance[2,2],3)*100)*"%)"))
metal_pca_sub0
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/pca_metals_sub0.pdf",
       metal_pca_sub0, device = "pdf", width = 4 * (30.6/22.6), height = 4)



#Substituting the detection limit
water_chemistry_subDL <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/waterchem_with_averagedsonde_labmetals_data.csv") %>%
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
water_chemistry_subDL <- as.data.frame(apply(water_chemistry_subDL, 2,
                                            FUN = function(x) ifelse(grepl("<",x), as.numeric(str_remove(x, "<")), ifelse(x == "", NA, x)))) %>%
  mutate_if(is_all_numeric,as.numeric)


pca_ready_metals_subDL <- water_chemistry_subDL %>%
  filter(!(Site.Code %in% c("TLMC", "COAC")),
         !(Drainage %in% c("Tanchachin","Claro"))) %>%
  dplyr::select(c(Date, Time, Site.Name, Drainage, Site.Code, Lon, Lat, matches("_t$")), month) %>%
  dplyr::select(-c(Te_t, Th_t, W_t, Zr_t, Bi_t, Sn_t, Tl_t, Ag_t, Be_t)) %>% 
  rename_with(function(x) gsub("_t", "", x)) %>% # Removing metals never observed or only observed once
  drop_na() %>%
  dplyr::select_if(function(x) (any(x != 0)))

wc.pca.metals_subDL <- prcomp(pca_ready_metals_subDL[,8:37], center = TRUE,scale. = TRUE)
summary(wc.pca.metals_subDL)
metals_contrib_subDL <- get_pca_var(wc.pca.metals_subDL)$contrib
plot(wc.pca.metals_subDL$sdev^2 / sum(wc.pca.metals_subDL$sdev^2), xlab="Principal component", ylab="Proportion of variance explained", type='b')

metal_pca_subDL <- fviz_pca_biplot(wc.pca.metals_subDL, axes = c(1,2), repel = T, label = "var", alpha.var = 0.20, mean.point = FALSE, col.var = "black", habillage = pca_ready_metals$Drainage, title = NULL)+
  theme(legend.position = "inside",
        legend.position.inside = c(.17, .83),
        legend.background = element_rect()) + scale_x_reverse() +
  scale_color_manual(values = section_colors) +
  scale_shape_manual(values = section_shapes) +
  labs(x = bquote("PC 1 ("*.(round(summary(wc.pca.metals_subDL)$importance[2,1],3)*100)*"%)"), y = bquote("PC 2 ("*.(round(summary(wc.pca.metals_subDL)$importance[2,2],3)*100)*"%)"))
metal_pca_subDL
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/pca_metals_subDL.pdf",
       metal_pca_subDL, device = "pdf", width = 4 * (30.6/22.6), height = 4)
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/pca_metals_subDL.png",
       metal_pca_subDL, device = "png", width = 4 * (30.6/22.6), height = 4)

vary_DLsub_PCAs <- plot_grid(metal_pca_sub0, metal_pca_subDL, nrow = 2, ncol = 1, labels = c("A","B"))
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/pca_metals_sub0_subDL.pdf",
       vary_DLsub_PCAs, device = "pdf", width = 6.5, height = 8)
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/pca_metals_sub0_subDL.png",
       vary_DLsub_PCAs, device = "png", width = 6.5, height = 8)



