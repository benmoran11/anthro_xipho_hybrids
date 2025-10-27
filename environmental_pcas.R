library(tidyverse)
library(ggbiplot)
library(ggforce)
library(scales)
library(cowplot)
library(factoextra)
library(car)
library(vegan)

# Colors to be used for different drainages - splitting the Calnali drainage upstream and downstream of town
section_colors <- c('#3C153B','#8B1E3F','#89BD9E','#EBB65C','#DB4C40')
names(section_colors) <- c("Calnali - Downstream","Calnali - Upstream", "Conzintla", "Huazalingo","Pochula")

# Colors to be used for different drainages - keeping Calnali together
river_colors <- c('#8B1E3F','#89BD9E','#EBB65C','#DB4C40')
names(river_colors) <- c("Calnali","Conzintla", "Huazalingo","Pochula")

# Shapes to be used for different drainages - splitting the Calnali drainage upstream and downstream of town
section_shapes <- c(17, 16, 16, 16, 16)
names(section_colors) <- c("Calnali - Downstream","Calnali - Upstream", "Conzintla", "Huazalingo","Pochula")

# Function to be used later in extracting legend from dummy plot for multi-panel figure
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

# Utility to check for numeric values - useful for parsing cases where numerics  
# not initially parsed successfully because of "<" signs present in ICP/MS data
is_all_numeric <- function(x) {
  !any(is.na(suppressWarnings(as.numeric(na.omit(x))))) & is.character(x)
}

# Loading geographic data
site_distances <- read.csv("site_riverkm_distances.csv") %>%
  dplyr::select(-Drainage)

# And water chemistry data
water_chemistry <- read.csv("waterchem_with_averagedsonde_labmetals_data.csv") %>%
  filter(Site.Code != "N/A", !Date %in% c("15-Jun-22","16-Jun-22")) %>%     # Removing measurements from within the onset of the rainy season - not representative
  mutate(day = sapply(Date, function(x) strsplit(x, split = "-")[[1]][1]),
         month = sapply(Date, function(x) strsplit(x, split = "-")[[1]][2]),
         year = sapply(Date, function(x) strsplit(x, split = "-")[[1]][3]),
         dist_from_hybridstart = ifelse(Drainage == "Calnali", Lon - -98.585, ifelse(Drainage == "Conzintla", Lon - -98.62412, ifelse(Drainage == "Pochula", Lon - -98.57946, ifelse(Drainage == "Huazalingo", Lon - -98.63713, NA)))),
         Turbidity = ifelse(Turbidity_instrument == "orion", Turbidity, NA),    # removing measurements from a different turbidity meter - not comparable to other measurements
         Drainage = ifelse(Site.Code %in% c("PEZM", "CAPS", "PLAZ","CALL", "TLCD"), "Calnali - Downstream", ifelse(Drainage == "Calnali", "Calnali - Upstream", Drainage))) %>%
  filter(year %in% c(22,23,24)) %>%  # removing early years with more limited measurements - large number of NAs limit PCA variables
  filter(!(month == "Sep")) %>%  # removing data from the rainy season in this analysis
  left_join(site_distances)
water_chemistry <- as.data.frame(apply(water_chemistry, 2,
                                       FUN = function(x) ifelse(grepl("<",x), as.numeric(str_remove(x, "<"))/2, ifelse(x == "", NA, x)))) %>%     # ICP/MS measurements below detection limit were coded as 1/2 the detection limit, by default
  mutate_if(is_all_numeric,as.numeric)

pca_ready_nonmetals <- water_chemistry %>%
  filter(!(Site.Code %in% c("TLMC", "COAC")),             # Removing measurements from other watersheds
         !(Drainage %in% c("Tanchachin","Claro"))) %>%.   # Removing measurements from other watersheds
  dplyr::select(-matches("_d$")) %>%
  dplyr::select(-matches("_t$"))%>%
  dplyr::select(-c(DO_sat, Temp, fDOM_RFU, DOC_lab, N_lab, NO3, Calcium_Hardness, Color_Apparent, Color_True, Sulfite, Sulfide,
                   Turbidity_instrument, Fe_colorimeter, Cu_free_colorimeter, 
                   Cu_total_colorimeter, Ni_colorimeter, Mn_colorimeter, Weather, Observations, day, month, year, dist_from_hybridstart, dist_from_first_site_km, dist_from_prev_site_km, elevation, graphing_elevation, Site.Upstream)) %>%
  rename("fDOM" = fDOM_QSU, "Oxygen" = DO_conc, "Hardness" = Total_Hardness) %>%
  drop_na() %>%
  select_if(function(x) (any(x != 0)))


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
ggsave("pca_nonmetals.pdf",
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
ggsave("pca_legend.pdf",
       leg_plot, device = "pdf", width = 2, height = 2)


pca_ready_metals <- water_chemistry %>%
  filter(!(Site.Code %in% c("TLMC", "COAC")),            # Removing measurements from other watersheds
         !(Drainage %in% c("Tanchachin","Claro"))) %>%   # Removing measurements from other watersheds
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
ggsave("pca_metals.pdf",
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
ggsave("pca_metals_highlight_outliers.pdf",
       metal_pca_supp, device = "pdf", width = 6.5 , height = 6.5* (24.8/30.1))
ggsave("pca_metals_highlight_outliers.png",
       metal_pca_supp, device = "png", width = 6.5 , height = 6.5* (24.8/30.1))


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
ggsave("all_siteavged_PCs.pdf",
       all_siteavged_PCs, width = 6.5, height = 6.5)
ggsave("all_siteavged_PCs.png",
       all_siteavged_PCs, width = 6.5, height = 6.5)


###### Redo with different substitutions for detection limit (Fig. S4)

# Substituting 0
water_chemistry_sub0 <- read.csv("waterchem_with_averagedsonde_labmetals_data.csv") %>%
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
ggsave("pca_metals_sub0.pdf",
       metal_pca_sub0, device = "pdf", width = 4 * (30.6/22.6), height = 4)



#Substituting the detection limit
water_chemistry_subDL <- read.csv("waterchem_with_averagedsonde_labmetals_data.csv") %>%
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
ggsave("pca_metals_subDL.pdf",
       metal_pca_subDL, device = "pdf", width = 4 * (30.6/22.6), height = 4)
ggsave("pca_metals_subDL.png",
       metal_pca_subDL, device = "png", width = 4 * (30.6/22.6), height = 4)

vary_DLsub_PCAs <- plot_grid(metal_pca_sub0, metal_pca_subDL, nrow = 2, ncol = 1, labels = c("A","B"))
ggsave("pca_metals_sub0_subDL.pdf",
       vary_DLsub_PCAs, device = "pdf", width = 6.5, height = 8)
ggsave("pca_metals_sub0_subDL.png",
       vary_DLsub_PCAs, device = "png", width = 6.5, height = 8)



