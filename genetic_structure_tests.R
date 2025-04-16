library(tidyverse)
library(diptest)
library(cowplot)
library(grid)
library(gridExtra) 
library(ggpubr)
library(moments)


##### Plots ######
malcol=rgb(0/255,0/255,139/255)
hetcol=rgb(65/255,105/255,225/255)
bircol=rgb(255/255,0/255,0/255)

site_distances <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/site_riverkm_distances_noduplicates.csv")

clines <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/huazalingo_pochula_calnali_conzintla_allhybridindex.csv") %>%
  mutate(Site_order = as.factor(Site_order)) %>%
  filter(Site.Code != "LESP") %>%
  left_join(site_distances) %>%
  filter(malcount + bircount > 10000, Collection_Year > 2004)
clines <- left_join(clines, clines %>% group_by(Site.Code, Drainage) %>% dplyr::summarize(Site_n = dplyr::n())) %>%
  group_by(Site.Code) %>%
  arrange(hybrid_index) %>%
  mutate(anc_index = row_number())
clines$ancgroup <- cut(clines$hybrid_index, breaks = seq(0,1, 0.04))
clines$Drainage <- factor(clines$Drainage, levels = c("Huazalingo", "Calnali", "Pochula", "Conzintla"))
hybridramp <- colorRampPalette(c(bircol,hetcol))(25)

sample_sizes <- clines %>% group_by(Drainage, Site.Code, Site_order) %>% summarize(n = n())
median(sample_sizes$n)

ancestry_diptests <- lapply(unique(clines$Site.Code), function(pop) {
  test <- dip.test(filter(clines, Site.Code == pop)$hybrid_index)
  return(c("Site.Code"=pop,"statistic"=test$statistic,"p.value"=test$p.value,"nobs"=test$nobs))
}) %>%
  bind_rows() %>%
  mutate(statistic.D = as.numeric(statistic.D),
         p.value = as.numeric(p.value),
         nobs = as.numeric(nobs),
         bonferroni.sig = p.value < 0.05/46) ### 46 = number of dip tests performed in study; 43 on whole sites (-La Esperanza) + 3 for CRZN annual subsets
cor.test(ancestry_diptests$nobs, ancestry_diptests$statistic.D)

# For use in later downstream stats, going to downsample to the lowest sample size in sites with signs of sympatry
overlap_sites <- c("CALM","AGZC", "PILO", "PLAZ", "CALL","CAPS", "PEZM", "XOCH","MXLA", "CRZN", "TULA", "TOTO")
minobs <- min(as.numeric(filter(ancestry_diptests, Site.Code %in% overlap_sites)$nobs))

set.seed(47)
downsampled_diptests <- lapply(overlap_sites, function(pop) {
  rep_dipstat_noreplace <- lapply(1:1000, function(r) {
    subset <- slice_sample(filter(clines, Site.Code == pop), n = minobs, replace = F)
    test <- dip.test(subset$hybrid_index)$statistic
    pval <- dip.test(subset$hybrid_index)$p
    list(test = test, pval = pval)
  })
  rep_dipstat_replace <- lapply(1:1000, function(r) {
    subset <- slice_sample(filter(clines, Site.Code == pop), n = minobs, replace = T)
    test <- dip.test(subset$hybrid_index)$statistic
    pval <- dip.test(subset$hybrid_index)$p
    list(test = test, pval = pval)
  })
  return(c("Site.Code"=pop,"mean_dipstat"=mean(unlist(do.call(rbind,rep_dipstat_noreplace)[,1])), "prop_sig" = sum(do.call(rbind,rep_dipstat_noreplace)[,2] < 0.05) / 1000, "prop_sig_bonf" = sum(do.call(rbind,rep_dipstat_noreplace)[,2] < (0.05 / 46)) / 1000, "mean_dipstat_replace" = mean(unlist(do.call(rbind,rep_dipstat_replace)[,1])), "prop_sig_replace" = sum(do.call(rbind,rep_dipstat_replace)[,2] < 0.05) / 1000))
}) %>%
  bind_rows() %>%
  mutate(mean_dipstat = as.numeric(mean_dipstat),
         mean_dipstat_replace = as.numeric(mean_dipstat_replace))
downsampled_diptests_wn <- downsampled_diptests %>% left_join(sample_sizes)


set.seed(47)
downsampled_diptests_33 <- lapply(overlap_sites, function(pop) {
  rep_dipstat_noreplace <- lapply(1:1000, function(r) {
    subset <- slice_sample(filter(clines, Site.Code == pop), n = 33, replace = F)
    test <- dip.test(subset$hybrid_index)$statistic
    pval <- dip.test(subset$hybrid_index)$p
    list(test = test, pval = pval)
  })
  rep_dipstat_replace <- lapply(1:1000, function(r) {
    subset <- slice_sample(filter(clines, Site.Code == pop), n = 33, replace = T)
    test <- dip.test(subset$hybrid_index)$statistic
    pval <- dip.test(subset$hybrid_index)$p
    list(test = test, pval = pval)
  })
  return(c("Site.Code"=pop,"mean_dipstat"=mean(unlist(do.call(rbind,rep_dipstat_noreplace)[,1])), "prop_sig" = sum(do.call(rbind,rep_dipstat_noreplace)[,2] < 0.05) / 1000, "mean_dipstat_replace" = mean(unlist(do.call(rbind,rep_dipstat_replace)[,1])), "prop_sig_replace" = sum(do.call(rbind,rep_dipstat_replace)[,2] < 0.05) / 1000))
}) %>%
  bind_rows() %>%
  mutate(mean_dipstat = as.numeric(mean_dipstat),
         mean_dipstat_replace = as.numeric(mean_dipstat_replace))
downsampled_diptests_33_wn <- downsampled_diptests %>% left_join(sample_sizes)



ancstats <- clines %>%
  group_by(Site.Code, Drainage) %>%
  summarize(mean_anc = mean(hybrid_index),
            sd_anc = sd(hybrid_index),
            var_anc = var(hybrid_index),
            skew_anc = skewness(hybrid_index),
            kurt_anc = kurtosis(hybrid_index),
            mean_het = mean(heterzygosity),
            sd_het = sd(heterzygosity),
            nobs_het = n())

structure_stats <- left_join(ancestry_diptests, ancstats) %>%
  left_join(distinct(dplyr::select(clines, Site.Code, Site_order, Drainage)))

ggplot(structure_stats, aes(x = mean_anc, y = mean_het, color = statistic.D)) +
  facet_wrap(~Drainage) +
  geom_point()


crzn <- filter(clines, Site.Code == "CRZN", malcount + bircount > 10000) 

crzn_diptests <- lapply(unique(crzn$Collection_Year), function(yr) {
  test <- dip.test(filter(crzn, Collection_Year == yr)$hybrid_index)
  return(c("Collection_Year"=yr,"statistic"=test$statistic,"p.value"=test$p.value,"nobs"=test$nobs))
}) %>%
  bind_rows() %>%
  mutate(bonferroni.sig = p.value < 0.05/46)

conzintla_minor <- clines %>%
  filter(Site.Code %in% c("XOCH", "MXLA")) %>%
  mutate(min_par_anc = ifelse(hybrid_index > 1 - hybrid_index, 1 - hybrid_index, hybrid_index)) %>%
  group_by(Site.Code) %>%
  summarize(mean_minor = mean(min_par_anc))

######## Genetics Plots #############

all_clines <- ggplot(filter(clines, malcount + bircount > 10000) %>% group_by(Drainage, Site_order) %>% mutate(mean = mean(hybrid_index)), aes(y = hybrid_index, fill = ancgroup)) +
  theme_cowplot()+
  coord_cartesian(y = c(0,1)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none") +
  facet_grid(Drainage ~ Site_order, scales = "free") +
  scale_fill_discrete(type = hybridramp) +
  geom_histogram() +
  geom_hline(aes(yintercept = mean))
all_clines

pochula_cline <- ggplot(filter(clines, malcount + bircount > 10000, Drainage == "Pochula"), aes(y = hybrid_index, fill = ancgroup)) +
  theme_cowplot()+
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        axis.line.x.bottom = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") +
  coord_cartesian(y = c(0, 1)) +
  facet_grid(Drainage ~ Site_order, scales = "free") +
  scale_fill_discrete(type = hybridramp, drop = F) +
  geom_histogram()
pochula_cline
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/Pochula_cline.png",
       pochula_cline, device = "png", width = 11, height = 2)

huazalingo_cline <- ggplot(filter(clines, malcount + bircount > 10000, Drainage == "Huazalingo"),
                           aes(y = hybrid_index, fill = ancgroup)) +
  theme_cowplot()+
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line.x.bottom = element_blank(),
        legend.position = "none") +
  facet_grid(Drainage ~ Site_order, scales = "free") +
  scale_fill_discrete(type = hybridramp, drop = F) +
  geom_histogram()
huazalingo_cline
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/Huazalingo_cline.png",
       huazalingo_cline, device = "png", width = 11, height = 2)

crzn_cline <- ggplot(filter(crzn, malcount + bircount > 10000),
                     aes(y = hybrid_index, fill = ancgroup)) +
  coord_cartesian(y = c(0,1)) +
  theme_cowplot()+
  labs(y = expression("Proportion"~italic("X. malinche")~"Ancestry")) +
  theme(legend.position = "none") +
  facet_grid(. ~ Collection_Year, scales = "free") +
  scale_fill_discrete(type = hybridramp, drop = F) +
  geom_histogram()
crzn_cline
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/Corazon_cline_byyear.png",
       crzn_cline, device = "png", width = 7, height = 4)

calnali_cline <- ggplot(filter(clines, malcount + bircount > 10000, Drainage == "Calnali"),
                        aes(y = hybrid_index, fill = ancgroup)) +
  theme_cowplot()+
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line.x.bottom = element_blank(),
        legend.position = "none") +
  facet_grid(.~ Site_order, scales = "free", drop = T) +
  scale_fill_discrete(type = hybridramp, drop = F) +
  geom_histogram()
calnali_cline
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/Calnali_cline.png",
       calnali_cline, device = "png", width = 11, height = 2)

conzintla_cline <- ggplot(filter(clines, malcount + bircount > 10000, Drainage == "Conzintla"),
                          aes(y = hybrid_index, fill = ancgroup)) +
  theme_cowplot()+
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x.bottom = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none") +
  facet_grid(. ~ Site_order, scales = "free", drop = T) +
  scale_fill_discrete(type = hybridramp, drop = F) +
  geom_histogram()
conzintla_cline
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/Conzintla_cline.png",
       conzintla_cline, device = "png", width = 11, height = 2)

clines_confluence_managed <- clines %>%
  mutate(Drainage = ifelse(Site.Code %in% c("PAPA","HZNP"), "Pochula", as.character(Drainage)))
clines_confluence_managed$Drainage <- factor(clines_confluence_managed$Drainage, levels = c("Huazalingo","Calnali","Pochula", "Conzintla"))

viridis_scale = "magma"

all_spectrum <- ggplot(arrange(filter(clines_confluence_managed, malcount + bircount > 10000, Site.Code != "CULH"), Drainage, Site_order, hybrid_index), 
                       aes(y = dist_from_first_site_km, x = 1 / (Site_n) * anc_index - 1/2 / (Site_n), width = 1/Site_n, height = .5, fill = hybrid_index)) +
  facet_grid( ~ Drainage) + 
  theme_cowplot()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(.8, .2))+
  lims(y = c(0, 50), x = c(0, 1)) + 
  #  facet_grid(Drainage ~ Site_order, scales = "free") +
  scale_fill_viridis_c(option = viridis_scale, direction = -1, limits = c(0,1)) +
  scale_y_reverse() +
  labs(x = "", y = "Distance from First Site (km)") +
  #  scale_fill_discrete(type = hybridramp, drop = F) +
  geom_tile()
all_spectrum
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/all_spectra_bydistance_vertical_draft.pdf",
       all_spectrum, device = "pdf", width = 6.5, height = 8)


all_spectrum_elev <- ggplot(arrange(filter(clines, malcount + bircount > 10000), Drainage, Site_order, hybrid_index), 
                       aes(x = graphing_elevation, y = 1 / (Site_n) * anc_index - 1/2 / (Site_n), height = 1/Site_n, width = 25, fill = hybrid_index)) +
  facet_grid(Drainage ~ .) + 
  theme_cowplot()+
  theme(strip.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(size = 24),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        legend.position.inside = c(.08,.8)) +
  scale_x_reverse() +
  coord_cartesian(y = c(0, 1)) +
  scale_fill_viridis_c(option = viridis_scale, direction = -1, limits = c(0,1)) +
  labs(y = "Hybrid Ancestry", x = "Elevation (m)") +
  #  scale_fill_discrete(type = hybridramp, drop = F) +
  geom_tile()
all_spectrum_elev
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/Spectrum_cline_elevation.png",
       all_spectrum_elev, device = "png", width = 11, height = 7)

all_spectrum_elev_vert <- ggplot(arrange(filter(clines_confluence_managed, malcount + bircount > 10000, Site.Code != "CULH"), Drainage, Site_order, hybrid_index), 
                            aes(y = graphing_elevation, x = 1 / (Site_n) * anc_index - 1/2 / (Site_n), width = 1/Site_n, height = 25, fill = hybrid_index)) +
  facet_grid( ~ Drainage) + 
  theme_cowplot()+
  theme(#strip.text.y = element_blank(),
    #strip.background = element_blank(),
    #strip.text = element_blank(),
    #axis.text.y = element_text(size = 24),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(.3,.2)) +
  coord_cartesian(x = c(0, 1)) +
  scale_fill_viridis_c(option = viridis_scale, direction = -1, limits = c(0,1), name = expression(atop("Proportion",italic("X.malinche")))) +
  labs(x = "", y = "Elevation (m)") +
  #  scale_fill_discrete(type = hybridramp, drop = F) + 
  geom_tile()
all_spectrum_elev_vert
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/Spectrum_cline_elevation_vertical_draft.pdf",
       all_spectrum_elev_vert, device = "pdf", width = 6.5, height = 7)

dummy <- ggplot(arrange(filter(clines, malcount + bircount > 10000), Drainage, Site_order, hybrid_index), 
                            aes(x = graphing_elevation, y = 1 / (Site_n) * anc_index - 1/2 / (Site_n), height = 1/Site_n, width = 25, fill = hybrid_index)) +
  theme_void() +
  theme(legend.position = 'left',
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent",
                                         colour = NA_character_),
        legend.box.background = element_rect(fill = "transparent",
                                             colour = NA_character_),
        legend.key = element_rect(fill = "transparent",
                                  colour = NA_character_)) +
  scale_x_reverse() +
  coord_cartesian(y = c(0, 1)) +
  scale_fill_viridis_c(option = viridis_scale, direction = -1, limits = c(0,1)) +
  labs(y = "Hybrid Ancestry", x = "Elevation (m)") +
  geom_tile()
dummy

just_legend <- get_legend(dummy)
as_ggplot(just_legend)
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/Spectrum_cline_elevation_just_legend.pdf",
        just_legend, device = "pdf", width = 1, height = 1.5)

pochula_spectrum <-  ggplot(arrange(filter(clines, malcount + bircount > 10000, Drainage == "Pochula"), Site_order, hybrid_index), 
                       aes(x = dist_from_first_site_km, y = 1 / (Site_n) * anc_index - 1/2 / (Site_n), height = 1/Site_n, width = 0.5, fill = hybrid_index)) +
  theme_cowplot()+
  coord_cartesian(y = c(0, 1)) +
  #  facet_grid(Drainage ~ Site_order, scales = "free") +
  scale_fill_viridis_c(option = viridis_scale, direction = -1, limits = c(0,1)) +
  labs(y = "Hybrid Ancestry", x = "Distancee from First Site (km)") +
  #  scale_fill_discrete(type = hybridramp, drop = F) +
  geom_tile()
pochula_spectrum


huazalingo_spectrum <- ggplot(arrange(filter(clines, malcount + bircount > 10000, Drainage == "Huazalingo"), Site_order, hybrid_index), 
                              aes(x = dist_from_first_site_km, y = 1 / (Site_n) * anc_index - 1/2 / (Site_n), height = 1/Site_n, width = 0.5, fill = hybrid_index)) +
  theme_cowplot()+
  coord_cartesian(y = c(0, 1)) +
  #  facet_grid(Drainage ~ Site_order, scales = "free") +
  labs(y = "Hybrid Ancestry", x = "Distancee from First Site (km)") +
  scale_fill_viridis_c(option = viridis_scale, direction = -1, values = c(0,1)) +
  #  scale_fill_discrete(type = hybridramp, drop = F) +
  geom_tile()
huazalingo_spectrum

calnali_spectrum <- ggplot(arrange(filter(clines, malcount + bircount > 10000, Drainage == "Calnali"), Site_order, hybrid_index), 
                           aes(x = dist_from_first_site_km, y = 1 / (Site_n) * anc_index - 1/2 / (Site_n), height = 1/Site_n, width = 0.5, fill = hybrid_index)) +
  theme_cowplot()+
  coord_cartesian(y = c(0, 1)) +
  #  facet_grid(Drainage ~ Site_order, scales = "free") +
  scale_fill_viridis_c(option = viridis_scale, direction = -1, limits = c(0,1)) +
  labs(y = "Hybrid Ancestry", x = "Distancee from First Site (km)") +
  #  scale_fill_discrete(type = hybridramp, drop = F) +
  geom_tile()
calnali_spectrum

conzintla_spectrum <- ggplot(arrange(filter(clines, malcount + bircount > 10000, Drainage == "Conzintla"), Site_order, hybrid_index), 
                             aes(x = dist_from_first_site_km, y = 1 / (Site_n) * anc_index - 1/2 / (Site_n), height = 1/Site_n, width = 0.5, fill = hybrid_index)) +
  theme_cowplot()+
  coord_cartesian(y = c(0, 1)) +
  #  facet_grid(Drainage ~ Site_order, scales = "free") +
  labs(y = "Hybrid Ancestry", x = "Distancee from First Site (km)") +
  scale_fill_viridis_c(option = viridis_scale, direction = -1, limits = c(0,1)) +
  geom_tile()
conzintla_spectrum



### Calnali Tributaries

tribs <- read.csv("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/calnali_tributaries_hybrid_index_lowcountsremoved.csv") %>%
  mutate(Site_order = as.factor(Site_order)) %>%
  filter(malcount + bircount > 10000)
levels(tribs$Site_order) <- c("CATU","CATM","CATD","CDDR")
tribs <- left_join(tribs, tribs %>% group_by(Site.Code) %>% dplyr::summarize(Site_n = dplyr::n())) %>%
  group_by(Site.Code) %>%
  arrange(hybrid_index) %>%
  mutate(anc_index = row_number())
tribs$ancgroup <- cut(tribs$hybrid_index, breaks = seq(0,1, 0.04))
hybridramp <- colorRampPalette(c(bircol,hetcol))(25)

trib_hist <- ggplot(tribs, aes(y = hybrid_index, fill = ancgroup)) +
  theme_cowplot()+
  theme(legend.position = "none") +
  labs(x="Count", y = expression(Proportion~italic(X.~malinche))) +
  coord_cartesian(y = c(0, 1)) +
  facet_grid(. ~ Site_order) +
  scale_fill_discrete(type = hybridramp, drop = F) +
  geom_histogram()
trib_hist
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/tributary_ancestries.png",
       trib_hist, width = 6.5, height = 4)
