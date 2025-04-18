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

site_distances <- read.csv("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/site_riverkm_distances_noduplicates.csv")

clines <- read.csv("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/huazalingo_pochula_calnali_conzintla_allhybridindex.csv") %>%
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

yearly_ancstats <- clines %>% 
  group_by(Site.Code, Drainage, Collection_Year) %>% 
  summarize(mean_anc = mean(hybrid_index),
            sd_anc = sd(hybrid_index),
            var_anc = var(hybrid_index),
            skew_anc = skewness(hybrid_index),
            kurt_anc = kurtosis(hybrid_index),
            mean_het = mean(heterzygosity),
            sd_het = sd(heterzygosity),
            n = n())


multiyear_sites <- unique((yearly_ancstats %>% group_by(Site.Code, Drainage) %>% filter(n() > 1))$Site.Code)

multiyear_site_clines <- filter(clines, Site.Code %in% multiyear_sites)
multiyear_ancstats <- filter(yearly_ancstats, Site.Code %in% multiyear_sites)

yearly_ancestry_diptests <- lapply(unique(multiyear_site_clines$Site.Code), function(pop) {
  lapply(unique(filter(multiyear_site_clines, Site.Code == pop)$Collection_Year), function(year) {
    test <- dip.test(filter(multiyear_site_clines, Site.Code == pop, Collection_Year == year)$hybrid_index)
    return(c("Site.Code"=pop,"Year"=year,"statistic"=test$statistic,"p.value"=test$p.value,"nobs"=test$nobs))
  })
}) %>%
  bind_rows() %>%
  mutate(statistic.D = as.numeric(statistic.D),
         p.value = as.numeric(p.value),
         corrected.signif = p.value < 0.05 / length(statistic.D),
         nobs = as.numeric(nobs))

all_KS_comparisons <- lapply(multiyear_sites, function(site) {
  site_subset <- filter(multiyear_site_clines, Site.Code == site)
  allcomps <- lapply(unique(site_subset$Collection_Year), function(Year1) {
    comps <- lapply(unique(site_subset$Collection_Year[site_subset$Collection_Year > Year1]), function(Year2) {
      mean_anc_diff <- mean(filter(site_subset, Collection_Year == Year2)$hybrid_index) - mean(filter(site_subset, Collection_Year == Year1)$hybrid_index)
      var_anc_diff <- var(filter(site_subset, Collection_Year == Year2)$hybrid_index) - var(filter(site_subset, Collection_Year == Year1)$hybrid_index)
      pval <- ks.test(filter(site_subset, Collection_Year == Year1)$hybrid_index, filter(site_subset, Collection_Year == Year2)$hybrid_index)$p
      return(c(Year2, pval, mean_anc_diff, var_anc_diff))
    })
    return(lapply(comps, function(x) c(site, Year1, x)))
  }) 
})
test <- lapply(all_KS_comparisons, function(x) {
  lapply(x, function(y) {
    if(!is.null(y)) do.call(rbind, y)
  })
})

collapsed = {}
for(element in 1:length(test)) {
  for(subelement in 1:length(test[[element]])) {
    if(!is.null(test[[element]][[subelement]])) {
      collapsed <- rbind(collapsed, test[[element]][[subelement]])
    }
  }
}
comparison_df <- as.data.frame(collapsed) %>%
  rename(`Site.Code` = V1, `Year1` = V2, `Year2` = V3, `Pval` = V4, `mean_anc_diff` = V5, `var_anc_diff` = V6) %>%
  mutate(Year1 = as.numeric(Year1),
         Year2 = as.numeric(Year2),
         Pval = as.numeric(Pval),
         Bonferroni.signif = Pval < 0.05 / length(Pval),
         BH.signif = Pval < rank(Pval) * 0.05 / length(Pval), # Not valid, tests are not independent
         mean_anc_diff = as.numeric(mean_anc_diff),
         abs_anc_diff = abs(as.numeric(mean_anc_diff)),
         var_anc_diff = as.numeric(var_anc_diff),
         yeardiff = Year2-Year1)

cor.test(comparison_df$yeardiff, comparison_df$Pval, method = "spearman")
cor.test(comparison_df$yeardiff, comparison_df$var_anc_diff, method = "spearman")
cor.test(comparison_df$yeardiff, comparison_df$abs_anc_diff, method = "spearman")
cor.test(comparison_df$yeardiff, comparison_df$mean_anc_diff, method = "spearman")

plot(comparison_df$yeardiff, comparison_df$Pval)
plot(comparison_df$yeardiff, comparison_df$mean_anc_diff)
plot(comparison_df$yeardiff, comparison_df$var_anc_diff)
plot(comparison_df$yeardiff, comparison_df$abs_anc_diff)


no_crzn <- comparison_df %>% filter(!Site.Code %in% c("CRZN"))
cor.test(no_crzn$yeardiff, no_crzn$Pval, method = "spearman")
cor.test(no_crzn$yeardiff, no_crzn$mean_anc_diff, method = "spearman")
cor.test(no_crzn$yeardiff, no_crzn$var_anc_diff, method = "spearman")
cor.test(no_crzn$yeardiff, no_crzn$abs_anc_diff, method = "spearman")

plot(no_crzn$yeardiff, no_crzn$Pval)
plot(no_crzn$yeardiff, no_crzn$mean_anc_diff)
plot(no_crzn$yeardiff, no_crzn$var_anc_diff)
plot(no_crzn$yeardiff, no_crzn$abs_anc_diff)


stable <- comparison_df %>%
  group_by(Site.Code) %>%
  summarize(ancdiff = any(Bonferroni.signif))

all_yearly <- map(multiyear_sites, function(site) {
  ggplot(dplyr::filter(multiyear_site_clines, Site.Code == site), aes(x = hybrid_index)) +
    coord_cartesian(xlim = c(0,1)) +
    theme_bw() +
    facet_grid(Collection_Year~.) +
    geom_histogram() +
    ggtitle(site)
})
all_yearly

crzn_cline <- ggplot(filter(clines, malcount + bircount > 10000, Site.Code == "CRZN"),
                     aes(y = hybrid_index, fill = ancgroup)) +
  coord_cartesian(y = c(0,1)) +
  theme_cowplot()+
  labs(y = expression("Proportion"~italic("X. malinche")~"Ancestry")) +
  theme(legend.position = "none") +
  facet_grid(. ~ Collection_Year, scales = "free") +
  scale_fill_discrete(type = hybridramp, drop = F) +
  geom_hline(data = filter(multiyear_ancstats, Site.Code == "CRZN"), aes(yintercept = mean_anc), lty = 2) +
  geom_histogram() +
  ggtitle("Corazon")
crzn_cline

toto_cline <- ggplot(filter(clines, malcount + bircount > 10000, Site.Code == "TOTO"),
                     aes(y = hybrid_index, fill = ancgroup)) +
  coord_cartesian(y = c(0,1)) +
  theme_cowplot()+
  labs(y = expression("Proportion"~italic("X. malinche")~"Ancestry")) +
  theme(legend.position = "none") +
  facet_grid(. ~ Collection_Year, scales = "free") +
  scale_fill_discrete(type = hybridramp, drop = F) +
  geom_hline(data = filter(multiyear_ancstats, Site.Code == "TOTO"), aes(yintercept = mean_anc), lty = 2) +
  geom_histogram() +
  ggtitle("Totonicapa")
toto_cline

xoch_cline <- ggplot(filter(clines, malcount + bircount > 10000, Site.Code == "XOCH"),
                     aes(y = hybrid_index, fill = ancgroup)) +
  coord_cartesian(y = c(0,1)) +
  theme_cowplot()+
  labs(y = expression("Proportion"~italic("X. malinche")~"Ancestry")) +
  theme(legend.position = "none") +
  facet_grid(. ~ Collection_Year, scales = "free") +
  scale_fill_discrete(type = hybridramp, drop = F) +
  geom_hline(data = filter(multiyear_ancstats, Site.Code == "XOCH"), aes(yintercept = mean_anc), lty = 2) +
  geom_histogram() +
  ggtitle("Xochicoatlán")
xoch_cline

call_cline <- ggplot(filter(clines, malcount + bircount > 10000, Site.Code == "CALL"),
                     aes(y = hybrid_index, fill = ancgroup)) +
  coord_cartesian(y = c(0,1)) +
  theme_cowplot()+
  labs(y = expression("Proportion"~italic("X. malinche")~"Ancestry")) +
  theme(legend.position = "none") +
  facet_grid(. ~ Collection_Year, scales = "free") +
  scale_fill_discrete(type = hybridramp, drop = F) +
  geom_hline(data = filter(multiyear_ancstats, Site.Code == "CALL"), aes(yintercept = mean_anc), lty = 2) +
  geom_histogram() +
  ggtitle("Calnali Low")
call_cline

#clines_temporal_variation_top <- plot_grid(crzn_cline, toto_cline, nrow = 1, labels = c("A", "B"))
#clines_temporal_variation_bottom <- plot_grid(xoch_cline, call_cline, nrow = 1, labels = c("C", "D"), rel_widths = c(1,2))
#clines_temporal_variation <- plot_grid(clines_temporal_variation_top, clines_temporal_variation_bottom, nrow = 2)
#clines_temporal_variation <- plot_grid(crzn_cline, xoch_cline, call_cline, nrow = 4, labels = c("A", "B", "C", "D"))
clines_temporal_variation_top <- plot_grid(crzn_cline, xoch_cline, nrow = 1, labels = c("A", "B"), rel_widths = c(4,3))
clines_temporal_variation_bottom <- plot_grid(call_cline, nrow = 1, labels = c("C"))
clines_temporal_variation <- plot_grid(clines_temporal_variation_top, clines_temporal_variation_bottom, nrow = 2)
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/clines_temporal_variation.pdf",
       clines_temporal_variation, width = 8, height = 8)
