library(bayestestR)
library(tidyverse)
library(cowplot)

# mode of beta distribution = (a - 1) / (a + b -2)
# want the mode to be the ancestry of the individual

beta = function(alpha, anc) (alpha - 1)/anc - alpha + 2
alpha = 2.55
ancestry = 0.1
plot(seq(0,1,by = 0.01), dbeta(seq(0,1,by = 0.01), alpha, beta(alpha, ancestry)))
dbeta(0.2, alpha, beta(alpha, ancestry)) / dbeta(ancestry, alpha, beta(alpha, ancestry))

maternal_offspring_index_CALL<-read.csv(file="~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/hybrid_index_CALL_mothers_embryos_Sept2023_match.txt",sep="\t",head=TRUE)
pop_index_CALL<-read.csv(file="~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/hybrid_index_CALL_adults_August2023",sep="\t",head=TRUE) %>%
  filter(malcount + bircount > 10000)
maternal_offspring_index_CALM<-read.csv(file="~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/hybrid_index_readcounts_CALM_07_08_09_mother_embryo_Feb2025.csv",head=TRUE) %>%
  filter(read1_counts > 200000 & read1_counts.1 & 200000, heterzygosity < 0.75 & heterzygosity.1 < 0.75,
         !(maternal_id %in% c("CALM-08-F-129.R1.fastq", "CALM-08-F-106.R1.fastq", "CALM-08-F-106.R1.fastq"))) %>% # Three individuals have heterozygosities that would be semi-plausible if it weren't for the fact that their embryos are near-pure malinche - must be contamination
  mutate(diff = hybrid_index.1 - hybrid_index)
pop_index_CALM<-read.csv(file="~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/huazalingo_pochula_calnali_conzintla_allhybridindex.csv",head=TRUE) %>%
  filter(malcount + bircount > 10000,
         Site.Code == "CALM")


ABC_nulls_assortative_mating_betadist <- function(mother_embryo_indices, population_indices, seed = 33, reps = 100000, sibling_sd = 0.002, maxalpha = 20) {
  females<-unique(mother_embryo_indices$maternal_id)
  mates<-c(na.omit(population_indices$hybrid_index))
  
  all_sims <- {}
  empirical <- {}
  set.seed(seed)
  for (j in 1:reps){
    
    #alpha = exp(runif(1, 0, log(maxalpha)))
    alpha = runif(1, 1, maxalpha)
    
    #generate random subset for multi embryo females
    
    idxsub<-{}
    for(k in 1:length(females)){      # Draw female-offspring pairs
      
      focal<-subset(mother_embryo_indices,mother_embryo_indices$maternal_id==females[k])
      rand<-focal[sample(nrow(focal),1),]
      idxsub<-rbind(idxsub,cbind(rand))
      
    }
    empirical<-rbind(empirical,cbind(mean(abs(idxsub$diff),na.rm=TRUE),var(abs(idxsub$diff),na.rm=TRUE)))
    
    offspring_index<-{}
    maternal_index<-{}
    for (x in 1:length(females)){
      maternal<-sample(idxsub$hybrid_index,1)
      b = (alpha - 1)/maternal - alpha + 2
      maxdbeta = dbeta(maternal, alpha, b)
      mate = sample(mates, 1, prob = dbeta(mates, alpha, b))
      hold=0
      off<-rnorm(1,mean=(maternal+mate)/2,sd=sibling_sd)
      #average sd for the plot
      offspring_index<-c(offspring_index,off)
      maternal_index<-c(maternal_index,maternal)
    }
    
    all_sims<-rbind(all_sims,cbind(mean(abs(maternal_index-offspring_index),na.rm=TRUE),var(abs(maternal_index-offspring_index),na.rm=TRUE),alpha))
    
  }
  return(list(all_sims, empirical))
}

CALM_r=0.01947892
CALM_v=0.0002921749


nulls_CALM_maxalpha20_reps3M_sd0.01 <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/nulls_betadist_maxalpha20_reps3M_sd0.01_combined_CALM.csv")
accepted_CALM_maxalpha20_sd0.01 <-filter(nulls_CALM_maxalpha20_reps3M_sd0.01,sim_diff > (CALM_r-CALM_r*0.05), sim_diff < (CALM_r+CALM_r*0.05), sim_var >(CALM_v-CALM_v*0.05), sim_var < (CALM_v+CALM_v*0.05))
accepted_CALM_20_MAP_sd0.01 <- map_estimate(accepted_CALM_maxalpha20_sd0.01[,'mating_prop'])$MAP_Estimate         # 6.65
ci(accepted_CALM_maxalpha20_sd0.01[,'mating_prop'], method = "ETI")                                   # 2.99 to 19.27
hist(accepted_CALM_maxalpha20_sd0.01[,'mating_prop'], xlim = c(1,20))

nulls_CALM_maxalpha20_reps3M_sd0.002 <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/nulls_betadist_maxalpha20_reps3M_sd0.002_combined_CALM.csv")
accepted_CALM_maxalpha20_sd0.002 <-filter(nulls_CALM_maxalpha20_reps3M_sd0.002,sim_diff > (CALM_r-CALM_r*0.05), sim_diff < (CALM_r+CALM_r*0.05), sim_var >(CALM_v-CALM_v*0.05), sim_var < (CALM_v+CALM_v*0.05))
accepted_CALM_20_MAP_sd0.002 <- map_estimate(accepted_CALM_maxalpha20_sd0.002[,'mating_prop'])$MAP_Estimate         # 5.53
ci(accepted_CALM_maxalpha20_sd0.002[,'mating_prop'], method = "ETI")                                   # 2.65 to 18.36
hist(accepted_CALM_maxalpha20_sd0.002[,'mating_prop'], xlim = c(1,20))


CALL_r=0.066
CALL_v=0.0046

nulls_CALL_maxalpha20_reps3M_sd0.01 <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/nulls_betadist_maxalpha20_reps3M_sd0.01_combined_CALL.csv")
accepted_CALL_maxalpha20_sd0.01 <-filter(nulls_CALL_maxalpha20_reps3M_sd0.01,sim_diff > (CALL_r-CALL_r*0.05), sim_diff < (CALL_r+CALL_r*0.05), sim_var >(CALL_v-CALL_v*0.05), sim_var < (CALL_v+CALL_v*0.05))
accepted_CALL_20_MAP_sd0.01 <- map_estimate(accepted_CALL_maxalpha20_sd0.01[,'mating_prop'])$MAP_Estimate         # 2.82
ci(accepted_CALL_maxalpha20_sd0.01[,'mating_prop'], method = "ETI")                                   # 1.76 to 4.97
hist(accepted_CALL_maxalpha20_sd0.01[,'mating_prop'], xlim = c(1,20))

nulls_CALL_maxalpha20_reps3M_sd0.002 <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/nulls_betadist_maxalpha20_reps3M_sd0.002_combined_CALL.csv")
accepted_CALL_maxalpha20_sd0.002 <-filter(nulls_CALL_maxalpha20_reps3M_sd0.002,sim_diff > (CALL_r-CALL_r*0.05), sim_diff < (CALL_r+CALL_r*0.05), sim_var >(CALL_v-CALL_v*0.05), sim_var < (CALL_v+CALL_v*0.05))
accepted_CALL_20_MAP_sd0.002 <- map_estimate(accepted_CALL_maxalpha20_sd0.002[,'mating_prop'])$MAP_Estimate         # 2.94
ci(accepted_CALL_maxalpha20_sd0.002[,'mating_prop'], method = "ETI")                                   # 1.74 to 4.84
hist(accepted_CALL_maxalpha20_sd0.002[,'mating_prop'], xlim = c(1,20))

mating_alpha_diff <- accepted_CALM_maxalpha20_sd0.01[1:min(nrow(accepted_CALM_maxalpha20_sd0.01),nrow(accepted_CALL_maxalpha20_sd0.01)),] %>%
  mutate(mating_diff = mating_prop - accepted_CALL_maxalpha20_sd0.01[1:min(nrow(accepted_CALM_maxalpha20_sd0.01),nrow(accepted_CALL_maxalpha20_sd0.01)),"mating_prop"])
accepted_alpha_diff_sd0.01 <- map_estimate(mating_alpha_diff$mating_diff)$MAP_Estimate         # 3.68
ci(mating_alpha_diff$mating_diff, method = "ETI")                                   # -0.49 to 16.31

mating_alpha_diff_sd0.002 <- accepted_CALM_maxalpha20_sd0.002[1:min(nrow(accepted_CALM_maxalpha20_sd0.002),nrow(accepted_CALL_maxalpha20_sd0.002)),] %>%
  mutate(mating_diff = mating_prop - accepted_CALL_maxalpha20_sd0.002[1:min(nrow(accepted_CALM_maxalpha20_sd0.002),nrow(accepted_CALL_maxalpha20_sd0.002)),"mating_prop"])
accepted_alpha_diff_sd0.002 <- map_estimate(mating_alpha_diff_sd0.002$mating_diff)$MAP_Estimate         # 2.94
ci(mating_alpha_diff_sd0.002$mating_diff, method = "ETI")                                   # -0.71 to 15.38

#nulls_maxalpha200_reps200k <- ABC_nulls_assortative_mating(seed = 33, reps = 200000, maxalpha = 200)
#microbenchmark::microbenchmark(nulls_maxalpha100_reps100k <- ABC_nulls_assortative_mating(reps = 1, maxalpha = 100))
nulls_maxalpha20_reps100k <- ABC_nulls_assortative_mating(seed = 33, reps = 100000, maxalpha = 20)
accepted_dbeta_100k <- base::subset(nulls_maxalpha20_reps100k[[1]],nulls_maxalpha20_reps100k[[1]][,1]>(r-r*0.05) & nulls_maxalpha20_reps100k[[1]][,1]< (r+r*0.05) & nulls_maxalpha20_reps100k[[1]][,2]>(v-v*0.05) & nulls_maxalpha20_reps100k[[1]][,2]< (v+v*0.05))
hist(accepted_dbeta_100k[,'alpha'])
accepted_dbeta_MAP <- map_estimate(accepted_dbeta_100k[,'alpha'])$MAP_Estimate         # 2.83
accepted_dbeta_MAP_CI <- ci(accepted_dbeta_100k[,'alpha'], method = "ETI")             # 1.71 to 4.82

ancestry = median(pop_index_CALL$hybrid_index)
plot(seq(0,1,by = 0.01), dbeta(seq(0,1,by = 0.01), ci(accepted_dbeta_100k[,'alpha'], method = "ETI")$CI_high, beta(ci(accepted_dbeta_100k[,'alpha'], method = "ETI")$CI_high, ancestry)), type = "l", lty = 2)
points(seq(0,1,by = 0.01), dbeta(seq(0,1,by = 0.01), accepted_dbeta_MAP, beta(accepted_dbeta_MAP, ancestry)), type = "l")
points(seq(0,1,by = 0.01), dbeta(seq(0,1,by = 0.01), ci(accepted_dbeta_100k[,'alpha'], method = "ETI")$CI_low, beta(ci(accepted_dbeta_100k[,'alpha'], method = "ETI")$CI_low, ancestry)), type = "l", lty = 3)
points(x = na.omit(pop_index_CALL$hybrid_index), y = runif(length(na.omit(pop_index_CALL$hybrid_index)), 0, 2))


posterior_alpha_dbeta <- ggplot(as.data.frame(accepted_dbeta_100k), aes(x = alpha)) +
  theme_bw() +
  labs(y = "Frequency", x = expression("Simulated"~alpha)) +
  geom_histogram(fill = "cornflowerblue", binwidth = 0.25) +
  geom_vline(xintercept = accepted_dbeta_MAP, lty = 2, col = "red4") +
  geom_text(label = paste("MAP =",round(accepted_dbeta_MAP,2)), col = "red4", x = accepted_dbeta_MAP - .75, y = Inf, vjust = 2)
posterior_alpha_dbeta
ggsave("~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/posterior_distribution_assortative_mating_ABC_CALL_betadist_alpha.pdf",
       posterior_alpha_dbeta, width = 3.5, height = 2, units = "in")

mating_pref_dens <- data.frame(x = seq(0, 1, by = .001), 
                              MAP = dbeta(seq(0, 1, by = .001), accepted_dbeta_MAP, beta(accepted_dbeta_MAP, median(pop_index_CALL$hybrid_index))),
                              lowCI = dbeta(seq(0, 1, by = .001), accepted_dbeta_MAP_CI$CI_low, beta(accepted_dbeta_MAP_CI$CI_low, median(pop_index_CALL$hybrid_index))),
                              highCI = dbeta(seq(0, 1, by = .001), accepted_dbeta_MAP_CI$CI_high, beta(accepted_dbeta_MAP_CI$CI_high, median(pop_index_CALL$hybrid_index))))

mating_ci_MAP = list(CI_low = qbeta(p = 0.025, accepted_dbeta_MAP, beta(accepted_dbeta_MAP, median(pop_index_CALL$hybrid_index))), CI_high = qbeta(p = 0.975, accepted_dbeta_MAP, beta(accepted_dbeta_MAP, median(pop_index_CALL$hybrid_index))))
mating_ci_lowCI = list(CI_low = qbeta(p = 0.025, accepted_dbeta_MAP_CI$CI_low, beta(accepted_dbeta_MAP_CI$CI_low, median(pop_index_CALL$hybrid_index))), CI_high = qbeta(p = 0.975, accepted_dbeta_MAP_CI$CI_low, beta(accepted_dbeta_MAP_CI$CI_low, median(pop_index_CALL$hybrid_index))))
mating_ci_highCI = list(CI_low = qbeta(p = 0.025, accepted_dbeta_MAP_CI$CI_high, beta(accepted_dbeta_MAP_CI$CI_high, median(pop_index_CALL$hybrid_index))), CI_high = qbeta(p = 0.975, accepted_dbeta_MAP_CI$CI_high, beta(accepted_dbeta_MAP_CI$CI_high, median(pop_index_CALL$hybrid_index))))

hypothetical_mating_dist <- ggplot(mating_pref_dens, aes(x = x * 100)) +
  labs(x = expression(atop("Mate Ancestry Fraction","(%"~italic("X. malinche")*")")), y = "Relative Preference") +
  theme_cowplot() +
  theme(axis.text.y = element_blank()) +
  #geom_density(data = . %>% filter(between(x,mating_ci_MAP$CI_low, mating_ci_MAP$CI_high)), stat = "identity", fill = "cornflowerblue", alpha = .5) +
  geom_density(aes(y = MAP), stat = "identity", fill = "cornflowerblue") +
  #geom_density(aes(y = lowCI), data = . %>% filter(between(x,mating_ci_lowCI$CI_low, mating_ci_lowCI$CI_high)), stat = "identity", fill = "yellow", alpha = .5, color = NA) +
  #geom_density(aes(y = lowCI), stat = "identity", lty = "dotted") +
  #geom_density(aes(y = highCI), data = . %>% filter(between(x,mating_ci_highCI$CI_low, mating_ci_highCI$CI_high)), stat = "identity", fill = "grey", alpha = .5, color = NA) +
  #geom_density(aes(y = highCI), stat = "identity", lty = "dotted") #+
  #geom_density(data = filter(clines, Site.Code == "CALL"), aes(x = hybrid_index))
  geom_vline(xintercept = median(pop_index_CALL$hybrid_index * 100), lty = "dotted")
hypothetical_mating_dist
ggsave("~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/MAP_assortative_mating_ABC_CALL_betadist.pdf",
       hypothetical_mating_dist, width = 3.25, height = 3.25, units = "in")


hypothetical_mating_dist_CI_medfemale <- ggplot(mating_pref_dens, aes(x = x * 100)) +
  labs(x = expression(atop("Mate Ancestry Fraction","(%"~italic("X. malinche")*")")), y = expression(atop("Relative", "Preference"))) +
  theme_cowplot() +
  theme(axis.text.y = element_blank()) +
  #geom_density(data = . %>% filter(between(x,mating_ci_MAP$CI_low, mating_ci_MAP$CI_high)), stat = "identity", fill = "cornflowerblue", alpha = .5) +
  geom_density(aes(y = MAP), stat = "identity") +
  #geom_density(aes(y = lowCI), data = . %>% filter(between(x,mating_ci_lowCI$CI_low, mating_ci_lowCI$CI_high)), stat = "identity", fill = "yellow", alpha = .5, color = NA) +
  geom_density(aes(y = lowCI), stat = "identity", lty = "dotted") +
  #geom_density(aes(y = highCI), data = . %>% filter(between(x,mating_ci_highCI$CI_low, mating_ci_highCI$CI_high)), stat = "identity", fill = "grey", alpha = .5, color = NA) +
  geom_density(aes(y = highCI), stat = "identity", lty = "dashed") #+
  #geom_density(data = filter(clines, Site.Code == "CALL"), aes(x = hybrid_index))
hypothetical_mating_dist_CI_medfemale

mating_pref_dens_halfHI <- data.frame(x = seq(0, 1, by = .001), 
                                     MAP = dbeta(seq(0, 1, by = .001), accepted_dbeta_MAP, beta(accepted_dbeta_MAP, 0.5)),
                                     lowCI = dbeta(seq(0, 1, by = .001), accepted_dbeta_MAP_CI$CI_low, beta(accepted_dbeta_MAP_CI$CI_low, 0.5)),
                                     highCI = dbeta(seq(0, 1, by = .001), accepted_dbeta_MAP_CI$CI_high, beta(accepted_dbeta_MAP_CI$CI_high, 0.5)))

hypothetical_mating_dist_CI_halffemale <- ggplot(mating_pref_dens_halfHI, aes(x = x * 100)) +
  labs(x = expression(atop("Mate Ancestry Fraction","(%"~italic("X. malinche")*")")), y = expression(atop("Relative", "Preference"))) +
  theme_cowplot() +
  theme(axis.text.y = element_blank()) +
  #geom_density(data = . %>% filter(between(x,mating_ci_MAP$CI_low, mating_ci_MAP$CI_high)), stat = "identity", fill = "cornflowerblue", alpha = .5) +
  geom_density(aes(y = MAP), stat = "identity") +
  #geom_density(aes(y = lowCI), data = . %>% filter(between(x,mating_ci_lowCI$CI_low, mating_ci_lowCI$CI_high)), stat = "identity", fill = "yellow", alpha = .5, color = NA) +
  geom_density(aes(y = lowCI), stat = "identity", lty = "dotted") +
  #geom_density(aes(y = highCI), data = . %>% filter(between(x,mating_ci_highCI$CI_low, mating_ci_highCI$CI_high)), stat = "identity", fill = "grey", alpha = .5, color = NA) +
  geom_density(aes(y = highCI), stat = "identity", lty = "dashed") #+
#geom_density(data = filter(clines, Site.Code == "CALL"), aes(x = hybrid_index))
hypothetical_mating_dist_CI_halffemale

mating_pref_dens_minHI <- data.frame(x = seq(0, 1, by = .001), 
                               MAP = dbeta(seq(0, 1, by = .001), accepted_dbeta_MAP, beta(accepted_dbeta_MAP, min(pop_index_CALL$hybrid_index))),
                               lowCI = dbeta(seq(0, 1, by = .001), accepted_dbeta_MAP_CI$CI_low, beta(accepted_dbeta_MAP_CI$CI_low, min(pop_index_CALL$hybrid_index))),
                               highCI = dbeta(seq(0, 1, by = .001), accepted_dbeta_MAP_CI$CI_high, beta(accepted_dbeta_MAP_CI$CI_high, min(pop_index_CALL$hybrid_index))))

hypothetical_mating_dist_CI_minfemale <- ggplot(mating_pref_dens_minHI, aes(x = x * 100)) +
  labs(x = expression(atop("Mate Ancestry Fraction","(%"~italic("X. malinche")*")")), y = expression(atop("Relative", "Preference"))) +
  theme_cowplot() +
  theme(axis.text.y = element_blank()) +
  #geom_density(data = . %>% filter(between(x,mating_ci_MAP$CI_low, mating_ci_MAP$CI_high)), stat = "identity", fill = "cornflowerblue", alpha = .5) +
  geom_density(aes(y = MAP), stat = "identity") +
  #geom_density(aes(y = lowCI), data = . %>% filter(between(x,mating_ci_lowCI$CI_low, mating_ci_lowCI$CI_high)), stat = "identity", fill = "yellow", alpha = .5, color = NA) +
  geom_density(aes(y = lowCI), stat = "identity", lty = "dotted") +
  #geom_density(aes(y = highCI), data = . %>% filter(between(x,mating_ci_highCI$CI_low, mating_ci_highCI$CI_high)), stat = "identity", fill = "grey", alpha = .5, color = NA) +
  geom_density(aes(y = highCI), stat = "identity", lty = "dashed") #+
#geom_density(data = filter(clines, Site.Code == "CALL"), aes(x = hybrid_index))
hypothetical_mating_dist_CI_minfemale

mating_pref_dens_maxHI <- data.frame(x = seq(0, 1, by = .001), 
                                     MAP = dbeta(seq(0, 1, by = .001), accepted_dbeta_MAP, beta(accepted_dbeta_MAP, max(pop_index_CALL$hybrid_index))),
                                     lowCI = dbeta(seq(0, 1, by = .001), accepted_dbeta_MAP_CI$CI_low, beta(accepted_dbeta_MAP_CI$CI_low, max(pop_index_CALL$hybrid_index))),
                                     highCI = dbeta(seq(0, 1, by = .001), accepted_dbeta_MAP_CI$CI_high, beta(accepted_dbeta_MAP_CI$CI_high, max(pop_index_CALL$hybrid_index))))

hypothetical_mating_dist_CI_maxfemale <- ggplot(mating_pref_dens_maxHI, aes(x = x * 100)) +
  labs(x = expression(atop("Mate Ancestry Fraction","(%"~italic("X. malinche")*")")), y = expression(atop("Relative", "Preference"))) +
  theme_cowplot() +
  theme(axis.text.y = element_blank()) +
  #geom_density(data = . %>% filter(between(x,mating_ci_MAP$CI_low, mating_ci_MAP$CI_high)), stat = "identity", fill = "cornflowerblue", alpha = .5) +
  geom_density(aes(y = MAP), stat = "identity") +
  #geom_density(aes(y = lowCI), data = . %>% filter(between(x,mating_ci_lowCI$CI_low, mating_ci_lowCI$CI_high)), stat = "identity", fill = "yellow", alpha = .5, color = NA) +
  geom_density(aes(y = lowCI), stat = "identity", lty = "dotted") +
  #geom_density(aes(y = highCI), data = . %>% filter(between(x,mating_ci_highCI$CI_low, mating_ci_highCI$CI_high)), stat = "identity", fill = "grey", alpha = .5, color = NA) +
  geom_density(aes(y = highCI), stat = "identity", lty = "dashed") #+
#geom_density(data = filter(clines, Site.Code == "CALL"), aes(x = hybrid_index))
hypothetical_mating_dist_CI_maxfemale

dummy_data <- data.frame(Estimate = rep(c("MAP","CI Lower Bound", "CI Upper Bound"), 2), x = c(0.5, 7.5, 3, 5,6,4), y = c(3, 3, 3,8,2,3))
dummy_data$Estimate <- factor(dummy_data$Estimate, levels = c("MAP","CI Lower Bound", "CI Upper Bound"))
dummy_plot <- ggplot(data = dummy_data,aes(x = x, y = y, linetype = Estimate)) +
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.justification.bottom = 0.5) +
  scale_linetype_manual(values = c("solid","dotted","dashed")) +
  geom_line()
dummy_plot
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

leg <- get_legend_35(dummy_plot)

female_pref_CI_plots <- plot_grid(hypothetical_mating_dist_CI_medfemale, hypothetical_mating_dist_CI_halffemale, hypothetical_mating_dist_CI_minfemale, hypothetical_mating_dist_CI_maxfemale,
        nrow = 2, labels = "AUTO")
female_pref_CI_plots_wleg <- plot_grid(female_pref_CI_plots, leg, nrow = 2, rel_heights = c(15,1))
female_pref_CI_plots_wleg
ggsave("~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/MAP_assortative_mating_ABC_CALL_betadist_examples.pdf",
       female_pref_CI_plots_wleg, width = 6.5, height = 6.5, units = "in")

