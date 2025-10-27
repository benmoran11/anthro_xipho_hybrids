library(bayestestR)
library(tidyverse)
library(cowplot)

# mode of beta distribution = (a - 1) / (a + b -2)
# want the mode to be the ancestry of the individual
beta = function(alpha, anc) (alpha - 1)/anc - alpha + 2
alpha = 3.68
ancestry = 0.97
plot(seq(0,1,by = 0.01), dbeta(seq(0,1,by = 0.01), alpha, beta(alpha, ancestry)))
dbeta(0.87, alpha, beta(alpha, ancestry)) / dbeta(ancestry, alpha, beta(alpha, ancestry))

maternal_offspring_index_CALL<-read.csv(file="hybrid_index_CALL_mother_embryo.csv",head=TRUE) # Pre-cleaned
maternal_offspring_index_CALM <-read.csv(file="hybrid_index_CALM_mother_embryo.csv",head=TRUE)  %>%
  filter(!(maternal_id %in% c("CALM-08-F-129.R1.fastq","CALM-08-F-106.R1.fastq"))) #  mother has plausible moderate heterozygosity but the embryos are pure malinche - must be contamination
# ^ Add '"CALM-08-F-129.R1.fastq",' to the above filter statement to remove the single intermediate individual (HI {0.3,0.7} with outsized effect,
# as discussed in Supp. Info. 7 
pop_index_CALM<-read.csv(file="hybrid_index_CALM_adults.csv",head=TRUE) %>%
  filter(malcount + bircount > 10000)
pop_index_CALL<-read.csv(file="hybrid_index_CALL_adults.csv",head=TRUE) %>%
  filter(malcount + bircount > 10000)

#generate random subset for multi embryo females

idxsub_CALM<-{}
females_CALM<-unique(maternal_offspring_index_CALM$maternal_id)
set.seed(14)
for(k in 1:length(females_CALM)){
  
  focal<-subset(maternal_offspring_index_CALM,maternal_offspring_index_CALM$maternal_id==females_CALM[k])
  rand<-focal[sample(nrow(focal),1),]
  idxsub_CALM<-rbind(idxsub_CALM,cbind(rand))
  
}
mean(idxsub_CALM$hybrid_index.1)
mean(abs(idxsub_CALM$diff)) # 0.02568919, # 0.01639041 without F129
var(abs(idxsub_CALM$diff)) #  0.002615013 # 0.0001779394 without F129

## Now CALL 
idxsub_CALL<-{}
females_CALL<-unique(maternal_offspring_index_CALL$maternal_id)
set.seed(14)
for(k in 1:length(females_CALL)){
  
  focal<-subset(maternal_offspring_index_CALL,maternal_offspring_index_CALL$maternal_id==females_CALL[k])
  rand<-focal[sample(nrow(focal),1),]
  idxsub_CALL<-rbind(idxsub_CALL,cbind(rand))
  
}
mean(idxsub_CALL$hybrid_index.1)
mean(abs(idxsub_CALL$diff)) # 0.06554228
var(abs(idxsub_CALL$diff)) #  0.004638745

### Function for running assortative mating simulations with continuous mate preference based on a beta distribution


ABC_nulls_assortative_mating_betadist <- function(mother_embryo_indices, population_indices, seed = 33, reps = 100000, sibling_sd = 0.01, maxalpha = 20) {
  females<-unique(mother_embryo_indices$maternal_id)
  mates<-c(na.omit(population_indices$hybrid_index))
  
  all_sims <- {}
  empirical <- {}
  set.seed(seed)
  for (j in 1:reps){
    
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


### Simulations were run in batches of 100000 simulations
### In the input files below, "maxalpha" refers to the maxalpha parameter,
### "sd" refers to the sibling_sd parameter, 
### and the random seed was set to the index of the simulation batch (from 1 to 30)
### for each population, "mother_embryo_indices" was the corresponding maternal_offspring_index_#### variable, and
### "population_indices" was the corresponding pop_index_#### variable

CALM_r=0.02568919
CALM_v=0.002615013

CALM_noF129_r=0.01639041
CALM_noF129_v=0.0001779394

CALL_r=0.066
CALL_v=0.0046

### Performing rejection sampling on simulations performed separately (parallelized on cluster)

nulls_CALM_maxalpha20_reps3M_sd0.002 <- read.csv("nulls_cluster_betadist_maxalpha20_reps3M_sd0.002_CALM.csv")
accepted_CALM_maxalpha20_sd0.002 <-filter(nulls_CALM_maxalpha20_reps3M_sd0.002,sim_diff > (CALM_r-CALM_r*0.05), sim_diff < (CALM_r+CALM_r*0.05), sim_var >(CALM_v-CALM_v*0.05), sim_var < (CALM_v+CALM_v*0.05))
accepted_CALM_20_MAP_sd0.002 <- map_estimate(accepted_CALM_maxalpha20_sd0.002[,'alpha'])$MAP_Estimate         # 3.45
ci(accepted_CALM_maxalpha20_sd0.002[,'alpha'], method = "ETI")                                   # 1.92 to 10.02
hist(accepted_CALM_maxalpha20_sd0.002[,'alpha'], xlim = c(1,20))

nulls_CALL_maxalpha20_reps3M_sd0.002 <- read.csv("nulls_cluster_betadist_maxalpha20_reps3M_sd0.002_CALL.csv")
accepted_CALL_maxalpha20_sd0.002 <-filter(nulls_CALL_maxalpha20_reps3M_sd0.002,sim_diff > (CALL_r-CALL_r*0.05), sim_diff < (CALL_r+CALL_r*0.05), sim_var >(CALL_v-CALL_v*0.05), sim_var < (CALL_v+CALL_v*0.05))
accepted_CALL_20_MAP_sd0.002 <- map_estimate(accepted_CALL_maxalpha20_sd0.002[,'mating_prop'])$MAP_Estimate         # 2.93
ci(accepted_CALL_maxalpha20_sd0.002[,'mating_prop'], method = "ETI")                                   # 1.74 to 4.84
hist(accepted_CALL_maxalpha20_sd0.002[,'mating_prop'], xlim = c(1,20))

# Repeating CALM with summary stats calculated without female 129, which had an outsized effect on inference
nulls_CALM_noF129_maxalpha20_reps3M_sd0.002 <- read.csv("nulls_cluster_betadist_maxalpha20_reps3M_sd0.002_CALM_noF129.csv")
accepted_CALM_noF129_maxalpha20_sd0.002 <-filter(nulls_CALM_noF129_maxalpha20_reps3M_sd0.002,sim_diff > (CALM_noF129_r-CALM_noF129_r*0.05), sim_diff < (CALM_noF129_r+CALM_noF129_r*0.05), sim_var >(CALM_noF129_v-CALM_noF129_v*0.05), sim_var < (CALM_noF129_v+CALM_noF129_v*0.05))
accepted_CALM_noF129_20_MAP_sd0.002 <- map_estimate(accepted_CALM_noF129_maxalpha20_sd0.002[,'alpha'])$MAP_Estimate         # 10.98
ci(accepted_CALM_noF129_maxalpha20_sd0.002[,'alpha'], method = "ETI")                                   # 3.59 to 19.33
hist(accepted_CALM_noF129_maxalpha20_sd0.002[,'alpha'], xlim = c(1,20))


mating_alpha_diff_sd0.002 <- accepted_CALM_maxalpha20_sd0.002[1:min(nrow(accepted_CALM_maxalpha20_sd0.002),nrow(accepted_CALL_maxalpha20_sd0.002)),] %>%
  mutate(mating_diff = alpha - accepted_CALL_maxalpha20_sd0.002[1:min(nrow(accepted_CALM_maxalpha20_sd0.002),nrow(accepted_CALL_maxalpha20_sd0.002)),"mating_prop"])
accepted_alpha_diff_sd0.002 <- map_estimate(mating_alpha_diff_sd0.002$mating_diff)$MAP_Estimate         # 0.73
ci(mating_alpha_diff_sd0.002$mating_diff, method = "ETI")                                   # -1.73 to 7.16

# Plotting
posterior_alpha_dbeta <- ggplot(accepted_CALL_maxalpha20_sd0.002, aes(x = mating_prop)) +
  theme_bw() +
  labs(y = "Frequency", x = expression("Simulated"~alpha)) +
  geom_histogram(fill = "cornflowerblue", binwidth = 0.25) +
  geom_vline(xintercept = accepted_CALL_20_MAP_sd0.002, lty = 2, col = "red4") +
  geom_text(label = paste("MAP =",round(accepted_CALL_20_MAP_sd0.002,2)), col = "red4", x = accepted_CALL_20_MAP_sd0.002 - .75, y = Inf, vjust = 2)
posterior_alpha_dbeta
ggsave("posterior_distribution_assortative_mating_ABC_CALL_betadist_alpha.pdf",
       posterior_alpha_dbeta, width = 3.5, height = 2, units = "in")


### Generating inferred mating preference curves for hypothetical females of varying ancestry

mating_pref_dens_medfemaleHI <- data.frame(x = seq(0, 1, by = .001), 
                                      CALM = dbeta(seq(0, 1, by = .001), accepted_CALM_20_MAP_sd0.002, beta(accepted_CALM_20_MAP_sd0.002, median(pop_index_CALL$hybrid_index))),
                                      CALL = dbeta(seq(0, 1, by = .001), accepted_CALL_20_MAP_sd0.002, beta(accepted_CALL_20_MAP_sd0.002, median(pop_index_CALL$hybrid_index))),
                                      CALM_noF129 = dbeta(seq(0, 1, by = .001), accepted_CALM_noF129_20_MAP_sd0.002, beta(accepted_CALM_noF129_20_MAP_sd0.002, median(pop_index_CALL$hybrid_index))))

hypothetical_mating_dist_CI_medfemale <- ggplot(mating_pref_dens_medfemaleHI, aes(x = x * 100)) +
  labs(x = expression(atop("Mate Ancestry Fraction","(%"~italic("X. malinche")*")")), y = expression(atop("Relative", "Preference"))) +
  theme_cowplot() +
  theme(axis.text.y = element_blank()) +
  geom_density(aes(y = CALL), stat = "identity") +
  geom_density(aes(y = CALM), stat = "identity", lty = "dotted")  +
  geom_density(aes(y = CALM_noF129), stat = "identity", lty = "dashed")
hypothetical_mating_dist_CI_medfemale

mating_pref_dens_halfHI <- data.frame(x = seq(0, 1, by = .001), 
                                      CALM = dbeta(seq(0, 1, by = .001), accepted_CALM_20_MAP_sd0.002, beta(accepted_CALM_20_MAP_sd0.002, 0.5)),
                                      CALL = dbeta(seq(0, 1, by = .001), accepted_CALL_20_MAP_sd0.002, beta(accepted_CALL_20_MAP_sd0.002, 0.5)),
                                      CALM_noF129 = dbeta(seq(0, 1, by = .001), accepted_CALM_noF129_20_MAP_sd0.002, beta(accepted_CALM_noF129_20_MAP_sd0.002, 0.5)))



hypothetical_mating_dist_CI_halffemale <- ggplot(mating_pref_dens_halfHI, aes(x = x * 100)) +
  labs(x = expression(atop("Mate Ancestry Fraction","(%"~italic("X. malinche")*")")), y = expression(atop("Relative", "Preference"))) +
  theme_cowplot() +
  theme(axis.text.y = element_blank()) +
  geom_density(aes(y = CALL), stat = "identity") +
  geom_density(aes(y = CALM), stat = "identity", lty = "dotted")  +
  geom_density(aes(y = CALM_noF129), stat = "identity", lty = "dashed")
hypothetical_mating_dist_CI_halffemale

mating_pref_dens_minHI <- data.frame(x = seq(0, 1, by = .001), 
                                     CALM = dbeta(seq(0, 1, by = .001), accepted_CALM_20_MAP_sd0.002, beta(accepted_CALM_20_MAP_sd0.002, min(pop_index_CALL$hybrid_index))),
                                     CALL = dbeta(seq(0, 1, by = .001), accepted_CALL_20_MAP_sd0.002, beta(accepted_CALL_20_MAP_sd0.002, min(pop_index_CALL$hybrid_index))),
                                     CALM_noF129 = dbeta(seq(0, 1, by = .001), accepted_CALM_noF129_20_MAP_sd0.002, beta(accepted_CALM_noF129_20_MAP_sd0.002, min(pop_index_CALL$hybrid_index))))

hypothetical_mating_dist_CI_minfemale <- ggplot(mating_pref_dens_minHI, aes(x = x * 100)) +
  labs(x = expression(atop("Mate Ancestry Fraction","(%"~italic("X. malinche")*")")), y = expression(atop("Relative", "Preference"))) +
  theme_cowplot() +
  theme(axis.text.y = element_blank()) +
  geom_density(aes(y = CALL), stat = "identity") +
  geom_density(aes(y = CALM), stat = "identity", lty = "dotted")  +
  geom_density(aes(y = CALM_noF129), stat = "identity", lty = "dashed")
hypothetical_mating_dist_CI_minfemale

mating_pref_dens_maxHI <- data.frame(x = seq(0, 1, by = .001), 
                                     CALM = dbeta(seq(0, 1, by = .001), accepted_CALM_20_MAP_sd0.002, beta(accepted_CALM_20_MAP_sd0.002, max(pop_index_CALL$hybrid_index))),
                                     CALL = dbeta(seq(0, 1, by = .001), accepted_CALL_20_MAP_sd0.002, beta(accepted_CALL_20_MAP_sd0.002, max(pop_index_CALL$hybrid_index))),
                                     CALM_noF129 = dbeta(seq(0, 1, by = .001), accepted_CALM_noF129_20_MAP_sd0.002, beta(accepted_CALM_noF129_20_MAP_sd0.002, max(pop_index_CALL$hybrid_index))))

hypothetical_mating_dist_CI_maxfemale <- ggplot(mating_pref_dens_maxHI, aes(x = x * 100)) +
  labs(x = expression(atop("Mate Ancestry Fraction","(%"~italic("X. malinche")*")")), y = expression(atop("Relative", "Preference"))) +
  theme_cowplot() +
  theme(axis.text.y = element_blank()) +
  geom_density(aes(y = CALL), stat = "identity") +
  geom_density(aes(y = CALM), stat = "identity", lty = "dotted")  +
  geom_density(aes(y = CALM_noF129), stat = "identity", lty = "dashed")
hypothetical_mating_dist_CI_maxfemale

dummy_data <- data.frame(Estimate = rep(c("CALL","CALM","CALM filtered"), 2), x = c(0.5, 7.5, 3, 5, 8, 1), y = c(3, 3, 3, 8, 8, 3))
dummy_data$Estimate <- factor(dummy_data$Estimate, levels = c("CALL","CALM","CALM filtered"))
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
ggsave("MAP_assortative_mating_ABC_CALL_CALM_betadist_examples.pdf",
       female_pref_CI_plots_wleg, width = 6.5, height = 6.5, units = "in")

