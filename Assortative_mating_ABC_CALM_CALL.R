library(bayestestR)
library(tidyverse)
library(cowplot)

maternal_offspring_index_CALL<-read.csv(file="hybrid_index_CALL_mother_embryo.csv",head=TRUE) # Pre-cleaned
maternal_offspring_index_CALM <-read.csv(file="hybrid_index_CALM_mother_embryo.csv",head=TRUE)  %>%
  filter(!(maternal_id %in% c("CALM-08-F-106.R1.fastq"))) #  mother has plausible moderate heterozygosity but the embryos are pure malinche - must be contamination
pop_index_CALM<-read.csv(file="hybrid_index_CALM_adults.csv",head=TRUE) %>%
  filter(malcount + bircount > 10000)
pop_index_CALL<-read.csv(file="hybrid_index_CALL_adults.csv",head=TRUE) %>%
  filter(malcount + bircount > 10000)

females_CALM<-unique(maternal_offspring_index_CALM$maternal_id)
mates_CALM<-c(na.omit(pop_index_CALM$hybrid_index))


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
mean(abs(idxsub_CALM$diff)) # 0.02568919
var(abs(idxsub_CALM$diff)) #  0.002615013

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

### Getting intra-brood variation ###

sibling_sd_CALM <- maternal_offspring_index_CALM %>%
  group_by(maternal_id) %>%
  summarize(sibling_sd = sd(hybrid_index.1))
median(sibling_sd_CALM$sibling_sd, na.rm = T)
quantile(sibling_sd_CALM$sibling_sd, probs = 0.05, na.rm = T) # ~= 0.002, same as CALL



### Function for running assortative mating simulations

ABC_nulls_assortative_mating_cluster <- function(mother_embryo_indices, population_indices, seed = 33, reps = 100000, assortative_cutoff = 0.05, sibling_sd = 0.002) {
  females<-unique(mother_embryo_indices$maternal_id)
  mates<-c(na.omit(population_indices$hybrid_index))
  
  all_sims<-{}
  empirical<-{}
  set.seed(seed)
  for (j in 1:reps){
    
    mating_prop=runif(1,0,1)
    
    
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
      mate<-sample(mates,1)
      
      hold=0
      while(hold==0){
        if(abs(maternal-mate)<assortative_cutoff){
          hold=1
        }
        if(abs(maternal-mate)>assortative_cutoff){
          pref=rbinom(1,1,1-mating_prop)
          if(pref == 1){
            hold=1
          } else{
            mate<-sample(mates,1)
          }
        }
      }
      
      off<-rnorm(1,mean=(maternal+mate)/2,sd=sibling_sd)
      #average sd for the plot
      offspring_index<-c(offspring_index,off)
      maternal_index<-c(maternal_index,maternal)
    }
    
    all_sims<-rbind(all_sims,cbind(mean(abs(maternal_index-offspring_index),na.rm=TRUE),var(abs(maternal_index-offspring_index),na.rm=TRUE),mating_prop))
    
  }
  return(list(all_sims, empirical))
}


### Simulations were run in batches of 100000 simulations
### In the input files below, "cutoff" refers to the assortative_cutoff parameter,
### "sd" refers to the sibling_sd parameter, 
### and the random seed was set to the index of the simulation batch (from 1 to 30)
### for each population, "mother_embryo_indices" was the corresponding maternal_offspring_index_#### variable, and
### "population_indices" was the corresponding pop_index_#### variable

### Performing rejection sampling on simulations performed separately (parallelized on cluster)

CALM_r=0.02568919
CALM_v=0.002615013

CALL_r=0.066
CALL_v=0.0046

# First for CALM (Calnali Mid)
# Loading simulation results file
nulls_CALM_cutoff05_reps3M_sd0.002 <- read.csv("nulls_cluster_cutoff0.05_reps3M_sd0.002_CALM.csv")
# performing rejection sampling
accepted_CALM_05_sd0.002 <-filter(nulls_CALM_cutoff05_reps3M_sd0.002,sim_diff > (CALM_r-CALM_r*0.05), sim_diff < (CALM_r+CALM_r*0.05), sim_var >(CALM_v-CALM_v*0.05), sim_var < (CALM_v+CALM_v*0.05))
accepted_CALM_05_MAP_sd0.002 <- map_estimate(accepted_CALM_05_sd0.002[,'mating_prop'])$MAP_Estimate         # 0.97
ci(accepted_CALM_05_sd0.002[,'mating_prop'], method = "ETI")                                   # 0.79 to 0.99

# Now for CALL (Calnali Low)
# Loading simulation results file
nulls_CALL_cutoff05_reps3M_sd0.002 <- read.csv("nulls_cluster_cutoff0.05_reps3M_sd0.002_CALL.csv")
# performing rejection sampling
accepted_CALL_05_sd0.002 <-filter(nulls_CALL_cutoff05_reps3M_sd0.002,sim_diff > (CALL_r-CALL_r*0.05), sim_diff < (CALL_r+CALL_r*0.05), sim_var >(CALL_v-CALL_v*0.05), sim_var < (CALL_v+CALL_v*0.05))
accepted_CALL_05_MAP_sd0.002 <- map_estimate(accepted_CALL_05_sd0.002[,'mating_prop'])$MAP_Estimate         # 0.68
ci(accepted_CALL_05_sd0.002[,'mating_prop'], method = "ETI")                                   # 0.20 to 0.82

# Getting difference in inferred strength of assortative mating from random pairs of simulations accepted for CALL & CALM
mating_prop_diff_sd0.002 <- accepted_CALM_05_sd0.002[1:min(nrow(accepted_CALM_05_sd0.002),nrow(accepted_CALL_05_sd0.002)),] %>%
  mutate(mating_diff = mating_prop - accepted_CALL_05_sd0.002[1:min(nrow(accepted_CALM_05_sd0.002),nrow(accepted_CALL_05_sd0.002)),"mating_prop"])
accepted_mating_diff_sd0.002 <- map_estimate(mating_prop_diff_sd0.002$mating_diff)$MAP_Estimate         # 0.273
ci(mating_prop_diff_sd0.002$mating_diff, method = "ETI")                                   # 0.08 to 0.72
hist(mating_prop_diff_sd0.002$mating_diff)
wilcox.test(x = mating_prop_diff_sd0.002$mating_diff)

# Plotting
posterior_matecutoff05_tol05 <- ggplot(as.data.frame(accepted_CALM_05_sd0.002), aes(x = mating_prop)) +
  theme_bw() +
  lims(x =c(0,1), y = c(0, 15)) +
  labs(y = "Density", x = expression(italic(P[assortative]))) +
  geom_vline(xintercept = accepted_CALM_05_MAP_sd0.002, lty = 2, col = "blue4") +
  geom_vline(xintercept = accepted_CALL_05_MAP_sd0.002, lty = 2, col = "red4") +
  geom_density(fill = "cornflowerblue", alpha = 0.75) +
  geom_density(data = as.data.frame(accepted_CALL_05_sd0.002), fill = "salmon", alpha = 0.75) +
  geom_text(label = round(accepted_CALM_05_MAP_sd0.002,2), col = "blue4", x = accepted_CALM_05_MAP_sd0.002 - .1, y = Inf, vjust = 2) +
  geom_text(label = round(accepted_CALL_05_MAP_sd0.002,2), col = "red4", x = accepted_CALL_05_MAP_sd0.002 - .1, y = 7.5, vjust = 2)
posterior_matecutoff05_tol05
ggsave("posterior_distribution_assortative_mating_ABC_CALL_CALM_cutoff05_siblingSD0.002_tolerance05.pdf",
       posterior_matecutoff05_tol05, width = 2.5, height = 2.5, units = "in")

### Now with stricter tolerance in the rejection sampling

accepted_CALM_05_sd0.002_tol0.025 <-filter(nulls_CALM_cutoff05_reps3M_sd0.002,sim_diff > (CALM_r-CALM_r*0.025), sim_diff < (CALM_r+CALM_r*0.025), sim_var >(CALM_v-CALM_v*0.025), sim_var < (CALM_v+CALM_v*0.025))
accepted_CALM_05_MAP_sd0.002_tol0.025 <- map_estimate(accepted_CALM_05_sd0.002_tol0.025[,'mating_prop'])$MAP_Estimate         # 0.98
ci(accepted_CALM_05_sd0.002_tol0.025[,'mating_prop'], method = "ETI")                                   # 0.80 to 0.99

accepted_CALL_05_sd0.002_tol0.025 <-filter(nulls_CALL_cutoff05_reps3M_sd0.002,sim_diff > (CALL_r-CALL_r*0.025), sim_diff < (CALL_r+CALL_r*0.025), sim_var >(CALL_v-CALL_v*0.025), sim_var < (CALL_v+CALL_v*0.025))
accepted_CALL_05_MAP_sd0.002_tol0.025 <- map_estimate(accepted_CALL_05_sd0.002_tol0.025[,'mating_prop'])$MAP_Estimate         # 0.68
ci(accepted_CALL_05_sd0.002_tol0.025[,'mating_prop'], method = "ETI")                                   # 0.22 to 0.81

mating_prop_diff_sd0.002_tol0.025 <- accepted_CALM_05_sd0.002_tol0.025[1:min(nrow(accepted_CALM_05_sd0.002_tol0.025),nrow(accepted_CALL_05_sd0.002_tol0.025)),] %>%
  mutate(mating_diff = mating_prop - accepted_CALL_05_sd0.002_tol0.025[1:min(nrow(accepted_CALM_05_sd0.002_tol0.025),nrow(accepted_CALL_05_sd0.002_tol0.025)),"mating_prop"])
accepted_mating_diff_sd0.002_tol0.025 <- map_estimate(mating_prop_diff_sd0.002_tol0.025$mating_diff)$MAP_Estimate         # 0.294
ci(mating_prop_diff_sd0.002_tol0.025$mating_diff, method = "ETI")                                   # 0.10 to 0.71
hist(mating_prop_diff_sd0.002_tol0.025$mating_diff)

posterior_matecutoff05_tol025 <- ggplot(as.data.frame(accepted_CALM_05_sd0.002_tol0.025), aes(x = mating_prop)) +
  theme_bw() +
  lims(x =c(0,1), y = c(0, 15)) +
  labs(y = "Density", x = expression(italic(P[assortative]))) +
  geom_vline(xintercept = accepted_CALM_05_MAP_sd0.002_tol0.025, lty = 2, col = "blue4") +
  geom_vline(xintercept = accepted_CALL_05_MAP_sd0.002_tol0.025, lty = 2, col = "red4") +
  geom_density(fill = "cornflowerblue", alpha = 0.75) +
  geom_density(data = as.data.frame(accepted_CALL_05_sd0.002_tol0.025), fill = "salmon", alpha = 0.75) +
  geom_text(label = round(accepted_CALM_05_MAP_sd0.002_tol0.025,2), col = "blue4", x = accepted_CALM_05_MAP_sd0.002_tol0.025 - .1, y = Inf, vjust = 2) +
  geom_text(label = round(accepted_CALL_05_MAP_sd0.002_tol0.025,2), col = "red4", x = accepted_CALL_05_MAP_sd0.002_tol0.025 - .1, y = 7.5, vjust = 2)
posterior_matecutoff05_tol025
ggsave("posterior_distribution_assortative_mating_ABC_CALL_CALM_cutoff05_siblingSD0.002_tolerance025.pdf",
       posterior_matecutoff05_tol025, width = 3.5, height = 2, units = "in")


### Now with laxer tolerance in the rejection sampling

accepted_CALM_05_sd0.002_tol0.10 <-filter(nulls_CALM_cutoff05_reps3M_sd0.002,sim_diff > (CALM_r-CALM_r*0.10), sim_diff < (CALM_r+CALM_r*0.10), sim_var >(CALM_v-CALM_v*0.10), sim_var < (CALM_v+CALM_v*0.10))
accepted_CALM_05_MAP_sd0.002_tol0.10 <- map_estimate(accepted_CALM_05_sd0.002_tol0.10[,'mating_prop'])$MAP_Estimate         # 0.98
ci(accepted_CALM_05_sd0.002_tol0.10[,'mating_prop'], method = "ETI")                                   # 0.79 to 0.99

accepted_CALL_05_sd0.002_tol0.10 <-filter(nulls_CALL_cutoff05_reps3M_sd0.002,sim_diff > (CALL_r-CALL_r*0.10), sim_diff < (CALL_r+CALL_r*0.10), sim_var >(CALL_v-CALL_v*0.10), sim_var < (CALL_v+CALL_v*0.10))
accepted_CALL_05_MAP_sd0.002_tol0.10 <- map_estimate(accepted_CALL_05_sd0.002_tol0.10[,'mating_prop'])$MAP_Estimate         # 0.73
ci(accepted_CALL_05_sd0.002_tol0.10[,'mating_prop'], method = "ETI")                                   # 0.18 to 0.85

mating_prop_diff_sd0.002_tol0.10 <- accepted_CALM_05_sd0.002_tol0.10[1:min(nrow(accepted_CALM_05_sd0.002_tol0.10),nrow(accepted_CALL_05_sd0.002_tol0.10)),] %>%
  mutate(mating_diff = mating_prop - accepted_CALL_05_sd0.002_tol0.10[1:min(nrow(accepted_CALM_05_sd0.002_tol0.10),nrow(accepted_CALL_05_sd0.002_tol0.10)),"mating_prop"])
accepted_mating_diff_sd0.002_tol0.10 <- map_estimate(mating_prop_diff_sd0.002_tol0.10$mating_diff)$MAP_Estimate         # 0.220
ci(mating_prop_diff_sd0.002_tol0.10$mating_diff, method = "ETI")                                   # 0.04 to 0.75
hist(mating_prop_diff_sd0.002_tol0.10$mating_diff)

posterior_matecutoff05_tol10 <- ggplot(as.data.frame(accepted_CALM_05_sd0.002_tol0.10), aes(x = mating_prop)) +
  theme_bw() +
  lims(x =c(0,1), y = c(0, 15)) +
  labs(y = "Density", x = expression(italic(P[assortative]))) +
  geom_vline(xintercept = accepted_CALM_05_MAP_sd0.002_tol0.10, lty = 2, col = "blue4") +
  geom_vline(xintercept = accepted_CALL_05_MAP_sd0.002_tol0.10, lty = 2, col = "red4") +
  geom_density(fill = "cornflowerblue", alpha = 0.75) +
  geom_density(data = as.data.frame(accepted_CALL_05_sd0.002_tol0.10), fill = "salmon", alpha = 0.75) +
  geom_text(label = round(accepted_CALM_05_MAP_sd0.002_tol0.10,2), col = "blue4", x = accepted_CALM_05_MAP_sd0.002_tol0.10 - .1, y = Inf, vjust = 2) +
  geom_text(label = round(accepted_CALL_05_MAP_sd0.002_tol0.10,2), col = "red4", x = accepted_CALL_05_MAP_sd0.002_tol0.10 - .1, y = 7.5, vjust = 2)
posterior_matecutoff05_tol10
ggsave("posterior_distribution_assortative_mating_ABC_CALL_CALM_cutoff05_siblingSD0.002_tolerance10.pdf",
       posterior_matecutoff05_tol10, width = 3.5, height = 2, units = "in")


### Now with broader cutoff for mating treated as assortative

nulls_CALM_cutoff10_reps3M_sd0.002 <- read.csv("nulls_cluster_cutoff0.1_reps3M_sd0.002_CALM.csv")
accepted_CALM_10_sd0.002 <-filter(nulls_CALM_cutoff10_reps3M_sd0.002,sim_diff > (CALM_r-CALM_r*0.05), sim_diff < (CALM_r+CALM_r*0.05), sim_var >(CALM_v-CALM_v*0.05), sim_var < (CALM_v+CALM_v*0.05))
accepted_CALM_10_MAP_sd0.002 <- map_estimate(accepted_CALM_10_sd0.002[,'mating_prop'])$MAP_Estimate         # 0.99
ci(accepted_CALM_10_sd0.002[,'mating_prop'], method = "ETI")                                   # 0.81 to 1.00

nulls_CALL_cutoff10_reps3M_sd0.002 <- read.csv("nulls_cluster_cutoff0.1_reps3M_sd0.002_CALL.csv")
accepted_CALL_10_sd0.002 <-filter(nulls_CALL_cutoff10_reps3M_sd0.002,sim_diff > (CALL_r-CALL_r*0.05), sim_diff < (CALL_r+CALL_r*0.05), sim_var >(CALL_v-CALL_v*0.05), sim_var < (CALL_v+CALL_v*0.05))
accepted_CALL_10_MAP_sd0.002 <- map_estimate(accepted_CALL_10_sd0.002[,'mating_prop'])$MAP_Estimate         # 0.69
ci(accepted_CALL_10_sd0.002[,'mating_prop'], method = "ETI")                                   # 0.30 to 0.83

mating_prop_diff_10_sd0.002 <- accepted_CALM_10_sd0.002[1:min(nrow(accepted_CALM_10_sd0.002),nrow(accepted_CALL_10_sd0.002)),] %>%
  mutate(mating_diff = mating_prop - accepted_CALL_10_sd0.002[1:min(nrow(accepted_CALM_10_sd0.002),nrow(accepted_CALL_10_sd0.002)),"mating_prop"])
accepted_mating_diff_cutoff10_sd0.002 <- map_estimate(mating_prop_diff_10_sd0.002$mating_diff)$MAP_Estimate         # 0.254
ci(mating_prop_diff_10_sd0.002$mating_diff, method = "ETI")                                   # 0.07 to 0.67
hist(mating_prop_diff_sd0.002$mating_diff)

posterior_matecutoff10_tol05 <- ggplot(as.data.frame(accepted_CALM_10_sd0.002), aes(x = mating_prop)) +
  theme_bw() +
  lims(x =c(0,1), y = c(0, 15)) +
  labs(y = "Density", x = expression(italic(P[assortative]))) +
  geom_vline(xintercept = accepted_CALM_10_MAP_sd0.002, lty = 2, col = "blue4") +
  geom_vline(xintercept = accepted_CALL_10_MAP_sd0.002, lty = 2, col = "red4") +
  geom_density(fill = "cornflowerblue", alpha = 0.75) +
  geom_density(data = as.data.frame(accepted_CALL_10_sd0.002), fill = "salmon", alpha = 0.75) +
  geom_text(label = round(accepted_CALM_10_MAP_sd0.002,2), col = "blue4", x = accepted_CALM_10_MAP_sd0.002 - .1, y = Inf, vjust = 2) +
  geom_text(label = round(accepted_CALL_10_MAP_sd0.002,2), col = "red4", x = accepted_CALL_10_MAP_sd0.002 - .1, y = 7.5, vjust = 2)
posterior_matecutoff10_tol05
ggsave("posterior_distribution_assortative_mating_ABC_CALL_CALM_cutoff10_siblingSD0.002_tolerance05.pdf",
       posterior_matecutoff10_tol05, width = 3.5, height = 2, units = "in")

### Now with tighter cutoff of mating treated as assortative

nulls_CALM_cutoff025_reps3M_sd0.002 <- read.csv("nulls_cluster_cutoff0.025_reps3M_sd0.002_CALM.csv")
accepted_CALM_025_sd0.002 <-filter(nulls_CALM_cutoff025_reps3M_sd0.002,sim_diff > (CALM_r-CALM_r*0.05), sim_diff < (CALM_r+CALM_r*0.05), sim_var >(CALM_v-CALM_v*0.05), sim_var < (CALM_v+CALM_v*0.05))
accepted_CALM_025_MAP_sd0.002 <- map_estimate(accepted_CALM_025_sd0.002[,'mating_prop'])$MAP_Estimate         # 0.95
ci(accepted_CALM_025_sd0.002[,'mating_prop'], method = "ETI")                                   # 0.77 to 0.99

nulls_CALL_cutoff025_reps3M_sd0.002 <- read.csv("nulls_cluster_cutoff0.025_reps3M_sd0.002_CALL.csv")
accepted_CALL_025_sd0.002 <-filter(nulls_CALL_cutoff025_reps3M_sd0.002,sim_diff > (CALL_r-CALL_r*0.05), sim_diff < (CALL_r+CALL_r*0.05), sim_var >(CALL_v-CALL_v*0.05), sim_var < (CALL_v+CALL_v*0.05))
accepted_CALL_025_MAP_sd0.002 <- map_estimate(accepted_CALL_025_sd0.002[,'mating_prop'])$MAP_Estimate         # 0.74
ci(accepted_CALL_025_sd0.002[,'mating_prop'], method = "ETI")                                   # 0.15 to 0.86

mating_prop_diff_025_sd0.002 <- accepted_CALM_025_sd0.002[1:min(nrow(accepted_CALM_025_sd0.002),nrow(accepted_CALL_025_sd0.002)),] %>%
  mutate(mating_diff = mating_prop - accepted_CALL_025_sd0.002[1:min(nrow(accepted_CALM_025_sd0.002),nrow(accepted_CALL_025_sd0.002)),"mating_prop"])
accepted_mating_diff_cutoff025_sd0.002 <- map_estimate(mating_prop_diff_025_sd0.002$mating_diff)$MAP_Estimate         # 0.237
ci(mating_prop_diff_025_sd0.002$mating_diff, method = "ETI")                                   # 0.05 to 0.79
hist(mating_prop_diff_sd0.002$mating_diff)

posterior_matecutoff025_tol05 <- ggplot(as.data.frame(accepted_CALM_025_sd0.002), aes(x = mating_prop)) +
  theme_bw() +
  lims(x =c(0,1), y = c(0, 15)) +
  labs(y = "Density", x = expression(italic(P[assortative]))) +
  geom_vline(xintercept = accepted_CALM_025_MAP_sd0.002, lty = 2, col = "blue4") +
  geom_vline(xintercept = accepted_CALL_025_MAP_sd0.002, lty = 2, col = "red4") +
  geom_density(fill = "cornflowerblue", alpha = 0.75) +
  geom_density(data = as.data.frame(accepted_CALL_025_sd0.002), fill = "salmon", alpha = 0.75) +
  geom_text(label = round(accepted_CALM_025_MAP_sd0.002,2), col = "blue4", x = accepted_CALM_025_MAP_sd0.002 - .1, y = Inf, vjust = 2) +
  geom_text(label = round(accepted_CALL_025_MAP_sd0.002,2), col = "red4", x = accepted_CALL_025_MAP_sd0.002 - .1, y = 7.5, vjust = 2)
posterior_matecutoff025_tol05
ggsave("posterior_distribution_assortative_mating_ABC_CALL_CALM_cutoff025_siblingSD0.002_tolerance05.pdf",
       posterior_matecutoff025_tol05, width = 3.5, height = 2, units = "in")


supplemental_ABC_distr <- plot_grid(posterior_matecutoff025_tol05, posterior_matecutoff10_tol05, posterior_matecutoff05_tol025, posterior_matecutoff05_tol10, align = "hv", nrow = 2, labels = "AUTO")
supplemental_ABC_distr
ggsave("posterior_distribution_assortative_mating_ABC_CALL_CALM_SuppFigs.pdf",
       supplemental_ABC_distr, width = 7, height = 4, units = "in")

