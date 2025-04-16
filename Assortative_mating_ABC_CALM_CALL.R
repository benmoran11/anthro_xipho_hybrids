library(bayestestR)
library(tidyverse)
library(cowplot)

maternal_offspring_index_CALL<-read.csv(file="~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/hybrid_index_CALL_mothers_embryos_Sept2023_match.txt",sep="\t",head=TRUE)
pop_index_CALL<-read.csv(file="~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/hybrid_index_CALL_adults_August2023",sep="\t",head=TRUE) %>%
  filter(malcount + bircount > 10000)
maternal_offspring_index_CALM <-read.csv(file="~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/hybrid_index_readcounts_CALM_07_08_09_mother_embryo_Feb2025.csv",head=TRUE) %>%
  filter(read1_counts > 200000, read1_counts.1 > 200000, heterzygosity < 0.75, heterzygosity.1 < 0.75,
         !(maternal_id %in% c("CALM-08-F-129.R1.fastq", "CALM-08-F-106.R1.fastq"))) %>% # Three individuals have heterozygosities that would be semi-plausible if it weren't for the fact that their embryos are near-pure malinche - must be contamination
  mutate(diff = hybrid_index.1 - hybrid_index)
pop_index_CALM<-read.csv(file="~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/huazalingo_pochula_calnali_conzintla_allhybridindex.csv",head=TRUE) %>%
  filter(malcount + bircount > 10000,
         Site.Code == "CALM")
females_CALM<-unique(maternal_offspring_index_CALM$maternal_id)
females_CALM_all<-unique(maternal_offspring_index_CALM_all$maternal_id)
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
mean(abs(idxsub_CALM$diff)) # 0.06062598
var(abs(idxsub_CALM$diff)) #  0.0001890117

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

mom_corrected_CALM <- lm(hybrid_index.1 ~ maternal_id, data = filter(maternal_offspring_index_CALM, !(maternal_id %in% c("CALM-08-F-101.R1.fastq", "CALM-08-F-107.R1.fastq", "CALM-18-07-F-12.R1.fastq","CALM-18-07-F-17.R1.fastq"))))
sd(mom_corrected_CALM$residuals)

### Function for running simulations

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

CALM_r=0.01947892
CALM_v=0.0002921749


#nulls_cutoff05 <- ABC_nulls_assortative_mating_cluster(seed = 33, assortative_cutoff = 0.05)
#write.csv(nulls_cutoff05_reps100K, "Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/nulls_cutoff_05_reps1M_cluster_seed33_CALM.csv", row.names = F)
#accepted_05 <- base::subset(nulls_cutoff05_reps100K[[1]],nulls_cutoff05_reps100K[[1]][,1]>(CALM_r-CALM_r*0.05) & nulls_cutoff05_reps100K[[1]][,1]< (CALM_r+CALM_r*0.05) & nulls_cutoff05_reps100K[[1]][,2]>(CALM_v-CALM_v*0.05) & nulls_cutoff05_reps100K[[1]][,2]< (CALM_v+CALM_v*0.05))
#accepted_05_MAP <- map_estimate(accepted_05[,'mating_prop'])$MAP_Estimate         # 0.65
#ci(accepted_05[,'mating_prop'], method = "ETI")                                   # 0.23 to 0.82

#nulls_CALM_cutoff05_reps100K_sd0.01 <- ABC_nulls_assortative_mating_cluster(mother_embryo_indices = maternal_offspring_index_CALM, 
#                                                                       population_indices = pop_index_CALM,
#                                                                       seed = 33, reps = 100000, assortative_cutoff = 0.05, sibling_sd = .01)
# #ggplot(as.data.frame(nulls_CALM_cutoff05_reps100K_sd0.01[[1]]), aes(x = V1, y = V2, color = mating_prop)) +
#  geom_point(alpha = 0.2) +
#  annotate(geom = "rect", xmax = CALM_r+CALM_r*0.05, xmin = CALM_r - CALM_r * 0.05, ymax = CALM_v+CALM_v*0.05, ymin = CALM_v - CALM_v * 0.05, color = "red")
#accepted_CALM_05_sd0.01 <- base::subset(nulls_CALM_cutoff05_reps100K_sd0.01[[1]],nulls_CALM_cutoff05_reps100K_sd0.01[[1]][,1]>(CALM_r-CALM_r*0.05) & nulls_CALM_cutoff05_reps100K_sd0.01[[1]][,1]< (CALM_r+CALM_r*0.05) & nulls_CALM_cutoff05_reps100K_sd0.01[[1]][,2]>(CALM_v-CALM_v*0.05) & nulls_CALM_cutoff05_reps100K_sd0.01[[1]][,2]< (CALM_v+CALM_v*0.05))
nulls_CALM_cutoff05_reps3M_sd0.01 <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/nulls_cluster_cutoff05_reps3M_sd0.01_combined_CALM.csv")
ggplot(sample_n(nulls_CALM_cutoff05_reps3M_sd0.01, size = 200000, replace = F), aes(x = sim_diff, y = sim_var, color = mating_prop)) +
  geom_point(alpha = 0.2) +
  annotate(geom = "rect", xmax = CALM_r+CALM_r*0.05, xmin = CALM_r - CALM_r * 0.05, ymax = CALM_v+CALM_v*0.05, ymin = CALM_v - CALM_v * 0.05, color = "red")
accepted_CALM_05_sd0.01 <-filter(nulls_CALM_cutoff05_reps3M_sd0.01,sim_diff > (CALM_r-CALM_r*0.05), sim_diff < (CALM_r+CALM_r*0.05), sim_var >(CALM_v-CALM_v*0.05), sim_var < (CALM_v+CALM_v*0.05))
accepted_CALM_05_MAP_sd0.01 <- map_estimate(accepted_CALM_05_sd0.01[,'mating_prop'])$MAP_Estimate         # 0.94
ci(accepted_CALM_05_sd0.01[,'mating_prop'], method = "ETI")                                   # 0.66 to 0.99
hist(accepted_CALM_05_sd0.01[,'mating_prop'], xlim = c(0,1))

nulls_CALM_cutoff05_reps3M_sd0.002 <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/nulls_cluster_cutoff05_reps3M_sd0.002_combined_CALM.csv")
ggplot(sample_n(nulls_CALM_cutoff05_reps3M_sd0.002, size = 200000, replace = F), aes(x = sim_diff, y = sim_var, color = mating_prop)) +
  geom_point(alpha = 0.2) +
  annotate(geom = "rect", xmax = CALM_r+CALM_r*0.05, xmin = CALM_r - CALM_r * 0.05, ymax = CALM_v+CALM_v*0.05, ymin = CALM_v - CALM_v * 0.05, color = "red")
accepted_CALM_05_sd0.002 <-filter(nulls_CALM_cutoff05_reps3M_sd0.002,sim_diff > (CALM_r-CALM_r*0.05), sim_diff < (CALM_r+CALM_r*0.05), sim_var >(CALM_v-CALM_v*0.05), sim_var < (CALM_v+CALM_v*0.05))
accepted_CALM_05_MAP_sd0.002 <- map_estimate(accepted_CALM_05_sd0.002[,'mating_prop'])$MAP_Estimate         # 0.86
ci(accepted_CALM_05_sd0.002[,'mating_prop'], method = "ETI")                                   # 0.51 to 0.97
hist(accepted_CALM_05_sd0.002[,'mating_prop'], xlim = c(0,1))


CALL_r=0.066
CALL_v=0.0046

#nulls_CALL_cutoff05_reps100K_sd0.002 <- ABC_nulls_assortative_mating_cluster(mother_embryo_indices = maternal_offspring_index_CALL, 
#                                                                            population_indices = pop_index_CALL,
#                                                                            seed = 45, reps = 100000, assortative_cutoff = 0.05, sibling_sd = .002)
nulls_CALL_cutoff05_reps3M_sd0.01 <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/nulls_cluster_cutoff05_reps3M_sd0.01_combined_CALL.csv")
ggplot(sample_n(nulls_CALL_cutoff05_reps3M_sd0.01, size = 100000, replace = F), aes(x = sim_diff, y = sim_var, color = mating_prop)) +
  geom_point(alpha = 0.2) +
  annotate(geom = "rect", xmax = CALL_r+CALL_r*0.05, xmin = CALL_r - CALL_r * 0.05, ymax = CALL_v+CALL_v*0.05, ymin = CALL_v - CALL_v * 0.05, color = "red")
accepted_CALL_05_sd0.01 <-filter(nulls_CALL_cutoff05_reps3M_sd0.01,sim_diff > (CALL_r-CALL_r*0.05), sim_diff < (CALL_r+CALL_r*0.05), sim_var >(CALL_v-CALL_v*0.05), sim_var < (CALL_v+CALL_v*0.05))
accepted_CALL_05_MAP_sd0.01 <- map_estimate(accepted_CALL_05_sd0.01[,'mating_prop'])$MAP_Estimate         # 0.70
ci(accepted_CALL_05_sd0.01[,'mating_prop'], method = "ETI")                                   # 0.25 to 0.84
hist(accepted_CALL_05_sd0.01[,'mating_prop'], xlim = c(0,1))

#nulls_CALL_cutoff05_reps3M_sd0.01 <- ABC_nulls_assortative_mating_cluster(mother_embryo_indices = maternal_offspring_index_CALL, 
#                                                                            population_indices = pop_index_CALL,
#                                                                            seed = 67, reps = 200000, assortative_cutoff = 0.05, sibling_sd = .01)
nulls_CALL_cutoff05_reps3M_sd0.002 <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/nulls_cluster_cutoff05_reps3M_sd0.002_combined_CALL.csv")
ggplot(sample_n(nulls_CALL_cutoff05_reps3M_sd0.002, 100000, replace = F), aes(x = sim_diff, y = sim_var, color = mating_prop)) +
  geom_point(alpha = 0.2) +
  annotate(geom = "rect", xmax = CALL_r+CALL_r*0.05, xmin = CALL_r - CALL_r * 0.05, ymax = CALL_v+CALL_v*0.05, ymin = CALL_v - CALL_v * 0.05, color = "red")
accepted_CALL_05_sd0.002 <-filter(nulls_CALL_cutoff05_reps3M_sd0.002,sim_diff > (CALL_r-CALL_r*0.05), sim_diff < (CALL_r+CALL_r*0.05), sim_var >(CALL_v-CALL_v*0.05), sim_var < (CALL_v+CALL_v*0.05))
accepted_CALL_05_MAP_sd0.002 <- map_estimate(accepted_CALL_05_sd0.002[,'mating_prop'])$MAP_Estimate         # 0.68
ci(accepted_CALL_05_sd0.002[,'mating_prop'], method = "ETI")                                   # 0.20 to 0.82
hist(accepted_CALL_05_sd0.002[,'mating_prop'], xlim = c(0,1))

nulls_CALL_cutoff05_reps3M_sd0.01 <- read.csv("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/nulls_cluster_cutoff05_reps3M_sd0.01_combined_CALL.csv")
accepted_CALL_05_sd0.01 <-filter(nulls_CALL_cutoff05_reps3M_sd0.01,sim_diff > (CALL_r-CALL_r*0.05), sim_diff < (CALL_r+CALL_r*0.05), sim_var >(CALL_v-CALL_v*0.05), sim_var < (CALL_v+CALL_v*0.05))
accepted_CALL_05_MAP_sd0.01 <- map_estimate(accepted_CALL_05_sd0.01[,'mating_prop'])$MAP_Estimate         # 0.70
ci(accepted_CALL_05_sd0.01[,'mating_prop'], method = "ETI")                                   # 0.25 to 0.84
hist(accepted_CALL_05_sd0.01[,'mating_prop'], xlim = c(0,1))

mating_prop_diff <- accepted_CALM_05_sd0.01[1:min(nrow(accepted_CALM_05_sd0.01),nrow(accepted_CALL_05_sd0.01)),] %>%
  mutate(mating_diff = mating_prop - accepted_CALL_05_sd0.01[1:min(nrow(accepted_CALM_05_sd0.01),nrow(accepted_CALL_05_sd0.01)),"mating_prop"])
accepted_mating_diff_sd0.01 <- map_estimate(mating_prop_diff$mating_diff)$MAP_Estimate         # 0.196
ci(mating_prop_diff$mating_diff, method = "ETI")                                   # -0.05 to 0.63
hist(mating_prop_diff$mating_diff)

mating_prop_diff_sd0.002 <- accepted_CALM_05_sd0.002[1:min(nrow(accepted_CALM_05_sd0.002),nrow(accepted_CALL_05_sd0.002)),] %>%
  mutate(mating_diff = mating_prop - accepted_CALL_05_sd0.002[1:min(nrow(accepted_CALM_05_sd0.002),nrow(accepted_CALL_05_sd0.002)),"mating_prop"])
accepted_mating_diff_sd0.002 <- map_estimate(mating_prop_diff_sd0.002$mating_diff)$MAP_Estimate         # 0.193
ci(mating_prop_diff_sd0.002$mating_diff, method = "ETI")                                   # -0.09 to 0.58
hist(mating_prop_diff_sd0.002$mating_diff)


nulls_CALL_cutoff05_reps100K_sd0.01 <- ABC_nulls_assortative_mating_cluster(mother_embryo_indices = maternal_offspring_index_CALL, 
                                                                            population_indices = pop_index_CALL,
                                                                            seed = 45, reps = 100000, assortative_cutoff = 0.05, sibling_sd = .01)


posterior_matecutoff05_tol05 <- ggplot(as.data.frame(accepted_05), aes(x = mating_prop)) +
  theme_bw() +
  xlim(c(0.75,1)) +
  labs(y = "Frequency", x = expression("Simulated"~italic(P[assortative]))) +
  geom_histogram(fill = "cornflowerblue") +
  #geom_histogram(fill = "cornflowerblue", binwidth = 0.05) +
  geom_vline(xintercept = accepted_05_MAP, lty = 2, col = "red4") +
  geom_text(label = paste("MAP =",round(accepted_05_MAP,2)), col = "red4", x = accepted_05_MAP - .25, y = Inf, vjust = 2)
posterior_matecutoff05_tol05
ggsave("~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/posterior_distribution_assortative_mating_ABC_CALL_matecutoff05_tolerance05.pdf",
       posterior_matecutoff05_tol05, width = 3.5, height = 2, units = "in")

#pdf("~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/posterior_distribution_assortative_mating_ABC_CALL_matecutoff05_tolerance05.pdf",
#    width = 3.5, height = 3) 
#hist(accepted_05[,3],main="",xlab=expression("Simulated"~italic(P[assortative])),xlim = c(0,1), col="cornflowerblue")
#abline(v=0.65,lty=2,lwd=2,col="red4")
#dev.off()
#write.table(accepted_05,file="~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/posterior_distribution_assortative_mating_ABC_CALL_Sept2023_matecutoff05.txt",sep="\t",row.names=FALSE,col.names=FALSE)


stringent_05 <- base::subset(nulls_cutoff05[[1]],nulls_cutoff05[[1]][,1]>(r-r*0.025) & nulls_cutoff05[[1]][,1]< (r+r*0.025) & nulls_cutoff05[[1]][,2]>(v-v*0.025) & nulls_cutoff05[[1]][,2]< (v+v*0.025))
stringent_05_MAP <- map_estimate(stringent_05[,'mating_prop'])$MAP_Estimate         # 0.67
ci(stringent_05[,'mating_prop'], method = "ETI")   # 0.22 to 0.78

posterior_matecutoff05_tol025 <- ggplot(as.data.frame(stringent_05), aes(x = mating_prop)) +
  theme_bw() +
  xlim(c(0,1)) +
  labs(y = "Frequency", x = expression("Simulated"~italic(P[assortative]))) +
  geom_histogram(fill = "cornflowerblue", binwidth = 0.05) +
  geom_vline(xintercept = stringent_05_MAP, lty = 2, col = "red4") +
  geom_text(label = paste("MAP =",round(stringent_05_MAP,2)), col = "red4", x = stringent_05_MAP - .25, y = Inf, vjust = 2)
posterior_matecutoff05_tol025
ggsave("~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/posterior_distribution_assortative_mating_ABC_CALL_matecutoff05_tolerance025.pdf",
       posterior_matecutoff05_tol025, width = 3.5, height = 2, units = "in")

#pdf("~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/posterior_distribution_assortative_mating_ABC_CALL_matecutoff05_tolerance025.pdf",
#    width = 3.5, height = 3) 
#hist(stringent_05[,3],main="",xlab=expression("Simulated"~italic(P[assortative])),xlim = c(0,1), col="cornflowerblue")
#abline(v=0.67,lty=2,lwd=2,col="red4")
#dev.off()
#write.table(stringent_05,file="~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/posterior_distribution_assortative_mating_ABC_CALL_Sept2023_matecutoff05_tolerance025.txt",sep="\t",row.names=FALSE,col.names=FALSE)

relaxed_05 <- base::subset(nulls_cutoff05_reps100k[[1]],nulls_cutoff05_reps100k[[1]][,1]>(r-r*0.1) & nulls_cutoff05_reps100k[[1]][,1]< (r+r*0.1) & nulls_cutoff05_reps100k[[1]][,2]>(v-v*0.1) & nulls_cutoff05_reps100k[[1]][,2]< (v+v*0.1))
relaxed_05_MAP <- map_estimate(relaxed_05[,'mating_prop'])$MAP_Estimate         # 0.91
ci(relaxed_05[,'mating_prop'], method = "ETI")    # 0.77 to 0.99

posterior_matecutoff05_tol10 <- ggplot(as.data.frame(relaxed_05), aes(x = mating_prop)) +
  theme_bw() +
  xlim(c(0,1)) +
  labs(y = "Frequency", x = expression("Simulated"~italic(P[assortative]))) +
  geom_histogram(fill = "cornflowerblue", binwidth = 0.05) +
  geom_vline(xintercept = relaxed_05_MAP, lty = 2, col = "red4") +
  geom_text(label = paste("MAP =",round(relaxed_05_MAP,2)), col = "red4", x = relaxed_05_MAP - .25, y = Inf, vjust = 2)
posterior_matecutoff05_tol10
ggsave("~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/posterior_distribution_assortative_mating_ABC_CALL_matecutoff05_tolerance10.pdf",
       posterior_matecutoff05_tol10, width = 3.5, height = 2, units = "in")


#pdf("~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/posterior_distribution_assortative_mating_ABC_CALL_matecutoff05_tolerance10.pdf",
#    width = 3.5, height = 3) 
#hist(relaxed_05[,3],main="",xlab=expression("Simulated"~italic(P[assortative])),xlim = c(0,1), col="cornflowerblue")
#abline(v=0.71,lty=2,lwd=2,col="red4")
#dev.off()
#write.table(relaxed_05,file="~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/posterior_distribution_assortative_mating_ABC_CALL_Sept2023_matecutoff05_tolerance10.txt",sep="\t",row.names=FALSE,col.names=FALSE)


# Trying again with a different mating threshold
#nulls_cutoff10 <- ABC_nulls_assortative_mating_cluster(seed = 33, assortative_cutoff = 0.10)
accepted_10 <- base::subset(nulls_cutoff10[[1]],nulls_cutoff10[[1]][,1]>(r-r*0.05) & nulls_cutoff10[[1]][,1]< (r+r*0.05) & nulls_cutoff10[[1]][,2]>(v-v*0.05) & nulls_cutoff10[[1]][,2]< (v+v*0.05))
accepted_10_MAP <- map_estimate(accepted_10[,'mating_prop'])$MAP_Estimate         # 0.71
ci(accepted_10[,'mating_prop'], method = "ETI")   # 0.30 to 0.84

posterior_matecutoff10_tol05 <- ggplot(as.data.frame(accepted_10), aes(x = mating_prop)) +
  theme_bw() +
  xlim(c(0,1)) +
  labs(y = "Frequency", x = expression("Simulated"~italic(P[assortative]))) +
  geom_histogram(fill = "cornflowerblue", binwidth = 0.05) +
  geom_vline(xintercept = accepted_10_MAP, lty = 2, col = "red4") +
  geom_text(label = paste("MAP =",round(accepted_10_MAP,2)), col = "red4", x = accepted_10_MAP - .25, y = Inf, vjust = 2)
posterior_matecutoff10_tol05
ggsave("~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/posterior_distribution_assortative_mating_ABC_CALL_matecutoff10_tolerance05.pdf",
       posterior_matecutoff10_tol05, width = 3.5, height = 2, units = "in")

#pdf("~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/posterior_distribution_assortative_mating_ABC_CALL_matecutoff10_tolerance05.pdf",
#    width = 3.5, height = 3) 
#hist(accepted_10[,3],main="",xlab=expression("Simulated"~italic(P[assortative])),xlim = c(0,1), col="cornflowerblue")
#abline(v=0.71,lty=2,lwd=2,col="red4")
#dev.off()
#write.table(accepted_10,file="~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/posterior_distribution_assortative_mating_ABC_CALL_Sept2023_matecutoff10.txt",sep="\t",row.names=FALSE,col.names=FALSE)


#nulls_cutoff025_reps200k <- ABC_nulls_assortative_mating_cluster(seed = 33, reps = 200000, assortative_cutoff = 0.025)
#nulls_cutoff025 <- ABC_nulls_assortative_mating_cluster(seed = 33, reps = 100000, assortative_cutoff = 0.025)
accepted_025 <- base::subset(nulls_cutoff025_reps200k[[1]],nulls_cutoff025_reps200k[[1]][,1]>(r-r*0.05) & nulls_cutoff025_reps200k[[1]][,1]< (r+r*0.05) & nulls_cutoff025_reps200k[[1]][,2]>(v-v*0.05) & nulls_cutoff025_reps200k[[1]][,2]< (v+v*0.05))
accepted_025_MAP <- map_estimate(accepted_025[,'mating_prop'])$MAP_Estimate         # 0.65
ci(accepted_025[,'mating_prop'], method = "ETI")   # 0.16 to 0.85

posterior_matecutoff025_tol05 <- ggplot(as.data.frame(accepted_025), aes(x = mating_prop)) +
  theme_bw() +
  xlim(c(0,1)) +
  labs(y = "Frequency", x = expression("Simulated"~italic(P[assortative]))) +
  geom_histogram(fill = "cornflowerblue", binwidth = 0.05) +
  geom_vline(xintercept = accepted_025_MAP, lty = 2, col = "red4") +
  geom_text(label = paste("MAP =",round(accepted_025_MAP,2)), col = "red4", x = accepted_025_MAP - .25, y = Inf, vjust = 2)
posterior_matecutoff025_tol05
ggsave("~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/posterior_distribution_assortative_mating_ABC_CALL_matecutoff025_tolerance05.pdf",
       posterior_matecutoff025_tol05, width = 3.5, height = 2, units = "in")

#pdf("~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/posterior_distribution_assortative_mating_ABC_CALL_matecutoff025_tolerance05.pdf",
#    width = 3.5, height = 3) 
#hist(accepted_025[,3],main="",xlab=expression("Simulated"~italic(P[assortative])),xlim = c(0,1), col="cornflowerblue")
#abline(v=0.65,lty=2,lwd=2,col="red4")
#dev.off()
#write.table(accepted_025,file="~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/posterior_distribution_assortative_mating_ABC_CALL_Sept2023_matecutoff025.txt",sep="\t",row.names=FALSE,col.names=FALSE)

supplemental_ABC_distr <- plot_grid(posterior_matecutoff025_tol05, posterior_matecutoff10_tol05, posterior_matecutoff05_tol025, posterior_matecutoff05_tol10, align = "hv", nrow = 2, labels = "AUTO")
supplemental_ABC_distr
ggsave("~/Swordtail\ Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/posterior_distribution_assortative_mating_ABC_CALL_SuppFigs.pdf",
       supplemental_ABC_distr, width = 7, height = 4, units = "in")

