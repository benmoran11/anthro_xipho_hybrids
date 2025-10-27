library(tidyverse)
library(ggbeeswarm)

# Loading and cleaning input files
maternal_offspring_index_CALL<-read.csv(file="hybrid_index_CALL_mother_embryo.csv",head=TRUE) # Pre-cleaned
maternal_offspring_index_CALM <-read.csv(file="hybrid_index_CALM_mother_embryo.csv",head=TRUE)  %>%
  filter(!(maternal_id %in% c("CALM-08-F-106.R1.fastq"))) #  mother has plausible moderate heterozygosity but the embryos are pure malinche - must be contamination
pop_index_CALM<-read.csv(file="hybrid_index_CALM_adults.csv",head=TRUE) %>%
  filter(malcount + bircount > 10000)
pop_index_CALL<-read.csv(file="hybrid_index_CALL_adults.csv",head=TRUE) %>%
  filter(malcount + bircount > 10000)

#generate random subset of mother-embryo pairs for multi embryo females
idxsub_CALL<-{}
females_CALL<-unique(maternal_offspring_index_CALL$maternal_id)
set.seed(14)
for(k in 1:length(females_CALL)){
  
  focal<-subset(maternal_offspring_index_CALL,maternal_offspring_index_CALL$maternal_id==females_CALL[k])
  rand<-focal[sample(nrow(focal),1),]
  idxsub_CALL<-rbind(idxsub_CALL,cbind(rand))
  
}
mean(idxsub_CALL$embryo_hybrid_index)
mean(abs(idxsub_CALL$diff)) #.0655
var(abs(idxsub_CALL$diff)) # .0046

idxsub_CALM<-{}
females_CALM<-unique(maternal_offspring_index_CALM$maternal_id)
set.seed(14)
for(k in 1:length(females_CALM)){
  
  focal<-subset(maternal_offspring_index_CALM,maternal_offspring_index_CALM$maternal_id==females_CALM[k])
  rand<-focal[sample(nrow(focal),1),]
  idxsub_CALM<-rbind(idxsub_CALM,cbind(rand))
  
}
mean(idxsub_CALM$embryo_hybrid_index)
mean(abs(idxsub_CALM$diff)) # 0.02568919
var(abs(idxsub_CALM$diff)) #  0.002615013

mean(abs(maternal_offspring_index_CALL$embryo_hybrid_index - maternal_offspring_index_CALL$hybrid_index))
mean(abs(maternal_offspring_index_CALM$embryo_hybrid_index - maternal_offspring_index_CALM$hybrid_index))


n_occur <- data.frame(table(maternal_offspring_index_CALL$maternal_id))
multiembryo_broods <- filter(maternal_offspring_index_CALL, maternal_id %in% filter(n_occur, Freq > 1)$Var1)

# P-value of observed CALL mother-embryo difference vs. the AGZC one from Schumer et al. 2017 using bootstrap

set.seed(14)
CALL_mother_offspring_resampled_means <- sapply(1:1000, function(i) {
  idxsub_CALL<-{}
  drawn_females = sample(females_CALL, size = length(females_CALL), replace = T)
  for(k in 1:length(drawn_females)){
    
    focal<-subset(maternal_offspring_index_CALL,maternal_offspring_index_CALL$maternal_id==drawn_females[k])
    rand<-focal[sample(nrow(focal),1),]
    idxsub_CALL<-rbind(idxsub_CALL,cbind(rand))
    
  }
  mean(abs(idxsub_CALL$diff))
})

CALL_empirical_p_motheroffspring <- ecdf(CALL_mother_offspring_resampled_means)
CALL_empirical_p_motheroffspring(0.022)

hist(CALL_mother_offspring_resampled_means, breaks = 50)



# Same for CALM

idxsub_CALM<-{}
females_CALM<-unique(maternal_offspring_index_CALM$maternal_id)
set.seed(14)
for(k in 1:length(females_CALM)){
  
  focal<-subset(maternal_offspring_index_CALM,maternal_offspring_index_CALM$maternal_id==females_CALM[k])
  rand<-focal[sample(nrow(focal),1),]
  idxsub_CALM<-rbind(idxsub_CALM,cbind(rand))
  
}
set.seed(14)
CALM_mother_offspring_resampled_means <- sapply(1:1000, function(i) {
  idxsub_CALM<-{}
  drawn_females = sample(females_CALM, size = length(females_CALM), replace = T)
  for(k in 1:length(drawn_females)){
    
    focal<-subset(maternal_offspring_index_CALM,maternal_offspring_index_CALM$maternal_id==drawn_females[k])
    rand<-focal[sample(nrow(focal),1),]
    idxsub_CALM<-rbind(idxsub_CALM,cbind(rand))
    
  }
  mean(abs(idxsub_CALM$diff))
})

CALM_empirical_p_motheroffspring <- ecdf(CALM_mother_offspring_resampled_means)
CALM_empirical_p_motheroffspring(0.022)

hist(CALM_mother_offspring_resampled_means, breaks = 50)

###########################################
#compare observed mother-embryo difference in CALL with random index expectations
###########################################


CALL_mates<-c(pop_index_CALL$hybrid_index)
CALL_offspring_index<-{}
CALL_maternal_index<-{}
set.seed(42)
for (x in 1:length(females_CALL)){
  maternal<-sample(idxsub_CALL$hybrid_index,1)
  mate<-sample(CALL_mates,1)
  off<-rnorm(1,mean=(maternal+mate)/2,sd=0.002)
  #average sd for the plot
  CALL_offspring_index<-c(CALL_offspring_index,off)
  CALL_maternal_index<-c(CALL_maternal_index,maternal)
}



par(mai=c(1.02,1.1,0.82,0.42))
plot(CALL_maternal_index*100,(CALL_maternal_index-CALL_offspring_index)*100,col=rgb(11/255,166/255,255/255,alpha=0.5),xlim=c(0.15*100,1*100),ylim=c(-0.5*100,0.5*100),cex=2,pch=20,xlab="Maternal index (% malinche)",ylab="Maternal index\n- offspring index",cex.lab=1.3,cex.axis=1.3)
points(idxsub_CALL$hybrid_index*100,(idxsub_CALL$diff)*100,pch=20,cex=2)
mean(idxsub_CALL$embryo_hybrid_index)


abline(h=0,lty=2)
legend("topleft", c("real data","random mating"), pch=c(20,20,17),cex=1.25,pt.cex=1.2,col=c("black",rgb(11/255,166/255,255/255,alpha=0.5)),box.col="white")


###############
# Same for CALM
###############

CALM_mates<-c(pop_index_CALM$hybrid_index)
CALM_offspring_index<-{}
CALM_maternal_index<-{}
set.seed(24)
for (x in 1:length(females_CALM)){
  maternal<-sample(idxsub_CALM$hybrid_index,1)
  mate<-sample(CALM_mates,1)
  off<-rnorm(1,mean=(maternal+mate)/2,sd=0.002)
  #average sd for the plot
  CALM_offspring_index<-c(CALM_offspring_index,off)
  CALM_maternal_index<-c(CALM_maternal_index,maternal)
}



par(mai=c(1.02,1.1,0.82,0.42))
plot(CALM_maternal_index*100,(CALM_maternal_index-CALM_offspring_index)*100,col=rgb(11/255,166/255,255/255,alpha=0.5),xlim=c(0.15*100,1*100),ylim=c(-0.5*100,0.5*100),cex=2,pch=20,xlab="Maternal index (% malinche)",ylab="Maternal index\n- offspring index",cex.lab=1.3,cex.axis=1.3)
points(idxsub_CALM$hybrid_index*100,(idxsub_CALM$diff)*100,pch=20,cex=2)
mean(idxsub_CALM$embryo_hybrid_index)


abline(h=0,lty=2)
legend("topleft", c("real data","random mating"), pch=c(20,20,17),cex=1.25,pt.cex=1.2,col=c("black",rgb(11/255,166/255,255/255,alpha=0.5)),box.col="white")





####################
#####generate distribution to compare to real data
####################

CALL_all_sims<-{}

set.seed(15)
for (j in 1:1000){
  
  
  #generate random subset for multi embryo females
  
  idxsub_CALL<-{}
  females_CALL<-unique(maternal_offspring_index_CALL$maternal_id)
  for(k in 1:length(females_CALL)){
    
    focal<-subset(maternal_offspring_index_CALL,maternal_offspring_index_CALL$maternal_id==females_CALL[k])
    rand<-focal[sample(nrow(focal),1),]
    idxsub_CALL<-rbind(idxsub_CALL,cbind(rand))
    
  }


  offspring_index<-{}
  maternal_index<-{}
  for (x in 1:length(females_CALL)){
    maternal<-sample(idxsub_CALL$hybrid_index,1)
    mate<-sample(CALL_mates,1)
    off<-rnorm(1,mean=(maternal+mate)/2,sd=0.002)
    #average sd for the plot
    offspring_index<-c(offspring_index,off)
    maternal_index<-c(maternal_index,maternal)
  }
  
  CALL_all_sims<-c(CALL_all_sims,mean(abs(maternal_index-offspring_index),na.rm=TRUE))
  
}

CALM_all_sims<-{}

set.seed(15)
for (j in 1:1000){
  
  
  #generate random subset for multi embryo females
  
  idxsub_CALM<-{}
  females_CALM<-unique(maternal_offspring_index_CALM$maternal_id)
  for(k in 1:length(females_CALM)){
    
    focal<-subset(maternal_offspring_index_CALM,maternal_offspring_index_CALM$maternal_id==females_CALM[k])
    rand<-focal[sample(nrow(focal),1),]
    idxsub_CALM<-rbind(idxsub_CALM,cbind(rand))
    
  }
  
  
  offspring_index<-{}
  maternal_index<-{}
  for (x in 1:length(females_CALM)){
    maternal<-sample(idxsub_CALM$hybrid_index,1)
    mate<-sample(CALM_mates,1)
    off<-rnorm(1,mean=(maternal+mate)/2,sd=0.002)
    #average sd for the plot
    offspring_index<-c(offspring_index,off)
    maternal_index<-c(maternal_index,maternal)
  }
  
  CALM_all_sims<-c(CALM_all_sims,mean(abs(maternal_index-offspring_index),na.rm=TRUE))
  
}



rand_diff <- ggplot(cbind(sim = 1:1000, diff = CALL_all_sims), aes(x = CALL_all_sims)) +
  theme_bw() +
  lims(x = c(0,.30)) +
  geom_density(fill = "salmon", alpha = 0.75) +
  geom_density(data = data.frame(sim = 1:1000, CALM_all_sims), aes(x = CALM_all_sims), fill = "cornflowerblue", alpha = 0.75) +
  geom_vline(xintercept = 0.06554228,lty=2,col="red4",lwd=1) +
  geom_vline(xintercept = 0.02568919,lty=2,col="blue4",lwd=1) +
  labs(x = "Average difference in\nmaternal-offspring ancestry", y = "Density")
  #geom_density(data = data.frame(CALL_mother_offspring_resampled_means), aes(x = CALL_mother_offspring_resampled_means), color = "cornflowerblue") +
  #geom_density(data = data.frame(CALM_mother_offspring_resampled_means), aes(x = CALM_mother_offspring_resampled_means), color = "salmon")

rand_diff
ggsave("Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/maternal_embryo_diff_random_true.pdf",
       rand_diff, width = 3, height = 3)




#### Miscellanea

# Exploring variance among siblings
brood_anc_variance <- maternal_offspring_index_CALL %>%
  group_by(maternal_id) %>%
  summarize(nobs = n(),
            brood_var = var(embryo_hybrid_index),
            brood_sd = sd(embryo_hybrid_index),
            brood_range = max(embryo_hybrid_index) - min(embryo_hybrid_index),
            maternal_hi = mean(hybrid_index))
median(brood_anc_variance$brood_sd, na.rm = T)
median(filter(brood_anc_variance, nobs > 2)$brood_sd, na.rm = T)

ggplot(multiembryo_broods, aes(x = maternal_id, y = embryo_hybrid_index)) + 
  geom_beeswarm() + 
stat_summary(geom = "pointrange", fun.data = "mean_sd")

# Recreating the estimate of full sibling variance from Schumer et al. 2017 PNAS
set.seed(77)
paired_anc_variance <- lapply(1:1000, function(i) {    # Note that the # iterations is  smaller than the # possible pairs (1020)
  e1 <- slice_sample(multiembryo_broods, n = 1)
  e2 <- slice_sample(filter(maternal_offspring_index_CALL, maternal_id == e1$maternal_id, offspring_id != e1$offspring_id), n = 1)
  return(var(c(e1$embryo_hybrid_index, e2$embryo_hybrid_index)))
}) %>%
  unlist()

sqrt(quantile(paired_anc_variance, probs = 0.5)) # The median of the ancestry variance distribution
sqrt(quantile(paired_anc_variance, probs = 0.05)) # The 5th percentile of the variance distribution (matching Schumer 2017)
# ^ used the above for setting sibling variance parameters in ABC

hist(paired_anc_variance, breaks = 50)

sib_var <- ecdf(paired_anc_variance)
plot(sib_var)



