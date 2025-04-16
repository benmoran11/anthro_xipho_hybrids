setwd("Documents/schumer_lab/")
library(tidyverse)

options(scipen=999)

infile <- "~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/genotypes_PLAZ.tsv_plot.txt"
ids<-"~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/PLAZ_id_list"
outfile<-"~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Data/ancestry_tracts_PLAZ.txt"

indivs<-read.table(file=ids,sep="\t",head=TRUE)
indivs<-indivs[,1]
indivs<-indivs[!grepl("PLAZ_61_F_read_1.fastq", indivs)]
indivs<-indivs[!grepl("PLAZ_64_F_read_1.fastq", indivs)]


data<-read.table(file=infile,sep="\t",head=TRUE,as.is=T)
data<-subset(data, select=-c(PLAZ_61_F_read_1.fastq,PLAZ_64_F_read_1.fastq))

data_diploid <- filter(data, !(chrom %in% ("chr-21-Y")))

alltracts <- {}


for(c in unique(data_diploid$chrom)) {
  

focalchr<-subset(data_diploid,data_diploid[,1]==c)


for(y in 1:length(indivs)){
  
  focalid=as.character(indivs[y])
  indiv <- na.omit(cbind(focalchr$pos,focalchr[focalid]))
  
  tract_lengths <- {}
  geno_curr = indiv[,2][1]
  start = indiv[,1][1]
  
  #write.table(cbind(geno_curr,start))
  
  for(x in 1:length(indiv[,1])){
    focal = indiv[,2][x]

    if(geno_curr != focal){
      stop = indiv[,1][x]-1
      tract_lengths <- rbind(tract_lengths,cbind(c, start,stop,geno_curr,focalid))
      start = indiv[,1][x]
      geno_curr = focal
    }
  }
  
  stop = indiv[,1][x]
  tract_lengths <- rbind(tract_lengths,cbind(c, start,stop,geno_curr,focalid))
  
  alltracts <- rbind(alltracts,tract_lengths)
  
}

}

write.table(alltracts,file=outfile,sep="\t",row.names=FALSE,col.names=c("chr","start","stop","ancestry","indiv"),quote=FALSE)


indiv_mean_genomewide <- as.data.frame(alltracts) %>% 
  mutate(start = as.numeric(start),
         stop = as.numeric(stop),
         geno_curr = as.numeric(geno_curr)) %>%
  group_by(focalid) %>%
  summarize(mean_anc_tract = weighted.mean(geno_curr, w = stop-start, na.rm = T) / 2)

AIMs <- read.csv("PLAZ_hybrid_index_readcounts.csv") %>%
  filter(!(X %in% c("PLAZ_61_F_read_1.fastq", "PLAZ_64_F_read_1.fastq"))) %>%
  rename(`focalid` = X) %>%
  left_join(indiv_mean_genomewide) %>%
  mutate(diff = mean_anc_tract - hybrid_index)

range(AIMs$diff)
mean(AIMs$diff)
mean(abs(AIMs$diff))


compare_ancsumms <- ggplot(AIMs, aes(x = hybrid_index, y = mean_anc_tract)) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme_bw() +
  labs(x = "Allele Count Method", y = "Ancestry Tract Method") +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  #geom_smooth(method = "lm", se = F) +
  geom_point()
compare_ancsumms
ggsave("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Population_structure_breakdown/Figures/compare_ancsumms.pdf",
       compare_ancsumms, width = 4, height = 4)
