library(tidyverse)
library(cowplot)
library(scales)
options(scipen=999)

#### Part 1: effects of methods for converting hard calls to hybrid index on estimated genome-wide ancestry



### Placeholder variables to be changed to match position of files locally
infile <- "genotypes_PLAZ.tsv_plot.txt"
outfile<-"ancestry_tracts_PLAZ.txt"

### Load local ancestry data from PLAZ
data<-read.table(file=infile,sep="\t",head=TRUE,as.is=T)
### Remove individuals that failed to sequence
data<-subset(data, select=-c(PLAZ_61_F_read_1.fastq,PLAZ_64_F_read_1.fastq))

data_diploid <- filter(data, !(chrom %in% ("chr-21-Y")))

indivs <- colnames(data)[-c(1,2)]

# Now calculate tracts based on genotypes at each ancestry-informative marker
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

# Summary stats of tracts
indiv_mean_genomewide <- as.data.frame(alltracts) %>% 
  mutate(start = as.numeric(start),
         stop = as.numeric(stop),
         geno_curr = as.numeric(geno_curr)) %>%
  group_by(focalid) %>%
  summarize(mean_anc_tract = weighted.mean(geno_curr, w = stop-start, na.rm = T) / 2)

# now hybrid indexes based on counts of AIMs - to compare to tract-based estimates
AIMs <- read.csv("PLAZ_hybrid_index_readcounts.csv") %>%
  filter(!(indiv %in% c("PLAZ_61_F_read_1.fastq", "PLAZ_64_F_read_1.fastq"))) %>%
  rename(`focalid` = indiv) %>%
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
ggsave("compare_ancsumms.pdf",
       compare_ancsumms, width = 4, height = 4)


#### Part 2: effects of sequencing depth and library type on genome-wide ancestry estimates

# Loading in hybrid index file with downsampled individuals from 150bp paired-end and 100bp single-read libraries
downsampled_ancestrycalls <- read_tsv("hybrid_index_downsampling_PE.tsv") %>%
  mutate(indiv = sapply(readfile, function(f) strsplit(f, "_")[[1]][1]) %>% str_remove(".R1.fastq"),
         readtype = "PE",
         reads = as.numeric(sapply(readfile, function(f) strsplit(f, "_")[[1]][2]) %>% str_remove("reads.R1.fastq"))) %>%
  group_by(indiv) %>%
  filter(reads <= 2000000) %>%
  mutate(HI_drift = abs(hybrid_index - 0.5) - min(abs(hybrid_index - 0.5), na.rm = T)) # stat of variation from each individual's estimate with maximum minor parent ancestry (proxy for truth)

downsampled_ancestrycalls_SE <- read_tsv("hybrid_index_downsampling_SE.tsv") %>%
mutate(indiv = sapply(readfile, function(f) strsplit(f, "_")[[1]][1]),
       readtype = "SE") %>%
  select(everything(),reads) %>%
       group_by(indiv) %>%
  filter(reads <= 2000000) %>%
  mutate(HI_drift = abs(hybrid_index - 0.5) - min(abs(hybrid_index - 0.5), na.rm = T)) # stat of variation from each individual's estimate with maximum minor parent ancestry (proxy for truth)

downsampled_ancestrycalls_all <- rbind(downsampled_ancestrycalls,downsampled_ancestrycalls_SE)
max(filter(downsampled_ancestrycalls, reads == 300000)$HI_drift) # getting maximum deviation from truth proxy among samples at the lowest coverage allowed in the study
max(filter(downsampled_ancestrycalls_SE, reads == 300000)$HI_drift) # getting maximum deviation from truth proxy among samples at the lowest coverage allowed in the study
max(filter(downsampled_ancestrycalls_all, reads == 300000)$HI_drift) # getting maximum deviation from truth proxy among samples at the lowest coverage allowed in the study


covg_anc_var <- ggplot(downsampled_ancestrycalls_all, aes(x = reads, y = hybrid_index, color = indiv, shape = readtype)) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 3e5) +
  scale_x_continuous(labels = label_number_auto())+
  labs(x = "Read Count", y = expression("Proportion"~italic("X. malinche")~"Ancestry"))+
  geom_point()
covg_anc_var
ggsave("../Figures/covg_downsampling_hybrid_index.pdf",
       covg_anc_var, width = 4, height = 4)


technical_checks <- plot_grid(compare_ancsumms, covg_anc_var, nrow = 2, labels = c("A","B"))
ggsave("../Figures/HI_calcmethod_covg_sensitivity.pdf",
       technical_checks, width = 5, height = 8)

