#!/bin/bash

cd /location/of/ancestryinfer/run/

# Convert ancestry probabilities to hard-called local genotypes with 90% cutoff
perl /location/of/parsetsv_to_genotypes_Dec2017_v2.pl ancestry-probs-par1_allchrs.tsv ancestry-probs-par2_allchrs.tsv genotypes.tsv
# Convert ancestry probabilities to local genotypes and sum parental alleles across loci 
perl /location/of/parsetsv_ancestry_v2.pl ancestry-probs-par1_allchrs.tsv ancestry-probs-par2_allchrs.tsv > hybrid_index.tsv