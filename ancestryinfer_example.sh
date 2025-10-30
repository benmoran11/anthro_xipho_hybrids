#!/bin/bash
#SBATCH --job-name=run_hmm
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH -p hns
#SBATCH --cpus-per-task=1
#SBATCH --mem=32000
#SBATCH --mail-user=benmoran@stanford.edu

# This script runs the Ancestry HMM from ancestryinfer (https://github.com/Schumerlab/ancestryinfer) based on a
# .cfg file (example provided as ancestryinfer_example.cfg)

module load python/3.9
module load boost/1.76.0
module load armadillo
module load biology
module load samtools
module load bcftools
module load bwa
module load gcc/10.1.0
module load gsl
module load R
module load java
export PATH="/home/location/of/Ancestry_HMM/src:$PATH"
export PATH="/home/location/of/shared_bin:$PATH"

perl /home/location/of/shared_bin/ancestryinfer_July2022/Ancestry_HMM_parallel_v7.pl ancestryinfer_example.cfg