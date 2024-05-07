## in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr/vcf_filt
# transfer 012 matrix indv to local machine to add 0-59 as first field - saved as stickles.filtered.012.indv.join

## if doing on Cedar 
#conda create -n rv4 -c conda-forge r-base=4.3.3
#conda activate rv4
#R

library(tidyverse)

## load in data to create genotype 012 matrix
df <- read_delim("~/Desktop/MRU_Faculty/Research/ucr_stickles/012_matrix/stickles.filtered.012", delim = '\t', col_names = "num")
indv <- read_delim("~/Desktop/MRU_Faculty/Research/ucr_stickles/012_matrix/stickles.filtered.012.indv.join", delim = '\t', col_names = c("num", "id"))

## write out 012 matrix for SNPrune
prune <- left_join(indv, df, by = 'num') %>%
  select(-c(num)) %>%
  filter(id != 'dedup/BAM_11.dedup.bam') %>%
  filter(id != 'dedup/OBBB_1.dedup.bam') %>%
  filter(id != 'dedup/OOB_1.dedup.bam')

# replace -1 from vcftools with -9 
prune[ prune == -1 ] <- -9
write_delim(prune, "~/Desktop/stickles.012matrix.txt", delim = ' ', )

