## in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr/vcf_filt

conda activate vcfFilt

vcftools --vcf stickles.filtered.recode.vcf --012 --out stickles.filtered
# transfer 012 files to local machine to add 0-59 as first field of indv - saved as stickles.filtered.012.indv.join

####### if doing on Cedar 
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
  filter(id != 'dedup/BAM_11.dedup.bam') %>% #low mappability
  filter(id != 'dedup/OBBB_1.dedup.bam') %>% # f1
  filter(id != 'dedup/OOB_1.dedup.bam')      # f1

# replace -1 from vcftools with -9 
prune[ prune == -1 ] <- -9
write_delim(prune, "~/Desktop/stickles.012matrix.txt", delim = ' ', )

####
# sftp stickles.012matrix.txt to Cedar in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr/vcf_filt

git clone https://github.com/mariocalus/SNPrune.git
cd SNPrune
chmod +x SNPrune
./SNPrune



##############################

## calculate IBD in PLINK for pedigree determination
#plink --vcf stickles.filtered.recode.vcf --allow-extra-chr --recode --out test
#plink --file test --allow-extra-chr --genome --out test2

#plink --vcf stickles.filtered.recode.vcf --allow-extra-chr --set-all-var-ids --recode --out stickles.plink
#plink --file stickles.plink --allow-extra-chr --indep-pairwise 100000 50000 0.5 --out stickles.plink




