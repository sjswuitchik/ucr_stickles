## in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr/vcf_filt
# transfer 012 matrix indv to local machine to add 0-59 as first field

## if doing on Cedar 
#conda create -n rv4 -c conda-forge r-base=4.3.3
#conda activate rv4
#R

library(sequoia) 
library(tidyverse)

df <- read_delim("stickles.filtered.012", delim = '\t', col_names = "num")
indv <- read_delim("stickles.filtered.012.indv.join", delim = '\t', col_names = c("num", "id"))
join <- left_join(indv, df, by = 'num') %>%
  select(-c(num)) %>%
  remove_rownames %>% 
  column_to_rownames(var="id")
# replace -1 from vcftools with -9 for sequoia 
join[ join == -1 ] <- -9

clean <- as.matrix(join)

# read in life history data - read_delim was being a bit weird so use read.delim instead
lh.df <- read.delim("~/Desktop/MRU_Faculty/Research/ucr_stickles/012_matrix/lh.df.tsv", sep = '\t') 

# run parental assignment 
par.out <- sequoia(GenoM = clean,
                  LifeHistData = lh.df,
                  Module = 'ped')
