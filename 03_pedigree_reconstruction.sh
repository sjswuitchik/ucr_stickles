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
clean <- left_join(indv, df, by = 'num') %>%
  select(-c(num))
clean[ clean == -1 ] <- -9
