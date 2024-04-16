## in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr/vcf_filt
# transfer 012 matrix indv to local machine to add 0-59 as first field then back to Cedar

conda activate rv4
conda update -c conda-forge r-base

R
library(tidyverse)
df <- read_delim("stickles.filtered.012", delim = '\t', col_names = "num")
indv <- read_delim("stickles.filtered.012.indv.join", delim = '\t', col_names = c("num", "id"))
clean <- left_join(df, indv, by = 'num')
