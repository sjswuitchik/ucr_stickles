## in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr/linkageMapping

conda activate vcfFilt
vcftools --vcf stickles.filtered.recode.vcf --remove-indv dedup/OBBB_1.dedup.bam --remove-indv dedup/OOB_1.dedup.bam --recode --recode-INFO-all --out stickles.filt.f2s

module load StdEnv/2023 java/21.0.1
java -cp lepmap/bin/ IBD vcfFile=../stickles.filt.f2s.recode.vcf > stickles.filt.f2s.ibd

# separate by tank to attempt IBD again 
vcftools --vcf stickles.filt.f2s.recode.vcf --keep tank2.indv --recode --recode-INFO-all --out stickles.filt.f2s.tank2
vcftools --vcf stickles.filt.f2s.recode.vcf --keep tank7.indv --remove-indv dedup/BAM_11.dedup.bam --recode --recode-INFO-all --out stickles.filt.f2s.tank7

java -cp lepmap/bin/ IBD vcfFile=stickles.filt.f2s.tank7.recode.vcf > stickles.filt.f2s.tank7.ibd
sort -n -r -k 3,3 stickles.filt.f2s.tank7.ibd | less # all comps above 0.5

# with F1s
vcftools --vcf stickles.filtered.recode.vcf --remove-indv dedup/BAM_11.dedup.bam --keep tank7.indv --recode --recode-INFO-all --out stickles.filt.tank7
java -cp lepmap/bin/ IBD vcfFile=stickles.filt.tank7.recode.vcf > stickles.filt.tank7.ibd
sort -n -r -k 3,3 stickles.filt.tank7.ibd | less # nope

java -cp lepmap/bin/ IBD vcfFile=stickles.filt.f2s.tank2.recode.vcf > stickles.filt.f2s.tank2.ibd
sort -n -r -k 3,3 stickles.filt.f2s.tank2.ibd | less # worked better, but only 3 comps below 0.5


#######################################
# stuff that hasn't worked but I'm not ready to throw away yet 

## calculate IBD in PLINK for pedigree determination
#plink --vcf stickles.filtered.recode.vcf --allow-extra-chr --recode --out test
#plink --file test --allow-extra-chr --genome --out test2

#plink --vcf stickles.filtered.recode.vcf --allow-extra-chr --set-all-var-ids --recode --out stickles.plink
#plink --file stickles.plink --allow-extra-chr --indep-pairwise 100000 50000 0.5 --out stickles.plink

####### sequoia 
library(sequoia)
library(tidyverse)

## load in data to create genotype 012 matrix
df <- read_delim("~/Desktop/MRU_Faculty/Research/ucr_stickles/012_matrix/012_matrix_pruned/stickles.filtered.012", delim = '\t', col_names = "num")
indv <- read_delim("~/Desktop/MRU_Faculty/Research/ucr_stickles/012_matrix/012_matrix_pruned/stickles.filtered.012.indv.join", delim = '\t', col_names = c("num", "id"))
join <- left_join(indv, df, by = 'num') %>%
  select(-c(num)) %>%
  filter(id != 'dedup/BAM_11.dedup.bam') #%>% # low mappability
  filter(id != 'dedup/OBBB_1.dedup.bam') %>% # include these two males in tank7 filter
  filter(id != 'dedup/OOB_1.dedup.bam') 

## split by tanks 
tank2 <- join %>%
  filter(id == 'dedup/BAM_16.dedup.bam' | id == 'dedup/BAM_18.dedup.bam' | id == 'dedup/BAM_20.dedup.bam' | id == 'dedup/BAM_21.dedup.bam' | id == 'dedup/BAM_39.dedup.bam' | id == 'dedup/BAM_40.dedup.bam' | id == 'dedup/BAM_42.dedup.bam' | id == 'dedup/BAM_43.dedup.bam' | id == 'dedup/BAM_44.dedup.bam' | id == 'dedup/BAM_45.dedup.bam' | id == 'dedup/BAM_46.dedup.bam' | id == 'dedup/BAM_47.dedup.bam' | id == 'dedup/BAM_48.dedup.bam' | id == 'dedup/BAM_49.dedup.bam' | id == 'dedup/BAM_50.dedup.bam' | id == 'dedup/BAM_51.dedup.bam' | id == 'dedup/BAM_52.dedup.bam' | id == 'dedup/BAM_53.dedup.bam' | id == 'dedup/BAM_54.dedup.bam' | id == 'dedup/BAM_55.dedup.bam' | id == 'dedup/BAM_56.dedup.bam' | id == 'dedup/BAM_57.dedup.bam' | id == 'dedup/BAM_16.dedup.bam' | id == 'dedup/BAM_19.dedup.bam' | id == 'dedup/BAM_20.dedup.bam' | id == 'dedup/BAM_40.dedup.bam' | id == 'dedup/BAM_42.dedup.bam' | id == 'dedup/BAM_47.dedup.bam' | id == 'dedup/BAM_50.dedup.bam' | id == 'dedup/BAM_52.dedup.bam') %>%
  remove_rownames %>% 
  column_to_rownames(var="id")
# replace -1 from vcftools with -9 for sequoia 
tank2[ tank2 == -1 ] <- -9

# convert to matrix
clean2 <- as.matrix(tank2)

tank7 <- join %>%
  filter(id == 'dedup/BAM_01.dedup.bam' | id == 'dedup/BAM_02.dedup.bam' | id == 'dedup/BAM_03.dedup.bam' | id == 'dedup/BAM_04.dedup.bam' | id == 'dedup/BAM_05.dedup.bam' | id == 'dedup/BAM_06.dedup.bam' | id == 'dedup/BAM_07.dedup.bam' | id == 'dedup/BAM_08.dedup.bam' | id == 'dedup/BAM_09.dedup.bam' | id == 'dedup/BAM_10.dedup.bam' | id == 'dedup/BAM_11.dedup.bam' | id == 'dedup/BAM_12.dedup.bam' | id == 'dedup/BAM_13.dedup.bam' | id == 'dedup/BAM_14.dedup.bam' | id == 'dedup/BAM_15.dedup.bam' | id == 'dedup/BAM_22.dedup.bam' | id == 'dedup/BAM_23.dedup.bam' | id == 'dedup/BAM_24.dedup.bam' | id == 'dedup/BAM_25.dedup.bam' | id == 'dedup/BAM_26.dedup.bam' | id == 'dedup/BAM_27.dedup.bam' | id == 'dedup/BAM_28.dedup.bam' | id == 'dedup/BAM_29.dedup.bam' | id == 'dedup/BAM_30.dedup.bam' | id == 'dedup/BAM_31.dedup.bam' | id == 'dedup/BAM_32.dedup.bam' | id == 'dedup/BAM_33.dedup.bam' | id == 'dedup/BAM_34.dedup.bam' | id == 'dedup/BAM_35.dedup.bam' | id == 'dedup/BAM_36.dedup.bam' | id == 'dedup/BAM_37.dedup.bam' | id == 'dedup/BAM_38.dedup.bam' | id == 'dedup/BAM_41.dedup.bam' | id == 'dedup/BAM_58.dedup.bam' | id == 'dedup/BAM_11.dedup.bam' | id == 'dedup/BAM_23.dedup.bam' | id == 'dedup/BAM_25.dedup.bam' | id == 'dedup/BAM_26.dedup.bam' | id == 'dedup/BAM_28.dedup.bam' | id == 'dedup/BAM_34.dedup.bam' | id == 'dedup/BAM_35.dedup.bam' | id == 'dedup/BAM_38.dedup.bam' | id == 'dedup/BAM_58.dedup.bam' | id == 'dedup/OBBB_1.dedup.bam' | id == 'dedup/OOB_1.dedup.bam') %>%
  remove_rownames %>% 
  column_to_rownames(var="id")
# replace -1 from vcftools with -9 for sequoia 
tank7[ tank7 == -1 ] <- -9

# convert to matrix
clean7 <- as.matrix(tank7)

# read in life history data - read_delim was being a bit weird so use read.delim instead
lh.df.2 <- read.delim("~/Desktop/MRU_Faculty/Research/ucr_stickles/012_matrix/012_matrix_pruned/lh.df.2.txt", sep = '\t')

lh.df.7 <- read.delim("~/Desktop/MRU_Faculty/Research/ucr_stickles/012_matrix/012_matrix_pruned/lh.df.7.txt", sep = '\t')

# run parental assignment 
GetMaybeRel(GenoM = clean2, LifeHistData = lh.df.2, Module = 'ped', Complex = 'simp')

ped.2.out <- sequoia(GenoM = clean2,
                   LifeHistData = lh.df.2,
                   Plot = T,
                   Module = 'ped',
                   Complex = 'mono')
SummarySeq(ped.2.out)

GetMaybeRel(GenoM = clean7, LifeHistData = lh.df.7, Module = 'ped', Complex = 'simp')

ped.7.out <- sequoia(GenoM = clean7,
                     LifeHistData = lh.df.7,
                     Plot = T,
                     Module = 'ped',
                     Complex = 'simp')
SummarySeq(ped.7.out)

##################################

# download COLONY for Linux from https://www.zsl.org/about-zsl/resources/software/colony and transfer files to working dir

conda activate vcfFilt
vcftools --vcf stickles.filtered.recode.vcf --maf 0.01 --max-missing 1 --recode --recode-INFO-all --out stickles.prune
## kept 60/60 indv, 2040/16466 sites
vcftools --vcf stickles.prune.recode.vcf --012 --out stickles.filtered
# transfer 012 files to local machine to add 0-59 as first field of indv - saved as stickles.filtered.012.indv.join





