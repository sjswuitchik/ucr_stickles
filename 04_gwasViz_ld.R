# libraries & set up 
library(qqman)
library(kableExtra)
library(tidyverse)

# reading in all continuous GWAS results, cleaning, and saving as unique objects
file_dir <- "~/Desktop/MRU_Faculty/Research/ucr_stickles/gwas_results"
file_paths <- list.files(file_dir, pattern="gasAcu.plink.*.glm.linear", full.names = TRUE)


for (file in file_paths) {
  
  name <- str_split(basename(file), "\\.", simplify = TRUE)[, 3]
  
  data <- read_table(file) %>%
    rename(CHR = "#CHROM") %>%
    select(CHR, POS, P) %>%
    filter((CHR <= 22 | is.na(CHR))) %>%
    mutate(snp = c(1:15861)) %>%
    filter(CHR != 1.1 | is.na(CHR)) %>%
    filter(CHR != 21.1 | is.na(CHR)) %>%
    mutate(CHR = replace_na(CHR, 23)) %>%
    rename(chr = CHR, pos = POS, p = P, SNP = snp)
  
  assign(name, data, envir = .GlobalEnv)  
}

# reading in binary GWAS results & cleaning
sf <- read_table("gwas_results/gasAcu.plink.sf.glm.logistic") %>%
  rename(CHR = "#CHROM") %>%
  select(CHR, POS, P) %>%
  filter((CHR <= 22 | is.na(CHR))) %>%
  mutate(snp = c(1:15861)) %>%
  filter(CHR != 1.1 | is.na(CHR)) %>%
  filter(CHR != 21.1 | is.na(CHR)) %>%
  mutate(CHR = replace_na(CHR, 23)) %>%
  rename(chr = CHR, pos = POS, p = P, SNP = snp)

# Manhattan plots & tables of SNPs with genome-wide significance (where relevant)
threshold <- 5E-8

manhattan(dist, main = "Dist", chr = "chr", bp = "pos", p = "p", snp = "SNP", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # no genome-wide sig (GWS)

manhattan(maxCranElev, main = "Max Cran Elev", chr = "chr", bp = "pos", p = "p", snp = "SNP", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # genome wide sig

maxCranElev %>%
  filter(p < threshold) %>%
  kbl() %>%
  kable_minimal()

manhattan(maxDecel, main = "Max Decel", chr = "chr", bp = "pos", p = "p", snp = "SNP", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # GWS

maxDecel %>%
  filter(p < threshold) %>%
  kbl() %>%
  kable_minimal()

manhattan(maxGape, main = "Max Gape", chr = "chr", bp = "pos", p = "p", snp = "SNP", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # no GWS

manhattan(maxHD, main = "Max HD", chr = "chr", bp = "pos", p = "p", snp = "SNP", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # no GWS

manhattan(maxJP, main = "Max JP", chr = "chr", bp = "pos", p = "p", snp = "SNP", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # no GWS

manhattan(PPD_MG, main = "PPD_MG", chr = "chr", bp = "pos", p = "p", snp = "SNP", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # GWS

PPD_MG %>%
  filter(p < threshold) %>%
  kbl() %>%
  kable_minimal()

manhattan(PPD_SI, main = "PPD SI", chr = "chr", bp = "pos", p = "p", snp = "SNP", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # no GWS

manhattan(ramSpeed, main = "Ram Speed", chr = "chr", bp = "pos", p = "p", snp = "SNP", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # GWS

ramSpeed %>%
  filter(p < threshold) %>%
  kbl() %>%
  kable_minimal()

manhattan(sf, main = "Success/Failure", chr = "chr", bp = "pos", p = "p", snp = "SNP", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # no GWS

manhattan(time_HDvMG, main = "Time - HD v MG", chr = "chr", bp = "pos", p = "p", snp = "SNP", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # GWS

time_HDvMG %>%
  filter(p < threshold) %>%
  kbl() %>%
  kable_minimal()

manhattan(ttpg, main = "TTPG", chr = "chr", bp = "pos", p = "p", snp = "SNP", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # GWS

ttpg %>%
  filter(p < threshold) %>%
  kbl() %>%
  kable_minimal()
