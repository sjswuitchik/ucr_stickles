library(qqman)
library(tidyverse)


# read in each continuous GWAS result & clean

setwd("~/Desktop/MRU_Faculty/Research/ucr_stickles/gwas_results")

file_dir <- "~/Desktop/MRU_Faculty/Research/ucr_stickles/gwas_results"
file_paths <- list.files(file_dir, pattern="gasAcu.plink.*.glm.linear", full.names = TRUE)


for (file in file_paths) {
  
  name <- str_split(basename(file), "\\.", simplify = TRUE)[, 3]
  
  data <- read.table(file, col.names = c("CHR",	"POS", "ID", "REF",	"ALT", "A1",	"TEST",	"OBS_CT",	"OR",	"LOG(OR)_SE",	"Z_STAT",	"P")) %>%
    dplyr::select(CHR, POS, ID, P) %>%
    filter((CHR <= 22 | is.na(CHR))) %>%
    filter(CHR != 1.1 | is.na(CHR)) %>%
    filter(CHR != 21.1 | is.na(CHR)) %>%
    mutate(CHR = replace_na(CHR, 23))
  
  assign(name, data, envir = .GlobalEnv)
}

# read in binary GWAS result

sf <- read.table("gasAcu.plink.sf.glm.logistic", col.names = c("CHR",	"POS", "ID", "REF",	"ALT", "A1",	"TEST",	"OBS_CT",	"OR",	"LOG(OR)_SE",	"Z_STAT",	"P")) %>%
  dplyr::select(CHR, POS, ID, P) %>%
  filter((CHR <= 22 | is.na(CHR))) %>%
  filter(CHR != 1.1 | is.na(CHR)) %>%
  filter(CHR != 21.1 | is.na(CHR)) %>%
  mutate(CHR = replace_na(CHR, 23))


# make manhattan plots for each phenotype & pull SNP IDs with genome-wide significance
## NB: blue (suggestive) line -log10(1E-5), red (genome-wide) line -log10(5E-8)

threshold <- 5E-8

manhattan(dist, main = "Dist", chr = "CHR", bp = "POS", p = "P", snp = "ID", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # no genome-wide sig (GWS)

manhattan(maxCranElev, main = "Max Cran Elev", chr = "CHR", bp = "POS", p = "P", snp = "ID", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # genome wide sig

mce_sig <- maxCranElev %>%
  filter(P < threshold)

manhattan(maxDecel, main = "Max Decel", chr = "CHR", bp = "POS", p = "P", snp = "ID", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # GWS

md_sig <- maxDecel %>%
  filter(P < threshold)

manhattan(maxGape, main = "Max Gape", chr = "CHR", bp = "POS", p = "P", snp = "ID", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # no GWS

manhattan(maxHD, main = "Max HD", chr = "CHR", bp = "POS", p = "P", snp = "ID", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # no GWS

manhattan(maxJP, main = "Max JP", chr = "CHR", bp = "POS", p = "P", snp = "ID", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # no GWS

manhattan(PPD_MG, main = "PPD_MG", chr = "CHR", bp = "POS", p = "P", snp = "ID", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # GWS

ppdmg_sig <- PPD_MG %>%
  filter(P < threshold)

manhattan(PPD_SI, main = "PPD SI", chr = "CHR", bp = "POS", p = "P", snp = "ID", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # no GWS

manhattan(ramSpeed, main = "Ram Speed", chr = "CHR", bp = "POS", p = "P", snp = "ID", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # GWS

rs_sig <- ramSpeed %>%
  filter(P < threshold)

manhattan(sf, main = "Success/Failure", chr = "CHR", bp = "POS", p = "P", snp = "ID", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # no GWS

manhattan(time_HDvMG, main = "Time - HD v MG", chr = "CHR", bp = "POS", p = "P", snp = "ID", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # GWS

thdvmg_sig <- time_HDvMG %>%
  filter(P < threshold)

manhattan(ttpg, main = "TTPG", chr = "CHR", bp = "POS", p = "P", snp = "ID", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) # GWS

ttpg_sig <- ttpg %>%
  filter(P < threshold)
