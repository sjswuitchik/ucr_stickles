#library(BiocManager)
#BiocManager::install("GENESIS", force = T)
library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(tidyverse)
library(qqman)
library(topr)
library(reshape2)
library(viridis)
library(readxl)

## loading phenotypes
sf.pheno <- read_delim("~/Desktop/MRU_Faculty/Research/stickles_ucr/gwas_results/phenotypes_gwas/phenos_successFail_plink2_caseControl.tsv") %>%
  rename(scanID = IID) %>%
  mutate(sf = if_else(sf == 1, 0, sf),
         sf = if_else(sf == 2, 1, sf))

ppdmg.pheno <- read_delim("~/Desktop/MRU_Faculty/Research/stickles_ucr/gwas_results/phenotypes_gwas/phenos_ppdmg.tsv") %>%
  rename(scanID = IID)

cont.phenos <- read_delim("~/Desktop/MRU_Faculty/Research/stickles_ucr/gwas_results/phenotypes_gwas/phenos_cont_plink2.tsv") %>%
  rename(scanID = IID)

## create scan annotation dfs
scanAnnot.sf <- ScanAnnotationDataFrame(sf.pheno)
scanAnnot.cont <- ScanAnnotationDataFrame(cont.phenos)
scanAnnot.ppdmg <- ScanAnnotationDataFrame(ppdmg.pheno)

# attempting to loop continuous GWAS
files <- c("dist", "maxCranElev", "maxDecel", "maxGape", "maxHD", "maxJP", "ppdmg", "PPD_SI", "ramSpeed", "time_HDvMG", "time_maxDecelvMG", "ttpg")

# make GDS function 
make_gds <- function(pheno) { 
  snpgdsBED2GDS(bed.fn = paste0("gasAcu.plink19.", pheno, ".bed"),
                bim.fn = paste0("gasAcu.plink19.", pheno, ".bim"),
                fam.fn = paste0("gasAcu.plink19.", pheno, ".fam"),
                out.gdsfn = paste0(pheno, ".gds"),
                cvt.chr = "char")
}

for (pheno in files) {
  make_gds(pheno)
}

# create KING matrices & geno data
create_geno <- function(pheno) {
  geno <- GdsGenotypeReader(filename = paste0(pheno, ".gds"))
  genoData <- GenotypeData(geno)
  assign(paste0("genoData.", pheno), genoData, envir = .GlobalEnv)
}

for (pheno in files) {
  create_geno(pheno)
}

create_kin <- function(pheno) {
  gds <- snpgdsOpen(paste0(pheno, ".gds"), readonly = F, allow.duplicate = T)
  kin <- snpgdsIBDKING(gds)
  kin.mat <- kingToMatrix(kin)
  assign(paste0("kin.mat.", pheno), kin.mat, envir = .GlobalEnv)
  snpgdsClose(gds)
}

for (pheno in files) {
  create_kin(pheno)
}

# scanAnnot.cont names i.e., no s/f, no ppdmg
cont.files <- c("dist", "maxCranElev", "maxDecel", "maxGape", "maxHD", "maxJP", "PPD_SI", "ramSpeed", "time_HDvMG", "time_maxDecelvMG", "ttpg")

# create null models 
for (pheno in cont.files) {
  kin_name <- paste0("kin.mat.", pheno)
  null_mod <- fitNullModel(
    scanAnnot.cont,
    outcome = pheno,
    cov.mat = get(kin_name),
    family = "gaussian"
  )

  assign(paste0("null.mod.", pheno), null_mod, envir = .GlobalEnv)
}

# run GWAS analyses for continous phenotypes, clean up results, write out to env
for (pheno in cont.files) {
  geno_data <- get(paste0("genoData.", pheno))
  genoIterator <- GenotypeBlockIterator(geno_data, snpBlock = 10000)
  assoc <- assocTestSingle(genoIterator, null.model = get(paste0("null.mod.", pheno)), BPPARAM = BiocParallel::SerialParam())
  
  assign(paste0("assoc.", pheno), assoc, envir = .GlobalEnv)
  
  assoc_clean <- assoc %>%
    mutate(chr = if_else(chr == "U", "23", chr),
           chr = as.numeric(chr)) %>%
    rename(chrom = chr, pos = pos, p = Score.pval)
  
  assign(paste0("assoc.clean.", pheno), assoc_clean, envir = .GlobalEnv)
  
  
}

## ppdmg GWAS 
snpgdsBED2GDS(bed.fn = "gasAcu.plink19.ppdmg.bed",
              bim.fn = "gasAcu.plink19.ppdmg.bim",
              fam.fn = "gasAcu.plink19.ppdmg.fam",
              out.gdsfn = "gasAcu.plink19.ppdmg.gds",
              cvt.chr="char")

geno.ppdmg <- GdsGenotypeReader(filename = "gasAcu.plink19.ppdmg.gds")
genoData.ppdmg <- GenotypeData(geno.ppdmg)

gds.ppdmg <- snpgdsOpen("gasAcu.plink19.ppdmg.gds", readonly = F, allow.duplicate = T)

kin.ppdmg <- snpgdsIBDKING(gds.ppdmg)

kin.mat.ppdmg <- kingToMatrix(kin.ppdmg)

snpgdsClose(gds.ppdmg)

null.mod.ppdmg <- fitNullModel(scanAnnot.ppdmg, outcome = "PPD_MG", cov.mat = kin.mat.ppdmg, family = "gaussian")

genoIterator.ppdmg <- GenotypeBlockIterator(genoData.ppdmg, snpBlock=10000)

assoc.ppdmg <- assocTestSingle(genoIterator.ppdmg, null.model = null.mod.ppdmg,
                            BPPARAM = BiocParallel::SerialParam())

assoc.ppdmg.clean <- assoc.ppdmg %>%
  mutate(chr = if_else(chr == "U", "23", chr),
         chr = as.numeric(chr)) %>%
  rename(chrom = chr, pos = pos, p = Score.pval)


## binary s/f GWAS
snpgdsBED2GDS(bed.fn = "gasAcu.plink19.sf.bed",
              bim.fn = "gasAcu.plink19.sf.bim",
              fam.fn = "gasAcu.plink19.sf.fam",
              out.gdsfn = "gasAcu.plink19.sf.gds",
              cvt.chr="char")

geno.sf <- GdsGenotypeReader(filename = "gasAcu.plink19.sf.gds")
genoData.sf <- GenotypeData(geno.sf)

gds.sf <- snpgdsOpen("gasAcu.plink19.sf.gds", readonly = F, allow.duplicate = T)

kin.sf <- snpgdsIBDKING(gds.sf)

kin.mat.sf <- kingToMatrix(kin.sf)

snpgdsClose(gds.sf)

null.mod.sf <- fitNullModel(scanAnnot.sf, outcome = "sf", cov.mat = kin.mat.sf, family = "binomial")

genoIterator.sf <- GenotypeBlockIterator(genoData.sf, snpBlock=10000)

assoc.sf <- assocTestSingle(genoIterator.sf, null.model = null.mod.sf,
                            BPPARAM = BiocParallel::SerialParam())

assoc.sf.clean <- assoc.sf %>%
  mutate(chr = if_else(chr == "U", "23", chr),
         chr = as.numeric(chr)) %>%
  rename(chrom = chr, pos = pos, p = Score.pval)


## Manhattan plots 

# s/f - no GWS
qqman::manhattan(assoc.sf.clean, main = "S/F", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.sf.clean, title = "S/F", annotate = 1e-5, ymin = 0, ymax = 10)

# ppdmg - no GWS
qqman::manhattan(assoc.ppdmg.clean, main = "PPD_MG", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.ppdmg.clean, title = "PPD_MG", ymin = 0, ymax = 10)

# dist - no GWS
qqman::manhattan(assoc.clean.dist, main = "Dist", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.dist, title = "Dist", annotate = 1e-5, ymin = 0, ymax = 10)

# maxCranElev - GWS!
qqman::manhattan(assoc.clean.maxCranElev, main = "Max Cran Elev", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.maxCranElev, title = "Max Cran Elev", ymin = 0, ymax = 10)

# maxDecel - GWS!
qqman::manhattan(assoc.clean.maxDecel, main = "Max Decel", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.maxDecel, title = "Max Decel", ymin = 0, ymax = 10)

# maxGape - no GWS
qqman::manhattan(assoc.clean.maxGape, main = "Max Gape", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.maxGape, title = "Max Gape", ymin = 0, ymax = 10)

# maxHD - no GWS
qqman::manhattan(assoc.clean.maxHD, main = "Max HD", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.maxHD, title = "Max HD", ymin = 0, ymax = 10)

# maxJP - no GWS
qqman::manhattan(assoc.clean.maxJP, main = "Max JP", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.maxJP, title = "Max JP", ymin = 0, ymax = 10)

# PPD_SI - no GWS
qqman::manhattan(assoc.clean.PPD_SI, main = "PPD_SI", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.PPD_SI, title = "PPD_SI", ymin = 0, ymax = 10)

# ramSpeed - no GWS
qqman::manhattan(assoc.clean.ramSpeed, main = "Ram Speed", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.ramSpeed, title = "Ram Speed", ymin = 0, ymax = 10)

# time_HDvMG - very close to GWS, worth investigating
qqman::manhattan(assoc.clean.time_HDvMG, main = "Time - HD v MG", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.time_HDvMG, title = "Time - HD v MG", ymin = 0, ymax = 10)

# time_maxDecelvMG - very close to GWS, worth investigating
qqman::manhattan(assoc.clean.time_maxDecelvMG, main = "Time - Max Decel v MG", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.time_maxDecelvMG, title = "Time - Max Decel v MG", ymin = 0, ymax = 10)

# ttpg - GWS!
qqman::manhattan(assoc.clean.ttpg, main = "TTPG", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.ttpg, title = "TTPG", ymin = 0, ymax = 10)

# sig or close to: maxCranElev, maxDecel, time_HDvMG, time_maxDecelvMG, ttpg

manhattan(list(assoc.clean.maxCranElev, assoc.clean.maxDecel, assoc.clean.time_HDvMG, assoc.clean.time_maxDecelvMG, assoc.clean.ttpg), ntop = 3, size = 1.5, legend_labels = c("Max Cran Elev", "Max Decel", "Time HD vMG", "Time MaxDecelvMG", "TTPG"), ymin = -10, ymax = 10)


# sig SNPs, or suggestive sig for the time variables 
threshold <- 5E-8

mce.sig <- assoc.clean.maxCranElev %>%
  filter(p < threshold) 

write_delim("gwas_grm_plots/mce.sig.csv", mce.sig, delim = ',')

md.sig <- assoc.clean.maxDecel %>%
  filter(p < threshold)

write_delim("gwas_grm_plots/md.sig.csv", md.sig, delim = ',')

time.hdvmg.sugg.sig <- assoc.clean.time_HDvMG %>%
  filter(p < 1E-5)

time.mdvmg.sugg.sig <- assoc.clean.time_maxDecelvMG %>%
  filter(p < 1E-5)

ttpg.sig <- assoc.clean.ttpg %>%
  filter(p < threshold)


##########################
## Morphological Traits ##
##########################

morph.phenos <- read_delim("~/Desktop/MRU_Faculty/Research/stickles_ucr/gwas_results/morphology/phenos/morph_phenos.tsv") %>%
  rename(scanID = IID)

scanAnnot.morph <- ScanAnnotationDataFrame(morph.phenos)

morph <- c("length", "height", "eye.dia", "caudal.area", "pec.length", "pec.area", "ray")

# make GDS function 
make_gds <- function(pheno) { 
  snpgdsBED2GDS(bed.fn = paste0("plink_out/gasAcu.plink19.", pheno, ".bed"),
                bim.fn = paste0("plink_out/gasAcu.plink19.", pheno, ".bim"),
                fam.fn = paste0("plink_out/gasAcu.plink19.", pheno, ".fam"),
                out.gdsfn = paste0("plink_out/", pheno, ".gds"),
                cvt.chr = "char")
}

for (pheno in morph) {
  make_gds(pheno)
}

# create KING matrices & geno data
create_geno <- function(pheno) {
  geno <- GdsGenotypeReader(filename = paste0("plink_out/", pheno, ".gds"))
  genoData <- GenotypeData(geno)
  assign(paste0("genoData.", pheno), genoData, envir = .GlobalEnv)
}

for (pheno in morph) {
  create_geno(pheno)
}

# create regularization function becuase otherwise nothing works (need positive definite matrix)
regularize_matrix <- function(K, lambda = 0.01) {
  (1 - lambda) * K + lambda * diag(nrow(K))
}

# create KING matrix
create_kin <- function(pheno) {
  gds <- snpgdsOpen(paste0("plink_out/", pheno, ".gds"), readonly = F, allow.duplicate = T)
  kin <- snpgdsIBDKING(gds)
  kin.mat <- kingToMatrix(kin)
  assign(paste0("kin.mat.", pheno), kin.mat, envir = .GlobalEnv)
  snpgdsClose(gds)
}

for (pheno in morph) {
  create_kin(pheno)
  kin.mat <- get(paste0("kin.mat.", pheno))
  kin.mat.reg <- regularize_matrix(kin.mat, lambda = 0.05)  # adjust lambda as needed
  assign(paste0("kin.mat.", pheno), kin.mat.reg, envir = .GlobalEnv)
}


# create null models 
for (pheno in morph) {
  kin_name <- paste0("kin.mat.", pheno)
  null_mod <- fitNullModel(
    scanAnnot.morph,
    outcome = pheno,
    cov.mat = get(kin_name),
    family = "gaussian"
  )
  
  assign(paste0("null.mod.", pheno), null_mod, envir = .GlobalEnv)
}

# run GWAS 
for (pheno in morph) {
  geno_data <- get(paste0("genoData.", pheno))
  genoIterator <- GenotypeBlockIterator(geno_data, snpBlock = 10000)
  assoc <- assocTestSingle(genoIterator, null.model = get(paste0("null.mod.", pheno)), BPPARAM = BiocParallel::SerialParam())
  
  assign(paste0("assoc.", pheno), assoc, envir = .GlobalEnv)
  
  assoc_clean <- assoc %>%
    mutate(chr = if_else(chr == "U", "23", chr),
           chr = as.numeric(chr)) %>%
    rename(chrom = chr, pos = pos, p = Score.pval)
  
  assign(paste0("assoc.clean.", pheno), assoc_clean, envir = .GlobalEnv)
  
  
}

# Manhattans
# length
qqman::manhattan(assoc.clean.length, main = "Standard Length", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.length, title = "Std Length", annotate = 1e-5, ymin = 0, ymax = 10)

# height 
qqman::manhattan(assoc.clean.height, main = "Height", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.height, title = "Height", annotate = 1e-5, ymin = 0, ymax = 10)

# eye diameter
qqman::manhattan(assoc.clean.eye.dia, main = "Eye Diameter", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.eye.dia, title = "Eye Diameter", annotate = 1e-5, ymin = 0, ymax = 10)

# caudal fin area
qqman::manhattan(assoc.clean.caudal.area, main = "Caudal Fin Area", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.caudal.area, title = "Caudal Fin Area", annotate = 1e-5, ymin = 0, ymax = 10)

# pec fin length
qqman::manhattan(assoc.clean.pec.length, main = "Pec Fin Length", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.pec.length, title = "Pec Fin Length", annotate = 1e-5, ymin = 0, ymax = 10)

# pec fin area
qqman::manhattan(assoc.clean.pec.area, main = "Pec Fin Area", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.pec.area, title = "Pec Fin Area", annotate = 1e-5, ymin = 0, ymax = 10)

# ray count - only pheno with GWS
qqman::manhattan(assoc.clean.ray, main = "Ray Count", chr = "chrom", bp = "pos", p = "p", snp = "variant.id", ylim = c(0, 10), chrlabs = c(1:21, "Y", "MT"), col = c("skyblue", "grey")) 

manhattan(assoc.clean.ray, title = "Ray Count", ymin = 0, ymax = 10)

# only dealing with ray count from here on out
threshold <- 5E-8

ray.sig <- assoc.clean.ray %>%
  filter(p < threshold) 

write.table(ray.sig, "gwas_out/ray.sig.csv", row.names = F, sep = ',')
