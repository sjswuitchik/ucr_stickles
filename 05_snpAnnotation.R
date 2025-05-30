library(tidyverse)

file_paths_pvar <- list.files(file_dir, pattern="gasAcu.plink.*\\.pvar$", full.names = TRUE)

for (file in file_paths_pvar) {
  
  parts <- str_split(basename(file), "\\.", simplify = TRUE)
  name_pvar <- paste0(parts[, 3], ".pvar")
  
  data_pvar <- read.table(file, header = F, col.names = c("CHR", "POS", "ID", "REF", "ALT", "QUAL", "INFO")) %>%
    filter((CHR <= 22 | is.na(CHR))) %>%
    filter(CHR != 1.1 | is.na(CHR)) %>%
    filter(CHR != 21.1 | is.na(CHR)) %>%
    mutate(CHR = replace_na(CHR, 23)) %>%
    dplyr::select(-c(QUAL, INFO))
  
  assign(name_pvar, data_pvar, envir = .GlobalEnv)
}

# associate sig SNPs with .pvar files, extract the significant SNPs, and write out

mce_sig_ann <- mce_sig %>% 
  inner_join(maxCranElev.pvar, by = c("CHR", "POS", "ID")) %>%
  dplyr::select(ID) %>%
  write_delim(., "mce_sig_snps.txt", delim = '\t')

md_sig_ann <- md_sig %>% 
  inner_join(maxDecel.pvar, by = c("CHR", "POS", "ID")) %>%
  dplyr::select(ID) %>%
  write_delim(., "md_sig_snps.txt", delim = '\t')

ppdmg_sig_ann <- ppdmg_sig %>% 
  inner_join(PPD_MG.pvar, by = c("CHR", "POS", "ID")) %>%
  dplyr::select(ID) %>%
  write_delim(., "ppdmg_sig_snps.txt", delim = '\t')

rs_sig_ann <- rs_sig %>% 
  inner_join(ramSpeed.pvar, by = c("CHR", "POS", "ID")) %>%
  dplyr::select(ID) %>%
  write_delim(., "rs_sig_snps.txt", delim = '\t')

thdvmg_sig_ann <- thdvmg_sig %>% 
  inner_join(time_HDvMG.pvar, by = c("CHR", "POS", "ID")) %>%
  dplyr::select(ID) %>%
  write_delim(., "thdvmg_sig_snps.txt", delim = '\t')

ttpg_sig_ann <- ttpg_sig %>% 
  inner_join(ttpg.pvar, by = c("CHR", "POS", "ID")) %>%
  dplyr::select(ID) %>%
  write_delim(., "ttpg_sig_snps.txt", delim = '\t')

##### make .BED file for gene association
mce_bed <- mce_sig %>% 
  inner_join(maxCranElev.pvar, by = c("chr", "pos", "ID")) %>%
  select(-c(REF, ALT, p)) %>%
  rename(end = pos) %>%
  mutate(start = end - 1,
         pheno = "mce") %>%
  select(chr, start, end, ID, pheno)

md_bed <- md_sig %>% 
  inner_join(maxDecel.pvar, by = c("chr", "pos", "ID")) %>%
  select(-c(REF, ALT, p)) %>%
  rename(end = pos) %>%
  mutate(start = end - 1,
         pheno = "md") %>%
  select(chr, start, end, ID, pheno)

ppdmg_bed <- ppdmg_sig %>% 
  inner_join(PPD_MG.pvar, by = c("chr", "pos", "ID")) %>%
  select(-c(REF, ALT, p)) %>%
  rename(end = pos) %>%
  mutate(start = end - 1,
         pheno = "ppdmg") %>%
  select(chr, start, end, ID, pheno)

rs_bed <- rs_sig %>% 
  inner_join(ramSpeed.pvar, by = c("chr", "pos", "ID")) %>%
  select(-c(REF, ALT, p)) %>%
  rename(end = pos) %>%
  mutate(start = end - 1,
         pheno = "rs") %>%
  select(chr, start, end, ID, pheno)

thdvmg_bed <- thdvmg_sig %>% 
  inner_join(time_HDvMG.pvar, by = c("chr", "pos", "ID")) %>%
  select(-c(REF, ALT, p)) %>%
  rename(end = pos) %>%
  mutate(start = end - 1,
         pheno = "thdvmg") %>%
  select(chr, start, end, ID, pheno)

ttpg_bed <- ttpg_sig %>% 
  inner_join(ttpg.pvar, by = c("chr", "pos", "ID")) %>%
  select(-c(REF, ALT, p)) %>%
  rename(end = pos) %>%
  mutate(start = end - 1,
         pheno = "ttpg") %>%
  select(chr, start, end, ID, pheno)

clean_bed <- bind_rows(mce_bed, md_bed, ppdmg_bed, rs_bed, thdvmg_bed, ttpg_bed)
write_delim(clean_bed, "geneAssoc/cleanGenes.bed", delim = '\t', col_names = F)
