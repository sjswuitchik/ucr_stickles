#### Linkage Disequilibrium Analysis in SNPs with genome-wide significance

# Read `.pvar` files in and save as new objects
files.pvar <- c("maxCranElev", "maxDecel", "time_HDvMG", "time_maxDecelvMG", "ttpg")

for (file in files.pvar) {
  df <- read.table(paste0("gasAcu.plink.", file, ".pvar"), 
                   header = F, 
                   col.names = c("CHR", "POS", "ID", "REF", "ALT", "QUAL", "INFO")) %>%
    filter((CHR <= 22 | is.na(CHR))) %>%
    filter(CHR != 1.1 | is.na(CHR)) %>%
    filter(CHR != 21.1 | is.na(CHR)) %>%
    mutate(CHR = replace_na(CHR, 23)) %>%
    dplyr::select(-c(QUAL, INFO)) %>%
    rename(chrom = CHR, pos = POS, variant.id = ID,)
  
  assign(paste0(file, ".pvar"), df, envir = .GlobalEnv)
}


# Associate significant SNPs with `.pvar` files, extract the significant SNPs, and write out

mce.sig.ann <- mce.sig %>% 
  inner_join(maxCranElev.pvar, by = c("chrom", "pos", "variant.id")) %>%
  dplyr::select(variant.id) %>%
  write_delim(., "mce_sig_snps.txt", delim = '\t')

md.sig.ann <- md.sig %>% 
  inner_join(maxDecel.pvar, by = c("chrom", "pos", "variant.id")) %>%
  dplyr::select(variant.id) %>%
  write_delim(., "md_sig_snps.txt", delim = '\t')

time.hdvmg.sig.ann <- time.hdvmg.sugg.sig %>% 
  inner_join(time_HDvMG.pvar, by = c("chrom", "pos", "variant.id")) %>%
  dplyr::select(variant.id) %>%
  write_delim(., "time_hdvmg_sugg_sig_snps.txt", delim = '\t')

time.mdvmg.sig.ann <- time.mdvmg.sugg.sig %>% 
  inner_join(time_HDvMG.pvar, by = c("chrom", "pos", "variant.id")) %>%
  dplyr::select(variant.id) %>%
  write_delim(., "time_mdvmg_sugg_sig_snps.txt", delim = '\t')

ttpg.sig.ann <- ttpg.sig %>% 
  inner_join(ttpg.pvar, by = c("chrom", "pos", "variant.id")) %>%
  dplyr::select(variant.id) %>%
  write_delim(., "ttpg_sig_snps.txt", delim = '\t')
