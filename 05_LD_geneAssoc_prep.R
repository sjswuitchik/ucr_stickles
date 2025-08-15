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

#### Make .BED for gene association
mce_bed <- mce.sig %>% 
  inner_join(maxCranElev.pvar, by = c("chrom", "pos", "variant.id")) %>%
  select(-c(REF, ALT, p)) %>%
  rename(end = pos) %>%
  mutate(start = end - 1,
         pheno = "mce") %>%
  select(chrom, start, end, variant.id, pheno)

md_bed <- md.sig %>% 
  inner_join(maxDecel.pvar, by = c("chrom", "pos", "variant.id")) %>%
  select(-c(REF, ALT, p)) %>%
  rename(end = pos) %>%
  mutate(start = end - 1,
         pheno = "md") %>%
  select(chrom, start, end, variant.id, pheno)

thdvmg_bed <- time.hdvmg.sig %>% 
  inner_join(time_HDvMG.pvar, by = c("chrom", "pos", "variant.id")) %>%
  select(-c(REF, ALT, p)) %>%
  rename(end = pos) %>%
  mutate(start = end - 1,
         pheno = "thdvmg") %>%
  select(chrom, start, end, variant.id, pheno)

tmdvmg_bed <- time.mdvmg.sugg.sig %>% 
  inner_join(time_HDvMG.pvar, by = c("chrom", "pos", "variant.id")) %>%
  select(-c(REF, ALT, p)) %>%
  rename(end = pos) %>%
  mutate(start = end - 1,
         pheno = "tmdvmg") %>%
  select(chrom, start, end, variant.id, pheno)

ttpg_bed <- ttpg.sig %>% 
  inner_join(ttpg.pvar, by = c("chrom", "pos", "variant.id")) %>%
  select(-c(REF, ALT, p)) %>%
  rename(end = pos) %>%
  mutate(start = end - 1,
         pheno = "ttpg") %>%
  select(chrom, start, end, variant.id, pheno)

clean_bed <- bind_rows(mce_bed, md_bed, thdvmg_bed, tmdvmg_bed, ttpg_bed)
write_delim(clean_bed, "geneAssoc/cleanGenes.bed", delim = '\t', col_names = F)

##########################
## Morphological Traits ##
##########################

ray.pvar <- read.table("plink_out/gasAcu.plink.ray.pvar", 
                 header = F, 
                 col.names = c("CHR", "POS", "ID", "REF", "ALT", "QUAL", "INFO")) %>%
  filter((CHR <= 22 | is.na(CHR))) %>%
  filter(CHR != 1.1 | is.na(CHR)) %>%
  filter(CHR != 21.1 | is.na(CHR)) %>%
  mutate(CHR = replace_na(CHR, 23)) %>%
  dplyr::select(-c(QUAL, INFO)) %>%
  rename(chrom = CHR, pos = POS, variant.id = ID,)

ray.sig.ann <- ray.sig %>% 
  inner_join(ray.pvar, by = c("chrom", "pos", "variant.id")) %>%
  dplyr::select(variant.id) %>%
  write_delim(., "gwas_out/ray_sig_snps.txt", delim = '\t')

# LD matrix
ray.ld.mat <- read.table("gwas_out/ray_ld_matrix.ld", header = F)
rownames(ray.ld.mat) <- ray.sig.ann$variant.id
colnames(ray.ld.mat) <- ray.sig.ann$variant.id

ray.ld.long <- melt(as.matrix(ray.ld.mat), varnames = c("SNP_A", "SNP_B"), value.name = "R2")

ggplot(ray.ld.long, aes(x = SNP_A, y = SNP_B, fill = R2)) +
  geom_tile() +
  scale_fill_viridis(option = "viridis") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank()
  ) +
  labs(title = "LD Heatmap (r²) Ray Count", x = "SNP", y = "SNP", fill = expression(r^2))

ray_bed <- ray.sig %>% 
  inner_join(ray.pvar, by = c("chrom", "pos", "variant.id")) %>%
  select(-c(REF, ALT, p)) %>%
  rename(end = pos) %>%
  mutate(start = end - 1,
         pheno = "ray") %>%
  select(chrom, start, end, variant.id, pheno)

write_delim(ray_bed, "gwas_out/rayGenes.bed", delim = '\t', col_names = F)

