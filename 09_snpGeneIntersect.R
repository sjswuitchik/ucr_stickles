library(tidyverse)

mce <- read.table("geneAssoc/mce_windows.bed", header = F) %>%
  select(-c(V5, V6, V7, V9, V10, V11, V13)) %>%
  rename(chr = V1, start = V2, end = V3, snp = V4, id = V8, feat = V12, notes = V14) %>%
  filter(feat == "CDS" | feat == "gene") %>%
  separate(col = notes, into = c(NA, "int"), sep = "gene=") %>%
  separate(col = int, into = c("gene", NA), sep = ";", extra = "drop")

mce_geneWindows <- mce %>%
  group_by(snp, start, end) %>%
  distinct(gene)

mce_int <- as_tibble(unique(mce_geneWindows$gene)) %>%
  rename(gene = value)

mce_locSub <- mce_int[grep('LOC*', mce_int$gene),]

mce_geneList <- anti_join(mce_int, mce_locSub)

mce_cleanGenes <- inner_join(mce_geneWindows, mce_geneList, by = "gene") %>%
  write_delim(., "geneAssoc/mce_cleanGenes.tsv", delim = '\t')

md <- read.table("geneAssoc/md_windows.bed", header = F) %>%
  select(-c(V5, V6, V7, V9, V10, V11, V13)) %>%
  rename(chr = V1, start = V2, end = V3, snp = V4, id = V8, feat = V12, notes = V14) %>%
  filter(feat == "CDS" | feat == "gene") %>%
  separate(col = notes, into = c(NA, "int"), sep = "gene=") %>%
  separate(col = int, into = c("gene", NA), sep = ";", extra = "drop")

md_geneWindows <- md %>%
  group_by(snp, start, end) %>%
  distinct(gene)

md_int <- as_tibble(unique(md_geneWindows$gene)) %>%
  rename(gene = value)

md_locSub <- md_int[grep('LOC*', md_int$gene),]

md_geneList <- anti_join(md_int, md_locSub)

md_cleanGenes <- inner_join(md_geneWindows, md_geneList, by = "gene") %>%
  write_delim(., "geneAssoc/md_cleanGenes.tsv", delim = '\t')

ppdmg <- read.table("geneAssoc/ppdmg_windows.bed", header = F) %>%
  select(-c(V5, V6, V7, V9, V10, V11, V13)) %>%
  rename(chr = V1, start = V2, end = V3, snp = V4, id = V8, feat = V12, notes = V14) %>%
  filter(feat == "CDS" | feat == "gene") %>%
  separate(col = notes, into = c(NA, "int"), sep = "gene=") %>%
  separate(col = int, into = c("gene", NA), sep = ";", extra = "drop")

ppdmg_geneWindows <- ppdmg %>%
  group_by(snp, start, end) %>%
  distinct(gene)

ppdmg_int <- as_tibble(unique(ppdmg_geneWindows$gene)) %>%
  rename(gene = value)

ppdmg_locSub <- ppdmg_int[grep('LOC*', ppdmg_int$gene),]

ppdmg_geneList <- anti_join(ppdmg_int, ppdmg_locSub)

ppdmg_cleanGenes <- inner_join(ppdmg_geneWindows, ppdmg_geneList, by = "gene") %>%
  write_delim(., "geneAssoc/ppdmg_cleanGenes.tsv", delim = '\t')

rs <- read.table("geneAssoc/rs_windows.bed", header = F) %>%
  select(-c(V5, V6, V7, V9, V10, V11, V13)) %>%
  rename(chr = V1, start = V2, end = V3, snp = V4, id = V8, feat = V12, notes = V14) %>%
  filter(feat == "CDS" | feat == "gene") %>%
  separate(col = notes, into = c(NA, "int"), sep = "gene=") %>%
  separate(col = int, into = c("gene", NA), sep = ";", extra = "drop")

rs_geneWindows <- rs %>%
  group_by(snp, start, end) %>%
  distinct(gene)

rs_int <- as_tibble(unique(rs_geneWindows$gene)) %>%
  rename(gene = value)

rs_locSub <- rs_int[grep('LOC*', rs_int$gene),]

rs_geneList <- anti_join(rs_int, rs_locSub)

rs_cleanGenes <- inner_join(rs_geneWindows, rs_geneList, by = "gene") %>%
  write_delim(., "geneAssoc/rs_cleanGenes.tsv", delim = '\t')

thdvmg <- read.table("geneAssoc/thdvmg_windows.bed", header = F) %>%
  select(-c(V5, V6, V7, V9, V10, V11, V13)) %>%
  rename(chr = V1, start = V2, end = V3, snp = V4, id = V8, feat = V12, notes = V14) %>%
  filter(feat == "CDS" | feat == "gene") %>%
  separate(col = notes, into = c(NA, "int"), sep = "gene=") %>%
  separate(col = int, into = c("gene", NA), sep = ";", extra = "drop")

thdvmg_geneWindows <- thdvmg %>%
  group_by(snp, start, end) %>%
  distinct(gene)

thdvmg_int <- as_tibble(unique(thdvmg_geneWindows$gene)) %>%
  rename(gene = value)

thdvmg_locSub <- thdvmg_int[grep('LOC*', thdvmg_int$gene),]

thdvmg_geneList <- anti_join(thdvmg_int, thdvmg_locSub)

thdvmg_cleanGenes <- inner_join(thdvmg_geneWindows, thdvmg_geneList, by = "gene") %>%
  write_delim(., "geneAssoc/thdvmg_cleanGenes.tsv", delim = '\t')

ttpg <- read.table("geneAssoc/ttpg_windows.bed", header = F) %>%
  select(-c(V5, V6, V7, V9, V10, V11, V13)) %>%
  rename(chr = V1, start = V2, end = V3, snp = V4, id = V8, feat = V12, notes = V14) %>%
  filter(feat == "CDS" | feat == "gene") %>%
  separate(col = notes, into = c(NA, "int"), sep = "gene=") %>%
  separate(col = int, into = c("gene", NA), sep = ";", extra = "drop")

ttpg_geneWindows <- ttpg %>%
  group_by(snp, start, end) %>%
  distinct(gene)

ttpg_int <- as_tibble(unique(ttpg_geneWindows$gene)) %>%
  rename(gene = value)

ttpg_locSub <- ttpg_int[grep('LOC*', ttpg_int$gene),]

ttpg_geneList <- anti_join(ttpg_int, ttpg_locSub)

ttpg_cleanGenes <- inner_join(ttpg_geneWindows, ttpg_geneList, by = "gene") %>%
  write_delim(., "geneAssoc/ttpg_cleanGenes.tsv", delim = '\t')
