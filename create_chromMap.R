library(tidyverse)

seq.rep <- read_delim("~/Desktop/MRU_Faculty/Research/ucr_stickles/sequence_report.tsv") %>%
  select("RefSeq seq accession", "Sequence name")
chrom <- read_delim("~/Desktop/MRU_Faculty/Research/ucr_stickles/gasAcu.chrom.map.txt", col_names = F) %>%
  select(X1)

seq.rep.chroms <- seq.rep[1:22,]
seq.rep.scaffs <- seq.rep[23:250,]

seq.rep.chroms.clean <- seq.rep.chroms %>%
  separate(col = "Sequence name", into = c("NA", "char"), sep = "chr") %>%
  select(-c("NA")) %>%
  mutate(chrom = c(1:22)) %>%
  select(-c("char"))

seq.rep.scaffs.clean <- seq.rep.scaffs %>%
  separate(col = "Sequence name", into = c("NA", "char"), sep = "Contig") %>%
  select(-c("NA")) %>%
  mutate(ex = "1") %>%
  unite("chrom", char:ex, sep = '.') %>%
  mutate(chrom = as.numeric(chrom))

clean.map <- bind_rows(seq.rep.chroms.clean, seq.rep.scaffs.clean)

write_delim(clean.map, "~/Desktop/MRU_Faculty/Research/ucr_stickles/gasAcu.chromMap")
