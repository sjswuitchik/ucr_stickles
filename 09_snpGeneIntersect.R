#### Create gene lists for associated SNPs

genes <- read.table("geneAssoc/gasAcu.genes.clean.bed", sep = '\t', 
                    col.names = c("chr", "start", "end", "genes", "snp", "pheno")) %>%
  arrange(pheno, chr, start, end) %>%
  select(snp, pheno, genes) 

gene_split <- as_tibble(str_split(genes$genes, " ", n = Inf, simplify = T))

gene_clean <- genes %>%
  bind_cols(., gene_split) %>%
  filter(snp != "12:16000501") %>% # gets rid of one leftover LOC*
  select(-c(genes))

write_delim(gene_clean, "geneAssoc/gasAcu.final.cleanGenes.bed", delim = "\t")
write_delim(gene_clean, "geneAssoc/gasAcu.final.cleanGenes.txt", delim = "\t")
write_delim(gene_clean, "geneAssoc/gasAcu.final.cleanGenes.csv", delim = ",")

snpList <- gene_clean %>%
  select(-c(snp, pheno)) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene)

write_delim(snpList, "geneAssoc/snpList", delim = "\t", col_names = F)
