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

## PANTHER interlude

list <- read_excel("geneAssoc/geneAssoc.xlsx")

panther <- read_delim("geneAssoc/pantherGeneList_full.txt", col_names = F) %>%
  rename(gene = "X2", sum = "X3", ref = "X6") %>%
  select(gene, sum, ref) %>%
  separate(sum, into = c("summary", NA), sep = ';', extra = "drop")

cleanList <- full_join(list, panther, by = "gene", relationship = "many-to-many") %>%
  mutate(ref = if_else(is.na(ref), "Homo sapiens", ref))
write_delim(cleanList, "geneAssoc/table_geneAssoc.tsv", delim = "\t")


##########################
## Morphological Traits ##
##########################
#### Create gene list for ray SNPs

genes <- read.table("gwas_out/gasAcu.rayGenes.clean.bed", sep = '\t', 
                    col.names = c("chr", "start", "end", "genes", "snp", "pheno")) %>%
  arrange(pheno, chr, start, end) %>%
  select(snp, pheno, genes) 

gene_split <- as_tibble(str_split(genes$genes, " ", n = Inf, simplify = T))

gene_clean <- genes %>%
  bind_cols(., gene_split) %>%
  select(-c(genes))

write_delim(gene_clean, "gwas_out/gasAcu.final.rayGenes.bed", delim = "\t")
write_delim(gene_clean, "gwas_out/gasAcu.final.rayGenes.txt", delim = "\t")
write_delim(gene_clean, "gwas_out/gasAcu.final.rayGenes.csv", delim = ",")

snpList <- gene_clean %>%
  select(-c(snp, pheno)) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene)

write_delim(snpList, "gwas_out/snpList", delim = "\t", col_names = F)

## PANTHER interlude

list <- read_excel("gwas_out/rayGeneAssoc.xlsx")

panther <- read_delim("gwas_out/pantherGeneList.txt", col_names = F) %>%
  rename(gene = "X2", sum = "X3", ref = "X6") %>%
  select(gene, sum, ref) %>%
  separate(sum, into = c("summary", NA), sep = ';', extra = "drop")

cleanList <- full_join(list, panther, by = "gene", relationship = "many-to-many") %>%
  mutate(ref = if_else(is.na(ref), "Homo sapiens", ref))
write_delim(cleanList, "gwas_out/table_rayGeneAssoc.tsv", delim = "\t")
