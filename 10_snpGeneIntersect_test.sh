# in projects/def-sjsmith/sjsmith/stickles_ucr/geneAssoc
conda activate vcfFilt

vcftools --vcf ../gasAcu.filter.chrRename.vcf --bed cleanGenes.bed --recode --recode-INFO-all --out gasAcu.snps
