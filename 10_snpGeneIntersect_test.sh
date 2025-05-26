# in projects/def-sjsmith/sjsmith/stickles_ucr/geneAssoc
conda activate vcfFilt

bedtools intersect -a cleanGenes.bed -b ../gasAcu.filter.chrRename.vcf -wa
