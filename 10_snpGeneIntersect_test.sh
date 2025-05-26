# in projects/def-sjsmith/sjsmith/stickles_ucr/geneAssoc
conda activate vcfFilt

vcftools --vcf ../gasAcu.filter.chrRename.vcf --bed cleanGenes.bed --recode --recode-INFO-all --out gasAcu.snps

# took this output to use in a trial of OmicsBox to run Genetic Variation > Variant Annotation 
# using GO version March 16, 2025, with reference genome and annotation from GAculeatus_UGA_version5 (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016920845.1/) 
