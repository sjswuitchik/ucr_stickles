# in projects/def-sjsmith/sjsmith/stickles_ucr/geneAssoc
conda activate vcfFilt

vcftools --vcf ../gasAcu.filter.chrRename.vcf --bed cleanGenes.bed --recode --recode-INFO-all --out gasAcu.snps

cp ../gasAcu.chromMap .

cat gasAcu.chromMap | awk '{print $2 "\t" $1}' > gasAcu.revChrMap

bcftools annotate gasAcu.snps.recode.vcf --rename-chrs gasAcu.revChrMap -o gasAcu.snps.chrRename.vcf -O v


# took this output to use in a trial of OmicsBox v.3.4.6 to run Genetic Variation > Variant Annotation 
# using GO version March 16, 2025, with reference genome and annotation from GAculeatus_UGA_version5 (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016920845.1/) 
