## just trying to bang this out so it's not going to be organized in any way at all or run as a standalone script until I redo the repo and actually organize things. You've been warned 

# in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr
conda activate gwas 

# remove individuals that don't have morph phenos
vcftools --vcf gasAcu.chrRename.final.recode.vcf --remove-indv dedup/BAM_18.dedup.bam --remove-indv dedup/BAM_59.dedup.bam --recode --recode-INFO-all --out gasAcu.chrRename.morph

conda deactivate 
conda activate plink2 

while IFS= read -r file
do
  plink2 --vcf gasAcu.chrRename.morph.recode.vcf --make-pgen --allow-extra-chr --set-all-var-ids @:# --snps-only --hwe 0.05 --pheno phenos.cont.plink2.tsv --pheno-name $file --out gasAcu.plink.$file
done < "cont.phenos"


