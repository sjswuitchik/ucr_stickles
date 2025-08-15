## just trying to bang this out so it's not going to be organized in any way at all or run as a standalone script until I redo the repo and actually organize things. You've been warned 

# in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr
conda activate vcfFilt

# remove individuals that don't have morph phenos
vcftools --vcf gasAcu.chrRename.final.recode.vcf --remove-indv dedup/BAM_18.dedup.bam --remove-indv dedup/BAM_59.dedup.bam --recode --recode-INFO-all --out gasAcu.chrRename.morph
# kept 56/58 invd, 16466/16466 sites

conda deactivate 
conda activate plink2 

while read file
do
  plink2 --vcf gasAcu.chrRename.morph.recode.vcf --make-pgen --allow-extra-chr --set-all-var-ids @:# --snps-only --hwe 0.05 --pheno morph_phenos.tsv --pheno-name $file --out gasAcu.plink.$file
done < morph.phenos

while read file
do 
  plink2 --pfile gasAcu.plink.$file --make-bed --allow-extra-chr --out gasAcu.plink19.$file
done < morph.phenos

mkdir -p gwas_results/morph
mv gasAcu.plink* gwas_results/morph


## transfer files to comp to use GENESIS in R 

## back for LD matrix
# in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr/gwas_results/morph
plink --bfile gasAcu.plink19.ray --extract ray_sig_snps.txt --r2 square --allow-extra-chr --out ray_ld_matrix

## transfer to comp, continue in R

## back for some BED stuff
# in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr/gwas_results/morph
conda deactivate 
conda activate bedtools 

cp ../../geneAssoc/gasAcu.windows.repl.bed .

# intersect clean genes BED with windows
bedtools sort -i rayGenes.bed | bedtools intersect -a gasAcu.windows.repl.bed -b - -wb | bedtools sort -i - | bedtools merge -i - -d 1 -c 4,8,9 -o distinct > gasAcu.rayGenes.intersect.bed

# clean up gene names field by removing unnamed LOCs
sed 's/LOC[0-9]\+,//g' gasAcu.rayGenes.intersect.bed | sed 's/','/ /g' > gasAcu.rayGenes.clean.bed 
