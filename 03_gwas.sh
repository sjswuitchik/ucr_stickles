# create env for gemma gwas
## NB: ended up not using GEMMA for GWAS but the bcf/vcftools still worked in this conda env
conda create -n gwas -c esgf -c bioconda -c ostrokach -c conda-forge python=3.10 plink vcftools bcftools gemma libgfortran5 
conda activate gwas

# in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr

## chromosomes need to be converted from alpahnumeric to numeric values. 
# create chromosome map
bcftools view -H stickles.filtered.recode.vcf | cut -f 1 | uniq | awk '{print $0}' > gasAcu.chrom.map.txt

## edit in R (because I don't have the capacity to loop this rn) gasAcu.chrom.map.txt using sequence report from NCBI (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016920845.1/) to assign numerical chr names, then save as gasAcu.chromMap

bcftools annotate stickles.filtered.recode.vcf --rename-chrs gasAcu.chromMap -o gasAcu.filter.chrRename.vcf -O v

## remove individuals with no phenotypes (not F2)
vcftools --vcf gasAcu.filter.chrRename.vcf --remove-indv dedup/OBBB_1.dedup.bam --remove-indv dedup/OOB_1.dedup.bam --recode --recode-INFO-all --out gasAcu.chrRename.final
# kept 16466 out of a possible 16466 Sites, 58/60 individuals

##### use the gasAcu.chrRename.final.recode.vcf for the S vs F GWAS 

# there are two individuals who only had failed trials and therefore have to be removed 
vcftools --vcf gasAcu.chrRename.final.recode.vcf --remove-indv dedup/BAM_19.dedup.bam --remove-indv dedup/BAM_59.dedup.bam --recode --recode-INFO-all --out gasAcu.chrRename.noFails
# kept 56 out of 58 Individuals

##### use the gasAcu.chrRename.noFails.recode.vcf for the continuous variable GWAS

conda deactivate 
conda create -n plink2 -c bioconda plink2
conda activate plink2

#success-failure
plink2 --vcf gasAcu.chrRename.final.recode.vcf --make-pgen --allow-extra-chr --set-all-var-ids @:# --snps-only --hwe 0.05 --pheno sf.caseControlphenos.tsv --out gasAcu.plink.sf

# continuous traits 
while IFS= read -r file
do
  plink2 --vcf gasAcu.chrRename.noFails.recode.vcf --make-pgen --allow-extra-chr --set-all-var-ids @:# --snps-only --hwe 0.05 --pheno phenos.cont.plink2.tsv --pheno-name $file --out gasAcu.plink.$file
done < "cont.phenos"

## old stuff

plink2 --pfile gasAcu.plink.sf --allow-extra-chr --glm --adjust --out gasAcu.plink
mv gasAcu.plink.log gasAcu.plink.sf.log

while IFS= read -r file
do
  plink2 --pfile gasAcu.plink.$file --allow-extra-chr --glm --adjust --out gasAcu.plink
  mv gasAcu.plink.log gasAcu.plink.$file.log
done < "cont.phenos"

mkdir gwas_results

while IFS= read -r file
do
  mv gasAcu.plink.$file.* gwas_results/
done < "cont.phenos"

mv gasAcu.plink.sf* gwas_results/

# convert plink 2.0 to plink 1.9 for LD analyses
cp cont.phenos gwas_results/

while read file
do 
  plink2 --pfile gasAcu.plink.$file --make-bed --allow-extra-chr --out gasAcu.plink19.$file
done < cont.phenos

plink2 --pfile gasAcu.plink.sf --make-bed --allow-extra-chr --out gasAcu.plink19.sf




### random new things that aren't working/correct but I'm not ready to get rid of yet

plink2 --pfile gasAcu.plink.sf --allow-extra-chr --make-rel --out gasAcu.sf.grm
plink2 --pfile gasAcu.plink.sf --allow-extra-chr --make-grm-bin --out gasAcu.sf.grm.bin
plink2 --pfile gasAcu.plink.sf --allow-extra-chr --glm --grm-bin gasAcu.sf.grm.bin --adjust --out gasAcu.plink.grm



plink2 --pfile gasAcu.plink.sf --allow-extra-chr --make-king-table --out gasAcu.sf.king



gcta64 --pfile gasAcu.plink.sf --make-grm --autosome --out gasAcu.plink.sf

gcta64 --grm <pheno.grm> --make-bK-sparse 0.05 --out <pheno>.sp.grm

## alternatively, can do: 
gcta64 --pfile gasAcu.plink.sf --make-grm --sparse-cutoff 0.05 --out gasAcu.plink.sf

gcta64 --pfile gasAcu.plink.sf --fastGWA-mlm --grm-sparse gasAcu.plink.sf.grm.sp --nofilter --out gasAcu.plink.sf.gcta
--pheno ../sf.caseControlphenos.tsv



gcta64 --pfile <pheno> --make-grm --autosome --out <pheno.grm>

gcta64 --grm <pheno.grm> --make-bK-sparse 0.05 --out <pheno>.sp.grm

## alternatively, can do: 
gcta64 --pfile <pheno> --make-grm --sparse-cutoff 0.05 --out <pheno.sp.grm>

gcta64 --pfile <pheno> --pheno <pheno> --fastGWA-mlm --grm-sparse <pheno.sp.grm> --nofilter -- --out <newresults>

