# create env for gemma gwas
conda create -n gwas -c esgf -c bioconda -c ostrokach -c conda-forge python=3.12 plink vcftools bcftools gemma libgfortran5 
conda activate gwas
conda install -c conda-forge libgfortran5 ## need to update to py3 before installing libgfortran5 


# get gemma
git clone https://github.com/johanzi/gwas_gemma.git

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


cp gasAcu.chrRename.final.recode.vcf gasAcu.chrRename.noFails.recode.vcf gwas_gemma/
cd gwas_gemma

bcftools query -l gasAcu.chrRename.final.recode.vcf > order_accession.txt



## use run_gwas_gemma.sh from required_files on this repo to preserve edits that are required to run GWAS on nonhuman genome
## submit run_gemma.sh from required_files in the gwas_gemma directory with the phenos.tsv and VCF in the same directory

sbatch --account=def-sjsmith run_gemma.sh
