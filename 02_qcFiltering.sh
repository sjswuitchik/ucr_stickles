## in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr/

# conda create -n vcfFilt -c bioconda plink vcftools htslib bcftools rename
conda activate vcfFilt

mkdir vcf_filt
cp seq_data/02_align/dedup/03_vcf/stickles_ucr.dedup.vcf vcf_filt
cd vcf_filt

bgzip stickles_ucr.dedup.vcf -c > stickles_ucr.dedup.vcf.gz
tabix -p vcf stickles_ucr.dedup.vcf.gz

vcftools --gzvcf stickles_ucr.dedup.vcf.gz --remove-indels --max-missing 0.5 --mac 3 --minQ 30 --recode --recode-INFO-all --out stickles.filt1
# kept 470374/944649
vcftools --vcf stickles.filt1.recode.vcf --minDP 3 --recode --recode-INFO-all --out stickles.filt2
# kept 470374/470374
vcftools --vcf stickles.filt2.recode.vcf --missing-indv
# no one to remove
vcftools --vcf stickles.filt2.recode.vcf --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 0.75 --minDP 3 --recode --recode-INFO-all --out stickles.filt3
# kept 470374/470374
vcftools --vcf stickles.filt3.recode.vcf --maf 0.01 --minGQ 10 --max-meanDP 200 --recode --recode-INFO-all --out stickles.filt4
# kept 453918/470374 sites
vcftools --vcf stickles.filt4.recode.vcf --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 0.75 --maf 0.01 --minDP 3 --minGQ 10 --recode --recode-INFO-all --out stickles.filtered
# final step; kept 16466/470374 sites

mkdir int_filt
mv out.* stickles.filt* int_filt
tar zcvf int_filt.tar.gz int_filt/
rm -r int_filt/

