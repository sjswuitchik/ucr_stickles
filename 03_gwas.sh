# create env for gemma gwas
conda conda create -n gwas -c esgf -c bioconda -c ostrokach python plink vcftools bcftools gemma

# get gemma
git clone https://github.com/johanzi/gwas_gemma.git

## chromosomes need to be converted from alpahnumeric to numeric values. 
# create chromosome map
bcftools view -H stickles.filtered.recode.vcf | cut -f 1 | uniq | awk '{print $0}' > gasAcu.chrom.map.txt

## edit in R (because I don't have the capacity to loop this rn) gasAcu.chrom.map.txt using sequence report from NCBI (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016920845.1/) to assign numerical chr names, then save as gasAcu.chromMap

bcftools annotate stickles.filtered.recode.vcf --rename-chrs gasAcu.chromMap -o gasAcu.filter.chrRename.vcf -O v

##### notes for myself for when I pick up tomorrow 
# - need to remove OBBB_1 and OOB_1 from VCF - can also refilter for --max-missing-count 3 and --min-ac 1 and doublecheck that min/max alleles is at 2, and indels are removed

cp gasAcu.filter.chrRename.vcf gwas_gemma/
cd gwas_gemma

## use run_gwas_gemma.sh from required_files on this repo to preserve edits that are required to run GWAS on nonhuman genome
## submit run_gemma.sh from required_files in the gwas_gemma directory with the phenos.tsv and VCF in the same directory

sbatch --account=def-sjsmith run_gemma.sh
