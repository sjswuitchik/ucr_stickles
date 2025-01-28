#!/bin/bash
#SBATCH -J gemma
#SBATCH -o out
#SBATCH -e err
#SBATCH -n 1
#SBATCH --mem 9000
#SBATCH -t 10000

conda init
conda activate gwas
bash run_gwas_gemma.sh phenos.tsv gasAcu.filter.chrRename.vcf
