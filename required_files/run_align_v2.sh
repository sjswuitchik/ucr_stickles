#!/bin/bash
#SBATCH -J align
#SBATCH -o 02_align/out
#SBATCH -e 02_align/err
#SBATCH -n 1
#SBATCH --mem 10000
#SBATCH -t 24:00:00

# submit from ~/projects/def-sjsmith/sjsmith/stickles_ucr/seq_data/runv2

module load bwa

while read file
do
  bwa mem -O 5 -B 3 -a -M -t 16 -R ../../reference/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna 01_demulti/trimmed/filtered/$file.1.1.fq 01_demulti/trimmed/filtered/$file.1.2.fq > 02_align/$file.sam 2> 02_align/$file.bwa.log
done < 01_demulti/samples
