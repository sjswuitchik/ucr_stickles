#!/bin/bash
#SBATCH -J align
#SBATCH -o 02_align/out
#SBATCH -e 02_align/err
#SBATCH -n 1
#SBATCH --mem 9000
#SBATCH -t 03:00:00

# submit from ~/projects/def-sjsmith/sjsmith/stickles_ucr/seq_data

conda activate stacks

while read file
do
  bwa mem -t 16 ../reference/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna 01_process_fastq/filtered/$file.1.1.fq 01_process_fastq/filtered/$file.1.2.fq > 02_align/$file.sam 2> 02_align/$file.bwa.log
done < 01_process_fastq/samples
