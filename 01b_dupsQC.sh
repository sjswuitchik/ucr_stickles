# in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr/seq_data/02_align
## removing dups with sambamba then generating VCF to compare with impact of downstream removal of dups from VCF (and/or impact of dups on variant calling)

conda activate fastq2vcf

mkdir dedup

## Remove duplicates
while read file
do
  sambamba markdup -r $file.sort.bam dedup/$file.dedup.bam
done < samples

# QC & filtering
while read file
do
  samtools flagstat dedup/$file.dedup.bam > $file.dedup.aln.stat
done < samples

## Call variants & generate VCF with BCFtools
ls *.dedup.bam > bamList 
cd ..
mkdir 03_vcf

bcftools mpileup -C 50 -E -I --max-depth 8000 -f ../reference/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna -b bamList -O u > 03_vcf/stickles_ucr.dedup.bcf
bcftools call -v -c -f GQ 03_vcf/stickles_ucr.dedup.bcf > 03_vcf/stickles_ucr.dedup.vcf
