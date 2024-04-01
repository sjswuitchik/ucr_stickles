# from /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr/seq_data/runv2
## testing if issues with initial processing

conda activate stacks
mkdir 01_demulti

process_radtags -1 ../00_raw_fastq/NS.LH00147_0009.006.B723---B503.THigham_20230515_plate1_R1.fastq.gz -2 ../00_raw_fastq/NS.LH00147_0009.006.B723---B503.THigham_20230515_plate1_R2.fastq.gz -o 01_demulti -i gzfastq -b barcodes.txt -e apeKI -y fastq -r 
# 178395276 total sequences
# 82910 ambiguous barcode drops (0.0%)
# 0 low quality read drops (0.0%)
# 1383657 ambiguous RAD-Tag drops (0.8%)
# 176928709 retained reads (99.2%)

rm *rem*

mkdir 01_demulti/fastqc
for file in 01_demulti/*.fq
do
  fastqc $file -o 01_demulti/fastqc
done
multiqc 01_demulti/fastqc

mv multiqc_data/ multiqc_data_demulti/
mv multiqc_report.html multiqc_report_demulti.html

cd 01_demulti
ls *.1.fq | sed '/\.1\.fq/s///' > samples
mkdir -p trimmed/fastP_out

while read file
do
  fastp --in1 $file.1.fq --in2 $file.2.fq --out1 trimmed/$file.R1.fq --out2 trimmed/$file.R2.fq -q 15 -u 50 -t 1 -T 1 -c --dedup -h fastP_out/$file.fp.html
done < samples

mkdir -p trimmed/fastqc
cd ..

for file in 01_demulti/trimmed/*.fq
do
  fastqc $file -o 01_demulti/trimmed/fastqc
done
multiqc 01_demulti/trimmed/fastqc

mv multiqc_data/ multiqc_data_filtered/
mv multiqc_report.html multiqc_report_filtered.html

# align reads to ref
mkdir 02_align

while read file
do
  bwa mem -O 5 -B 3 -a -M -t 16 ../../reference/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna 01_demulti/trimmed/$file.R1.fq 01_demulti/trimmed/$file.R2.fq > 02_align/$file.sam 2> 02_align/$file.bwa.log
done < 01_demulti/samples

# check alignment stats 
cd 02_align
ls *.sam | sed '/\.sam/s///' > samples

while read file
do
  samtools view -Sbt ../../../reference/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna 02_align/$file.sam | samtools flagstat - > $file.aln.log
done < samples

# convert SAM to BAM, sort, & index  

while read file 
do
  samtools view -q 20 -b -S $file.sam > $file.bam
  samtools sort $file.bam -o $file.sort.bam
  samtools index $file.bam -o $file.sort.idx
done < samples

# mark dups
while read file
do
  sambamba markdup $file.sort.bam $file.markdup.bam
done < samples

# QC & filtering
while read file
do
  samtools coverage --output $file.cov.txt $file.sort.bam
  samtools flagstat -O tsv $file.markdup.bam > $file.markdup.aln.stat
done < samples










