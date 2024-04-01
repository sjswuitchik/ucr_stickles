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
mkdir -p trimmed/filtered
mkdir -p trimmed/fastP_out

while read file
do
  fastp --in1 $file.1.fq --in2 $file.2.fq --out1 trimmed/$file.R1.fq --out2 trimmed/$file.R2.fq -q 15 -u 50 -t 1 -T 1 -c --dedup -h fastP_out/$file.fp.html
done < samples

while read file
do
  clone_filter -1 trimmed/$file.R1.fq -2 trimmed/$file.R2.fq -i fastq -D -o trimmed/filtered >> $file.log
done < samples 

mkdir -p trimmed/filtered/fastqc
for file in trimmed/filtered/*.fq.gz
do
  fastqc $file -o trimmed/filtered/fastqc
done
multiqc trimmed/filtered/fastqc

mv multiqc_data/ multiqc_data_filtered/
mv multiqc_report.html multiqc_report_filtered.html

cd ..
mkdir 02_align

## NTS: check the file extensions after clone_filter before running. 
while read file
do
  bwa mem -O 5 -B 3 -a -M -R ../../reference/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna 01_demulti/trimmed/filtered/$file.1.1.fq 01_demulti/trimmed/filtered/$file.1.2.fq > 02_align/$file.sam 2> 02_align/$file.bwa.log
done < 01_demulti/samples

# check alignment stats 
cd 02_align
ls *.sam | sed '/\.sam/s///' > samples

while read file
do
  samtools -view Sbt ../../reference/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna 02_align/$file.sam | samtools flagstat - 
done < 02_align/samples

# convert SAM to BAM, & sort 

while read file 
do
  samtools view -q 20 -b -S $file.sam > $file.bam
  samtools sort $file.bam -o $file.sort.bam
  samtools index $file.bam -o $file.sort.index
done < samples

# calculate coverage 








