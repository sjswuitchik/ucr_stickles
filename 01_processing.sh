# in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr/seq_data

conda create -n stacks -c bioconda stacks fastqc multiqc bwa samtools bcftools

# process radtags
mkdir 01_process_fastq

process_radtags -1 00_raw_fastq/NS.LH00147_0009.006.B723---B503.THigham_20230515_plate1_R1.fastq.gz -2 00_raw_fastq/NS.LH00147_0009.006.B723---B503.THigham_20230515_plate1_R2.fastq.gz -o 01_process_fastq -b barcodes.txt -e apeKI -y fastq -r -c -q

# 178395276 total sequences
# 82910 barcode not found drops (0.0%)
# 995028 low quality read drops (0.6%)
# 1383657 RAD cutsite not found drops (0.8%)
# 175933681 retained reads (98.6%)

stacks-dist-extract --pretty 01_process_fastq/process_radtags.log per_barcode_raw_read_counts > 01_process_fastq/process.sum

rm 01_process_fastq/*.rem.*

# extract fastqc reports
mkdir -p 01_process_fastq/fastqc

for file in 01_process_fastq/*.fq
do
  fastqc $file -o 01_process_fastq/fastqc
done

multiqc 01_process_fastq/fastqc

# remove PCR duplicates
cd 01_process_fastq
ls *.1.fq | sed '/\.1\.fq/s///' > samples

while read file
do
  clone_filter -1 $file.1.fq -2 $file.1.fq -i fastq -D -o filtered >> $file.log
done < samples 

cd ..

# align to reference genome
mkdir 02_align

while read file
do
  bwa mem -t 16 ../reference/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna 01_process_fastq/filtered/$file.1.1.fq 01_process_fastq/filtered/$file.1.2.fq > 02_align/$file.sam 2> 02_align/$file.bwa.log
done < 01_process_fastq/samples

# remove ambiguous reads & created sorted BAMs
cd 02_align
ls *.sam | sed '/\.sam/s///' > samples

while read file 
do
  samtools view -q 20 -b -S $file > $file.bam
  samtools sort $file.bam - o $file.sort.bam
done < samples
