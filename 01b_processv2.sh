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

cd 01_demulti
ls *.1.fq | sed '/\.1\.fq/s///' > samples
mkdir filtered 

while read file
do
  clone_filter -1 $file.1.fq -2 $file.2.fq -i fastq -D -o filtered >> $file.log
done < samples 

## fastp

cd ..
mkdir 02_align
sbatch --account=def-sjsmith run_align_v2.sh

