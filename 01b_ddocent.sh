# from /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr/seq_data/ddocent
## trying the dDocent pipeline to test if issues with initial processing

conda create -n ddocent -c bioconda ddocent

conda activate stacks
mkdir 01_demulti

process_radtags -1 ../00_raw_fastq/NS.LH00147_0009.006.B723---B503.THigham_20230515_plate1_R1.fastq.gz -2 ../00_raw_fastq/NS.LH00147_0009.006.B723---B503.THigham_20230515_plate1_R2.fastq.gz -o 01_demulti -i gzfastq -b ../barcodes.txt -e apeKI -y fastq -r 
# 178395276 total sequences
# 82910 ambiguous barcode drops (0.0%)
# 0 low quality read drops (0.0%)
# 1383657 ambiguous RAD-Tag drops (0.8%)
# 176928709 retained reads (99.2%)

