# in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr/fastq_dl

module load StdEnv/2020 stacks/2.64

process_radtags -1 NS.LH00147_0009.006.B723---B503.THigham_20230515_plate1_R1.fastq.gz -2 NS.LH00147_0009.006.B723---B503.THigham_20230515_plate1_R2.fastq.gz -o ./process_fastq -b barcodes.txt -e apeKI -r -c -q
