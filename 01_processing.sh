# in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr/fastq_dl
## would be better with a conda env 

conda create -n stacks -c bioconda stacks fastqc

module load StdEnv/2020 stacks/2.64 

# process radtags
mkdir process_fastq
process_radtags -1 NS.LH00147_0009.006.B723---B503.THigham_20230515_plate1_R1.fastq.gz -2 NS.LH00147_0009.006.B723---B503.THigham_20230515_plate1_R2.fastq.gz -o process_fastq -b barcodes.txt -e apeKI -y fastq -r -c -q
rm process_fastq/*.mem.*

# extract fastqc reports
module load StdEnv/2023 fastqc/0.12.1

for file in process_fastq/*.1.fq
do
  fastqc $file.1.fq -o demulti/fastqc
  fastqc $file.2.fq -o demulti/fastqc
done
