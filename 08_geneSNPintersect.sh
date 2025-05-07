# in ~/projects/def-sjsmith/sjsmith/stickles_ucr/reference/GCF_016920845.1

conda activate bedtools 
module load perl # some conda issues with CC

chmod +x replace_chrs.pl

convert2bed --input=gff --output=bed < gasAcu.gff > gasAcu.sorted.bed

