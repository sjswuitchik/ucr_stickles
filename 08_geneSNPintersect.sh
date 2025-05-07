# in ~/projects/def-sjsmith/sjsmith/stickles_ucr/reference/GCF_016920845.1

conda activate bedtools 
module load perl # some conda issues with CC

gff2bed < gasAcu.gff > gasAcu_sorted.bed

