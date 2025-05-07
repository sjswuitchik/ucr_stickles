# in ~/projects/def-sjsmith/sjsmith/stickles_ucr/reference/GCF_016920845.1

conda activate bedtools 
module load perl # some conda issues with CC

chmod +x replace_chrs.pl
./replace_chrs.pl gasAcu_acckey gasAcu.gff > gasAcu.repl.gff
# manually add gff3 header back into gasAcu.repl.gff because I don't have time to re-write replace_chrs.pl right now

convert2bed --input=gff --output=bed < gasAcu.repl.gff > gasAcu_sorted.bed

