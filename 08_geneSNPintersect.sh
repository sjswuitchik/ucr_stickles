# in ~/projects/def-sjsmith/sjsmith/stickles_ucr/

conda activate bedtools 

mkdir geneAssoc
cp reference/GCF_016920845.1/gasAcu.gff eference/GCF_016920845.1/gasAcu_acckey geneAssoc/
cd geneAssoc/

# convert GFF to BED
convert2bed --input=gff --output=bed < gasAcu.gff > gasAcu.bed

# GFF has whitespaces in the notes, which doesn't work with current version of replace_chrs.pl and I don't have time to rewrite that, so we're doing this instead
sed 's/ //g' gasAcu.bed > gasAcu.nws.bed
./replace_chrs.pl gasAcu_acckey gasAcu.nws.bed > gasAcu.nws.repl.bed

bedtools sort -i gasAcu.nws.repl.bed > gasAcu.sorted.bed

bedtools sort -i mce_sig_forBED.bed | bedtools intersect -a - -b gasAcu.sorted.bed -wb > mce_windows.bed
bedtools sort -i md_sig_forBED.bed | bedtools intersect -a - -b gasAcu.sorted.bed -wb > md_windows.bed
bedtools sort -i ppdmg_sig_forBED.bed | bedtools intersect -a - -b gasAcu.sorted.bed -wb > ppdmg_windows.bed
bedtools sort -i rs_sig_forBED.bed | bedtools intersect -a - -b gasAcu.sorted.bed -wb > rs_windows.bed
bedtools sort -i thdvmg_sig_forBED.bed | bedtools intersect -a - -b gasAcu.sorted.bed -wb > thdvmg_windows.bed
bedtools sort -i ttpg_sig_forBED.bed | bedtools intersect -a - -b gasAcu.sorted.bed -wb > ttpg_windows.bed



