# in ~/projects/def-sjsmith/sjsmith/stickles_ucr/reference/GCF_016920845.1

conda activate bedtools 

# convert GFF to BED
convert2bed --input=gff --output=bed < gasAcu.gff > gasAcu.sorted.bed

# GFF has whitespaces in the notes, which doesn't work with current version of replace_chrs.pl and I don't have time to rewrite that, so we're doing this instead
join gasAcu_acckey gasAcu.sorted.bed > gasAcu.sorted.repl.int.bed
awk -F" " 'BEGIN{OFS=" "}{for(i=2;i<NF;i++){printf "%s%s",$i,OFS} print $NF} END{OFS="\t"}' gasAcu.sorted.repl.int.bed > gasAcu.sorted.repl.bed
######## this didn't work, try after conference

mkdir ../../geneAssoc
mv gasAcu.sorted.* ../../geneAssoc
cd ../../geneAssoc

bedtools sort mce_sig_forBED.bed | bedtools intersect -a - -b gasAcu.sorted.repl.bed | less


