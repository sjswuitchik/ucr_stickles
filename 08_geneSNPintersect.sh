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
awk '$8 == "gene"' gasAcu.sorted.bed > gasAcu.gene.bed

## make chrom sizes file for windows  
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod +x ./faToTwoBit

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
chmod +x ./twoBitInfo

 ./faToTwoBit ../reference/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna gasAcu.fa.2bit
 ./twoBitInfo gasAcu.fa.2bit stdout | sort -k2rn > gasAcu.chrom.sizes

bedtools slop -i gasAcu.gene.bed -g gasAcu.chrom.sizes -b 100000 | cut -f1,2,3,10 > gasAcu.windows.bed

bedtools merge -i gasAcu.windows.bed -d -1 -c 4 -o distinct | bedtools sort -i - | less -S  

##### this needs to be dealt with tomorrow to make sure the SNP is being captured within the window
##### intersect, merge, perhaps slop?
### some more dev: 
awk '$3 == "CDS"' gasAcu.gff > test.onlyCDS.gff
awk -f gff2bed.awk test.onlyCDS.gff > test.onlyCDS.bed 
cat test.onlyCDS.bed | python genenames.py > test.onlyCDS.genes.bed
bedtools slop -i test.onlyCDS.genes.bed -g gasAcu.chrom.sizes -b 100000 | bedtools sort -i - > test.windows.bed
sed 's/NA/MT/g' gasAcu_acckey > test.acckey
./replace_chrs.pl test.acckey test.windows.bed > test.windows.repl.bed

##### the stuff between the hashes works, clean this up post-conference 

bedtools sort -i mce_sig_forBED.bed | bedtools intersect -a - -b gasAcu.sorted.bed -wb > mce_windows.bed
bedtools sort -i md_sig_forBED.bed | bedtools intersect -a - -b gasAcu.sorted.bed -wb > md_windows.bed
bedtools sort -i ppdmg_sig_forBED.bed | bedtools intersect -a - -b gasAcu.sorted.bed -wb > ppdmg_windows.bed
bedtools sort -i rs_sig_forBED.bed | bedtools intersect -a - -b gasAcu.sorted.bed -wb > rs_windows.bed
bedtools sort -i thdvmg_sig_forBED.bed | bedtools intersect -a - -b gasAcu.sorted.bed -wb > thdvmg_windows.bed
bedtools sort -i ttpg_sig_forBED.bed | bedtools intersect -a - -b gasAcu.sorted.bed -wb > ttpg_windows.bed



