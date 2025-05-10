# in ~/projects/def-sjsmith/sjsmith/stickles_ucr/

conda activate bedtools 

mkdir geneAssoc
cp reference/GCF_016920845.1/gasAcu.gff eference/GCF_016920845.1/gasAcu_acckey geneAssoc/
cd geneAssoc/

# convert GFF to BED & pull out only CDS regions
awk '$3 == "CDS"' gasAcu.gff > gasAcu.onlyCDS.gff
awk -f gff2bed.awk gasAcu.onlyCDS.gff > gasAcu.onlyCDS.bed 
cat gasAcu.onlyCDS.bed | python genenames.py > gasAcu.onlyCDS.genes.bed

## make chrom sizes file for window creation  
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod +x ./faToTwoBit

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
chmod +x ./twoBitInfo

./faToTwoBit ../reference/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna gasAcu.fa.2bit
./twoBitInfo gasAcu.fa.2bit stdout | sort -k2rn > gasAcu.chrom.sizes

# create windows
bedtools slop -i gasAcu.onlyCDS.genes.bed -g gasAcu.chrom.sizes -b 100000 | bedtools sort -i - > gasAcu.windows.bed

# replace NCBI chromosome names with linkage group numbers (and MT instead of NA)
sed 's/NA/MT/g' gasAcu_acckey > gasAcu.acckey
./replace_chrs.pl gasAcu.acckey gasAcu.windows.bed > gasAcu.windows.repl.bed

# intersect clean genes BED with windows
bedtools sort -i cleanGenes.bed | bedtools intersect -a gasAcu.windows.repl.bed -b - -wb | bedtools sort -i - | bedtools merge -i - -d 1 -c 4,8,9 -o distinct > gasAcu.genes.intersect.bed

# clean up gene names field by removing unnamed LOCs
sed 's/LOC[0-9]\+,//g' gasAcu.genes.intersect.bed > gasAcu.genes.clean.bed 
