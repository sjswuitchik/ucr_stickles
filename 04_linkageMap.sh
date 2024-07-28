# from /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr

module load StdEnv/2023 java/21.0.1

## Download LepMap3 from SourceForge
mkdir -p linkageMapping/lepmap
cd linkageMapping/lepmap 
wget https://sourceforge.net/projects/lep-map3/files/latest/download
unzip download

## Download SnpEff
# NTS: may not use snpEff? TBD
cd ..
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip
mkdir -p snpEff/data/gasAcu

cp ../vcf_filt/stickles_ucr.dedup.vcf .

conda activate vcfFilt
vcftools --vcf stickles_ucr.dedup.vcf --remove-indv dedup/OBBB_1.dedup.bam --remove-indv dedup/OOB_1.dedup.bam --012 --out stickles.dedup.f2s



