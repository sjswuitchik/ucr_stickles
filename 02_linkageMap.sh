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

# create a PLINK env to generate IBD stats for manaul clustering 
conda create -n plink -c bioconda plink
conda activate plink

cp seq_data/02_align/dedup/03_vcf/stickles_ucr.dedup.vcf linkageMapping/



