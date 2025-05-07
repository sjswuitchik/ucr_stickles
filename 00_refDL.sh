# in projects/def-sjsmith/sjsmith/stickles_ucr/reference

conda create -n ncbi -c conda-forge ncbi-datasets-cli
conda activate ncbi

datasets download genome accession GCF_016920845.1 --include gff3,genome,seq-report
unzip ncbi_dataset.zip
mv ncbi_dataset/data/* .
rm -r *.zip ncbi_dataset *.json*
mv GCF_016920845.1/genomic.gff GCF_016920845.1/gasAcu.gff

conda deactivate
#conda create -n stacks -c bioconda stacks fastqc multiqc bwa samtools bcftools (also outlined in 00_setupEnv.sh)
conda activate stacks

bwa index GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna
