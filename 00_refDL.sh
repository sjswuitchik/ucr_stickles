# in projects/def-sjsmith/sjsmith/stickles_ucr/reference

conda create -n ncbi -c conda-forge ncbi-datasets-cli
conda activate ncbi

datasets download genome accession GCF_016920845.1 --include gff3,genome,seq-report
unzip ncbi_dataset.zip
mv ncbi_dataset/data/* .
rm -r *.zip ncbi_dataset *.json*

conda deactivate
conda activate stacks
conda update samtools --force-reinstall # issue with libcrypto

samtools faidx GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna
