# different conda envs used in different scripts
conda create -n fastq2vcf -c bioconda stacks fastqc multiqc bwa bcftools samtools fastp sambamba java-jdk
conda activate fastq2vcf
conda install -c bioconda openssl=1.0 # for bcftools issue

conda create -n vcfFilt -c bioconda plink vcftools htslib bcftools rename

conda create -n gwas -c esgf -c bioconda -c ostrokach -c conda-forge python=3.10 plink vcftools bcftools gemma libgfortran5 

conda create -n plink2 -c bioconda plink2 gcta

conda create -n bedtools -c bioconda -c conda-forge bedtools bedops
