conda create -n fastq2vcf -c bioconda stacks fastqc multiqc bwa bcftools samtools fastp sambamba java-jdk
conda activate fastq2vcf
conda install -c bioconda openssl=1.0 # for bcftools issue
