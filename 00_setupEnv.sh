conda create -n fastq2vcf -c bioconda stacks fastqc multiqc bwa bcftools samtools fastp sambamba java-jdk
conda install -c bioconda openssl=1.0 # for bcftools issue

conda create -n ddocent -c bioconda ddocent # in case I decide to use dDocent

conda create -n gatk -c bioconda gatk4 picard sambamba java-jdk
