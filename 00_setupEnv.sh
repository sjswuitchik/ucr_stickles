conda create -n stacks -c bioconda stacks fastqc multiqc bwa bcftools samtools fastp
conda install -c bioconda openssl=1.0 # for bcftools issue

conda create -n ddocent -c bioconda ddocent # in case I decide to use dDocent