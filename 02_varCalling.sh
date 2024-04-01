# in /home/sjsmith/projects/def-sjsmith/sjsmith/stickles_ucr/seq_data/runv2
conda activate gatk

mkdir -p 03_vcf

gatk HaplotypeCaller -R ../../reference/GCF_016920845.1/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna -I 02_align/BAM_01.markdup.bam -O 03_gatk/gvcfs/BAM_01.gvcf --emit-ref-confidence GVCF --min-pruning 1 --min-dangling-branch-length 1 
