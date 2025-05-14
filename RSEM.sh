#Make RSEM file (Gene expression)

#1. Prepare RSEM reference
rsem-prepare-reference \
  --gtf Homo_sapiens.GRCh38.109.gtf \
  --star \
  --star-path /home/fff0996/STAR/bin/Linux_x86_64_static \
  Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  rsem_ref_GRCh38
rsem-extract-reference-transcripts rsem_ref_GRCh38 0 Homo_sapiens.GRCh38.109.gtf None 0 Homo_sapiens.GRCh38.dna.primary_assembly.fa

#2. Make RSEM file (Calculate gene expression)
 rsem-calculate-expression   --star   --star-path /home/fff0996/STAR/bin/Linux_x86_64_static
--star-gzipped-read-file   -p 8   --paired-end S01-002TR220408_1.fastq.gz S01-002TR220408_2.fastq.gz
/home/fff0996/Reference/rsem_ref_GRCh38   S01-002
