fn_bed=./genome_snp.bed
fn_bam=aln.sorted.bam

samtools bedcov $fn_bed $fn_bam|awk '{sum += $4}; END {print sum}'
