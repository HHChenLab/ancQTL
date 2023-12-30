phaser.py \
        --bam $RNAseq.waspfilter.bam \
        --vcf $subjectID.vcf.gz \
        --sample $subjectID \
        --o $output \
        --mapq 255 --baseq 0 --paired_end  1 \
        --blacklist hg38_hla.chr.bed.gz \
        --haplo_count_blacklist hg38_haplo_count_blacklist.chr.bed.gz \
        --threads 4 \
        --gw_phase_vcf 1
phaser_gene_ae.py \
        --haplotypic_counts $output.haplotypic_counts.txt \
        --features gencode.v26.GRCh38.genes.chr.bed \
        --o $output.gene_ae.txt
