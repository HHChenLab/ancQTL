samtools view -h $STAR_output.Aligned.out.bam | \
        grep -v "vW:i:[2-7]" | \
        samtools view -h1 | samtools sort -@ 3 -T $output/temp --output-fmt bam > $output.waspfilter.bam
samtools index -@ 3 $output.waspfilter.bam
