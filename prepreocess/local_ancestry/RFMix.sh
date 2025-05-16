rfmix -f phased_target_individuals.vcf.gz \
        -r phased_references.vcf.gz \
        -m reference_ancestry.txt \
        -g genetic_map.gmap \
        -o output \
        --chromosome=chr$i
