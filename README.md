# ancQTL
local ancestry-based eQTL mapping

## Tools
* [STAR](https://github.com/alexdobin/STAR)
* [bcftools](https://samtools.github.io/bcftools/bcftools.html)
* [phASER](https://github.com/secastel/phaser)
* [RNA-SeQC2](https://github.com/getzlab/rnaseqc)
* [mixQTL](https://github.com/hakyimlab/mixqtl)
* [RFMix2] 
## Packages
### R
* DESeq2
* PEER
* GAP
* data.table

## Local ancestry inference
Step 1:
RFMix2
Step 2:
Converting to gene level

## Haplotypic expression estimation 
Step 1:  
STAR alignment  
Step 2:  
WASP filter  
Step 3:  
phaser count  
Step 4:  
rnaseqc count  
step5:  
PEER and library size  

