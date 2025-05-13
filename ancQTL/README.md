# ancQTL

`ancQTL` is an R script adapted from [mixQTL](https://github.com/hakyimlab/mixqtl) for eQTL  mapping. This version extends `mixQTL` by including **local ancestry as an additional covariate**, making it particularly useful for admixed populations where ancestry effects are of interest.

## Features

* Based on the established `mixQTL` framework
* Adds support for **local ancestry covariate**
* Compatible with mixQTL input formats
* Suitable for **local ancestry based (LAB) eQTL** and **sQTL** analyses

## Requirements

* R (>= 4.0)
* Data inputs compatible with `mixQTL` (see below)
* Local ancestry estimates per individual per gene

## Installation

In R:

```r
source("path/to/ancQTL/ancQTL.R")
```

## Quick Start

The usage closely mirrors `mixQTL`, with one key addition: an input for local ancestry.

### Input Data

* **Genotype matrix (`genotypes`)**: Individuals × SNPs
* **Expression vector (`expression`)**: Gene expression values for individuals
* **Covariates matrix (`covariates`)**: Standard covariates (e.g., PCs, batch effects)
* **Local ancestry vector (`local_ancestry`)**: Same dimension as `genotypes`, with local ancestry dosages

### Example Usage

```r
result <- ancqtl(
  geno1 = genotype1, 
  geno2 = genotype2,
  y1 = haplotype_specific_count1,
  y2 = haplotype_specific_count2,
  anc1 = local_ancestry1,
  anc2 = local_ancestry2,
  ytotal = total_read_count,
  lib_size = seq_lib_size,
  cov_offset = covariate_offset_matrix
)
```

The function returns a table of association results including **local ancrstry based** effect sizes and p-values.

## Output

The output is a data frame with:

* SNP ID
* Effect size
* Standard error
* p-value
* Other relevant statistics

## Differences from mixQTL

| Feature             | mixQTL | ancQTL |
| ------------------- | ------ | ------ |
| Standard covariates | ✔      | ✔      |
| Local ancestry      | ✖      | ✔      |

