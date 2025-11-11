---
title: Using SummarizedExperiment (with JMF amplicon projects)
author: Bela Hausmann (Joint Microbiome Facility)
output: 
  html_document: 
    toc: true
    keep_md: true
---



``` r
library(conflicted)
library(jmf)
quiet()
```


``` r
library(tidyverse)
suppressPackageStartupMessages(library(mia))

conflicts_prefer(dplyr::first)
```

```
## [conflicted] Will prefer dplyr::first over any other package.
```


# Load


``` r
se <-
  fs::dir_ls(glob = "*.qs2", recurse = TRUE) |>
  first() |>
  print() |>
  qs2::qs_read()
```

```
## data/JMF-1906-4_16S_rRNA_V4_samples_refined_gte_1000_SE.qs2
```


# Inspect


``` r
se
```

```
## class: SingleCellExperiment 
## dim: 3464 26 
## metadata(4): version taxonomy_ranks filtered_features clr_pseudocount
## assays(3): counts relabundance clr
## rownames(3464): ASV_104_sx7 ASV_10f_39x ... ASV_zpu_k9u ASV_zzb_gea
## rowData names(11): Domain Phylum ... sequence_length decontam_p_value
## colnames(26): JMF-1906-4-0001 JMF-1906-4-0002 ... JMF-1906-4-0025
##   JMF-1906-4-0027
## colData names(41): JMF_sample_ID User_sample_ID ...
##   .alpha_diversity_at_2173_chao1_se .alpha_diversity_at_2173_shannon
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
```

``` r
assay(se) |> head(c(5L, 3L))
```

```
##             JMF-1906-4-0001 JMF-1906-4-0002 JMF-1906-4-0003
## ASV_104_sx7               2               0               1
## ASV_10f_39x               0               0               0
## ASV_10z_d3k               0               0               2
## ASV_116_8mi               0               0               0
## ASV_11y_e80               0               0               0
```

``` r
(assay(se, "relabundance") * 100.0) |> head(c(5L, 3L))
```

```
##             JMF-1906-4-0001 JMF-1906-4-0002 JMF-1906-4-0003
## ASV_104_sx7      0.01124164               0     0.007674008
## ASV_10f_39x      0.00000000               0     0.000000000
## ASV_10z_d3k      0.00000000               0     0.015348016
## ASV_116_8mi      0.00000000               0     0.000000000
## ASV_11y_e80      0.00000000               0     0.000000000
```

``` r
assay(se, "clr") |> head(c(5L, 3L))
```

```
##             JMF-1906-4-0001 JMF-1906-4-0002 JMF-1906-4-0003
## ASV_104_sx7       0.6986293      -0.5936341       0.1118882
## ASV_10f_39x      -0.9108086      -0.5936341      -0.9867241
## ASV_10z_d3k      -0.9108086      -0.5936341       0.6227138
## ASV_116_8mi      -0.9108086      -0.5936341      -0.9867241
## ASV_11y_e80      -0.9108086      -0.5936341      -0.9867241
```


``` r
colData(se) |> head(c(5L, 3L))
```

```
## DataFrame with 5 rows and 3 columns
##                   JMF_sample_ID User_sample_ID BioSample_accession
##                     <character>    <character>         <character>
## JMF-1906-4-0001 JMF-1906-4-0001     A1 mulched                  NA
## JMF-1906-4-0002 JMF-1906-4-0002     A2 mulched                  NA
## JMF-1906-4-0003 JMF-1906-4-0003     A3 mulched                  NA
## JMF-1906-4-0004 JMF-1906-4-0004   A4 grassland                  NA
## JMF-1906-4-0005 JMF-1906-4-0005   A5 grassland                  NA
```


``` r
taxonomyRanks(se)
```

```
## [1] "Domain"  "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
```

``` r
rowData(se) |> head(c(5L, 3L))
```

```
## DataFrame with 5 rows and 3 columns
##                  Domain          Phylum               Class
##             <character>     <character>         <character>
## ASV_104_sx7    Bacteria  Actinomycetota     Thermoleophilia
## ASV_10f_39x    Bacteria Planctomycetota      Planctomycetes
## ASV_10z_d3k    Bacteria  Pseudomonadota Gammaproteobacteria
## ASV_116_8mi    Bacteria  Pseudomonadota Alphaproteobacteria
## ASV_11y_e80    Bacteria              NA                  NA
```


# Further reading:

* https://microbiome.github.io/OMA/
