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
library(tidyverse)
```

* `mia` loads many functions that mask other packages:
  * To hide the warnings we use `suppressPackageStartupMessages`.
  * To avoid unpredictable results, we use `conflicted::conflicts_prefer` to choose which functions we actually want.
  * All function with duplicated names that are not specified with `conflicts_prefer(...)` will result in an error.


``` r
suppressPackageStartupMessages(library(mia))
conflicts_prefer(dplyr::first, .quiet = TRUE)
```

For packages what will be cited in the paper's methods:


``` r
packageVersion("mia")
```

```
## [1] '1.18.0'
```

``` r
citation("mia")
```

```
## To cite package 'mia' in publications use:
## 
##   Borman T, Ernst F, Shetty S, Lahti L (2025). _mia: Microbiome
##   analysis_. doi:10.18129/B9.bioc.mia
##   <https://doi.org/10.18129/B9.bioc.mia>, R package version 1.18.0,
##   <https://bioconductor.org/packages/mia>.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {mia: Microbiome analysis},
##     author = {Tuomas Borman and Felix G.M. Ernst and Sudarshan A. Shetty and Leo Lahti},
##     year = {2025},
##     note = {R package version 1.18.0},
##     url = {https://bioconductor.org/packages/mia},
##     doi = {10.18129/B9.bioc.mia},
##   }
```


# Load

* Make liberal use of R's native pipe (`|>`).
* Instead of loading every package with `library(...)` we can also call it directly with `package::function`.
* `qs2` is a modern version of `rds` (`?readRDS`).


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

* `head(c(number_or_rows, number_of_cols))` keeps the example output small. **Do not use when further analyzing any of the output.**
* A `L` trailing a number specified the number is an integer.


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
# raw counts
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
# relative abundance (%)
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
# center log-ratios
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
metadata(se)[["clr_pseudocount"]]
```

```
## [1] 0.5
```

## Samples: BioSample/Library metadata


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

## Features: ASV/taxa taxonomy and metadata


``` r
# standardized ranks
taxonomyRanks(se)
```

```
## [1] "Domain"  "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
```

``` r
# JMF-style ranks
metadata(se)[["taxonomy_ranks"]][["current"]]
```

```
## [1] "Domain"  "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
## [8] "ASV_ID"
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

## more metadata

* `?str` keeps the example output small. **Do not use when further analyzing any of the output.**


``` r
# Bioconductor-style metadata (generated by the JMF)
metadata(se) |> str()
```

```
## List of 4
##  $ version          : int 1
##  $ taxonomy_ranks   :List of 2
##   ..$ initial: chr [1:8] "Domain" "Phylum" "Class" "Order" ...
##   ..$ current: chr [1:8] "Domain" "Phylum" "Class" "Order" ...
##  $ filtered_features:List of 5
##   ..$ decontam_p_value: tibble [3 × 2] (S3: tbl_df/tbl/data.frame)
##   .. ..$ Feature_ID      : chr [1:3] "ASV_5pr_40t" "ASV_lwc_fxo" "ASV_pu2_753"
##   .. ..$ decontam_p_value: Factor w/ 10 levels "[0,0.1)","[0.1,0.2)",..: 1 1 1
##   ..$ Domain          : tibble [2 × 2] (S3: tbl_df/tbl/data.frame)
##   .. ..$ Feature_ID: chr [1:2] "ASV_j8t_qvz" "ASV_qs0_9qr"
##   .. ..$ Domain    : chr [1:2] NA "Eukaryota"
##   ..$ Order           : tibble [23 × 2] (S3: tbl_df/tbl/data.frame)
##   .. ..$ Feature_ID: chr [1:23] "ASV_3wo_c77" "ASV_529_cyh" "ASV_5ie_w3h" "ASV_6mm_8ja" ...
##   .. ..$ Order     : chr [1:23] "Chloroplast" "Chloroplast" "Chloroplast" "Chloroplast" ...
##   ..$ Family          : tibble [137 × 2] (S3: tbl_df/tbl/data.frame)
##   .. ..$ Feature_ID: chr [1:137] "ASV_10n_49z" "ASV_150_v2b" "ASV_1hf_jy2" "ASV_1vo_m1p" ...
##   .. ..$ Family    : chr [1:137] "Mitochondria" "Mitochondria" "Mitochondria" "Mitochondria" ...
##   ..$ ASV_ID          : tibble [0 × 2] (S3: tbl_df/tbl/data.frame)
##   .. ..$ Feature_ID: chr(0) 
##   .. ..$ ASV_ID    : chr(0) 
##  $ clr_pseudocount  : num 0.5
```

``` r
# provenance metadata generated by the JMF
se |>
  attr("provenance", exact = TRUE) |>
  str()
```

```
## List of 6
##  $ project      : chr "JMF-1906-4"
##  $ gene         : chr "16S rRNA V4 (515F/806R)"
##  $ resolution   : chr "samples"
##  $ state        : chr "refined"
##  $ rank         : chr "ASV_ID"
##  $ sample filter:List of 1
##   ..$ ≥: int 1000
```


# Further reading:

* https://microbiome.github.io/OMA/
