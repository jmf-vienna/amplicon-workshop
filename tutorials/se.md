---
title: Using SummarizedExperiment (with JMF amplicon projects)
author: Bela Hausmann (Joint Microbiome Facility)
output: 
  html_document: 
    toc: true
    keep_md: true
---



``` r
options(width = 120L)
```


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
conflicts_prefer(dplyr::filter, .quiet = TRUE)
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
##   Borman T, Ernst F, Shetty S, Lahti L (2025). _mia: Microbiome analysis_. doi:10.18129/B9.bioc.mia
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
## colnames(26): JMF-1906-4-0001 JMF-1906-4-0002 ... JMF-1906-4-0025 JMF-1906-4-0027
## colData names(41): JMF_sample_ID User_sample_ID ... .alpha_diversity_at_2173_chao1_se
##   .alpha_diversity_at_2173_shannon
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
```

``` r
assayNames(se)
```

```
## [1] "counts"       "relabundance" "clr"
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

Columns are always samples (biosamples, libraries, merged meta-samples, etc).


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
# sample metadata variable names:
se |>
  colData() |>
  names() |>
  cat(sep = "\n")
```

```
## JMF_sample_ID
## User_sample_ID
## BioSample_accession
## Sample_description
## Dataset
## JMF_project_ID
## BioProject_accessions
## Group
## Submit_to_SRA
## Environment_ID
## MIxS_environmental_package
## env_broad_scale
## env_local_scale
## env_medium
## organism
## collection_date
## geo_loc_name
## ph
## d_13C_per_12C
## d_15N_per_14N
## dw_per_ww
## Group4
## Location
## mg_C_per_g_DWS
## mg_N_per_g_DWS
## mg_TDNC_per_g_DWS
## mg_TDOC_per_g_DWS
## Mulched
## Sample_code
## Soil_type
## Soil_type_2
## TC_per_TN
## TDOC_per_TDN
## .Number_of_libraries
## .Libraries
## .alpha_diversity_chao1
## .alpha_diversity_chao1_se
## .alpha_diversity_shannon
## .alpha_diversity_at_2173_chao1
## .alpha_diversity_at_2173_chao1_se
## .alpha_diversity_at_2173_shannon
```

``` r
# access single variables
head(se$Group)
```

```
## [1] "mulched_A"   "mulched_A"   "mulched_A"   "grassland_A" "grassland_A" "grassland_A"
```

``` r
# the more formal way:
colData(se)[, "Group"] |> head()
```

```
## [1] "mulched_A"   "mulched_A"   "mulched_A"   "grassland_A" "grassland_A" "grassland_A"
```

* RStudio's table viewer can not use Bioconductor constructs directly, so conversion to data.frames or tibbles is needed.

```r
se |> colData() |> as.data.frame() |> View()
```



## Features: ASV/taxa taxonomy and metadata

Rows are always features, i.e., ASVs, OTUs, genera, etc.


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
## [1] "Domain"  "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species" "ASV_ID"
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

``` r
# Bioconductor/tidyverse mixed usage:
se |>
  rowData() |>
  as_tibble() |>
  select(Domain:ASV_ID) |>
  head()
```

```
## # A tibble: 6 × 8
##   Domain   Phylum          Class               Order               Family            Genus        Species ASV_ID     
##   <chr>    <chr>           <chr>               <chr>               <chr>             <chr>        <chr>   <chr>      
## 1 Bacteria Actinomycetota  Thermoleophilia     Solirubrobacterales 67-14             <NA>         <NA>    ASV_104_sx7
## 2 Bacteria Planctomycetota Planctomycetes      Pirellulales        Pirellulaceae     Pir4 lineage <NA>    ASV_10f_39x
## 3 Bacteria Pseudomonadota  Gammaproteobacteria Burkholderiales     <NA>              <NA>         <NA>    ASV_10z_d3k
## 4 Bacteria Pseudomonadota  Alphaproteobacteria <NA>                <NA>              <NA>         <NA>    ASV_116_8mi
## 5 Bacteria <NA>            <NA>                <NA>                <NA>              <NA>         <NA>    ASV_11y_e80
## 6 Bacteria Gemmatimonadota Gemmatimonadia      Gemmatimonadales    Gemmatimonadaceae Gemmatimonas <NA>    ASV_122_kp1
```

Extra feature metadata added by JMF:


``` r
rowData(se)$Lineage |> first()
```

```
## [1] "Bacteria ‣ Actinomycetota ‣ Thermoleophilia ‣ Solirubrobacterales ‣ 67-14 ‣ ‣ ‣ ASV_104_sx7"
```

``` r
rowData(se)$sequence_length |> first()
```

```
## [1] 253
```

``` r
# decontam P-values are only usable if negative controls were correctly assigned and used (see user sample sheet and config)
rowData(se)$decontam_p_value |> first()
```

```
## [1] 0.999018
```


## everything


``` r
mse <-
  se |>
  meltSE(assay.type = "relabundance", add.row = TRUE, add.col = TRUE) |>
  mutate(relabundance := relabundance * 100L)

dim(mse)
```

```
## [1] 90064    55
```


``` r
# Arbitrary example:
mse |>
  filter(Family == "Nitrososphaeraceae") |>
  group_by(across(c(User_sample_ID, Group, Family:Genus))) |>
  summarise(relabundance = mean(relabundance), .groups = "drop") |>
  arrange(-relabundance) |>
  head(4L)
```

```
## # A tibble: 4 × 5
##   User_sample_ID Group       Family             Genus                      relabundance
##   <chr>          <chr>       <chr>              <chr>                             <dbl>
## 1 T6 grassland   grassland_T Nitrososphaeraceae <NA>                              2.35 
## 2 T7 forested    forested_T  Nitrososphaeraceae <NA>                              1.15 
## 3 A8 forested    forested_A  Nitrososphaeraceae <NA>                              0.775
## 4 A8 forested    forested_A  Nitrososphaeraceae Candidatus Nitrosocosmicus        0.774
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
  attr("provenance") |>
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

* JMF workflow functions can also be quickly (re-)used like this for convenience:


``` r
source("https://raw.githubusercontent.com/jmf-vienna/amplicon-analysis/refs/heads/main/R/provenance.R")

se |>
  get_provenance() |>
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

``` r
provenance_as_short_title(se)
```

```
## [1] "samples/refined/ASV_ID/1000"
```

``` r
provenance_as_tibble(se)
```

```
## # A tibble: 1 × 6
##   project    gene                    resolution state   rank   `sample filter ≥`
##   <chr>      <chr>                   <chr>      <chr>   <chr>              <int>
## 1 JMF-1906-4 16S rRNA V4 (515F/806R) samples    refined ASV_ID              1000
```

## final tips


``` r
# number of features:
nrow(se)
```

```
## [1] 3464
```

``` r
# number of samples:
ncol(se)
```

```
## [1] 26
```

``` r
# internal feature IDs
rownames(se) |> head()
```

```
## [1] "ASV_104_sx7" "ASV_10f_39x" "ASV_10z_d3k" "ASV_116_8mi" "ASV_11y_e80" "ASV_122_kp1"
```

``` r
# internal sample IDs
colnames(se) |> head()
```

```
## [1] "JMF-1906-4-0001" "JMF-1906-4-0002" "JMF-1906-4-0003" "JMF-1906-4-0004" "JMF-1906-4-0005" "JMF-1906-4-0006"
```


# Amend metadata

## Changing existing sample metadata variables


``` r
# Converting a nominal variable to an ordinal variable (factor)
colData(se)$Group <- colData(se)$Group |> fct_infreq()
```

## Creating new sample metadata variables from existing ones

* Prefer `%in%` over `==` in vector comparison when `NA`s should NOT be preserved.


``` r
colData(se)$Replicate_number <- colData(se)$User_sample_ID |> str_extract("[0-9]+")
colData(se)$is_grassland <- colData(se)$Soil_type %in% "grassland soil"
colData(se)$is_grassland_with_NAs <- colData(se)$Soil_type == "grassland soil" # DIFFERENT RESULT!
```

## Added new data from external tables


``` r
extra_metadata <-
  fs::path("data", "soils.tsv") |>
  read_tsv() |>
  print()
```

```
## # A tibble: 3 × 2
##   Soil_type Color
##   <chr>     <chr>
## 1 grassland green
## 2 forested  brown
## 3 mulched   brown
```

* Use `left_join` to make sure no samples are missing after merging tables.
* Specify `by` when merging tables to catch errors earlier.


``` r
new_metadata <-
  se |>
  colData() |>
  as.data.frame() |>
  rownames_to_column() |>
  as_tibble() |>
  dplyr::left_join(extra_metadata, by = "Soil_type") |>
  column_to_rownames("rowname") |>
  as("DataFrame")

# ASSERTION: make sure nothing unexpected happened to number and order of samples!
new_metadata |>
  rownames() |>
  identical(colnames(se)) |>
  stopifnot()

# debug issues with `waldo`:
waldo::compare(colData(se), new_metadata)
```

```
## `old@listData` is length 44
## `new@listData` is length 45
## 
## `names(old@listData)[42:44]`: "Replicate_number" "is_grassland" "is_grassland_with_NAs"        
## `names(new@listData)[42:45]`: "Replicate_number" "is_grassland" "is_grassland_with_NAs" "Color"
## 
## `old@listData$Color` is absent
## `new@listData$Color` is a character vector ('brown', 'brown', 'brown', 'green', 'green', ...)
```

``` r
colData(se) <- new_metadata
```


# Subsetting

## Subset samples


``` r
se[, se$Soil_type %in% "forested"]
```

```
## class: SingleCellExperiment 
## dim: 3464 7 
## metadata(4): version taxonomy_ranks filtered_features clr_pseudocount
## assays(3): counts relabundance clr
## rownames(3464): ASV_104_sx7 ASV_10f_39x ... ASV_zpu_k9u ASV_zzb_gea
## rowData names(11): Domain Phylum ... sequence_length decontam_p_value
## colnames(7): JMF-1906-4-0007 JMF-1906-4-0008 ... JMF-1906-4-0025 JMF-1906-4-0027
## colData names(45): JMF_sample_ID User_sample_ID ... is_grassland_with_NAs Color
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
```

## Subset features


``` r
se[rowData(se)$Family %in% "Nitrososphaeraceae", ]
```

```
## class: SingleCellExperiment 
## dim: 19 26 
## metadata(4): version taxonomy_ranks filtered_features clr_pseudocount
## assays(3): counts relabundance clr
## rownames(19): ASV_1m6_gv5 ASV_1uo_lwx ... ASV_qbb_py4 ASV_qh9_8jg
## rowData names(11): Domain Phylum ... sequence_length decontam_p_value
## colnames(26): JMF-1906-4-0001 JMF-1906-4-0002 ... JMF-1906-4-0025 JMF-1906-4-0027
## colData names(45): JMF_sample_ID User_sample_ID ... is_grassland_with_NAs Color
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
```

``` r
se[!is.na(rowData(se)$decontam_p_value) & rowData(se)$decontam_p_value < 0.5, ]
```

```
## class: SingleCellExperiment 
## dim: 16 26 
## metadata(4): version taxonomy_ranks filtered_features clr_pseudocount
## assays(3): counts relabundance clr
## rownames(16): ASV_1rx_gqy ASV_2al_wxi ... ASV_s67_kdf ASV_snr_mtn
## rowData names(11): Domain Phylum ... sequence_length decontam_p_value
## colnames(26): JMF-1906-4-0001 JMF-1906-4-0002 ... JMF-1906-4-0025 JMF-1906-4-0027
## colData names(45): JMF_sample_ID User_sample_ID ... is_grassland_with_NAs Color
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
```


# DRY

https://en.wikipedia.org/wiki/Don't_repeat_yourself


``` r
# Oversimplified, there are better ways to test ASV abundance against sample variables. This is for education only!
my_analysis <- function(se, var, asv) {
  variable <-
    se[[var]] |>
    base::as.factor() |>
    as.numeric()

  abundance <- assay(se, "clr")[asv, ]

  cor.test(abundance, variable)
}

my_analysis(se, "Group", "ASV_1uo_lwx")
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  abundance and variable
## t = -2.0982, df = 24, p-value = 0.04659
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.67770755 -0.00750048
## sample estimates:
##        cor 
## -0.3937086
```


``` r
variables_of_interrest <- c("Group", "Soil_type", "Location")

map(variables_of_interrest, \(x) my_analysis(se, x, "ASV_1uo_lwx"))
```

```
## [[1]]
## 
## 	Pearson's product-moment correlation
## 
## data:  abundance and variable
## t = -2.0982, df = 24, p-value = 0.04659
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.67770755 -0.00750048
## sample estimates:
##        cor 
## -0.3937086 
## 
## 
## [[2]]
## 
## 	Pearson's product-moment correlation
## 
## data:  abundance and variable
## t = -6.1324, df = 23, p-value = 2.957e-06
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.9020769 -0.5700036
## sample estimates:
##        cor 
## -0.7877196 
## 
## 
## [[3]]
## 
## 	Pearson's product-moment correlation
## 
## data:  abundance and variable
## t = 0.012293, df = 24, p-value = 0.9903
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.3852172  0.3894826
## sample estimates:
##        cor 
## 0.00250919
```


# Further reading:

* Amplicon Analysis with Bioconductor https://microbiome.github.io/OMA/
* R for Data Science https://r4ds.hadley.nz/
