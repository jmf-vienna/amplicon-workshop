#!/usr/bin/env Rscript

if (!require("pak", quietly = TRUE)) {
  install.packages("pak")
}

pak::pkg_install(
  c(
    "tidyverse",
    "conflicted",
    "patchwork",
    "styler",
    "lintr",
    "usethis",
    "seqinr",
    "vegan",
    "ggrepel",
    "ggbeeswarm",
    "svglite",
    "rstatix",
    "ggpubr",
    "targets",
    "config",
    "waldo",
    # Bioconductor:
    "BiocVersion",
    "BiocManager",
    "dada2",
    "phyloseq",
    "ALDEx2",
    "DESeq2",
    "decontam",
    "mia",
    "miaViz",
    # GitHub:
    "david-barnett/microViz",
    "mikemc/speedyseq",
    "MadsAlbertsen/ampvis2@*release",
    "jmf-vienna/jmf-r@*release"
  ),
  upgrade = TRUE
)

if (packageVersion("BiocVersion") < "3.22") {
  BiocManager::install(version = "3.22")
}
