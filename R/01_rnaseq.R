# 30 enero 2024

# Prerquisites
## For installing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

## Install required packages
install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(
  c(
    "usethis", ## Utilities
    "here",
    "biocthis",
    "lobstr",
    "postcards",
    "sessioninfo",
    "SummarizedExperiment", ## Main containers / vis
    "iSEE",
    "edgeR", ## RNA-seq
    "ExploreModelMatrix",
    "limma",
    "recount3",
    "pheatmap", ## Visualization
    "ggplot2",
    "patchwork",
    "RColorBrewer",
    "ComplexHeatmap",
    "spatialLIBD" ## Advanced
  )
)

# Files
usethis::create_project('rnaseq_2024_notes')

usethis::use_r('01_rnaseq.R')

usethis::use_r("02-visualizar-mtcars.R")
