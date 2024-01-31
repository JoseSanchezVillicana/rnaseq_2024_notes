## Para hacer
## https://libd.shinyapps.io/SRP009615/

## Primero necesitamos configurar RStudio con
## shinyapps.io. Para eso necesitaremos:
# install.packages("rsconnect")

## También necesitamos verificar que tengamos todos
## los paquetes en versiones nuevas. Eso
## lo podemos hacer con:
# BiocManager::valid()

## Después necesitamos copiar y pegar la información
## de nuestra cuenta (numéro y token de acceso)

## Ahora si ya podemos continuar
options(repos = BiocManager::repositories())

library("recount3")


## URL de recount3
# options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")

## ----'quick_example'--------------------------------
## Revisemos todos los proyectos con datos de humano en recount3
# human_projects <- available_projects()
#
# ## Encuentra tu proyecto de interés. Aquí usaremos
# ## SRP009615 de ejemplo
# proj_info <- subset(
#     human_projects,
#     project == "SRP009615" & project_type == "data_sources"
# )
# dput(proj_info)

## Crear la info del proyecto de forma manual
proj_info <- structure(list(
  project = "SRP009615", organism = "human", file_source = "sra",
  project_home = "data_sources/sra", project_type = "data_sources",
  n_samples = 12L
), row.names = 1838L, class = "data.frame")

## Crea un objeto de tipo RangedSummarizedExperiment (RSE)
## con la información a nivel de genes
rse_gene_SRP009615 <- create_rse(proj_info)

## ----"tranform_counts"------------------------------
## Convirtamos las cuentas por nucleotido a cuentas por lectura
## usando compute_read_counts().
## Para otras transformaciones como RPKM y TPM, revisa transform_counts().
assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)


## ----"expand_attributes"----------------------------
## Para este estudio en específico, hagamos más fácil de usar la
## información del experimento
rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)

## Crear el sitio web interactivo
iSEE::iSEE(rse_gene_SRP009615)
