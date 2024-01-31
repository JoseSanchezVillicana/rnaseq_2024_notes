library(recount3)

## Revisemos todos los proyectos con datos de humano en recount3
human_projects <- available_projects()

## Encuentra tu proyecto de interés. Aquí usaremos
## SRP009615 de ejemplo
proj_info <- subset(
  human_projects,
  project == "SRP009615" & project_type == "data_sources"
)

## Crea un objeto de tipo RangedSummarizedExperiment (RSE)
## con la información a nivel de genes
rse_gene_SRP009615 <- create_rse(proj_info)

metadata(rse_gene_SRP009615)

'''
El assay de raw_count no es realmente los conteos tradicionales
Sino más bien son en base a los exones con los que hizo match
La fxn compute_read_counts transforma esos conteos a conteos tradicionales (# de lecturas)
'''

## Convirtamos las cuentas por nucleotido a cuentas por lectura
## usando compute_read_counts().
## Para otras transformaciones como RPKM y TPM, revisa transform_counts().
assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)


## Para este estudio en específico, hagamos más fácil de usar la
## información del experimento
rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)

iSEE::iSEE(rse_gene_SRP009615)
