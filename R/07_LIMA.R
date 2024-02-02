## ----download_SRP045638----------------
library("recount3")

human_projects <- available_projects()

rse_gene_SRP045638 <- create_rse(
  subset(
    human_projects,
    project == "SRP045638" & project_type == "data_sources"
  )
)
assay(rse_gene_SRP045638, "counts") <- compute_read_counts(rse_gene_SRP045638)


## ----describe_issue--------------------
rse_gene_SRP045638$sra.sample_attributes[1:3]

## ----solve_issue-----------------------
rse_gene_SRP045638$sra.sample_attributes <- gsub("dev_stage;;Fetal\\|", "", rse_gene_SRP045638$sra.sample_attributes)
rse_gene_SRP045638$sra.sample_attributes[1:3]

## ----attributes------------------------
rse_gene_SRP045638 <- expand_sra_attributes(rse_gene_SRP045638)

colData(rse_gene_SRP045638)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP045638)))
]

## ----re_cast---------------------------
## Pasar de character a numeric o factor
rse_gene_SRP045638$sra_attribute.age <- as.numeric(rse_gene_SRP045638$sra_attribute.age)
rse_gene_SRP045638$sra_attribute.disease <- factor(tolower(rse_gene_SRP045638$sra_attribute.disease))
rse_gene_SRP045638$sra_attribute.RIN <- as.numeric(rse_gene_SRP045638$sra_attribute.RIN)
rse_gene_SRP045638$sra_attribute.sex <- factor(rse_gene_SRP045638$sra_attribute.sex)

## Resumen de las variables de interés
summary(as.data.frame(colData(rse_gene_SRP045638)[
  ,
  grepl("^sra_attribute.[age|disease|RIN|sex]", colnames(colData(rse_gene_SRP045638)))
]))

## Encontraremos diferencias entre muestra prenatalas vs postnatales
rse_gene_SRP045638$prenatal <- factor(ifelse(rse_gene_SRP045638$sra_attribute.age < 0, "prenatal", "postnatal"))
table(rse_gene_SRP045638$prenatal)

## http://rna.recount.bio/docs/quality-check-fields.html
rse_gene_SRP045638$assigned_gene_prop <- rse_gene_SRP045638$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP045638$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP045638$assigned_gene_prop)

with(colData(rse_gene_SRP045638), plot(assigned_gene_prop, sra_attribute.RIN))

## Hm... veamos si hay una diferencia entre los grupos
with(colData(rse_gene_SRP045638), tapply(assigned_gene_prop, prenatal, summary))

## Guardemos nuestro objeto entero por si luego cambiamos de opinión
rse_gene_SRP045638_unfiltered <- rse_gene_SRP045638

## Eliminemos a muestras malas
hist(rse_gene_SRP045638$assigned_gene_prop)

rse_gene_SRP045638 <- rse_gene_SRP045638[, rse_gene_SRP045638$assigned_gene_prop > 0.3]

## Calculemos los niveles medios de expresión de los genes en nuestras
## muestras.
## Ojo: en un análisis real probablemente haríamos esto con los RPKMs o CPMs
## en vez de las cuentas.
## En realidad usariamos:
# edgeR::filterByExpr() https://bioconductor.org/packages/edgeR/ https://rdrr.io/bioc/edgeR/man/filterByExpr.html
# genefilter::genefilter() https://bioconductor.org/packages/genefilter/ https://rdrr.io/bioc/genefilter/man/genefilter.html
# jaffelab::expression_cutoff() http://research.libd.org/jaffelab/reference/expression_cutoff.html
#
gene_means <- rowMeans(assay(rse_gene_SRP045638, "counts"))
summary(gene_means)

## Eliminamos genes con bajos niveles de expresión
rse_gene_SRP045638 <- rse_gene_SRP045638[gene_means > 0.1, ]

## Dimensiones finales
dim(rse_gene_SRP045638)

## Porcentaje de genes que retuvimos
round(nrow(rse_gene_SRP045638) / nrow(rse_gene_SRP045638_unfiltered) * 100, 2)

    # Normalización de los datos

library("edgeR") # BiocManager::install("edgeR", update = FALSE)
dge <- DGEList(
  counts = assay(rse_gene_SRP045638, "counts"),
  genes = rowData(rse_gene_SRP045638)
)
dge <- calcNormFactors(dge)

    # Expresión Diferencial

library("ggplot2")
ggplot(as.data.frame(colData(rse_gene_SRP045638)), aes(y = assigned_gene_prop, x = prenatal)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Age Group")

mod <- model.matrix(~ prenatal + sra_attribute.RIN + sra_attribute.sex + assigned_gene_prop,
                    data = colData(rse_gene_SRP045638)
)
colnames(mod)
    # Un valor de -2 en prenatal sería un gen que tendría mayor expresión en prenatal que en postnatal
    # prenatalprenatal -> nombre variable, no referencia

# Podemos usar limma para realizar el análisis de expresión diferencial como tal

library("limma")
vGene <- voom(dge, mod, plot = TRUE)

eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef = 2, #'prenatalprenatal'
  number = nrow(rse_gene_SRP045638),
  sort.by = "none"
)
dim(de_results)
head(de_results)
  # logfoldchange logFC es el que nos habla de niveles de expresión

## Genes diferencialmente expresados entre pre y post natal con FDR < 5% (False Discovery Rate)
table(de_results$adj.P.Val < 0.05)
## Visualicemos los resultados estadísticos
plotMA(eb_results, coef = 2)


# Eje y tiene valores p
volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)

de_results[de_results$gene_name %in% c("ZSCAN2", "VASH2", "KIAA0922"), ]


    # Visualizando genes

## Extraer valores de los genes de interés
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ] # Mayor señal de expresión

## Creemos una tabla con información de las muestras
## y con nombres de columnas más amigables
df <- as.data.frame(colData(rse_gene_SRP045638)[, c("prenatal", "sra_attribute.RIN", "sra_attribute.sex")])
colnames(df) <- c("AgeGroup", "RIN", "Sex")

## Hagamos un heatmap
library("pheatmap")
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = df
)

## Para colores
library("RColorBrewer")

## Conviertiendo los grupos de edad a colores
col.group <- df$AgeGroup
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")

col.group <- as.character(col.group)

## MDS por grupos de edad
plotMDS(vGene$E, labels = df$AgeGroup, col = col.group)

pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE, # Esto sólo pone nombres de los ENS0
  show_colnames = FALSE,
  annotation_col = df
)

## Tenemos que usar gene_id y gene_name
rowRanges(rse_gene_SRP045638)

## Guardemos los IDs de nuestros 50 genes
nombres_originales <- rownames(exprs_heatmap)

## Con match() podemos encontrar cual es cual
rownames(exprs_heatmap) <- rowRanges(rse_gene_SRP045638)$gene_name[
  match(rownames(exprs_heatmap), rowRanges(rse_gene_SRP045638)$gene_id)
]

## Vean que tambien podriamos haber usado rank()
identical(
  which(rank(de_results$adj.P.Val) <= 50),
  match(nombres_originales, rowRanges(rse_gene_SRP045638)$gene_id)
)

## Por último podemos cambiar el valor de show_rownames de FALSE a TRUE
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE, # Ahora sí ya porque cambiamos los IDs por los nombres en el df que ocupa
  show_colnames = FALSE,
  annotation_col = df
)


# Ejercicio de Revisión
'''
¿Debemos explorar las relaciones entre nuestras variables con información de nuestras muestras previo a hacer un análisis de expresión diferencial?
¿Por qué usamos el paquete edgeR?
¿Por qué es importante el argumento sort.by en topTable()?
¿Por qué es importante el argumento coef en topTable()?
'''

speaqeasy_data <- file.path(tempdir(), "rse_speaqeasy.RData")
download.file("https://github.com/LieberInstitute/SPEAQeasy-example/blob/master/rse_speaqeasy.RData?raw=true", speaqeasy_data, mode = "wb")
library("SummarizedExperiment")
load(speaqeasy_data, verbose = TRUE)
rse_gene
# 60609 genes y 40 muestras

## Exploremos la variable de PrimaryDx
table(rse_gene$PrimaryDx)

## Eliminemos el diagnosis "Other" porque no tiene información
rse_gene$PrimaryDx <- droplevels(rse_gene$PrimaryDx)
table(rse_gene$PrimaryDx)

## Exploremos numéricamente diferencias entre grupos de diagnosis para
## varias variables
with(colData(rse_gene), tapply(totalAssignedGene, PrimaryDx, summary))

with(colData(rse_gene), tapply(mitoRate, PrimaryDx, summary))

## Podemos hacer lo mismo para otras variables
with(colData(rse_gene), tapply(mitoRate, BrainRegion, summary))

## Podemos resolver la primeras preguntas con iSEE
if (interactive()) iSEE::iSEE(rse_gene)
# Column data plot 1:
# Y -> totalAssigned Gene ; X -> PrimaryDx

## O hacer graficas nosotros mismos. Aquí les muestro una posible respuesta
## con ggplot2
library("ggplot2")
ggplot(
  as.data.frame(colData(rse_gene)),
  aes(y = totalAssignedGene, group = PrimaryDx, x = PrimaryDx)
) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  xlab("Diagnosis")

ggplot(
  as.data.frame(colData(rse_gene)),
  aes(y = totalAssignedGene, group = paste0(PrimaryDx, "", BrainRegion), x = paste0(PrimaryDx, "", BrainRegion))
) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  xlab("Diagnosis")

ggplot(
  as.data.frame(colData(rse_gene)),
  aes(y = mitoRate, group = PrimaryDx, x = PrimaryDx)
) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  xlab("Diagnosis")

# No hay diferencia significativa (visualmente)

## Encontremos el gene SNAP25
rowRanges(rse_gene)

## En este objeto los nombres de los genes vienen en la variable "Symbol"
i <- which(rowRanges(rse_gene)$Symbol == "SNAP25")
i

## Para graficar con ggplot2, hagamos un pequeño data.frame
df <- data.frame(
  expression = assay(rse_gene)[i, ],
  Dx = rse_gene$PrimaryDx
)

## Ya teniendo el pequeño data.frame, podemos hacer la gráfica
ggplot(df, aes(y = log2(expression + 0.5), group = Dx, x = Dx)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  xlab("Diagnosis") +
  ylab("SNAP25: log2(x + 0.5)")

## https://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html#3_Visualizing_expression_values
scater::plotExpression(
  as(rse_gene, "SingleCellExperiment"),
  features = rownames(rse_gene)[i],
  x = "PrimaryDx",
  exprs_values = "counts",
  colour_by = "BrainRegion",
  xlab = "Diagnosis"
)

if (requireNamespace("plotly", quietly = TRUE)) {
  ## Lo pueden instalar con
  # install.packages("plotly")

  ## Guardemos el resultado de plotExpression()
  p <- scater::plotExpression(
    as(rse_gene, "SingleCellExperiment"),
    features = rownames(rse_gene)[i],
    x = "PrimaryDx",
    exprs_values = "counts",
    colour_by = "BrainRegion",
    xlab = "Diagnosis"
  )
  ## scater::plotExpression() regresa un objeto de clase ggplot
  class(p)

  ## así que podemos usar plotly para crear una versión
  ## interactiva
  plotly::ggplotly(p)
}


## Para el modelo estadístico exploremos la información de las muestras
colnames(colData(rse_gene))

## Podemos usar región del cerebro porque tenemos suficientes datos
table(rse_gene$BrainRegion)

## Pero no podemos usar "Race" porque son solo de 1 tipo
table(rse_gene$Race)

## Ojo! Acá es importante que hayamos usado droplevels(rse_gene$PrimaryDx)
## si no, vamos a tener un modelo que no sea full rank
mod <- with(
  colData(rse_gene),
  model.matrix(~ PrimaryDx + totalAssignedGene + mitoRate + rRNA_rate + BrainRegion + Sex + AgeDeath)
)

## Exploremos el modelo de forma interactiva
if (interactive()) {
  ## Tenemos que eliminar columnas que tienen NAs.
  info_no_NAs <- colData(rse_gene)[, c(
    "PrimaryDx", "totalAssignedGene", "rRNA_rate", "BrainRegion", "Sex",
    "AgeDeath", "mitoRate", "Race"
  )]
  ExploreModelMatrix::ExploreModelMatrix(
    info_no_NAs,
    ~ PrimaryDx + totalAssignedGene + mitoRate + rRNA_rate + BrainRegion + Sex + AgeDeath
  )

  ## Veamos un modelo más sencillo sin las variables numéricas (continuas) porque
  ## ExploreModelMatrix nos las muestra como si fueran factors (categoricas)
  ## en vez de continuas
  ExploreModelMatrix::ExploreModelMatrix(
    info_no_NAs,
    ~ PrimaryDx + BrainRegion + Sex
  )

  ## Si agregamos + Race nos da errores porque Race solo tiene 1 opción
  # ExploreModelMatrix::ExploreModelMatrix(
  #     info_no_NAs,
  #     ~ PrimaryDx + BrainRegion + Sex + Race
  # )
}
