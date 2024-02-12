library("recount3")

    ## ----download_SRP140558--------------
## Sólo tiene dos parámetros

human_projects <- available_projects()

rse_gene <- create_rse(
  subset(
    human_projects,
    project == "SRP140558" & project_type == "data_sources"
  )
)
assay(rse_gene, "counts") <- compute_read_counts(rse_gene)
rse_gene

## ----attributes------------------------

rse_gene <- expand_sra_attributes(rse_gene)
colData(rse_gene)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene)))
]

## ----re_cast---------------------------
## Pasar de character a numeric o factor
rse_gene$sra_attribute.visit <- as.factor(rse_gene$sra_attribute.visit)
rse_gene$sra_attribute.group <- factor(tolower(rse_gene$sra_attribute.group))

colData(rse_gene)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene)))
]

## Resumen de las variables de interés
summary(as.data.frame(colData(rse_gene)[
  ,
  grepl("^sra_attribute.[group|visit]", colnames(colData(rse_gene)))
]))

# Proporción de expresiones
rse_gene$assigned_gene_prop <- rse_gene$recount_qc.gene_fc_count_all.assigned / rse_gene$recount_qc.gene_fc_count_all.total
summary(rse_gene$assigned_gene_prop)

with(colData(rse_gene), plot(assigned_gene_prop, sra_attribute.group))

with(colData(rse_gene), plot(assigned_gene_prop, sra_attribute.visit))

## Eliminemos a muestras malas
hist(rse_gene_unfiltered$assigned_gene_prop)
rse_gene_unfiltered <- rse_gene

rse_gene <- rse_gene[, rse_gene$assigned_gene_prop > 0.3]

gene_means <- rowMeans(assay(rse_gene, "counts"))
summary(gene_means)
## Eliminamos genes con bajos niveles de expresión
rse_gene <- rse_gene[gene_means > 0.1, ]

dim(rse_gene_unfiltered)
dim(rse_gene)

## Porcentaje de genes que retuvimos
round(nrow(rse_gene) / nrow(rse_gene_unfiltered) * 100, 2)
  # Retuvimos 60.87%

  # Normalización de los datos

library("edgeR") # BiocManager::install("edgeR", update = FALSE)
dge <- DGEList(
  counts = assay(rse_gene, "counts"),
  genes = rowData(rse_gene)
)
dge <- calcNormFactors(dge)

# Expresión Diferencial

library("ggplot2")
ggplot(as.data.frame(colData(rse_gene)), aes(y = assigned_gene_prop, x = sra_attribute.group)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Age Group")

  # Están casi iguales

mod <- model.matrix(~ sra_attribute.group + sra_attribute.visit + assigned_gene_prop,
                        data = colData(rse_gene)
)
colnames(mod)

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


## Genes diferencialmente expresados entre infants y niños con FDR < 5% (False Discovery Rate)
table(de_results$adj.P.Val < 0.05)
## Visualicemos los resultados estadísticos
plotMA(eb_results, coef = 2)


# Eje y tiene valores p
volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)

de_results[de_results$gene_name %in% c("MID2", "FAM153B", "IGF2BP3"), ]

    # Visualizando genes

## Extraer valores de los genes de interés
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ] # Mayor señal de expresión

## Creemos una tabla con información de las muestras
## y con nombres de columnas más amigables
df <- as.data.frame(colData(rse_gene)[, c("sra_attribute.group", "sra_attribute.visit")])
colnames(df) <- c("AgeGroup", "Visit")
df

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









