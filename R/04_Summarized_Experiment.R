'''
Summarized Experiment

Objeto Summarized Experiment
  rowRanges: genes, chr, start, id, etc...
  colData: muestras, diferentes variables
  assay(s): gene(i), muestras

                colData
  rowRanges   assay(s)


'''

## Lets build our first SummarizedExperiment object
library("SummarizedExperiment")
## ?SummarizedExperiment

## De los ejemplos en la ayuda oficial

## Creamos los datos para nuestro objeto de tipo SummarizedExperiment
## para 200 genes a lo largo de 6 muestras
nrows <- 200
ncols <- 6
## Números al azar de cuentas
set.seed(20210223)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
## Información de nuestros genes, del paquete GenomicRanges
rowRanges <- GenomicRanges::GRanges(
  rep(c("chr1", "chr2"), c(50, 150)),
  IRanges::IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
  strand = sample(c("+", "-"), 200, TRUE),
  feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))
## Información de nuestras muestras
colData <- DataFrame(
  Treatment = rep(c("ChIP", "Input"), 3),
  row.names = LETTERS[1:6]
)
## Juntamos ahora toda la información en un solo objeto de R
rse <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  rowRanges = rowRanges,
  colData = colData
)
  # Es Ranged Summarize Experiment porque tiene las coordenadas de los genes

## Exploremos el objeto resultante
rse
dim(rse)
dimnames(rse)
assayNames(rse)
assay(rse, 1)

rowRanges(rse)
rowData(rse)
rowData(rse)$feature_id[3]
colData(rse)$Treatment[3]
rse$Treatment[3]

# Ejercicio

rse[1:2, ]
  # Recupera esa parte del objeto, no los datos como tal
  # Nos está mostrando los primeros dos genes y todas las condiciones

rse[, c("A", "D", "F")]
  # Recupera esa parte del objeto, no los datos como tal
  # Todos los genes para las condiciones especificadas


# iSEE

## Explora el objeto rse de forma interactiva
library("iSEE")
iSEE::iSEE(rse)

# Ejercicio con spatialLIB

## Descarguemos unos datos de spatialLIBD
sce_layer <- spatialLIBD::fetch_data("sce_layer")
sce_layer

## Revisemos el tamaño de este objeto
lobstr::obj_size(sce_layer)

iSEE::iSEE(sce_layer)

