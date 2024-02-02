## ?model.matrix
mat <- with(trees, model.matrix(log(Volume) ~ log(Height) + log(Girth)))

class(mat)
dim(mat)
dim(trees)
head(mat)

# Resumen de toda la regresión lineal
summary(lm(log(Volume) ~ log(Height) + log(Girth), data = trees))
  # Coefficients: estadísticos de cada una de las columnas del modelo
    # El estimado de log(Height)(b1): cambio promedio en el logaritmo del volumen (y)

# Corremos el modelo de regresión para los valores de expresión de todos los genes en el ensayo
  # log(Height) sería el diagnóstico (expresión de los genes dado el diagnóstico)
  # log(Girth) puede ser otra variable como la edad de los pacientes
  # log(Volume) es la expresión de los genes
# lm(log(Volume) ~ log(Height) + log(Girth), data = trees)
#     y                b0             b1

# ExploreModelMatrix

# Ejemplo 1

(sampleData <- data.frame(genotype = rep(c("A", "B"), each = 4),
                          treatment = rep(c("ctrl", "trt"), 4)))

library(ExploreModelMatrix)
vd <- VisualizeDesign(sampleData = sampleData,
                      designFormula = ~ genotype + treatment,
                      textSizeFitted = 4)
cowplot::plot_grid(plotlist = vd$plotlist)
# Nos sirve para ver la combinación de coeficientes que necesitamos sumar y restar para obtener la interpretación del genotipo que queremos

app <- ExploreModelMatrix(sampleData = sampleData,
                          designFormula = ~ genotype + treatment)
if (interactive()) shiny::runApp(app)

# Ejemplo 2

# Algunos pacientes son resistentes y otros son senibles al tratamiento
# Estudiados antes y después del tratamiento
# Reordenados dependiendo del grupo de respuesta al que pertenecen

(sampleData <- data.frame(
  Response = rep(c("Resistant", "Sensitive"), c(12, 18)),
  Patient = factor(rep(c(1:6, 8, 11:18), each = 2)),
  Treatment = factor(rep(c("pre","post"), 15)),
  ind.n = factor(rep(c(1:6, 2, 5:12), each = 2))))

vd <- VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ Response + Response:ind.n + Response:Treatment,
  textSizeFitted = 3
)
cowplot::plot_grid(plotlist = vd$plotlist, ncol = 1)

app <- ExploreModelMatrix(
  sampleData = sampleData,
  designFormula = ~ Response + Response:ind.n + Response:Treatment
)
if (interactive()) shiny::runApp(app)
# Dice que hay muchas cosas que no puede estimar porque hay combinaciones lineales

vd <- VisualizeDesign(sampleData = sampleData,
                      designFormula = ~ Treatment + Response,
                      textSizeFitted = 4)
cowplot::plot_grid(plotlist = vd$plotlist, ncol = 1)

# Dos opciones para obtener el coeficiente de ResponseResistant.Treatmentpre
# E[Y | treatment = pre ; tipo de respuesta= resistant] - E[Y | treatment = post ; tipo de respuesta= resistant]
# E[Y | treatment = pre ; tipo de respuesta= sensitive] - E[Y | treatment = post ; tipo de respuesta= sensitive]


'''
El 0 al inicio de la fórmula en el ejercicio 3:

No estimamos el valor promedio para todas las muestras
  Estimamos el valor promedio para cada batch
Le estamos diciendo a R que no queremos el intercepto
'''


library(recount3)

human_projects <- available_projects()
head(human_projects)
table(
  human_projects$file_source[human_projects$project_type == "data_sources"]
)
head(subset(human_projects, file_source == "gtex"), 30)
proj_info <- subset(
  human_projects,
  project == "BLADDER" & project_type == "data_sources"
)
proj_info
rse_gene_SRP092465 <- create_rse(proj_info)
rse_gene_SRP092465
metadata(rse_gene_SRP092465)
rowData(rse_gene_SRP092465)
colData(rse_gene_SRP092465)
