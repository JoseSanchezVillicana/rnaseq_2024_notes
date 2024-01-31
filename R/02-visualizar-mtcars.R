library("sessioninfo")
library("here")
library("ggplot2")

## Hello world
print("Soy Leo")

## Directorios
dir_plots <- here::here("figuras")
dir_rdata <- here::here("processed-data")

## Crear directorio para las figuras y archivos
dir.create(dir_plots, showWarnings = FALSE)
dir.create(dir_rdata, showWarnings = FALSE)

## Hacer una imagen de ejemplo
library(ggplot2)
pdf(file.path(dir_plots, "mtcars_gear_vs_mpg.pdf"),
    useDingbats = FALSE
)
ggplot(mtcars, aes(group = gear, y = mpg)) +
  geom_boxplot()
dev.off()

## Para reproducir mi código
options(width = 120)
sessioninfo::session_info()


# Repositorio

usethis::create_github_token()

gitcreds::gitcreds_set()

usethis::edit_r_environ()

  ## Reinicié R

usethis::edit_git_config()

usethis::use_git()

usethis::use_github()


# Material del curso

## Opción más nueva:
library("gert")
repo <- git_clone(
  "https://github.com/lcolladotor/rnaseq_LCG-UNAM_2024",
  "~/rnaseq_LCG-UNAM_2024"
)
setwd(repo)
