## code to prepare `crispat_schraivogel` dataset goes here

## git clone https://github.com/velten-group/crispat ~/src/crispat

library(Matrix)

crispat_schraivogel <- read.csv(
  '~/src/crispat/example_data/Schraivogel/gRNA_counts.csv'
)

rownames(crispat_schraivogel) <- crispat_schraivogel[,1]
crispat_schraivogel <- crispat_schraivogel[,-1]

crispat_schraivogel <- as.matrix(crispat_schraivogel)
crispat_schraivogel <- as(crispat_schraivogel, 'CsparseMatrix')

usethis::use_data(crispat_schraivogel, overwrite = TRUE)
