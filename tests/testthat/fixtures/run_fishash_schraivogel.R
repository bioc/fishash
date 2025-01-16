library(fishash)

data(crispat_schraivogel)

res_fishash <- fishash(crispat_schraivogel)

saveRDS(
  res_fishash,
  'tests/testthat/fixtures/fishash_schraivogel.rds'
)
