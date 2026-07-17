library(fishash)

library(SingleCellExperiment)

data(tapseq_diffex)

res_fishash <- fishash(counts(altExp(tapseq_diffex)))

saveRDS(
    res_fishash,
    'tests/testthat/fixtures/fishash_schraivogel.rds'
)
