library(SummarizedExperiment)
library(Matrix)

test_that("Fishash results on test dataset haven't changed", {
  data(crispat_schraivogel)

  expected_result <- readRDS(test_path(
    'fixtures', 'fishash_schraivogel.rds'
  ))

  new_result <- fishash(crispat_schraivogel, exclude_empty=FALSE)

  abs_diff <- abs(
    assay(expected_result, 'log_pval') - assay(new_result, 'log_pval')
  )

  expect_equal(max(abs_diff), 0)

  expect_equal(metadata(expected_result)$log_pval_cutoff,
               metadata(new_result)$log_pval_cutoff)

  expect_identical(colData(expected_result)$demux_type,
                   colData(new_result)$demux_type)

  ## sensitive to order of comma-delimited string
  #expect_identical(colData(expected_result)$assignment,
  #                 colData(new_result)$assignment)

  expect_identical(assay(expected_result, 'assigned'),
                   assay(new_result, 'assigned'))
})
