library(SummarizedExperiment)
library(Matrix)

test_that("Test simulation results haven't changed", {
    set.seed(12345)
    res <- simulate_guidebender2(
        n_guides = 1000,
        n_cells = 1000,
        moi = .3,
        hurdle_prob = .1,
        snr = 4,
        count_per_cell = 100,
        frac_noise_endo = .75,
        return_sparse_only = TRUE
    )

    res2 <- readRDS(test_path(
        'fixtures',
        'example_simulate.rds'
    ))

    expect_equal(
        assay(res, 'ground_truth'),
        assay(res2, 'ground_truth')
    )

    expect_equal(
        assay(res, 'counts'),
        assay(res2, 'counts')
    )

    expect_equal(
        assay(res, 'counts_signal'),
        assay(res2, 'counts_signal')
    )
})
