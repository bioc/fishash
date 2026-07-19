devtools::load_all()
set.seed(12345)

res <- simulate_guidebender2(
    n_guides = 1000,
    n_cells = 1000,
    moi = .3,
    hurdle_prob = .1,
    snr = 4,
    count_per_cell = 100,
    frac_noise_endo = .75,
    #chunk_cells=100,
    return_sparse_only = TRUE
)

saveRDS(res, "tests/testthat/fixtures/example_simulate.rds")
