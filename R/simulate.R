#' Simulate guide count matrix based on a cellbender-like model
#'
#' @param n_guides Number of guides. Must be strictly greater than 1.
#' @param n_cells Number of cells. Must be strictly greater than 1.
#' @param moi The pre-selection multiplicity of infection. The
#'     Poisson-rate in the hurdle Poisson for the number of guides per
#'     cell.
#' @param hurdle_prob The probability a cell has zero guides; the
#'     hurdle in the hurdle-Poisson for the number of guides per cell.
#' @param guide_infection_alpha The guide infection relative
#'     frequencies (`p_g`) is Dirichlet-distributed with shape
#'     `guide_infection_alpha`.
#' @param d_sigma_guide The guide expression size factor
#'     (`d_g^{guide}`) is lognormal with mean 0 and SD `d_sigma_guide`
#' @param endo_shape_sum The endogenouse noise guide frequencies
#'     `\chi_g^a` is Dirichlet distributed; its shape parameters sum
#'     to `endo_shape_sum * n_guides`.
#' @param endo_shape_flat If this is 1, then the endogenous noise
#'     guide frequencies `\chi_g^a` are Dirichlet distributed with
#'     constant shape parameter equal to `endo_shape_sum`. If this is
#'     0, then the Dirichlet shape parameters are instead proportional
#'     to `d_g^{guide} * p_g`, and sum to `endo_shape_sum * n_guides`.
#'     If this is between 0 and 1, it interpolates between these 2
#'     options.
#' @param d_mu_drop The droplet ambient size factor `d_n^{drop}` is
#'     lognormal with location `d_mu_drop`.
#' @param d_sigma_drop The droplet ambient size factor `d_n^{drop}` is
#'     lognormal with scale `d_sigma_drop`.
#' @param d_mu_cell The cell size factor `d_n^{cell}` is lognormal
#'     with location `d_mu_cell`.
#' @param d_sigma_cell The cell size factor `d_n^{cell}` is lognormal
#'     with scale `d_sigma_cell`.
#' @param rho_alpha Exogeneous noise frequency (i.e. from PCR
#'     chimeras) is `Beta(rho_alpha, rho_beta)`.
#' @param rho_beta Exogeneous noise frequency (i.e. from PCR chimeras)
#'     is `Beta(rho_alpha, rho_beta)`.
#' @param eps_alpha Droplet-specificy capture efficiency param
#'     `\epsilon_n` is Gamma with mean 1 and shape `eps_alpha`
#' @param Phi_cell The negative binomial overdispersion of the signal
#'     counts
#' @param Phi_noise The negative binomial overdispersion of the noise
#'     counts
#' @param return_sparse_only Whether to return only the sparse
#'     matrices of counts and assignments, or whether to also return
#'     the dense matrices of Poisson rates and latent odds-ratios.
#' @param median_per_cell If non-null, the total number of counts is
#'     rescaled so that this is the median per cell.
#' @param chunk_cells Simulate in chunks, with this many cells per
#'     chunk.  This option can reduce memory usage, when used in
#'     combination with `return_sparse_only=TRUE`. Must perfectly
#'     divide `n_cells`, and be strictly greater than 1.
#'
#' @returns A SummarizedExperiment, with the following assays:
#' \describe{
#'   \item{ground_truth}{A logical matrix indicating which cells were infected
#'                       with which guides}
#'   \item{counts}{The observed guide counts}
#'   \item{counts_signal}{The counts coming from the true signal (i.e.
#'                        with noise counts subtracted).}
#' }
#' When `return_sparse_only` is `FALSE`, the following assays are also added:
#' \describe{
#'  \item{lambd_signal}{The Poisson rate of the signal counts}
#'  \item{lambd_noise}{The Poisson rate of the noise counts}
#'  \item{latent_odds_ratio}{The true odds ratio of each cell-guide pair.}
#' }
#'
#' @export
#' @examples
#' set.seed(123)
#' sim <- simulate_guidebender(n_guides = 5, n_cells = 20,
#'                             moi = 0.2, hurdle_prob = 0.1)
#' @importFrom extraDistr rtpois rdirichlet
#' @importFrom stats rbeta rbinom rgamma rmultinom rnorm rpois median
#' @importFrom SummarizedExperiment SummarizedExperiment cbind
simulate_guidebender <- function(
    n_guides, n_cells,
    moi, hurdle_prob,
    guide_infection_alpha=1,
    d_sigma_guide=0.5,
    endo_shape_sum=1,
    endo_shape_flat=0,
    d_mu_drop=log(10), d_sigma_drop=0.5,
    d_mu_cell=log(60), d_sigma_cell=0.5,
    rho_alpha=0.5, rho_beta=9.5,
    eps_alpha=50,
    Phi_cell=1,
    Phi_noise=0,
    return_sparse_only=FALSE,
    median_per_cell=NULL,
    chunk_cells=NULL
) {
    if (n_cells <= 1) {
        stop("n_cells must be strictly greater than 1.")
    }

    if (n_guides <= 1) {
        stop("n_guides must be strictly greater than 1.")
    }

    # sample number of guides per cell
    infections_per_cell <- rtpois(n_cells, moi, a=0)
    infections_per_cell[as.logical(rbinom(n_cells, 1, hurdle_prob))] <- 0

    # sample the guide frequencies
    guide_infection_freqs <- as.vector(rdirichlet(
        1, guide_infection_alpha * rep(1, n_guides)))

    # guide size factors for reads per observation
    d_g_guide <- exp(rnorm(n_guides, mean=0, sd=d_sigma_guide))

    # ambient background guide counts (\chi_g^a)
    chi_g_a <- guide_infection_freqs * d_g_guide
    chi_g_a <- chi_g_a / sum(chi_g_a)
    chi_g_a <- chi_g_a * n_guides * endo_shape_sum
    chi_g_a <- (1 - endo_shape_flat) * chi_g_a + endo_shape_flat * rep(endo_shape_sum, n_guides)
    chi_g_a <- as.vector(rdirichlet(1, chi_g_a))

    # fraction of exogenous reads
    rho_n <- rbeta(n_cells, rho_alpha, rho_beta)
    # droplet capture efficiency, mean of 1
    eps_n <- rgamma(n_cells, eps_alpha, rate=eps_alpha)

    # cell and droplet size factor
    d_n_cell <- exp(rnorm(n_cells, mean=d_mu_cell, sd=d_sigma_cell))
    d_n_drop <- exp(rnorm(n_cells, mean=d_mu_drop, sd=d_sigma_drop))

    cell_names <- paste0("cell_", seq_len(n_cells))

    if (is.null(chunk_cells)) {
        simulate_guidebender_helper(
            guide_infection_freqs=guide_infection_freqs,
            infections_per_cell=infections_per_cell,
            chi_g_a=chi_g_a,
            d_g_guide=d_g_guide,
            rho_n=rho_n,
            eps_n=eps_n,
            d_n_cell=d_n_cell,
            d_n_drop=d_n_drop,
            cell_names=cell_names,
            Phi_cell=Phi_cell,
            Phi_noise=Phi_noise,
            median_per_cell=median_per_cell,
            return_sparse_only=return_sparse_only
        )
    } else {
        if (chunk_cells <= 1) {
            stop("chunk_cells must be strictly greater than 1.")
        }

        if (n_cells %% chunk_cells != 0) {
            stop(sprintf(
                "chunk_cells must be a proper divisor of n_cells (%d %% %d == %d != 0)",
                n_cells, chunk_cells, n_cells %% chunk_cells
            ))
        }

        n_chunks <- n_cells / chunk_cells
        list_idxs <- lapply(
            seq_len(n_chunks),
            function(i) {
                seq(
                    from=(i-1) * chunk_cells + 1,
                    to=i * chunk_cells,
                    by=1
                )
            }
        )

        res <- lapply(
            list_idxs,
            function(idxs) {
                simulate_guidebender_helper(
                    guide_infection_freqs=guide_infection_freqs,
                    chi_g_a=chi_g_a,
                    d_g_guide=d_g_guide,
                    infections_per_cell=infections_per_cell[idxs],
                    rho_n=rho_n[idxs],
                    eps_n=eps_n[idxs],
                    d_n_cell=d_n_cell[idxs],
                    d_n_drop=d_n_drop[idxs],
                    cell_names=cell_names[idxs],
                    Phi_cell=Phi_cell,
                    Phi_noise=Phi_noise,
                    median_per_cell=median_per_cell,
                    return_sparse_only=return_sparse_only
                )
            }
        )

        do.call(cbind, res)
    }
}

simulate_guidebender_helper <- function(
    guide_infection_freqs,
    chi_g_a,
    d_g_guide,
    infections_per_cell,
    rho_n,
    eps_n,
    d_n_cell,
    d_n_drop,
    cell_names,
    Phi_cell,
    Phi_noise,
    median_per_cell,
    return_sparse_only
) {
    # sample the true infections
    mat_n_infections <- vapply(
        infections_per_cell,
        function(i) as.vector(rmultinom(1, i, guide_infection_freqs)),
        FUN.VALUE = integer(length(guide_infection_freqs))
    )

    stopifnot(colSums(mat_n_infections) == infections_per_cell)

    mat_truth <- mat_n_infections > 0

    # expected counts from the signal
    mu_ng_cell <- mat_n_infections
    mu_ng_cell <- sweep(mu_ng_cell, 1, d_g_guide, '*')
    mu_ng_cell <- sweep(mu_ng_cell, 2, (1-rho_n)* eps_n * d_n_cell, '*')

    # gamma-poisson with shape alpha and scale theta
    # is NB with r=alpha and p=1/(1+theta),
    # mu = r(1-p)/p = alpha * theta
    # var = mu + mu^2 / alpha
    # so, alpha=1/Phi

    # sample from gamma for overdispersion
    if (Phi_cell == 0) {
        lambd_ng_cell <- mu_ng_cell
    } else {
        lambd_ng_cell <- matrix(rgamma(
            length(mu_ng_cell),
            shape=1/Phi_cell,
            scale=Phi_cell*as.vector(mu_ng_cell)
        ), nrow=nrow(mu_ng_cell))
    }

    # expected endogenous noise
    lambd_ng_noise_endo <- chi_g_a %o% ((1-rho_n) * eps_n * d_n_drop)

    # bulk profile \bar{\chi}_g
    chi_bar_g <- rowSums(lambd_ng_cell + lambd_ng_noise_endo)
    chi_bar_g <- chi_bar_g / sum(chi_bar_g)

    # endogenous signal+noise
    lambd_ng_endo <- lambd_ng_cell + lambd_ng_noise_endo

    # expected exogenous noise
    lambd_ng_noise_exo <- sweep(
        rowSums(lambd_ng_endo) %o% colSums(lambd_ng_endo)
        / sum(lambd_ng_endo),
        2, rho_n / (1-rho_n), '*'
    )

    lambd_ng_noise <- lambd_ng_noise_endo + lambd_ng_noise_exo

    if (Phi_noise != 0) {
        lambd_ng_noise <- matrix(rgamma(
            length(lambd_ng_noise),
            shape=1/Phi_noise,
            scale=Phi_noise*as.vector(lambd_ng_noise)
        ), nrow=nrow(lambd_ng_noise))
    }

    if (!is.null(median_per_cell)) {
        orig_med <- median(colSums(lambd_ng_cell) + colSums(lambd_ng_noise))
        rescale_factor <- median_per_cell / orig_med
        lambd_ng_cell <- rescale_factor * lambd_ng_cell
        lambd_ng_noise <- rescale_factor * lambd_ng_noise
    }

    counts_cell <- matrix(rpois(
        length(lambd_ng_cell), as.vector(lambd_ng_cell)
    ), nrow=nrow(lambd_ng_cell))
    
    counts_noise <- matrix(rpois(
        length(lambd_ng_noise), as.vector(lambd_ng_noise)
    ), nrow=nrow(lambd_ng_noise))

    lambd_ng_tot <- lambd_ng_noise + lambd_ng_cell

    quo1 <- -sweep(lambd_ng_tot, 1, rowSums(lambd_ng_tot), "-")
    quo2 <- -sweep(lambd_ng_tot, 2, colSums(lambd_ng_tot), "-")

    or_latent_mat <- lambd_ng_tot / quo1 / quo2 * (
        sum(lambd_ng_tot) - quo1 - quo2 + lambd_ng_tot)

    if (return_sparse_only) {
        res <- SummarizedExperiment(
            assays=list(
                ground_truth = as(mat_truth, 'CsparseMatrix'),
                counts = as(counts_cell + counts_noise, 'CsparseMatrix'),
                counts_signal = as(counts_cell, 'CsparseMatrix')
            )
        )
    } else {
        res <- SummarizedExperiment(
            assays=list(
                ground_truth = mat_truth,
                counts = counts_cell + counts_noise,
                counts_signal = counts_cell,
                lambd_signal = lambd_ng_cell,
                lambd_noise = lambd_ng_noise,
                latent_odds_ratio = or_latent_mat
            )
        )
    }

    colnames(res) <- cell_names
    rownames(res) <- paste0("feature_", seq_len(nrow(res)))

    return(res)
}

#' Alternative interface to `simulate_guidebender` that tries to have
#' more interpretable input parameters, in particular: Signal to Noise
#' Ratio, count per cell, and fraction of noise that is endogeneous vs
#' exogeneous. These replace the following parameters from
#' `simulate_guidebender` which should not be used with this function:
#' `d_mu_drop`, `d_mu_cell`, `rho_alpha`, `rho_beta`.
#'
#' @param n_guides Number of guides
#' @param n_cells Number of cells
#' @param moi The pre-selection multiplicity of infection. The
#'     Poisson-rate in the hurdle Poisson for the number of guides per
#'     cell.
#' @param hurdle_prob The probability a cell has zero guides; the
#'     hurdle in the hurdle-Poisson for the number of guides per cell.
#' @param snr Ratio of counts from signal to counts from noise
#' @param count_per_cell Median or mean counts per cell. `use_median`
#'     determines whether this specifies the median or mean.
#' @param frac_noise_endo Fraction of the noise that is endogeneous
#'     (ie background extracellular debris) as opposed to exogeneous
#'     (ie PCR chimeras).
#' @param rho_sum Sum of `rho_alpha` and `rho_beta` from
#'     `simulate_guidebender`. Smaller values imply greater
#'     variability between cells for the fraction of their counts
#'     coming from exogenous noise (PCR chimeras).
#' @param d_sigma_guide The guide expression size factor
#'     (`d_g^{guide}`) is lognormal with mean 0 and SD `d_sigma_guide`
#' @param d_sigma_drop The droplet ambient size factor `d_n^{drop}` is
#'     lognormal with scale `d_sigma_drop`.
#' @param d_sigma_cell The cell size factor `d_n^{cell}` is lognormal
#'     with scale `d_sigma_cell`.
#' @param use_median Whether `counts_per_cell` refers to the median
#'     (TRUE) or mean (FALSE).
#' @param ... Additional arguments to pass to
#'     `simulate_guidebender`. Note the following arguments must NOT
#'     be specified: `d_mu_drop`, `d_mu_cell`, `rho_alpha`, `rho_beta`.
#'
#' @return A SummarizedExperiment as returned by `simulate_guidebender`.
#' @export
#' @examples
#' set.seed(123)
#' sim <- simulate_guidebender2(n_guides = 5, n_cells = 20,
#'                              moi = 0.2, hurdle_prob = 0.1,
#'                              snr = 5, count_per_cell = 100,
#'                              frac_noise_endo = 0.5)
simulate_guidebender2 <- function(
    n_guides, n_cells,
    moi, hurdle_prob,
    snr,
    count_per_cell,
    frac_noise_endo,
    rho_sum=10,
    d_sigma_drop=.5, d_sigma_cell=.5, d_sigma_guide=.5,
    use_median = TRUE,
    ...
) {
    avg_infections <- moi / (1-exp(-moi)) * (1 - hurdle_prob)

    frac_signal <- snr / (1 + snr)
    frac_endo <- frac_noise_endo * (1 - frac_signal)
    frac_exo <- 1 - frac_signal - frac_endo

    rho_alpha <- frac_exo * rho_sum
    rho_beta <- rho_sum - rho_alpha

    # expectation of lognormal is exp(mu + sigma^2 / 2)
    # frac_signal / frac_endo = avg_infections * exp(mu_c + sigma_c^2/2 + sigma_g^2/2) / exp(mu_d + sigma_d^2 / 2)
    d_mu_cell_drop_diff <- (
        log(frac_signal / frac_endo / avg_infections)
        - .5 * d_sigma_cell^2
        - .5 * d_sigma_guide^2
        + .5 * d_sigma_drop^2
    )

    # n = avg_infections * exp(mu_c+sigma_c^2/2+sigma_g^2/) + exp(mu_d+sigma_d^2/2)
    # n / exp(mu_d) = avg_infections * exp(diff)*exp(sigma_c^2/2+sigma_g^2/2) + exp(sigma_d^2/2)
    # exp(mu_d) = n / (avg_infections * exp(diff)*exp(sigma_c^2/2+sigma_g^2/2) + exp(sigma_d^2/2))
    d_mu_drop <- (
        log(count_per_cell) - log(
            avg_infections * exp(d_mu_cell_drop_diff + d_sigma_cell^2 / 2 + d_sigma_guide^2 / 2) +
                exp(d_sigma_drop^2)
        )
    )

    d_mu_cell <- d_mu_drop + d_mu_cell_drop_diff

    if (use_median) {
        median_per_cell <- count_per_cell
    } else {
        median_per_cell <- NULL
    }

    simulate_guidebender(
        n_guides = n_guides,
        n_cells = n_cells,
        moi = moi,
        hurdle_prob = hurdle_prob,
        d_sigma_guide = d_sigma_guide,
        d_mu_drop = d_mu_drop,
        d_sigma_drop = d_sigma_drop,
        d_mu_cell = d_mu_cell,
        d_sigma_cell = d_sigma_cell,
        rho_alpha = rho_alpha,
        rho_beta = rho_beta,
        median_per_cell = median_per_cell,
        ...
    )
}
