#' Simulate guide count matrix based on a cellbender-like model
#'
#' @param n_cells Number of cells
#' @param n_guides Number of guides
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
#' @param endo_shape_sum The guide frequencies in the
#'     "endogeneous noise" `\chi_g^a` is Dirichlet distributed; the
#'     shape parameter is proportional to `d_g^{guide} * p_g`, and
#'     sums to `endo_shape_sum * n_guides`.
#' @param d_mu_drop The droplet ambient size factor `d_n^{drop}` is
#'     lognormal with mean `d_mu_drop`.
#' @param d_sigma_drop The droplet ambient size factor `d_n^{drop}` is
#'     lognormal with SD `d_sigma_drop`.
#' @param d_mu_cell The cell size factor `d_n^{cell}` is lognormal
#'     with mean `d_mu_cell`.
#' @param d_sigma_cell The cell size factor `d_n^{cell}` is lognormal
#'     with SD `d_sigma_cell`.
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
#'
#' @returns A SummarizedExperiment, with the following assays:
#'     `"ground_truth"` a logical matrix indicating which cells were
#'     infected with which guides; `"counts"` the observed guide
#'     counts; `"counts_signal"` the counts coming from the true
#'     signal (i.e. with noise counts subtracted); `"lambd_signal"`
#'     the Poisson rate of the signal counts; `"lambd_noise"` the
#'     Poisson rate of the noise counts; `"latent_odds_ratio"` the
#'     true odds ratio of each cell-guide pair.
#' 
#' @export
#' @importFrom extraDistr rtpois rdirichlet
#' @importFrom stats rbeta rbinom rgamma rmultinom rnorm rpois
#' @importFrom SummarizedExperiment SummarizedExperiment
simulate_guidebender <- function(
    # Rows and columns of the count matrix
    n_guides, n_cells,
    # Hurdle-Poisson parameters for number of infections per cell
    moi, hurdle_prob,
    # probability of sampling over guides is dirichlet
    guide_infection_alpha=1,
    # guide size factor is lognormal(0, d_sigma^guide)
    d_sigma_guide=1,
    # endogenous noise is Dirichlet with shape proportional to the
    # guide counts in signal and summing to n_guides * endo_shape_sum
    endo_shape_sum=1,
    # droplet size factor. Roughly the expected ambient droplet counts per cell
    # d_n^drop ~ lognormal(d_mu^drop, d_sigma^drop)
    d_mu_drop=log(10), d_sigma_drop=1,
    # cell size factor. Roughly the expected "signal" counts per cell per guide
    # note it slightly differs from cellbender defn (counts per cell only).
    # d_n^cell ~ lognormal(d_mu^cell, d_sigma^cell)
    d_mu_cell=log(10), d_sigma_cell=1,
    # "exogenous" prob (index hopping, chimeras)
    # rho_n ~ Beta(rho_alpha, rho_beta)
    rho_alpha=1.5, rho_beta=50,
    # droplet-specific capture efficiency param, close to 1
    # epsilon_n ~ Gamma(shape=epsilon_alpha, mean=1)
    eps_alpha=50,
    # Overdispersion of the signal counts
    Phi_cell=1,
    # Overdispersion of the noise counts
    Phi_noise=0
) {
    # sample number of guides per cell
    infections_per_cell <- rtpois(n_cells, moi, a=0)
    infections_per_cell[as.logical(rbinom(n_cells, 1, hurdle_prob))] <- 0

    # sample the guide frequencies
    guide_infection_freqs <- as.vector(rdirichlet(
        1, guide_infection_alpha * rep(1, n_guides)))

    mat_n_infections <- sapply(infections_per_cell,
                               function(i) rmultinom(1, i, guide_infection_freqs))

    stopifnot(colSums(mat_n_infections) == infections_per_cell)

    mat_truth <- mat_n_infections > 0

    # guide size factors for reads per observation
    d_g_guide <- exp(rnorm(n_guides, mean=0, sd=d_sigma_guide))

    # ambient background guide counts (\chi_g^a)
    chi_g_a <- guide_infection_freqs * d_g_guide
    chi_g_a <- chi_g_a / sum(chi_g_a)
    chi_g_a <- as.vector(rdirichlet(
        1, chi_g_a * n_guides * endo_shape_sum))

    # fraction of exogenous reads
    rho_n <- rbeta(n_cells, rho_alpha, rho_beta)
    # droplet capture efficiency, mean of 1
    eps_n <- rgamma(n_cells, eps_alpha, rate=eps_alpha)

    # cell and droplet size factor
    d_n_cell <- exp(rnorm(n_cells, mean=d_mu_cell, sd=d_sigma_cell))
    d_n_drop <- exp(rnorm(n_cells, mean=d_mu_drop, sd=d_sigma_drop))

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

    return(SummarizedExperiment(
        assays=list(
            ground_truth = mat_truth,
            counts = counts_cell + counts_noise,
            counts_signal = counts_cell,
            lambd_signal = lambd_ng_cell,
            lambd_noise = lambd_ng_noise,
            latent_odds_ratio = or_latent_mat
        )
    ))
}
