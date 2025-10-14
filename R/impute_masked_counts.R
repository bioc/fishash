#' For a matrix of counts, impute masked entries via an alternating algorithm
#'
#' @param counts Matrix of counts
#' @param mask Matrix of 0/1 or TRUE/FALSE indicating which entries to mask
#' @param eps Stops when the relative change in colSums and rowSums less than this
#' @param max_iter Maximum iterations for the alternating algorithm
#' @param verbose Whether to produce verbose output
#'
#' @export
#' @importFrom Matrix colSums rowSums Diagonal
impute_masked_counts <- function(counts, mask,
                                 eps=1e-4, max_iter=10, verbose=FALSE) {
    # HACK logical indexing fails for large matrices
    counts0 <- counts - counts * mask
    counts <- counts0

    for (i in 1:max_iter) {
        guide_freqs <- rowSums(counts)
        guide_freqs <- guide_freqs / sum(guide_freqs)

        cell_sizes <- colSums(counts)
        if (i == 1) {
            sz_fac_correction <- Diagonal(x=guide_freqs) %*% mask
            sz_fac_correction <- 1 - colSums(sz_fac_correction)
            sz_fac_correction[sz_fac_correction == 0] <- 1

            cell_sizes <- cell_sizes / sz_fac_correction
        }

        mask_imputed <- Diagonal(x=guide_freqs) %*% (
            mask %*% Diagonal(x=cell_sizes))

        counts <- counts0 + mask_imputed

        if (i > 1) {
            max_guide_diff <- max(abs(log(guide_freqs) - log(prev_freqs)),
                                  na.rm=T)
            max_size_diff <- max(abs(log(cell_sizes) - log(prev_sizes)),
                                 na.rm=T)

            if (verbose) {
                message(sprintf(
                    "Iter %d: max_guide_diff = %f, max_size_diff = %f",
                    i, max_guide_diff, max_size_diff
                ))
            }

            if (max_guide_diff < eps & max_size_diff < eps) {
                break
            }
        }

        prev_freqs <- guide_freqs
        prev_sizes <- cell_sizes
    }

    counts
}
