#' For a matrix of counts, impute masked entries via an alternating algorithm
#'
#' More specifically, we perform rank-1 matrix completion on the masked counts
#' optimizing the Poisson likelihood of the unmasked values.
#'
#' @param counts Matrix of counts
#' @param mask Matrix of 0/1 or TRUE/FALSE indicating which entries to mask
#' @param eps Stops when the relative change in column and row factors less than this
#' @param max_iter Maximum iterations for the alternating algorithm
#' @param verbose Whether to produce verbose output
#'
#' @return A matrix of the same dimensions as `counts` with masked entries
#'   replaced by their rank-1 Poisson imputed values.
#' @export
#' @examples
#' counts <- Matrix::Matrix(c(1, 0, 2, 3, 0, 1, 0, 4, 2),
#'                          nrow = 3, sparse = TRUE)
#' mask <- Matrix::Matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 1),
#'                        nrow = 3, sparse = TRUE)
#' impute_masked_counts(counts, mask)
#' @importFrom Matrix colSums rowSums Diagonal
impute_masked_counts <- function(
    counts,
    mask,
    eps = 1e-4,
    max_iter = 10,
    verbose = FALSE
) {
    # HACK logical indexing fails for large matrices
    counts0 <- counts - counts * mask
    counts <- counts0

    skip_row <- rowSums(counts0) == 0
    skip_col <- colSums(counts0) == 0

    cell_sizes <- rep(1, ncol(counts))

    for (i in seq_len(max_iter)) {
        guide_freqs <- sum(cell_sizes) -
            rowSums(mask %*% Diagonal(x = cell_sizes))
        guide_freqs <- rowSums(counts0) / guide_freqs
        guide_freqs[skip_row] <- 0
        guide_freqs <- guide_freqs / sum(guide_freqs)

        cell_sizes <- 1 - colSums(Diagonal(x = guide_freqs) %*% mask)
        cell_sizes <- colSums(counts0) / cell_sizes
        cell_sizes[skip_col] <- 0

        mask_imputed <- Diagonal(x = guide_freqs) %*%
            (mask %*% Diagonal(x = cell_sizes))

        counts <- counts0 + mask_imputed

        if (i > 1) {
            max_guide_diff <- max(
                abs(log(guide_freqs) - log(prev_freqs)),
                na.rm = TRUE
            )
            max_size_diff <- max(
                abs(log(cell_sizes) - log(prev_sizes)),
                na.rm = TRUE
            )

            if (verbose) {
                message(sprintf(
                    "Iter %d: max_guide_diff = %f, max_size_diff = %f",
                    i,
                    max_guide_diff,
                    max_size_diff
                ))
            }

            if (max_guide_diff < eps && max_size_diff < eps) {
                break
            }
        }

        prev_freqs <- guide_freqs
        prev_sizes <- cell_sizes
    }

    counts
}
