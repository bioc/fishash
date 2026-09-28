#' Plot histogram and weighted histogram of UMI counts
#'
#' @param counts Matrix of counts, columns are cells and rows are
#'     features.
#' @param bins Number of bins for histogram
#' @param assigned Optional matrix of assigned entries, or
#'     SummarizedExperiment returned by `fishash()`. If non-NULL, the
#'     histogram is colored by guide assignment status.
#'
#' @returns A patchwork object.
#'
#' @export
#' @examples
#' data(tapseq_diffex)
#' library(SingleCellExperiment)
#' nonzero_histogram_with_weighted(counts(altExp(tapseq_diffex)))
#' @importFrom methods as is
#' @importFrom Matrix summary
#' @importFrom rlang .data
#' @importFrom ggplot2 aes ggplot geom_histogram scale_x_log10 xlab ylab theme_bw scale_y_reverse
#' @importFrom patchwork wrap_plots plot_layout
nonzero_histogram_with_weighted <- function(counts, bins = 50, assigned = NULL) {
    counts <- as(counts, 'CsparseMatrix')
    nz <- counts > 0

    nz_long <- data.frame(
        count = counts[nz]
    )

    if (is.null(assigned)) {
        p1 <- ggplot(nz_long, aes(x = .data$count))
        p2 <- ggplot(nz_long, aes(x = .data$count, weight = .data$count))
    } else {
        if (is(assigned, "SummarizedExperiment")) {
            assigned <- assay(assigned, 'assigned')
        }

        assigned <- as(assigned, 'CsparseMatrix')
        nz_long$assigned <- assigned[nz]

        p1 <- ggplot(nz_long, aes(x = .data$count, fill = .data$assigned))
        p2 <- ggplot(nz_long, aes(x = .data$count,
                                  weight = .data$count,
                                  fill = .data$assigned))
    }

    p1 <- p1 +
        geom_histogram(bins = bins) +
        scale_x_log10() +
        xlab("UMIs per matrix entry") +
        ylab("Frequency (matrix entries)") +
        theme_bw(base_size = 16)

    p2 <- p2 +
        geom_histogram(bins = bins) +
        scale_x_log10() +
        xlab("UMIs per matrix entry") +
        ylab("Weighted frequency (UMIs)") +
        scale_y_reverse() +
        theme_bw(base_size = 16)

    wrap_plots(p1, p2, ncol = 1) + plot_layout(axes = 'collect_x', guides = 'collect')
}
