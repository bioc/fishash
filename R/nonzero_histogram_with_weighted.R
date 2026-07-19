#' Plot histogram and weighted histogram of UMI counts
#'
#' @param counts Matrix of counts, columns are cells and rows are
#'     features.
#' @param bins Number of bins for histogram
#'
#' @returns A patchwork object.
#'
#' @export
#' @examples
#' data(tapseq_diffex)
#' library(SingleCellExperiment)
#' nonzero_histogram_with_weighted(counts(altExp(tapseq_diffex)))
#' @importFrom methods as
#' @importFrom Matrix summary
#' @importFrom rlang .data
#' @importFrom ggplot2 aes ggplot geom_histogram scale_x_log10 xlab ylab theme_bw scale_y_reverse
#' @importFrom patchwork wrap_plots plot_layout
nonzero_histogram_with_weighted <- function(counts, bins = 50) {
    counts <- as(counts, 'CsparseMatrix')
    nz_long <- as.data.frame(Matrix::summary(counts))
    nz_long <- nz_long[nz_long$x > 0, ]

    rename_columns <- c(
        i = "feature",
        j = "cell",
        x = "count"
    )
    colnames(nz_long) <- rename_columns[colnames(nz_long)]

    p1 <- ggplot(nz_long, aes(x = .data$count)) +
        geom_histogram(bins = bins) +
        scale_x_log10() +
        xlab("UMIs per matrix entry") +
        ylab("Frequency (matrix entries)") +
        theme_bw(base_size = 16)

    p2 <- ggplot(nz_long, aes(x = .data$count, weight = .data$count)) +
        geom_histogram(bins = bins) +
        scale_x_log10() +
        xlab("UMIs per matrix entry") +
        ylab("Weighted frequency (UMIs)") +
        scale_y_reverse() +
        theme_bw(base_size = 16)

    wrap_plots(p1, p2, ncol = 1) + plot_layout(axes = 'collect_x')
}
