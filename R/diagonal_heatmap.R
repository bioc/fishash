#' Plot guide count matrix as an approximately diagonal heatmap
#'
#' This function plots a heatmap after permuting the matrix columns to
#' "look diagonal" by sorting them according to their largest
#' element. The function is geared towards sparse count matrices (it
#' only draws the nonzero elements).
#'
#' @param counts Matrix of counts, cells are columns and features
#'     (guides) are rows
#' @param subsample_cells If smaller than `ncol(counts)`, subsamples
#'     the columns
#' @param rownames_size Scales the font size of the row (feature) names
#' @param plot If TRUE return a plot from ggplot2; if FALSE return the
#'     dataframe used to produce the plot.
#'
#' @returns Either a ggplot2 object or the dataframe used to produce
#'     it, depending on the `plot` argument.
#'
#' @export
#' @examples
#' data(tapseq_diffex)
#' library(SingleCellExperiment)
#' diagonal_heatmap(counts(altExp(tapseq_diffex)))
#' @importFrom methods as
#' @importFrom Matrix summary
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn scale_x_discrete scale_y_discrete theme_classic theme element_blank element_text
#' @seealso Uses [sort_columns_by_top_row()] to obtain the permuted
#'     matrix. Originally proposed by [pertpy](https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/guide_rna_assignment.html).
diagonal_heatmap <- function(
    counts,
    subsample_cells = 1000,
    rownames_size = 8,
    plot = TRUE
) {
    if (isTRUE(subsample_cells < ncol(counts))) {
        counts <- counts[, sample(seq_len(ncol(counts)), subsample_cells)]
    }

    counts <- sort_columns_by_top_row(counts)

    counts <- as(counts, 'CsparseMatrix')
    nz_long <- as.data.frame(Matrix::summary(counts))
    nz_long <- nz_long[nz_long$x > 0, ]

    rename_columns <- c(
        i = "feature",
        j = "cell",
        x = "count"
    )
    colnames(nz_long) <- rename_columns[colnames(nz_long)]

    feature_levels <- rownames(counts)
    if (is.null(feature_levels)) {
        feature_levels <- sprintf("feature%d", seq_len(nrow(counts)))
    }

    nz_long$feature <- factor(
        feature_levels[nz_long$feature],
        levels = feature_levels
    )

    cell_levels <- colnames(counts)
    if (is.null(cell_levels)) {
        cell_levels <- sprintf("cell%d", seq_len(ncol(counts)))
    }
    nz_long$cell <- factor(
        cell_levels[nz_long$cell],
        levels = cell_levels
    )

    nz_long <- nz_long[sample(seq_len(nrow(nz_long))), ]

    if (isTRUE(plot)) {
        ggplot(
            nz_long,
            aes(x = .data$cell, y = .data$feature, fill = .data$count)
        ) +
            geom_tile() +
            scale_fill_gradientn(
                colors = c('gold', 'red', 'black'),
                trans = 'log10'
            ) +
            scale_x_discrete(drop = FALSE) +
            scale_y_discrete(drop = FALSE) +
            theme_classic(base_size = 16) +
            theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = rownames_size)
            )
    } else {
        nz_long
    }
}

#' Sort matrix columns to place maxes near diagonal
#'
#' This re-orders the columns of a matrix according to the row-index
#' of its top element (breaking ties randomly). It is useful for
#' visualizing a guide count matrix as an approximately diagonal
#' heatmap.
#'
#' @param counts Matrix of counts, columns are cells and rows are
#'     features
#'
#' @returns A matrix of the same dimension as `counts`
#'
#' @export
#' @examples
#' data(tapseq_diffex)
#' library(SingleCellExperiment)
#' sort_columns_by_top_row(counts(altExp(tapseq_diffex)))
#' @importFrom nnet which.is.max
#' @seealso Used by [diagonal_heatmap()]. Originally proposed by [pertpy](https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/guide_rna_assignment.html).
sort_columns_by_top_row <- function(counts) {
    counts[, order(apply(counts, 2, which.is.max))]
}
