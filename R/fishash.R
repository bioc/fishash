#' Cell hashing by one-sided Fisher test
#'
#' @param counts Matrix of counts. Typically, columns are cells and
#'     rows are features (e.g. gRNAs).
#' @param padj_cutoff Threshold for false discovery rate
#' @param padj_method FDR correction type. "BH" for Benjamini-Hochberg
#'     "BY" for Benjamini-Yekutueli, "GS" for Guo & Sarkar 2020
#'     assuming independence across cells (and arbitrary dependence
#'     within).
#' @param min_count Minimum number of counts to call a feature
#'     present.
#' @param min_frac Minimum fraction of counts within a cell to call a
#'     feature present.
#' @param refit Whether to iteratively refit the model, re-estimating
#'     the background vs signal counts from the previous iteration. If
#'     0, no refitting is done; if a positive integer, reruns with
#'     that many iterations (or until convergence, whichever is
#'     first).  The significant entries from the previous iteration
#'     are subtracted from the counts to get a matrix of
#'     "background noise" counts, which are then used for non-cell
#'     entries of the 2x2 table.  This may mitigate effects from
#'     Simpson's paradox.
#' @param exclude_empty If TRUE, rows and columns that have 0 counts
#'     are ignored for multiple testing correction.
#' @param background Optional matrix of background noise counts,
#'     e.g. estimated from demuxEM, which are used for non-cell
#'     entries of the 2x2 table.
#'
#' @returns A SummarizedExperiment. The colData has columns for
#'   whether the cell is singlet, doublet, or unassigned
#'   ("demux_type"), and the assignment as a comma-delimited string
#'   ("assignment"). The metadata contains an entry for the p-value
#'   cutoff at the given FDR level ("log_pval_cutoff").  The
#'   SummarizedExperiment contains assays for the assignment matrix
#'   ("assigned"), the log-p-value ("log_pval"), the odds-ratio
#'   ("odds_ratio"), and a regularized version of the odds-ratio that
#'   avoids dividing by 0 by adding a pseudocount to the off-diagonal
#'   entries of the 2x2 contingency table
#'   ("odds_ratio_regularized"). Note the "odds_ratio_regularized"
#'   assay is EXPERIMENTAL and subject to change.
#'
#' @export
#' @examples
#' data(tapseq_diffex)
#' library(SingleCellExperiment)
#' result <- fishash(counts(altExp(tapseq_diffex)))
#' @importFrom SummarizedExperiment SummarizedExperiment assay `assay<-`
#' @importFrom S4Vectors metadata `metadata<-`
#' @importFrom Matrix colSums rowSums Diagonal sparseMatrix drop0
#' @importFrom sparseMatrixStats colMins
#' @importFrom dplyr case_when
#' @importFrom methods as
#' @importFrom stats p.adjust phyper
fishash <- function(
    counts,
    padj_cutoff = .05,
    padj_method = c("GS", "BY", "BH"),
    min_count = 2,
    min_frac = 0,
    refit = 10,
    exclude_empty = TRUE,
    background = NULL
) {
    padj_method <- match.arg(padj_method)

    counts <- as(counts, 'CsparseMatrix')

    if (!is.null(background) && refit > 0) {
        # the refitting procedure assumes the background is derived
        # from the counts, so it's not allowed to set a separate
        # background
        stop("Non-null background with refitting not allowed")
    }

    for (i in seq_len(refit + 1)) {
        if (is.null(background) & i == 1) {
            background <- counts
        } else if (i > 1) {
            background <- impute_masked_counts(counts, mask)
        }

        res <- fishash_internal(
            counts = counts,
            padj_cutoff = padj_cutoff,
            padj_method = padj_method,
            min_count = min_count,
            min_frac = min_frac,
            background = background,
            exclude_empty = exclude_empty
        )

        if (i > 3) {
            # use | to mask entries that were assigned in previous
            # iterations but not the most recent one. Often, these are
            # borderline significant, so worth masking. This also
            # guarantees convergence and prevents "alternating"
            # behavior where a borderline entry flips between assigned
            # or not after every iteration
            mask <- mask | assay(res, 'assigned')
        } else {
            # However, on the first couple iterations, the assignments
            # are less stable, and in particular Simpson's paradox may
            # cause substantial false positives on the first
            # iteration. So, we don't use | on the first couple
            # iterations.
            mask <- assay(res, 'assigned')
        }

        # if converged (no change in assignments), end iterations
        if (i > 1 && sum(abs(prev - assay(res, 'assigned'))) == 0) {
            break
        }

        prev <- assay(res, 'assigned')
    }
    metadata(res)$num_iter <- i

    assay(res, 'assigned') <- drop0(assay(res, 'assigned'))

    res
}

fishash_internal <- function(
    counts,
    padj_cutoff,
    padj_method,
    min_count,
    min_frac,
    background,
    exclude_empty
) {
    tot <- sum(counts)

    col_sums <- colSums(counts)
    row_sums <- rowSums(counts)

    if (exclude_empty) {
        n_rows <- sum(row_sums > 0)
        n_cols <- sum(col_sums > 0)
    } else {
        n_rows <- nrow(counts)
        n_cols <- ncol(counts)
    }
    n_entries <- as.numeric(n_rows) * as.numeric(n_cols)

    counts <- as(counts, 'CsparseMatrix')
    fracs <- counts %*% Diagonal(x = 1 / pmax(colSums(counts), 1))

    counts <- as(counts, 'TsparseMatrix')

    df <- data.frame(
        count = counts@x,
        row_idx = counts@i + 1,
        col_idx = counts@j + 1
    )

    col_sums_bg <- colSums(background)
    row_sums_bg <- rowSums(background)
    tot_bg <- sum(background)

    df$col_sum <- col_sums[df$col_idx]

    #df$row_sum <- row_sums[df$row_idx]
    df$row_sum <- (row_sums_bg[df$row_idx] -
        background[cbind(df$row_idx, df$col_idx)] +
        df$count)
    df$tot <- tot_bg - col_sums_bg[df$col_idx] + col_sums[df$col_idx]

    df$log_pval <- phyper(
        df$count - 1,
        df$row_sum,
        df$tot - df$row_sum,
        df$col_sum,
        lower.tail = FALSE,
        log.p = TRUE
    )

    df$odds_ratio <- (df$count *
        (df$tot - df$row_sum - df$col_sum + df$count) /
        (df$row_sum - df$count) /
        (df$col_sum - df$count))

    df$odds_ratio_regularized <- (df$count *
        (df$tot - df$row_sum - df$col_sum + df$count) /
        (df$row_sum - df$count + 1) /
        (df$col_sum - df$count + 1))

    mat_logpval <- sparseMatrix(
        i = df$row_idx,
        j = df$col_idx,
        x = df$log_pval,
        dims = dim(counts),
        dimnames = dimnames(counts),
        index1 = TRUE
    )

    if (padj_method == "BY") {
        cm <- log(n_entries) + 1 / (2 * n_entries) + .57721
    } else {
        cm <- 1
    }

    if (padj_method %in% c("BH", "BY")) {
        df <- df[order(df$log_pval), ]
        df$rank <- seq_len(nrow(df))

        df$padj <- pmin(n_entries * cm * exp(df$log_pval) / df$rank, 1)

        df <- df[nrow(df):1, ]
        df$padj <- cummin(df$padj)

        n_signif <- sum(df$padj <= padj_cutoff)

        logpval_cutoff <- (log(padj_cutoff) -
            log(cm) -
            log(n_entries) +
            log(n_signif))

        stopifnot(sum(df$log_pval <= logpval_cutoff) == n_signif)
    } else if (padj_method == "GS") {
        colmin_logpval <- colMins(mat_logpval)

        block_padj <- p.adjust(
            pmin(exp(colmin_logpval) * n_rows, 1),
            method = 'BH'
        )

        B <- sum(block_padj <= padj_cutoff)

        logpval_cutoff <- log(padj_cutoff) - log(n_entries) + log(B)
    } else {
        stop(sprintf("Unrecognized padj_method %s", padj_method))
    }

    mat_assigned <- mat_logpval <= logpval_cutoff

    # only apply threshold when nonzero to avoid dense matrix op
    if (min_count > 0) {
        mat_assigned <- mat_assigned & (counts >= min_count)
    }

    # only apply threshold when nonzero to avoid dense matrix op
    if (min_frac > 0) {
        mat_assigned <- mat_assigned & (fracs >= min_frac)
    }

    n_assigned <- colSums(mat_assigned)
    demux_type <- case_when(
        n_assigned == 1 ~ 'singlet',
        n_assigned > 1 ~ 'doublet',
        TRUE ~ 'unknown'
    )

    df$assigned <- mat_assigned[cbind(df$row_idx, df$col_idx)]

    df_assigned <- df[df$assigned, ]
    assignment <- rep("", ncol(counts))
    if (nrow(df_assigned) > 0) {
        assign_map <- tapply(
            rownames(counts)[df_assigned$row_idx],
            df_assigned$col_idx,
            paste,
            collapse = ','
        )
        assignment[as.integer(names(assign_map))] <- assign_map
    }

    SummarizedExperiment(
        assays = list(
            assigned = mat_assigned,
            log_pval = mat_logpval,
            odds_ratio = sparseMatrix(
                i = df$row_idx,
                j = df$col_idx,
                x = df$odds_ratio,
                dims = dim(counts),
                dimnames = dimnames(counts),
                index1 = TRUE
            ),
            odds_ratio_regularized = sparseMatrix(
                i = df$row_idx,
                j = df$col_idx,
                x = df$odds_ratio_regularized,
                dims = dim(counts),
                dimnames = dimnames(counts),
                index1 = TRUE
            )
        ),
        colData = data.frame(
            demux_type = demux_type,
            assignment = assignment
        ),
        metadata = list(log_pval_cutoff = logpval_cutoff)
    )
}
