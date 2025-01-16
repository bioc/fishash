#' Cell hashing by one-sided Fisher test
#'
#' @param counts Matrix of counts. Typically, columns are cells and
#'   rows are features (e.g. gRNAs).
#' @param padj_cutoff Threshold for false discovery rate
#' @param padj_method FDR correction type. "BH" for Benjamini-Hochberg
#'   "BY" for Benjamini-Yekutueli, "GS" for Guo & Sarkar 2020 assuming
#'   independence across cells (and arbitrary dependence within).
#' @param min_count Minimum number of counts to call a feature
#'   present. Note this threshold is applied after the FDR correction.
#' @param min_frac Minimum fraction of counts within a cell to call a
#'   feature present. Note this threshold is applied after FDR
#'   correction.
#'
#' @returns A SummarizedExperiment, containing assays for the
#'   assignment matrix ("assigned"), and log-p-value ("log_pval"). The
#'   colData has columns for whether the cell is singlet, doublet, or
#'   unassigned ("demux_type"), and the assignment as a
#'   comma-delimited string ("assignment"). The metadata contains an
#'   entry for the p-value cutoff at the given FDR level
#'   ("log_pval_cutoff").
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom Matrix colSums rowSums Diagonal
#' @importFrom sparseMatrixStats colMins
fishash <- function(counts, padj_cutoff=.05, padj_method=c("GS", "BY", "BH"),
                    min_count=2, min_frac=0) {
  padj_method <- match.arg(padj_method)

  tot <- sum(counts)

  col_sums <- colSums(counts)
  row_sums <- rowSums(counts)

  n_entries <- as.numeric(nrow(counts)) * as.numeric(ncol(counts))

  counts <- as(counts, 'CsparseMatrix')
  fracs <- counts %*% Diagonal(x=1/pmax(colSums(counts), 1))

  counts <- as(counts, 'TsparseMatrix')

  df <- data.frame(
    count=counts@x,
    row_idx=counts@i+1,
    col_idx=counts@j+1
  )

  df$col_sum <- col_sums[df$col_idx]
  df$row_sum <- row_sums[df$row_idx]

  df$log_pval <- phyper(
    df$count - 1,
    df$row_sum,
    tot-df$row_sum,
    df$col_sum,
    lower.tail=FALSE,
    log.p=TRUE
  )

  mat_logpval <- sparseMatrix(
    i=df$row_idx, j=df$col_idx,
    x=df$log_pval, dims=dim(counts), dimnames=dimnames(counts),
    index1=TRUE
  )

  if (padj_method == "BY") {
    cm <- log(n_entries) + 1/(2*n_entries) + .57721
  } else {
    cm <- 1
  }

  if (padj_method %in% c("BH", "BY")) {
    df <- df[order(df$log_pval),]
    df$rank <- 1:nrow(df)

    df$padj <- pmin(n_entries * cm * exp(df$log_pval) / df$rank, 1)

    df <- df[nrow(df):1,]
    df$padj <- cummin(df$padj)

    n_signif <- sum(df$padj <= padj_cutoff)
    logpval_cutoff <- log(padj_cutoff) - log(cm) - log(n_entries) + log(n_signif)
    stopifnot(sum(df$log_pval <= logpval_cutoff) == n_signif)
  } else if (padj_method == "GS") {
    colmin_logpval <- colMins(mat_logpval)
    block_padj <- p.adjust(exp(colmin_logpval) * nrow(counts), method='BH')
    B <- sum(block_padj <= padj_cutoff)

    logpval_cutoff <- log(padj_cutoff) - log(n_entries) + log(B)
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

  SummarizedExperiment(
    assays=list(assigned=mat_assigned, log_pval=mat_logpval),
    metadata=list(log_pval_cutoff=logpval_cutoff)
  )
}
