#' Cell hashing by one-sided Fisher test
#'
#' @param counts Matrix of counts. Typically, columns are cells and
#'   rows are features (e.g. gRNAs).
#' @param fdr_cutoff Threshold for false discovery rate
#' @param fdr_method FDR correction type -- either "BH" for
#'   Benjamini-Hochberg or "BY" for Benjamini-Yekutueli
#'
#' @returns A SummarizedExperiment, containing assays for the
#'   assignment ("assigned"), log-p-value ("log_pval"), 1 minus
#'   BH-corrected FDR ("padj_bh_1m"), 1 minus BY-corrected FDR
#'   ("padj_by_1m"). The colData has columns for whether the cell is
#'   singlet, doublet, or unassigned ("demux_type"), and the
#'   assignment ("assignment"). The metadata contains an entry for the
#'   p-value cutoff at the given FDR level ("log_pval_cutoff"), as
#'   well as cutoffs under both BH and BY ("log_pval_cutoff_bh",
#'   "log_pval_cutoff_by").
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom Matrix colSums rowSums
fishash <- function(counts, fdr_cutoff=.05, fdr_method=c("Holm", "BY", "BH")) {
  fdr_method <- match.arg(fdr_method)

  tot <- sum(counts)

  col_sums <- colSums(counts)
  row_sums <- rowSums(counts)

  n_entries <- as.numeric(nrow(counts)) * as.numeric(ncol(counts))

  cm <- log(n_entries) + 1/(2*n_entries) + .57721

  counts <- as(counts, 'CsparseMatrix')
  counts <- as(counts, 'TsparseMatrix')

  df <- data.frame(
    count=counts@x,
    row_idx=counts@i+1,
    col_idx=counts@j+1
  )

  df$col_sum <- col_sums[df$col_idx]
  df$row_sum <- col_sums[df$row_idx]

  df$log_pval <- phyper(
    df$count - 1,
    df$row_sum,
    tot-df$row_sum,
    df$col_sum,
    lower.tail=FALSE,
    log.p=TRUE
  )

  df <- df[order(df$log_pval),]
  df$rank <- 1:nrow(df)
  df$padj_by <- pmin(n_entries * cm * exp(df$log_pval) / df$rank, 1)
  df$padj_bh <- pmin(n_entries * 1 * exp(df$log_pval) / df$rank, 1)
  df$padj_holm <- exp(df$log_pval) * (n_entries - df$rank + 1)
  df$padj_holm <- pmin(cummax(df$padj_holm), 1)

  df <- df[nrow(df):1,]
  df$padj_by <- cummin(df$padj_by)
  df$padj_bh <- cummin(df$padj_bh)

  df$padj_bh_1m <- 1-df$padj_bh
  df$padj_by_1m <- 1-df$padj_by

  # TODO add assigned, odds_ratio
  assay_names <- c('log_pval', 'padj_by_1m', 'padj_bh_1m')
  names(assay_names) <- assay_names

  ret <- lapply(
    assay_names,
    function(x) {
      sparseMatrix(
        i=df$row_idx, j=df$col_idx,
        x=df[,x], dims=dim(counts), dimnames=dimnames(counts),
        index1=TRUE
      )
    }
  )

  ret <- SummarizedExperiment(ret)

  metadata(ret)$log_pval_cutoff_by <- max(df[df$padj_by <= fdr_cutoff,'log_pval'])
  metadata(ret)$log_pval_cutoff_bh <- max(df[df$padj_bh <= fdr_cutoff,'log_pval'])
  metadata(ret)$log_pval_cutoff_holm <- max(df[df$padj_holm <= fdr_cutoff,'log_pval'])

  if (fdr_method == 'Holm') {
    metadata(ret)$log_pval_cutoff <- metadata(ret)$log_pval_cutoff_holm
  } else if (fdr_method == 'BH') {
    metadata(ret)$log_pval_cutoff <- metadata(ret)$log_pval_cutoff_bh
  } else if (fdr_method == 'BY') {
    metadata(ret)$log_pval_cutoff <- metadata(ret)$log_pval_cutoff_by
  } else {
    stop(sprintf("Unrecognized FDR method %s", fdr_method))
  }

  # TODO add colData

  ret
}
