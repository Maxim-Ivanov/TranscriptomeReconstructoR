#' Call TSS or PAS sites as clusters of reads from 5'- or 3'-tag sequencing data
#'
#' @param data \code{GRangesList} object with all elements having the "score" mcols.
#' @param min_support Integer >= 1 (but not exceeding the number of replicates in \code{data})
#' @param min_tpm Non-negative numeric.
#' @param max_gap Non-negative integer.
#' @param q_trim Numeric in the range (0, 1), or NULL.
#' @param min_score Non-negative numeric.
#' @section Details:
#' The input \code{data} is returned by the \code{load_BAM_files()} function.\cr
#' Only genomic positions with non-zero signal in at least \code{min_support} replicates are considered for TSS/PAS calling.\cr
#' After merging the replicates, only genomic positions with at least \code{min_tpm} signal are considered as positions with valid signal.\cr
#' Genomic positions with valid signal separated by not more than \code{max_gap} bases are merged into tag cluster (TC) candidates.\cr
#' TC candidates are trimmed to retain only genomic positions containing \code{1 - q_trim} fraction of the signal.\cr
#' Only TCs with TPM (averaged among replicates) exceeding \code{min_score} are considered valid TSS/PAS.
#' @return \code{GRanges} object containing coordinates of called tag clusters.
#' @export
call_TCs <- function(data, min_support = length(data), min_tpm = 0.05, max_gap = 50, q_trim = 0.05, min_score = 0.25) {
  stopifnot(BiocGenerics::grepl("GRangesList", class(data)))
  stopifnot(BiocGenerics::lapply(data, function(x) { !is.null(BiocGenerics::score(x)) }) %>% BiocGenerics::unlist() %>% all())
  stopifnot(is.numeric(min_support) && length(min_support) == 1 && min_support %% 1 == 0 && min_support > 0 && min_support <= length(data))
  stopifnot(is.numeric(min_tpm) && length(min_tpm) == 1 && min_tpm >= 0)
  stopifnot(is.numeric(max_gap) && length(max_gap) == 1 && max_gap %% 1 == 0 && max_gap >= 0)
  stopifnot(is.null(q_trim) || (is.numeric(q_trim) && length(q_trim) == 1 && q_trim > 0 && q_trim < 1))
  stopifnot(is.numeric(min_score) && length(min_score) == 1 && min_score >= 0)
  message("Merging replicates..."); flush.console()
  cov_fw_full <- BiocGenerics::lapply(data, function(x) { GenomicRanges::coverage(x[BiocGenerics::strand(x) == "+"], weight = "score") }) %>% BiocGenerics::Reduce(`+`, .)
  cov_rev_full <- BiocGenerics::lapply(data, function(x) { GenomicRanges::coverage(x[BiocGenerics::strand(x) == "-"], weight = "score") }) %>% BiocGenerics::Reduce(`+`, .)
  cov_fw <- cov_fw_full
  cov_rev <- cov_rev_full
  message("Filtering by support (", min_support, " samples)..."); flush.console()
  if (min_support > 1) {
    sup_fw <- BiocGenerics::lapply(data, function(x) { GenomicRanges::coverage(x[BiocGenerics::strand(x) == "+"]) }) %>% BiocGenerics::Reduce(`+`, .)
    cov_fw[sup_fw > 0 & sup_fw < min_support] <- 0
    sup_rev <- BiocGenerics::lapply(data, function(x) { GenomicRanges::coverage(x[BiocGenerics::strand(x) == "-"]) }) %>% BiocGenerics::Reduce(`+`, .)
    cov_rev[sup_rev > 0 & sup_rev < min_support] <- 0
  }
  total_tags <- BiocGenerics::lapply(data, function(x) { sum(BiocGenerics::score(x) * BiocGenerics::width(x)) }) %>% BiocGenerics::unlist() %>% sum()
  min_tags <- ceiling(min_tpm * total_tags / 1e06)
  message("Applying static cutoff (", min_tpm, " TPM = ", min_tags, " tags)...")
  cov_fw[cov_fw < min_tags] <- 0
  cov_rev[cov_rev < min_tags] <- 0
  # Convert the merged coverage to GRanges:
  gr_fw <- GenomicRanges::bindAsGRanges(score = cov_fw)
  BiocGenerics::strand(gr_fw) <- "+"
  gr_rev <- GenomicRanges::bindAsGRanges(score = cov_rev)
  BiocGenerics::strand(gr_rev) <- "-"
  gr <- c(gr_fw, gr_rev) %>% BiocGenerics::sort()
  gr <- gr[BiocGenerics::score(gr) > 0]
  message("Calling TCs..."); flush.console()
  tc <- GenomicRanges::reduce(gr, min.gapwidth = max_gap)
  names(tc) <- 1:length(tc)
  tc_fw <- tc[BiocGenerics::strand(tc) == "+"]
  tc_rev <- tc[BiocGenerics::strand(tc) == "-"]
  if (!is.null(q_trim)) {
    # Quantile trimming:
    message("Trimming TCs to ", round((1 - q_trim) * 100), "% tags:"); flush.console()
    message("\tForward strand..."); flush.console()
    rel_fw <- BiocGenerics::lapply(cov_fw_full[tc_fw], quantile_trim, q = q_trim)
    seqnames_fw_rep <- rep(GenomeInfoDb::seqnames(tc_fw), times = S4Vectors::elementNROWS(rel_fw))
    start_fw_rep <- rep(BiocGenerics::start(tc_fw), times = S4Vectors::elementNROWS(rel_fw))
    tc_fw_new <- GenomicRanges::GRanges(seqnames = seqnames_fw_rep, IRanges::IRanges(start = start_fw_rep + BiocGenerics::unlist(rel_fw) - 1, width = 1), strand = "+", seqinfo = GenomeInfoDb::seqinfo(tc))
    message("\tReverse strand..."); flush.console()
    rel_rev <- BiocGenerics::lapply(cov_rev_full[tc_rev], quantile_trim, q = q_trim)
    seqnames_rev_rep <- rep(GenomeInfoDb::seqnames(tc_rev), times = S4Vectors::elementNROWS(rel_rev))
    start_rev_rep <- rep(BiocGenerics::start(tc_rev), times = S4Vectors::elementNROWS(rel_rev))
    tc_rev_new <- GenomicRanges::GRanges(seqnames = seqnames_rev_rep, IRanges::IRanges(start = start_rev_rep + BiocGenerics::unlist(rel_rev) - 1, width = 1), strand = "-", seqinfo = GenomeInfoDb::seqinfo(tc))
    tc <- c(tc_fw_new, tc_rev_new) %>% BiocGenerics::sort() %>% GenomicRanges::reduce(min.gapwidth = max_gap)
    names(tc) <- 1:length(tc)
    tc_fw <- tc[BiocGenerics::strand(tc) == "+"]
    tc_rev <- tc[BiocGenerics::strand(tc) == "-"]
  }
  message("Counting overlapping tags..."); flush.console()
  results <- vector("list", length(data))
  for (i in seq_along(data)) {
    message("\t", names(data)[[i]]); flush.console()
    curr <- data[[i]]
    curr_fw <- curr[BiocGenerics::strand(curr) == "+"]
    curr_rev <- curr[BiocGenerics::strand(curr) == "-"]
    tags_fw <- GenomicRanges::coverage(curr_fw, weight = "score") %>% `[`(tc_fw) %>% as("NumericList") %>% BiocGenerics::lapply(sum) %>% BiocGenerics::unlist()
    names(tags_fw) <- names(tc_fw)
    tags_rev <- GenomicRanges::coverage(curr_rev, weight = "score") %>% `[`(tc_rev) %>% as("NumericList") %>% BiocGenerics::lapply(sum) %>% BiocGenerics::unlist()
    names(tags_rev) <- names(tc_rev)
    tags <- c(tags_fw, tags_rev)
    tags <- tags[BiocGenerics::order(as.numeric(names(tags)))]
    results[[i]] <- tags
  }
  mat <- BiocGenerics::do.call(BiocGenerics::cbind, results)
  colnames(mat) <- names(data)
  # Calculate TC score (= average TPM):
  message("Calculating TC scores..."); flush.console()
  avg_tpm <- DelayedArray::apply(mat, 2, function(col) { col / sum(col) * 1e06 }) %>% Matrix::rowMeans() %>% round(5) %>% unname()
  # Find TC summits (= position with the highest expression):
  message("Finding TC summits on forward strand..."); flush.console()
  rel_pos_fw <- cov_fw_full[tc_fw] %>% BiocGenerics::lapply(BiocGenerics::which.max) %>% BiocGenerics::unlist()
  names(rel_pos_fw) <- names(tc_fw)
  message("Finding TC summits on reverse strand..."); flush.console()
  rel_pos_rev <- BiocGenerics::lapply(cov_rev_full[tc_rev], BiocGenerics::which.max) %>% BiocGenerics::unlist()
  names(rel_pos_rev) <- names(tc_rev)
  message("Constructing the output GRanges..."); flush.console()
  rel_pos <- c(rel_pos_fw, rel_pos_rev)
  rel_pos <- rel_pos[BiocGenerics::order(as.numeric(names(rel_pos)))]
  abs_pos <- BiocGenerics::start(tc) + rel_pos - 1
  summit <- GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(tc), IRanges::IRanges(start = abs_pos, width = 1), strand = BiocGenerics::strand(tc), seqinfo = GenomeInfoDb::seqinfo(tc))
  S4Vectors::mcols(tc) <- BiocGenerics::cbind(data.frame("score" = avg_tpm), BiocGenerics::as.data.frame(mat))
  S4Vectors::mcols(tc)$thick <- summit
  names(tc) <- NULL
  # Filter by min_score:
  if (min_score > 0) {
    tc <- tc[BiocGenerics::score(tc) >= min_score]
  }
  message(length(tc), " TCs returned;"); flush.console()
  return(tc)
}

# ----------------------------------------------------------------------------------------------------------------

quantile_trim <- function(x, q, value = "idx") {
  stopifnot(value %in% c("idx", "lgl"))
  idx <- BiocGenerics::order(x, decreasing = TRUE)
  x2 <- x[idx]
  last_good <- BiocGenerics::which(cumsum(x2) / sum(x2) >= (1 - q))[[1]]
  if (last_good == length(x)) {
    out <- 1:length(x)
  } else {
    out <- BiocGenerics::which(x >= x2[last_good])
  }
  if (value == "idx") {
    return(out)
  } else {
    lgl <- rep(FALSE, times = length(x))
    lgl[out] <- TRUE
    return(lgl)
  }
}
