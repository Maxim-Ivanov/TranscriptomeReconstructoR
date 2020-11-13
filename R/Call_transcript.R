#' Call transcript and gene models from corrected full-length RNA-seq reads
#'
#' @param long_reads \code{GRangesList} object.
#' @param skip_minor_tx Numeric in the range (0, 1), or NULL.
#' @param max_overlap_called Numeric in the range [0, 1).
#' @param min_read_width Positive integer.
#' @param min_overlap_fusion Numeric in the range (0, 1].
#' @param clust_threshold Numeric in the range (0, 1].
#' @section Details:
#' The input \code{GRangesList} object is returned by the \code{detect_alignment_errors()} function.\cr
#' Long reads either marked as truncated by \code{extend_long_reads_to_TSS_and_PAS()}, or containing a misaligned exon (as revealed by \code{detect_alignment_errors()}), are skipped from the transcript calling procedure.
#' The remaining long reads are collapsed into transcripts. The transcripts are classified into high confidence (HC), medium confidence (MC) and low confidence (LC) groups:
#' \itemize{
#'     \item HC transcripts are called from reads which start in a TSS and end in a PAS;
#'     \item MC transcripts are called from TSS-only or PAS-only reads which do not overlap with any HC transcript by more than \code{max_overlap_called} fraction of either read or transcript length;
#'     \item LC transcripts are called from reads which neither start in a TSS nor end in a PAS, and do not overlap with any HC or MC transcript by mode than \code{max_overlap_called}.
#' }
#' This iterative procedure of transcript calling ensures that highly expressed HC loci are not contaminated with less reliable MC or LC transcripts.
#' The MC/LC transcripts are not guaranteed to be full-length. To decrease the risk of picking up products of partial RNA degradation, MC and LC transcripts
#' can be called only from reads longer than \code{min_read_length} bp.\cr
#' The called HC, MC and LC transcripts are clustered into HC, MC and LC genes, respectively. A pair of transcripts of the same type having overlap (intersect/union) above the \code{clust_threshold}
#' are considered belonging to the same gene.\cr
#' Within each gene, the minor transcripts (collectively representing up to \code{skip_minor_tx} fraction of the reads) are skipped from further consideration. To suppress this behavior, set \code{skip_minor_tx = NULL}.\cr
#' Finally, transcripts which overlap at least two other disjoint transcripts by at least \code{min_overlap_fusion} fraction of their lengths, are considered fusion transcripts.
#' @return List of length 5:
#' \enumerate{
#'     \item \code{GRanges} object (HC, MC and LC genes);
#'     \item \code{GRangesList} object (HC, MC and LC transcripts);
#'     \item \code{GRanges} object (fusion genes);
#'     \item \code{GRangesList} object (fusion transcripts);
#'     \item \code{GRangesList} object (unused long reads outside of the called genes and transcripts).
#' }
#' @export
call_transcripts_and_genes <- function(long_reads, skip_minor_tx = 0.01, max_overlap_called = 0.1, min_read_width = 1000, min_overlap_fusion = 0.5, clust_threshold = 0.8) {
  stopifnot(BiocGenerics::grepl("GRangesList", class(long_reads)))
  stopifnot(is.numeric(skip_minor_tx) && length(skip_minor_tx) == 1 && skip_minor_tx > 0 && skip_minor_tx < 1)
  stopifnot(is.numeric(max_overlap_called) && length(max_overlap_called) == 1 && max_overlap_called >= 0 && max_overlap_called < 1)
  stopifnot(is.numeric(min_read_width) && length(min_read_width) == 1 && min_read_width %% 1 == 0 && min_read_width > 0)
  stopifnot(is.numeric(min_overlap_fusion) && length(min_overlap_fusion) == 1 && min_overlap_fusion > 0 && min_overlap_fusion <= 1)
  stopifnot(is.numeric(clust_threshold) && length(clust_threshold) == 1 && clust_threshold > 0 && clust_threshold <= 1)
  all_exons <- BiocGenerics::unlist(long_reads, use.names = FALSE)
  # Skip truncated reads, as well as reads containing alignment errors:
  mc <- S4Vectors::mcols(all_exons)
  bad_id <- mc$read_id[!mc$complete] %>% BiocGenerics::unique()
  message(length(bad_id), " truncated reads;")
  if (!is.null(mc$aln_error)) {
    bad_id_2 <- mc$read_id[mc$aln_error] %>% BiocGenerics::unique()
    message(length(bad_id_2), " reads with alignment errors;")
    bad_id <- c(bad_id, bad_id_2) %>% BiocGenerics::unique()
  }
  message(length(bad_id), " reads skipped from transcript calling;")
  to_skip <- mc$read_id %in% bad_id
  exons_unused <- all_exons[to_skip]
  reads_unused <- S4Vectors::split(exons_unused, S4Vectors::mcols(exons_unused)$read_id)
  good_exons <- all_exons[!to_skip]
  # Call HC transcripts (from TSS+PAS reads):
  lgl <- S4Vectors::mcols(good_exons)$over_tc == "both"
  res_hc <- good_exons[lgl] %>% call_tx_and_find_fusion(prefix = "HC", min_overlap_fusion = min_overlap_fusion, clust_threshold = clust_threshold, skip = TRUE)
  hc_tx <- res_hc[[1]]
  fusion_tx <- res_hc[[2]]
  #good_exons <- good_exons[!lgl] %>% update_free_reads(c(hc_tx, fusion_tx), max_overlap_called = max_overlap_called, min_read_width = min_read_width, clust_threshold = clust_threshold)
  good_exons <- good_exons[!lgl] %>% update_free_reads(hc_tx, max_overlap_called = max_overlap_called, min_read_width = min_read_width, clust_threshold = clust_threshold)
  # Call MC transcripts (from TSS-only or PAS-only reads):
  lgl <- S4Vectors::mcols(good_exons)$over_tc != "no"
  res_mc <- good_exons[lgl] %>% call_tx_and_find_fusion(prefix = "MC", known_tx = hc_tx, known_fusion = fusion_tx, min_overlap_fusion = min_overlap_fusion, clust_threshold = clust_threshold, skip = TRUE)
  hc_mc_tx <- res_mc[[1]]
  fusion_tx <- res_mc[[2]]
  #good_exons <- good_exons[!lgl] %>% update_free_reads(c(hc_mc_tx, fusion_tx), max_overlap_called = max_overlap_called, min_read_width = min_read_width, clust_threshold = clust_threshold)
  good_exons <- good_exons[!lgl] %>% update_free_reads(hc_mc_tx, max_overlap_called = max_overlap_called, min_read_width = min_read_width, clust_threshold = clust_threshold)
  # Call LC transcripts (from noTSS-noPAS reads):
  res_lc <- good_exons %>% call_tx_and_find_fusion(prefix = "LC", known_tx = hc_mc_tx, known_fusion = fusion_tx, min_overlap_fusion = min_overlap_fusion, clust_threshold = clust_threshold, skip = TRUE)
  hml_tx <- res_lc[[1]]
  fusion_tx <- res_lc[[2]]
  #exons_free <- exons_unused %>% update_free_reads(c(hml_tx, fusion_tx), max_overlap_called = max_overlap_called, min_read_width = min_read_width, clust_threshold = clust_threshold)
  exons_free <- exons_unused %>% update_free_reads(hml_tx, max_overlap_called = max_overlap_called, min_read_width = min_read_width, clust_threshold = clust_threshold)
  reads_free <- S4Vectors::split(exons_free, S4Vectors::mcols(exons_free)$read_id)
  # Split hml_tx into HC, MC and LC transcripts:
  tx_spl <- split_tx_by_type(hml_tx)
  hc_tx <- tx_spl[[1]]
  mc_tx <- tx_spl[[2]]
  lc_tx <- tx_spl[[3]]
  # Call HC, MC and LC genes:
  hc_genes <- call_genes(hc_tx, clust_threshold = clust_threshold)
  mc_genes <- call_genes(mc_tx, clust_threshold = clust_threshold) # observe that some MC genes may behave as HC or LC (i.e. appear as TSS+PAS or noTSS-noPAS)
  lc_genes <- call_genes(lc_tx, clust_threshold = clust_threshold)
  # Skip minor transcripts within each gene:
  if (!is.null(skip_minor_tx)) {
    message("Skipping minor transcripts...")
    out_hc <- skip_minor_transcripts(hc_genes, hc_tx, q = skip_minor_tx)
    hc_tx <- out_hc[[1]]
    message("\t", length(out_hc[[2]]), " minor HC transcripts skipped;")
    out_mc <- skip_minor_transcripts(mc_genes, mc_tx, q = skip_minor_tx)
    mc_tx <- out_mc[[1]]
    message("\t", length(out_mc[[2]]), " minor MC transcripts skipped;")
    out_lc <- skip_minor_transcripts(lc_genes, lc_tx, q = skip_minor_tx)
    lc_tx <- out_lc[[1]]
    message("\t", length(out_lc[[2]]), " minor LC transcripts skipped;")
  }
  # Recalculate genes (to update revmap) and generate unique gene IDs:
  hml_genes <- c(hc_tx, mc_tx, lc_tx) %>% sort_grl() %>% call_HC_MC_LC_genes(clust_threshold = clust_threshold, generate_ids = TRUE)
  # Generate unique transcript IDs:
  hc_tx <- hml_genes[S4Vectors::mcols(hml_genes)$type == "HC"] %>% generate_tx_id(hc_tx, .)
  mc_tx <- hml_genes[S4Vectors::mcols(hml_genes)$type == "MC"] %>% generate_tx_id(mc_tx, .)
  lc_tx <- hml_genes[S4Vectors::mcols(hml_genes)$type == "LC"] %>% generate_tx_id(lc_tx, .)
  # Combine HC, MC and LC transcripts:
  message("Final stats:")
  message("\t", length(hc_tx), " HC transcripts in ", length(hc_genes), " HC genes;")
  message("\t", length(mc_tx), " MC transcripts in ", length(mc_genes), " MC genes;")
  message("\t", length(lc_tx), " LC transcripts in ", length(lc_genes), " LC genes;")
  hml_tx <- c(hc_tx, mc_tx, lc_tx) %>% sort_grl()
  # Find more fusion reads among unused reads by overlap with called HC, MC and LC genes:
  fusion_idx <- reads_unused %>% range() %>% BiocGenerics::unlist() %>% find_fusion_tx(hml_genes, skip = FALSE)
  reads_fusion <- reads_unused %>% `[`(fusion_idx) %>% BiocGenerics::unlist(use.names = FALSE) %>% `[`(S4Vectors::mcols(.)$complete) %>% S4Vectors::split(S4Vectors::mcols(.)$read_id) %>% sort_grl()
  # Convert them to new fusion transcripts:
  reads_fusion_dedup <- deduplicate_grl(reads_fusion)
  reads_fusion_unl <- BiocGenerics::unlist(reads_fusion_dedup, use.names = FALSE)
  S4Vectors::mcols(reads_fusion_unl)$type <- "Fusion"
  reads_fusion_dedup <- BiocGenerics::relist(reads_fusion_unl, reads_fusion_dedup)
  fusion_tx <- c(fusion_tx, reads_fusion_dedup) %>% unname() %>% sort_grl()
  fusion_tx_unl <- unlist(fusion_tx, use.names = FALSE)
  S4Vectors::mcols(fusion_tx_unl)$grp <- rep(1:length(fusion_tx), times = S4Vectors::elementNROWS(fusion_tx))
  fusion_tx <- BiocGenerics::relist(fusion_tx_unl, fusion_tx)
  # Group fusion transcripts into "fusion TUs":
  fusion_tu <- fusion_tx %>% call_tu(clust_threshold = clust_threshold) %>% skip_nested()
  message("\t", length(fusion_tx), " fusion transcripts in ", length(fusion_tu), " fusion genes;")
  S4Vectors::mcols(fusion_tu)$type <- "Fusion"
  fusion_tu_id <- sprintf(paste0("%0", nchar(length(fusion_tu)), "i"), 1:length(fusion_tu))
  S4Vectors::mcols(fusion_tu)$name <- BiocGenerics::paste(S4Vectors::mcols(fusion_tu)$type, "TU", fusion_tu_id, sep = "_")
  # Generate unique IDs for fusion transcripts:
  fusion_tx <- generate_tx_id(fusion_tx, fusion_tu)
  S4Vectors::mcols(hml_genes)$revmap <- NULL
  S4Vectors::mcols(fusion_tu)$revmap <- NULL
  out <- list("hml_genes" = hml_genes, "hml_tx" = hml_tx, "fusion_genes" = fusion_tu, "fusion_tx" = fusion_tx, "reads_free" = reads_free)
  return(out)
}

# --------------------------------------------------------------------------------------------------------------------

call_transcripts <- function(all_exons, prefix) {
  all_exons <- all_exons[BiocGenerics::order(S4Vectors::mcols(all_exons)$exon_id)] # in fact, they are assumed to be already ordered...
  reads <- S4Vectors::split(all_exons, S4Vectors::mcols(all_exons)$read_id)
  tx <- deduplicate_grl(reads)
  tx_unl <- BiocGenerics::unlist(tx, use.names = FALSE)
  S4Vectors::mcols(tx_unl)$type <- prefix
  tx <- BiocGenerics::relist(tx_unl, tx)
  return(tx)
}

#------------------------------------------------------------------------------------------------------------------

find_free_reads <- function(exons, genes, max_overlap_called, min_read_width) {
  exons <- exons[BiocGenerics::order(S4Vectors::mcols(exons)$exon_id)]
  reads <- S4Vectors::split(exons, S4Vectors::mcols(exons)$read_id)
  gr <- reads %>% range() %>% BiocGenerics::unlist()
  over <- gr %over% genes
  hits <- GenomicRanges::findOverlaps(gr, genes)
  par1 <- gr[S4Vectors::queryHits(hits)]
  par2 <- genes[S4Vectors::subjectHits(hits)]
  overlap <- GenomicRanges::pintersect(par1, par2) %>% BiocGenerics::width()
  best <- BiocGenerics::tapply(overlap, S4Vectors::queryHits(hits), find_uniq_max, simplify = FALSE) %>% BiocGenerics::unlist() %>% unname()
  best_overlap <- overlap[best]
  hits <- hits[best]
  par1 <- par1[best]
  par2 <- par2[best]
  weak <- best_overlap / BiocGenerics::width(par1) <= max_overlap_called & best_overlap / BiocGenerics::width(par2) <= max_overlap_called & BiocGenerics::width(par1) >= min_read_width
  idx <- c(BiocGenerics::which(!over), S4Vectors::queryHits(hits)[weak]) %>% BiocGenerics::sort()
  message("\t", length(idx), " reads remain free;")
  out <- reads[idx] %>% BiocGenerics::unlist(use.names = FALSE)
  return(out)
}

#------------------------------------------------------------------------------------------------------------------

find_fusion_tx <- function(gr1, gr2 = NULL, min_overlap_fusion = 0.5, skip = TRUE) {
  if (is.null(gr2)) {
    gr2 <- gr1
    hits <- GenomicRanges::findOverlaps(gr1, gr2)
    hits <- hits[S4Vectors::queryHits(hits) != S4Vectors::subjectHits(hits)]
  } else {
    hits <- GenomicRanges::findOverlaps(gr1, gr2)
  }
  if (length(hits) == 0) {
    return(integer(0))
  }
  # Skip weak overlaps:
  par1 <- gr1[S4Vectors::queryHits(hits)] # fusion candidates
  par2 <- gr2[S4Vectors::subjectHits(hits)]
  valid <- BiocGenerics::width(GenomicRanges::pintersect(par1, par2)) / BiocGenerics::width(par2) >= min_overlap_fusion
  hits <- hits[valid]
  if (length(hits) == 0) {
    return(integer(0))
  }
  # Check that the fusion candidate covers at least 2 non-overlapping TUs:
  grl <- gr2[S4Vectors::subjectHits(hits)] %>% S4Vectors::split(S4Vectors::queryHits(hits))
  group <- rep(1:length(grl), times = S4Vectors::elementNROWS(grl))
  good <- grl %>% GenomicRanges::reduce() %>% S4Vectors::elementNROWS() %>% `>=`(2) %>% `[`(group)
  hits <- hits[good]
  if (length(hits) == 0) {
    return(integer(0))
  }
  # Skip fusion candidates which are transcribed higher than the overlapped TUs:
  if (isTRUE(skip)) {
    stopifnot(!is.null(BiocGenerics::score(gr1)) && !is.null(BiocGenerics::score(gr2)))
    par1 <- gr1[S4Vectors::queryHits(hits)]
    par2 <- gr2[S4Vectors::subjectHits(hits)]
    grl <- S4Vectors::split(par2, S4Vectors::queryHits(hits))
    group <- rep(1:length(grl), times = S4Vectors::elementNROWS(grl))
    score_sum <- BiocGenerics::score(par2) %>% S4Vectors::split(S4Vectors::queryHits(hits)) %>% BiocGenerics::lapply(sum) %>% BiocGenerics::unlist() %>% unname() %>% `[`(group)
    bad <- BiocGenerics::score(par1) > score_sum
    hits <- hits[!bad]
    if (length(hits) == 0) {
      return(integer(0))
    }
  }
  fusion_idx <- hits %>% S4Vectors::queryHits() %>% BiocGenerics::unique()
  return(fusion_idx)
}

#------------------------------------------------------------------------------------------------------------------

call_tu <- function(tx, clust_threshold) {
  gr <- tx %>% range() %>% BiocGenerics::unlist()
  hits <- GenomicRanges::findOverlaps(gr, gr)
  hits <- hits[S4Vectors::queryHits(hits) < S4Vectors::subjectHits(hits)]
  par1 <- gr[S4Vectors::queryHits(hits)]
  par2 <- gr[S4Vectors::subjectHits(hits)]
  valid <- BiocGenerics::width(GenomicRanges::pintersect(par1, par2)) / BiocGenerics::width(GenomicRanges::punion(par1, par2)) >= clust_threshold # nested transcripts remain unmerged
  hits <- hits[valid]
  used <- c(S4Vectors::queryHits(hits), S4Vectors::subjectHits(hits)) %>% BiocGenerics::unique()
  out_1 <- gr[-used] %>% unname()
  if (length(out_1) > 0) {
    S4Vectors::mcols(out_1)$revmap <- seq(1, length(gr))[-used] %>% S4Vectors::split(1:length(.)) %>% unname() %>% as("IntegerList")
  }
  # Merge tx with at least one strong overlap:
  mat <- hits %>% as.matrix() %>% unname()
  clust <- cluster_indexes(mat)
  group <- rep(1:length(clust), times = S4Vectors::elementNROWS(clust))
  out_2 <- gr[BiocGenerics::unlist(clust)] %>% S4Vectors::split(group) %>% GenomicRanges::reduce() %>% BiocGenerics::unlist() %>% unname()
  if (length(out_2) > 0) {
    S4Vectors::mcols(out_2)$revmap <- clust %>% as("IntegerList")
  }
  out <- c(out_1, out_2) %>% BiocGenerics::sort()
  # Calculate TU support (summary score of the merged transcripts):
  if (!is.null(BiocGenerics::score(BiocGenerics::unlist(tx)))) {
    tx_score <- BiocGenerics::unlist(tx) %>% BiocGenerics::score() %>% BiocGenerics::tapply(rep(1:length(tx), times = S4Vectors::elementNROWS(tx)), BiocGenerics::unique) %>% BiocGenerics::unlist() %>% unname()
    revmap <- S4Vectors::mcols(out)$revmap
    BiocGenerics::score(out) <- tx_score[BiocGenerics::unlist(revmap)] %>% S4Vectors::split(rep(1:length(revmap), times = S4Vectors::elementNROWS(revmap))) %>% BiocGenerics::lapply(sum) %>% BiocGenerics::unlist() %>% unname()
  }
  return(out)
}


# -----------------------------------------------------------------------------------------------------------------

skip_nested <- function(gr) {
  hits <- GenomicRanges::findOverlaps(gr, gr, type = "within") %>% as.matrix()
  c1 <- hits[, 1]
  c2 <- hits[, 2]
  good <- c1 <= c2
  mat <- BiocGenerics::cbind(ifelse(good, c1, c2), ifelse(good, c2, c1))
  mat <- mat[mat[, 1] < mat[, 2], , drop = FALSE]
  if (nrow(mat) == 0) {
    return(gr)
  }
  used <- c(mat[, 1], mat[, 2]) %>% BiocGenerics::unique()
  out <- gr[-used]
  if (length(out) > 0) {
    if (!is.null(S4Vectors::mcols(gr)$revmap)) {
      S4Vectors::mcols(out) <- BiocGenerics::subset(S4Vectors::mcols(out), select = c("revmap", "score"))
    } else {
      S4Vectors::mcols(out) <- BiocGenerics::subset(S4Vectors::mcols(out), select = "score")
    }
  } else {
    S4Vectors::mcols(out) <- NULL
  }
  if (nrow(mat) > 0) {
    clust <- cluster_indexes(mat)
    group <- rep(1:length(clust), times = S4Vectors::elementNROWS(clust))
    gr_2 <- gr[BiocGenerics::unlist(clust)]
    out_2 <- gr_2 %>% S4Vectors::split(group) %>% GenomicRanges::reduce() %>% BiocGenerics::unlist() %>% unname()
    if (!is.null(S4Vectors::mcols(gr)$revmap)) {
      rm_2 <- S4Vectors::mcols(gr_2)$revmap %>% BiocGenerics::tapply(group, function(x) { BiocGenerics::Reduce(c, x) }) %>% as("IntegerList") %>% unname()
      S4Vectors::mcols(out_2)$revmap <- rm_2
    }
    sc_2 <- BiocGenerics::score(gr_2) %>% BiocGenerics::tapply(group, sum) %>% BiocGenerics::unlist() %>% unname()
    BiocGenerics::score(out_2) <- sc_2
    out <- c(out, out_2)
  }
  return(out)
}

# ------------------------------------------------------------------------------------------------------------------

generate_tx_id <- function(tx, genes) {
  stopifnot(!is.null(S4Vectors::mcols(genes)$revmap))
  stopifnot(tx %>% BiocGenerics::unlist() %>% S4Vectors::mcols() %>% .$grp %>% is.null() %>% `!`)
  revmap <- S4Vectors::mcols(genes)$revmap
  tx <- BiocGenerics::unlist(tx, use.names = FALSE) %>% `[`(BiocGenerics::order(S4Vectors::mcols(.)$grp)) %>% S4Vectors::split(S4Vectors::mcols(.)$grp) %>% `[`(BiocGenerics::unlist(revmap))
  max_num <- S4Vectors::elementNROWS(tx) %>% max() %>% nchar()
  tx_num <- BiocGenerics::lapply(S4Vectors::elementNROWS(revmap), function(x) { sprintf(paste0("%0", max_num, "i"), seq(1, x)) }) %>% BiocGenerics::unlist()
  names(tx) <- S4Vectors::mcols(genes)$name %>% rep(times = S4Vectors::elementNROWS(revmap)) %>% BiocGenerics::paste("tx", tx_num, sep = "_")
  return(tx)
}


# ------------------------------------------------------------------------------------------------------------------

# <skip> fusion candidates which are transcribed higher than the overlapped TUs
call_tx_and_find_fusion <- function(input_exons, prefix, known_tx = GenomicRanges::GRangesList(), known_fusion = GenomicRanges::GRangesList(), min_overlap_fusion, clust_threshold, skip) {
  message("Calling ", prefix, " transcripts...")
  new_tx <- call_transcripts(input_exons, prefix = prefix) %>% BiocGenerics::sort() %>% sort_grl()
  message("\t", length(new_tx), " ", prefix, " transcripts called;"); flush.console()
  # Combine known and new transcripts:
  all_tx <- c(known_tx, new_tx) %>% sort_grl()
  # Call all TUs:
  all_tu <- all_tx %>% call_tu(clust_threshold = clust_threshold)
  # Find fusion transcripts on the combined data:
  idx1 <- find_fusion_tx(all_tu, min_overlap_fusion = min_overlap_fusion, skip = skip)
  idx2 <- S4Vectors::mcols(all_tu)$revmap %>% `[`(idx1) %>% BiocGenerics::unlist()
  message("\t", length(idx2), " fusion transcripts detected;"); flush.console()
  # Append novel fusion transcripts to the existing ones:
  if (length(idx2) > 0) {
    new_fusion <- all_tx[idx2]
    new_unl <- BiocGenerics::unlist(new_fusion, use.names = FALSE)
    S4Vectors::mcols(new_unl)$type <- "Fusion"
    new_fusion <- BiocGenerics::relist(new_unl, new_fusion)
    all_tx <- all_tx[-idx2]
    known_fusion <- c(known_fusion, new_fusion)
  }
  out <- list("all_tx" = all_tx, "known_fusion" = known_fusion)
  return(out)
}

# ------------------------------------------------------------------------------------------------------------------

update_unused_reads <- function(all_exons, used_id) {
  unused_exons <- all_exons[!S4Vectors::mcols(all_exons)$read_id %in% used_id]
  message("\t", S4Vectors::mcols(unused_exons)$read_id %>% BiocGenerics::unique() %>% length(), " reads remain unused;"); flush.console()
  return(unused_exons)
}

# --------------------------------------------------------------------------------------------------------------------

update_free_reads <- function(exons, known_tx, max_overlap_called, min_read_width, clust_threshold) {
  known_genes <- known_tx %>% call_tu(clust_threshold = clust_threshold) %>% skip_nested()
  out <- find_free_reads(exons, known_genes, max_overlap_called = max_overlap_called, min_read_width = min_read_width)
  return(out)
}

# -------------------------------------------------------------------------------------------------------------------

skip_minor_transcripts <- function(genes, tx, q) {
  stopifnot(!is.null(BiocGenerics::score(BiocGenerics::unlist(tx))))
  revmap <- S4Vectors::mcols(genes)$revmap
  revmap_unl <- BiocGenerics::unlist(revmap)
  grp1 <- rep(1:length(revmap), times = S4Vectors::elementNROWS(revmap))
  grp2 <- rep(1:length(tx), times = S4Vectors::elementNROWS(tx))
  score_tx <- tx %>% BiocGenerics::unlist() %>% BiocGenerics::score() %>% BiocGenerics::tapply(grp2, BiocGenerics::unique) %>% BiocGenerics::unlist() %>% unname() %>% `[`(revmap_unl)
  good <- score_tx %>% S4Vectors::split(grp1) %>% BiocGenerics::lapply(quantile_trim, q = q, value = "lgl") %>% BiocGenerics::unlist() %>% unname()
  skip_idx <- revmap_unl[!good]
  #message(length(skip_idx), " transcripts were skipped due to low expression;")
  if (length(skip_idx) > 0) {
    out <- list("good_tx" = tx[-skip_idx], "bad_tx" = tx[skip_idx])
  } else {
    out <- list("good_tx" = tx, "bad_tx" = GenomicRanges::GRangesList())
  }
  return(out)
}

# ------------------------------------------------------------------------------------------------------------

call_genes <- function(tx, clust_threshold, type = NULL) {
  genes <- tx %>% call_tu(clust_threshold = clust_threshold) %>% skip_nested()
  if (!is.null(type) && class(type) == "character") {
    S4Vectors::mcols(genes)$type <- type
  }
  return(genes)
}

# ------------------------------------------------------------------------------------------------------------

split_tx_by_type <- function(hml_tx) {
  hml_tx_unl <- BiocGenerics::unlist(hml_tx, use.names = FALSE)
  S4Vectors::mcols(hml_tx_unl)$grp <- rep(1:length(hml_tx), times = S4Vectors::elementNROWS(hml_tx))
  hc_tx <- hml_tx_unl %>% `[`(S4Vectors::mcols(.)$type == "HC") %>% S4Vectors::split(S4Vectors::mcols(.)$grp)
  mc_tx <- hml_tx_unl %>% `[`(S4Vectors::mcols(.)$type == "MC") %>% S4Vectors::split(S4Vectors::mcols(.)$grp)
  lc_tx <- hml_tx_unl %>% `[`(S4Vectors::mcols(.)$type == "LC") %>% S4Vectors::split(S4Vectors::mcols(.)$grp)
  out <- list("hc_tx" = hc_tx, "mc_tx" = mc_tx, "lc_tx" = lc_tx)
  return(out)
}

# -----------------------------------------------------------------------------------------------------------

call_HC_MC_LC_genes <- function(hml_tx, clust_threshold, generate_ids = FALSE) {
  # Split hml_tx into HC, MC and LC transcripts:
  tx_spl <- split_tx_by_type(hml_tx)
  # Call the genes:
  hc_genes <- tx_spl[[1]] %>% call_genes(clust_threshold = clust_threshold, type = "HC")
  mc_genes <- tx_spl[[2]] %>% call_genes(clust_threshold = clust_threshold, type = "MC")
  lc_genes <- tx_spl[[3]] %>% call_genes(clust_threshold = clust_threshold, type = "LC")
  hml_genes <- c(hc_genes, mc_genes, lc_genes) %>% BiocGenerics::sort()
  if (isTRUE(generate_ids)) {
    # Generate unique gene IDs:
    gene_id <- sprintf(paste0("%0", nchar(length(hml_genes)), "i"), 1:length(hml_genes))
    S4Vectors::mcols(hml_genes)$name <- BiocGenerics::paste(S4Vectors::mcols(hml_genes)$type, "gene", gene_id, sep = "_")
  }
  return(hml_genes)
}

