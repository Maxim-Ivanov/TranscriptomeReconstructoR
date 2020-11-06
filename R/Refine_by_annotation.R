#' Refine called transcripts by an existing transcript model
#'
#' @param hml_tx,annot_exons,fusion_tx \code{GRangesList} objects.
#' @param tss,pas \code{GRanges} objects.
#' @param max_exon_diff Positive integer.
#' @param tx_flanks_up,tx_flanks_down Integer vectors of length 2.
#' @param min_score_2 Non-negative numeric.
#' @param min_tx_cov Numeric in the range (0, 1].
#' @param clust_threshold Numeric in the range (0, 1].
#' @param min_overlap_fusion Numeric in the range (0, 1].
#' @section Details:
#' \code{hml_tx} is the called transcript model (the second element in the list returned by the \code{call_transcripts_and_genes()} function).
#' \code{fusion_tx} is the called set of fusion transcripts (the fourth element in the list returned by the \code{call_transcripts_and_genes()} function).
#' \code{annot_exons} is a known transcript model (returned by e.g. \code{exonsBy(txdb, by = "tx")}, where \code{txdb} is a \code{GenomicFeatures} object).\cr
#' The function aims to adjust the called transcripts by the annotated transcripts:
#' \itemize{
#'     \item 5'- and 3'-borders of called exons are adjusted to the most similar border of an annotated exon (by not more than \code{max_exon_diff} bp);
#'     \item Annotated transcripts are classified into valid and non-valid. A valid known transcript must overlap with called TSS and PAS (both having scores above \code{min_score_2})
#'     within \code{tx_flanks_up} and \code{tx_flanks_down} bp windows around its start and end, respectively;
#'     \item 5'- and/or 3'-borders of called MC and LC transcripts lacking overlap with TSS and/or PAS are adjusted to the borders of the most similar mate
#'     among the valid annotated transcripts (given that at least \code{min_tx_cov} fraction of the called transcript is covered by the annotated mate);
#'     \item Valid annotated transcripts which do not overlap with any called transcript are copied from the annotation to the called HC transcript set;
#' }
#' In addition, the set of fusion transcripts is updated by finding called transcripts which overlap at least two valid annotated transcripts
#' (or an annotated and a called transcript) by at least \code{min_overlap_fusion} fraction of their lengths.
#' @return List of length 4:
#' \enumerate{
#'     \item \code{GRanges} object (updated HC, MC and LC genes);
#'     \item \code{GRangesList} object (refined HC, MC and LC transcripts);
#'     \item \code{GRanges} object (updated fusion genes);
#'     \item \code{GRangesList} object (updated fusion transcripts);
#' }
#' @export
refine_transcripts_by_annotation <- function(hml_tx, annot_exons, tss, pas, fusion_tx = GenomicRanges::GRangesList(), max_exon_diff = 10,
                                             tx_flanks_up = c(-100, 100), tx_flanks_down = c(-100, 100), min_score_2 = 5,
                                             min_tx_cov = 0.95, clust_threshold = 0.8, min_overlap_fusion = 0.5) {
  stopifnot(BiocGenerics::grepl("GRangesList", class(hml_tx)))
  stopifnot(BiocGenerics::grepl("GRangesList", class(annot_exons)))
  stopifnot(BiocGenerics::grepl("GRangesList", class(fusion_tx)))
  stopifnot(is.numeric(max_exon_diff) && length(max_exon_diff) == 1 && max_exon_diff %% 1 == 0 && max_exon_diff > 0)
  stopifnot(is.numeric(tx_flanks_up) && length(tx_flanks_up) == 2 && tx_flanks_up[[1]] %% 1 == 0 && tx_flanks_up[[2]] %% 1 == 0)
  stopifnot(is.numeric(tx_flanks_down) && length(tx_flanks_down) == 2 && tx_flanks_down[[1]] %% 1 == 0 && tx_flanks_down[[2]] %% 1 == 0)
  stopifnot(is.numeric(min_score_2) && length(min_score_2) == 1 && min_score_2 >= 0)
  stopifnot(is.numeric(min_tx_cov) && length(min_tx_cov) == 1 && min_tx_cov > 0 && min_tx_cov <= 1)
  stopifnot(is.numeric(clust_threshold) && length(clust_threshold) == 1 && clust_threshold > 0 && clust_threshold <= 1)
  stopifnot(is.numeric(min_overlap_fusion) && length(min_overlap_fusion) == 1 && min_overlap_fusion > 0 && min_overlap_fusion <= 1)
  # Adjust called exon borders by the annotated exons:
  message("Adjusting exon borders...")
  hml_tx <- adjust_exons_of_called_transcripts(hml_tx, annot_exons, max_exon_diff = max_exon_diff)
  # Adjust borders of annotated transcripts to the called TSS and PAS:
  annot_exons <- extend_long_reads_to_TSS_and_PAS(annot_exons, tss, pas, read_flanks_up = tx_flanks_up, read_flanks_down = tx_flanks_down, mode = "transcripts", verbose = FALSE)
  # Extend called transcripts along annotated ones:
  message("Extending transcript borders along annotated transcripts...")
  hml_tx <- extend_called_transcripts(hml_tx, annot_exons, tss, mode = "start", min_score_2 = min_score_2, min_tx_cov = min_tx_cov)
  hml_tx <- extend_called_transcripts(hml_tx, annot_exons, pas, mode = "end", min_score_2 = min_score_2, min_tx_cov = min_tx_cov)
  # Update HC/MC/LC classification:
  hml_range <- hml_tx %>% range() %>% BiocGenerics::unlist(use.names = FALSE)
  start_in_tss <- GenomicRanges::resize(hml_range, 1, "start") %over% tss
  end_in_pas <- GenomicRanges::resize(hml_range, 1, "end") %over% pas
  type <- ifelse(start_in_tss & end_in_pas, "HC", ifelse(start_in_tss | end_in_pas, "MC", "LC"))
  hml_unl <- BiocGenerics::unlist(hml_tx, use.names = FALSE)
  S4Vectors::mcols(hml_unl)$type <- rep(type, times = S4Vectors::elementNROWS(hml_tx))
  hml_tx <- BiocGenerics::relist(hml_unl, hml_tx)
  # Find annotated transcripts with both TSS and PAS:
  annot_range <- annot_exons %>% range() %>% BiocGenerics::unlist(use.names = FALSE)
  annot_start <- GenomicRanges::resize(annot_range, 1, "start")
  annot_end <- GenomicRanges::resize(annot_range, 1, "end")
  valid_tss <- tss[BiocGenerics::score(tss) >= min_score_2]
  valid_pas <- pas[BiocGenerics::score(tss) >= min_score_2]
  annot_hc <- annot_exons[annot_start %over% valid_tss & annot_end %over% valid_pas]
  # Try to detect more fusion transcripts:
  hml_tu <- hml_tx %>% call_tu(clust_threshold = clust_threshold)
  BiocGenerics::score(hml_tu) <- NULL
  S4Vectors::mcols(hml_tu)$source <- "hml"
  annot_tu <- annot_hc %>% call_tu(clust_threshold = clust_threshold)
  S4Vectors::mcols(annot_tu)$source <- "ann"
  all_tu <- c(hml_tu, annot_tu) %>% BiocGenerics::sort()
  idx1 <- find_fusion_tx(all_tu, min_overlap_fusion = min_overlap_fusion, skip = FALSE)
  idx2 <- all_tu[idx1] %>% S4Vectors::mcols() %>% dplyr::as_tibble() %>% dplyr::filter(source == "hml") %>% .$revmap %>% BiocGenerics::unlist()
  if (length(idx2) > 0) {
    message("Found ", length(idx2), " new fusion transcripts;")
    new_fusion <- hml_tx[idx2]
    new_unl <- BiocGenerics::unlist(new_fusion, use.names = FALSE)
    S4Vectors::mcols(new_unl)$type <- "Fusion"
    new_fusion <- BiocGenerics::relist(new_unl, new_fusion)
    fusion_tx <- c(fusion_tx, new_fusion)
    fusion_tu <- fusion_tx %>% call_tu(clust_threshold = clust_threshold) %>% skip_nested()
    hml_tx <- hml_tx[-idx2]
    hml_range <- hml_tx %>% range() %>% BiocGenerics::unlist(use.names = FALSE)
  } else {
    message("Did not find new fusion transcripts;")
    if (length(fusion_tx) > 0) {
      fusion_tu <- fusion_tx %>% call_tu(clust_threshold = clust_threshold) %>% skip_nested()
    } else {
      fusion_tu <- GenomicRanges::GRangesList()
    }
  }
  # Find annotated transcripts which do not overlap with called ones:
  annot_range <- annot_hc %>% range() %>% BiocGenerics::unlist(use.names = FALSE)
  annot_good <- annot_hc[annot_range %outside% hml_range]
  if (length(annot_good) > 0) {
    message("Found ", length(annot_good), " annotated transcripts with TSS and PAS which were absent from the called set;")
    annot_good_unl <- BiocGenerics::unlist(annot_good, use.names = FALSE)
    S4Vectors::mcols(annot_good_unl) <- NULL
    S4Vectors::mcols(annot_good_unl)$score <- 0
    S4Vectors::mcols(annot_good_unl)$type <- "HC"
    max_grp <- hml_tx %>% BiocGenerics::unlist(use.names = FALSE) %>% S4Vectors::mcols() %>% .$grp %>% max()
    S4Vectors::mcols(annot_good_unl)$grp <- seq(max_grp + 1, max_grp + length(annot_good)) %>% rep(times = S4Vectors::elementNROWS(annot_good))
    annot_good <- BiocGenerics::relist(annot_good_unl, annot_good)
    hml_tx <- c(hml_tx, annot_good)
  }
  hml_tx <- sort_grl(hml_tx)
  # Recalculate genes, generate unique gene IDs:
  hml_genes <- call_HC_MC_LC_genes(hml_tx, clust_threshold = clust_threshold, generate_ids = TRUE)
  # Update unique transcript IDs:
  hml_tx <- update_tx_id(hml_tx, hml_genes)
  out <- list("hml_genes" = hml_genes, "hml_tx" = hml_tx, "fusion_tu" = fusion_tu, "fusion_tx" = fusion_tx)
  return(out)
}

# ---------------------------------------------------------------------------------------------------------

update_tx_id <- function(hml_tx, hml_genes) {
  tx_spl <- split_tx_by_type(hml_tx)
  hc_tx <- tx_spl[[1]]
  mc_tx <- tx_spl[[2]]
  lc_tx <- tx_spl[[3]]
  hc_genes <- hml_genes[S4Vectors::mcols(hml_genes)$type == "HC"]
  mc_genes <- hml_genes[S4Vectors::mcols(hml_genes)$type == "MC"]
  lc_genes <- hml_genes[S4Vectors::mcols(hml_genes)$type == "LC"]
  hc_tx <- generate_tx_id(hc_tx, hc_genes)
  mc_tx <- generate_tx_id(mc_tx, mc_genes)
  lc_tx <- generate_tx_id(lc_tx, lc_genes)
  out <- c(hc_tx, mc_tx, lc_tx) %>% sort_grl()
  return(out)
}

# ---------------------------------------------------------------------------------------------------------

summarize_grl_width <- function(grl) {
  stopifnot(BiocGenerics::grepl("GRangesList", class(grl)))
  out <- grl %>% BiocGenerics::width() %>% BiocGenerics::lapply(sum, simplify = FALSE) %>% BiocGenerics::unlist() %>% unname()
  return(out)
}

# ---------------------------------------------------------------------------------------------------------

extend_called_transcripts <- function(hml_tx, annot_exons, tc, mode, min_score_2, min_tx_cov) {
  stopifnot(mode %in% c("start", "end"))
  tc_type <- ifelse(mode == "start", "TSS", "PAS")
  stopifnot(class(tc) == "GRanges" && !is.null(score(tc)))
  # (known transcripts are expected to be already extended to summits of nearby TSS and PAS)
  # Skip known transcripts which do not have a valid TC:
  annot_tx_border <- annot_exons %>% range() %>% BiocGenerics::unlist() %>% GenomicRanges::resize(1, mode)
  valid_tc <- tc[BiocGenerics::score(tc) >= min_score_2]
  annot_exons_good <- annot_exons[annot_tx_border %over% valid_tc]
  # Enumerate the called transcripts:
  old_names <- names(hml_tx)
  names(hml_tx) <- 1:length(hml_tx)
  # Find called transcripts which lack TC:
  called_tx_border <- hml_tx %>% range() %>% BiocGenerics::unlist() %>% GenomicRanges::resize(1, mode)
  over_tc <- called_tx_border %over% tc
  out_1 <- hml_tx[over_tc]
  other <- hml_tx[!over_tc] ####################################################### there can be some HC guys... Check where do they come from!
  if (length(other) == 0) {
    message("\tNo transcripts were extended towards ", tc_type, ";")
    names(hml_tx) <- old_names
    return(hml_tx)
  }
  # Among them, find ones which overlap with known transcripts:
  over_ann <- other %over% annot_exons_good
  out_2 <- other[!over_ann]
  other <- other[over_ann]
  if (length(other) == 0) {
    message("\tNo transcripts were extended towards ", tc_type, ";")
    names(hml_tx) <- old_names
    return(hml_tx)
  }
  # Find all overlaps:
  hits <- GenomicRanges::findOverlaps(other, annot_exons_good)
  if (length(hits) == 0) {
    message("\tNo transcripts were extended towards ", tc_type, ";")
    names(hml_tx) <- old_names
    return(hml_tx)
  }
  par1 <- other[S4Vectors::queryHits(hits)]
  par2 <- annot_exons_good[S4Vectors::subjectHits(hits)]
  # Calculate percent overlaps at the exonic level:
  wi <- pintersect_punion_on_grl(par1, par2, mode = "pintersect") %>% summarize_grl_width()
  w1 <- summarize_grl_width(par1)
  # Consider only called tx with at least one valid overlap:
  valid_over <- wi / w1 >= min_tx_cov
  hits <- hits[valid_over]
  if (length(hits) == 0) {
    message("\tNo transcripts were extended towards ", tc_type, ";")
    names(hml_tx) <- old_names
    return(hml_tx)
  }
  idx <- S4Vectors::queryHits(hits) %>% BiocGenerics::unique()
  out_3 <- other[-idx]
  other <- other[idx]
  par1 <- par1[valid_over]
  par2 <- par2[valid_over]
  wi <- wi[valid_over]
  # For each called tx, choose unique annotated partner with the strongest mutual overlap:
  wu <- pintersect_punion_on_grl(par1, par2, mode = "punion") %>% summarize_grl_width()
  mutual_over <- wi / wu
  best_mutual <- BiocGenerics::tapply(mutual_over, S4Vectors::queryHits(hits), find_uniq_max, simplify = FALSE) %>% BiocGenerics::unlist() %>% unname()
  hits <- hits[best_mutual]
  par1 <- par1[best_mutual] # == other
  par2 <- par2[best_mutual]
  # Check that the annotated transcript extends the called one, not shrinks it:
  par1_border <- par1 %>% range() %>% BiocGenerics::unlist(use.names = FALSE) %>% GenomicRanges::resize(1, mode)
  par2_border <- par2 %>% range() %>% BiocGenerics::unlist(use.names = FALSE) %>% GenomicRanges::resize(1, mode)
  decision <- check_gr_up_down(par1_border, par2_border)
  if (mode == "start") {
    good <- decision == "up" # annotated border (par2) is upstream from the called border (par1)
  } else {
    good <- decision == "down"
  }
  if (!any(good)) {
    message("\tNo transcripts were extended towards ", tc_type, ";")
    names(hml_tx) <- old_names
    return(hml_tx)
  }
  out_4 <- par1[!good]
  par1 <- par1[good]
  par2 <- par2[good]
  par1_border <- par1_border[good]
  par2_border <- par2_border[good]
  # Finally, merge called and annotated transcripts in each pair:
  new_range <- GenomicRanges::pgap(par1_border, par2_border) %>% GenomicRanges::resize(width(.) + 2, "center")
  new_exons <- GenomicRanges::pintersect(par2, new_range) # observe zero-width exons...
  mc_uniq <- par1 %>% BiocGenerics::unlist(use.names = FALSE) %>% S4Vectors::mcols() %>% BiocGenerics::unique()
  rep_idx <- rep(1:length(new_exons), times = S4Vectors::elementNROWS(new_exons))
  new_exons_unl <- BiocGenerics::unlist(new_exons, use.names = FALSE)
  S4Vectors::mcols(new_exons_unl) <- mc_uniq[rep_idx, ]
  new_exons_unl <- new_exons_unl %>% `[`(BiocGenerics::width(.) > 0)
  merged <- par1 %>% BiocGenerics::unlist(use.names = FALSE) %>% c(new_exons_unl) %>% S4Vectors::split(S4Vectors::mcols(.)$grp) %>% GenomicRanges::reduce()
  names(merged) <- names(par1)
  merged_unl <- BiocGenerics::unlist(merged, use.names = FALSE)
  rep_idx2 <- rep(1:length(merged), times = S4Vectors::elementNROWS(merged))
  S4Vectors::mcols(merged_unl) <- mc_uniq[rep_idx2, ]
  merged <- BiocGenerics::relist(merged_unl, merged) %>% BiocGenerics::sort()
  message("\t", length(merged), " transcripts were extended towards ", tc_type, ";")
  out <- c(out_1, out_2, out_3, out_4, merged)
  out <- out[BiocGenerics::order(as.integer(names(out)))]
  names(out) <- old_names
  return(out)
}


# ---------------------------------------------------------------------------------------------------------

adjust_exons_of_called_transcripts <- function(hml_tx, annot_exons, max_exon_diff) {
  old_names <- names(hml_tx)
  # Unlist all called exons:
  tx_unl <- BiocGenerics::unlist(hml_tx, use.names = FALSE)
  names(tx_unl) <- 1:length(tx_unl) # to quickly recover the original order
  stopifnot(!is.null(S4Vectors::mcols(tx_unl)$type) && !is.null(S4Vectors::mcols(tx_unl)$grp))
  message("\t", length(tx_unl), " called exons;")
  # Enumerate the called exons:
  idx_fw <- BiocGenerics::lapply(elementNROWS(hml_tx), function(x) { 1:x }) %>% BiocGenerics::unlist()
  idx_rev <- BiocGenerics::lapply(elementNROWS(hml_tx), function(x) { x:1 }) %>% BiocGenerics::unlist()
  S4Vectors::mcols(tx_unl)$pos_fw <- ifelse(BiocGenerics::strand(tx_unl) == "+", idx_fw, idx_rev) # position of exon relative to the transcript start
  S4Vectors::mcols(tx_unl)$pos_rev <- ifelse(BiocGenerics::strand(tx_unl) == "+", idx_rev, idx_fw) # position of exon relative to the transcript end
  # Get unique annotated exons:
  annot_unl <- BiocGenerics::unlist(annot_exons, use.names = FALSE) %>% BiocGenerics::unique()
  # Adjust exon starts:
  tx_unl <- adjust_exon_borders_by_ann(tx_unl, annot_unl, mode = "start", max_exon_diff = max_exon_diff)
  # Adjust exon ends:
  tx_unl <- adjust_exon_borders_by_ann(tx_unl, annot_unl, mode = "end", max_exon_diff = max_exon_diff)
  # Prepare the output:
  S4Vectors::mcols(tx_unl)$pos_fw <- NULL
  S4Vectors::mcols(tx_unl)$pos_rev <- NULL
  names(tx_unl) <- NULL
  out <- S4Vectors::split(tx_unl, S4Vectors::mcols(tx_unl)$grp)
  names(out) <- old_names
  return(out)
}

# -----------------------------------------------------------------------------------------------------

adjust_exon_borders_by_ann <- function(tx_unl, annot_unl, mode, max_exon_diff) {
  stopifnot(!is.null(S4Vectors::mcols(tx_unl)$pos_fw) && !is.null(S4Vectors::mcols(tx_unl)$pos_rev))
  stopifnot(mode %in% c("start", "end"))
  if (mode == "start") {
    exclude <- S4Vectors::mcols(tx_unl)$pos_fw == 1
    opp <- "end"
  } else {
    exclude <- S4Vectors::mcols(tx_unl)$pos_rev == 1
    opp <- "start"
  }
  out_1 <- tx_unl[exclude]
  other <- tx_unl[!exclude]
  over_ann <- other %over% annot_unl
  out_2 <- other[!over_ann]
  other <- other[over_ann]
  hits <- GenomicRanges::findOverlaps(other, annot_unl)
  par1 <- other[S4Vectors::queryHits(hits)]
  par2 <- annot_unl[S4Vectors::subjectHits(hits)]
  # In case of multiple hits, choose the best match:
  overlap <- GenomicRanges::pintersect(par1, par2) %>% BiocGenerics::width()
  best <- BiocGenerics::tapply(overlap, S4Vectors::queryHits(hits), find_uniq_max, simplify = FALSE) %>% BiocGenerics::unlist() %>% unname()
  par1 <- par1[best]
  par2 <- par2[best]
  b1 <- GenomicRanges::resize(par1, 0, mode)
  b2 <- GenomicRanges::resize(par2, 0, mode)
  diff <- GenomicRanges::pgap(b1, b2) %>% BiocGenerics::width()
  good_diff <- diff > 0 & diff <= max_exon_diff
  out_3 <- other[!good_diff]
  out_4 <- other[good_diff]
  if (length(out_4) > 0) {
    par2 <- par2[good_diff]
    IRanges::ranges(out_4) <- GenomicRanges::pgap(GenomicRanges::resize(out_4, 0, opp), GenomicRanges::resize(par2, 0, mode)) %>% IRanges::ranges()
    message("\t", length(out_4), " exon ", mode, "s were adjusted by the annotation;")
  }
  out <- c(out_1, out_2, out_3, out_4)
  out <- out[BiocGenerics::order(as.integer(names(out)))]
  return(out)
}

