#' Adjust outer borders of full-length RNA-seq reads to the nearby TSS and PAS
#'
#' @param grl \code{GRangesList} object.
#' @param tss,pas \code{GRanges} objects.
#' @param read_flanks_up,read_flanks_down Integer vectors of length 2.
#' @param trunc_limit Positive integer.
#' @param mode Character. Must be either "reads" (default) or "transcripts".
#' @param verbose Logical.
#' @section Details:
#' If \code{mode == "reads"}, then \code{grl} is expected to be full-length RNA-seq reads returned by \code{load_BAM_files(mode = "long_read")}.\cr
#' If \code{mode == "transcripts"}, then \code{grl} is expected to be known transcripts returned by \code{exonsBy(txdb, by = "tx")} (where \code{txdb} is a \code{GenomicFeatures} object).\cr
#' \code{tss} and \code{pas} are returned by \code{call_TCs()} on 5'- and 3'-tag sequencing data, respectively.\cr
#' The function searches for called TSS within \code{read_flanks_up} bp around 5' ends of long reads (or known transcripts), and for PAS within \code{read_flanks_down} around their 3' ends.\cr
#' When neighboring TSS/PAS are found, the long reads are extended towards their summits.\cr
#' If multiple TSS or PAS are found within the search window, then the one with the highest score is chosen.\cr
#' The long read may be truncated, if the TSS/PAS summit is located inside of its terminal exon (however truncation to a length smaller than \code{trunc_limit} bp is disallowed).\cr
#' @return \code{GRangesList} object containing the adjusted exons + some extra mcols.
#' @export
extend_long_reads_to_TSS_and_PAS <- function(grl, tss, pas, read_flanks_up = c(-100, 10), read_flanks_down = c(-50, 50), trunc_limit = 20, mode = "reads", verbose = isTRUE(mode == "reads")) {
  stopifnot(BiocGenerics::grepl("GRangesList", class(grl)))
  stopifnot(class(tss) == "GRanges" && class(pas) == "GRanges")
  stopifnot(is.numeric(read_flanks_up) && length(read_flanks_up) == 2 && read_flanks_up[[1]] %% 1 == 0 && read_flanks_up[[2]] %% 1 == 0)
  stopifnot(is.numeric(read_flanks_down) && length(read_flanks_down) == 2 && read_flanks_down[[1]] %% 1 == 0 && read_flanks_down[[2]] %% 1 == 0)
  stopifnot(mode %in% c("reads", "transcripts"))
  stopifnot(is.numeric(trunc_limit) && length(trunc_limit) == 1 && trunc_limit %% 1 == 0 && trunc_limit > 0)
  stopifnot(is.logical(verbose) && length(verbose) == 1)
  if (mode == "transcripts") {
    old_mcols <- grl %>% BiocGenerics::unlist(use.names = FALSE) %>% S4Vectors::mcols()
    old_names <- names(grl)
  }
  # Annotate read/tx exons with their indexes and relative positions:
  all_exons <- grl %>% add_mcols_to_grl() %>% BiocGenerics::unlist(use.names = FALSE)
  # Analyze intersections of reads/tx with TSS and PAS:
  if (isTRUE(verbose)) {
    message(length(grl), " input ", mode, ";")
    raw_range <- grl %>% range() %>% BiocGenerics::unlist()
    over_tss <- GenomicRanges::resize(raw_range, 1, "start") %over% tss
    message("\t", sum(over_tss), " of them start in TSS;")
    over_pas <- GenomicRanges::resize(raw_range, 1, "end") %over% pas
    message("\t", sum(over_pas), " end in PAS;")
    message("\t", sum(over_pas & over_tss), " start in TSS and end in PAS;")
  }
  # Extend reads/tx to summits of nearby TSS and PAS (if any):
  all_exons <- all_exons %>% extend_reads(tc = tss, mode = "tss", flank_out = abs(read_flanks_up[[1]]), flank_in = abs(read_flanks_up[[2]]), trunc_limit = trunc_limit)
  if (isTRUE(verbose)) {
    w0 <- BiocGenerics::width(raw_range)
    w1 <- all_exons %>% S4Vectors::split(S4Vectors::mcols(.)$read_id) %>% range() %>% BiocGenerics::unlist() %>% BiocGenerics::width()
    message(sum(w1 != w0), " ", mode, " were extended to summits of nearby TSS;")
  }
  all_exons <- all_exons %>% extend_reads(tc = pas, mode = "pas", flank_out = abs(read_flanks_down[[2]]), flank_in = abs(read_flanks_down[[1]]), trunc_limit = trunc_limit)
  grl <- all_exons %>% S4Vectors::split(S4Vectors::mcols(.)$read_id)
  new_range <- grl %>% range() %>% BiocGenerics::unlist()
  if (isTRUE(verbose)) {
    w2 <- BiocGenerics::width(new_range)
    message(sum(w2 != w1), " ", mode, " were extended to summits of nearby PAS;")
    # Re-analyze intersections of extended long reads with TSS and PAS:
    message("After extension:")
    over_tss <- GenomicRanges::resize(new_range, 1, "start") %over% tss
    message("\t", sum(over_tss), " ", mode, " start in TSS;")
    over_pas <- GenomicRanges::resize(new_range, 1, "end") %over% pas
    message("\t", sum(over_pas), " end in PAS;")
    message("\t", sum(over_pas & over_tss), " start in TSS and end in PAS;")
  }
  # Find overlapping TSS and PAS, if any:
  idx_tss <- rep(NA, times = length(grl))
  idx_pas <- rep(NA, times = length(grl))
  h1 <- GenomicRanges::resize(new_range, 1, "start") %>% GenomicRanges::findOverlaps(tss)
  h2 <- GenomicRanges::resize(new_range, 1, "end") %>% GenomicRanges::findOverlaps(pas)
  idx_tss[S4Vectors::queryHits(h1)] <- S4Vectors::subjectHits(h1)
  idx_pas[S4Vectors::queryHits(h2)] <- S4Vectors::subjectHits(h2)
  S4Vectors::mcols(all_exons)$tss_id <- rep(idx_tss, times = S4Vectors::elementNROWS(grl))
  S4Vectors::mcols(all_exons)$pas_id <- rep(idx_pas, times = S4Vectors::elementNROWS(grl))
  # Classify reads/tx by overlap with TSS and PAS:
  S4Vectors::mcols(all_exons)$over_tc <- "no"
  has_tss <- !is.na(S4Vectors::mcols(all_exons)$tss_id)
  has_pas <- !is.na(S4Vectors::mcols(all_exons)$pas_id)
  S4Vectors::mcols(all_exons)[has_tss & has_pas, "over_tc"] <- "both"
  S4Vectors::mcols(all_exons)[has_tss & !has_pas, "over_tc"] <- "tss_only"
  S4Vectors::mcols(all_exons)[!has_tss & has_pas, "over_tc"] <- "pas_only"
  # Also add read/tx width:
  S4Vectors::mcols(all_exons)$read_width <- rep(BiocGenerics::width(new_range), times = S4Vectors::elementNROWS(grl))
  # Extend terminal exons of TSS-only and PAS-only reads/tx:
  if (isTRUE(verbose)) {
    message("Trying to extend terminal exons of TSS-only and PAS-only ", mode, "...")
  }
  fl_tbl <- S4Vectors::mcols(all_exons) %>% dplyr::as_tibble() %>% dplyr::filter(over_tc == "both")
  tss_in_use <-  fl_tbl$tss_id %>% BiocGenerics::unique()
  pas_in_use <- fl_tbl$pas_id %>% BiocGenerics::unique()
  tss_only <- S4Vectors::mcols(all_exons)$over_tc == "tss_only" & !S4Vectors::mcols(all_exons)$tss_id %in% tss_in_use
  pas_only <- S4Vectors::mcols(all_exons)$over_tc == "pas_only" & !S4Vectors::mcols(all_exons)$pas_id %in% pas_in_use
  tss_only_adj <- all_exons[tss_only] %>% extend_terminal_exons(mode = "tss")
  pas_only_adj <- all_exons[pas_only] %>% extend_terminal_exons(mode = "pas")
  other_exons <- all_exons[!tss_only & !pas_only]
  if (length(other_exons) > 0) {
    S4Vectors::mcols(other_exons)$term_adj <- NA
  }
  all_exons_2 <- c(other_exons, tss_only_adj, pas_only_adj) %>% `[`(BiocGenerics::order(S4Vectors::mcols(.)$exon_id))
  # TSS_only and PAS_only reads/tx with extended terminal exons can be regarded as "complete":
  adj_read <- BiocGenerics::tapply(S4Vectors::mcols(all_exons_2)$term_adj, S4Vectors::mcols(all_exons_2)$read_id, any, simplify = FALSE) %>% BiocGenerics::unlist() %>% unname() %>% ifelse(is.na(.), FALSE, .)
  S4Vectors::mcols(all_exons_2)$complete <- S4Vectors::mcols(all_exons_2)$over_tc == "both" | rep(adj_read, times = S4Vectors::elementNROWS(grl))
  if (isTRUE(verbose)) {
    message("\t", sum(adj_read), " TSS-only or PAS-only ", mode, " were re-classified as complete;")
    # Extend terminal exons of noTSS-noPAS reads/tx:
    message("Trying to extend terminal exons of noTSS-noPAS ", mode, "...")
  }
  trunc <- S4Vectors::mcols(all_exons_2)$over_tc == "no"
  other_1 <- all_exons_2[!trunc]
  if (length(other_1) > 0) {
    S4Vectors::mcols(other_1)$last_adj <- NA
    S4Vectors::mcols(other_1)$first_adj <- NA
  }
  exons_trunc <- all_exons_2[trunc]
  reads_trunc <- S4Vectors::split(exons_trunc, S4Vectors::mcols(exons_trunc)$read_id)
  range_trunc <- reads_trunc %>% range() %>% BiocGenerics::unlist()
  # Group by overlap, not by TSS/PAS indexes:
  groups <- GenomicRanges::reduce(range_trunc)
  support <- GenomicRanges::countOverlaps(groups, range_trunc)
  groups <- groups[support >= 2]
  hits <- GenomicRanges::findOverlaps(exons_trunc, groups, type = "within")
  other_2 <- exons_trunc[-S4Vectors::queryHits(hits)]
  if (length(other_2) > 0) {
    S4Vectors::mcols(other_2)$last_adj <- NA
    S4Vectors::mcols(other_2)$first_adj <- NA
  }
  exons_trunc <- exons_trunc[S4Vectors::queryHits(hits)]
  S4Vectors::mcols(exons_trunc)$grp_id <- S4Vectors::subjectHits(hits)
  exons_trunc <- exons_trunc %>% extend_terminal_exons(mode = "tss", colname_grp = "grp_id", colname_out = "last_adj")
  exons_trunc <- exons_trunc %>% extend_terminal_exons(mode = "pas", colname_grp = "grp_id", colname_out = "first_adj")
  all_exons_3 <- c(other_1, other_2, exons_trunc) %>% `[`(BiocGenerics::order(S4Vectors::mcols(.)$exon_id))
  # NoTSS-noPAS reads/tx with both terminal exons extended are also "complete":
  adj_last <- BiocGenerics::tapply(S4Vectors::mcols(all_exons_3)$last_adj, S4Vectors::mcols(all_exons_3)$read_id, any, simplify = FALSE) %>% BiocGenerics::unlist() %>% unname() %>% ifelse(is.na(.), FALSE, .)
  adj_first <- BiocGenerics::tapply(S4Vectors::mcols(all_exons_3)$first_adj, S4Vectors::mcols(all_exons_3)$read_id, any, simplify = FALSE) %>% BiocGenerics::unlist() %>% unname() %>% ifelse(is.na(.), FALSE, .)
  adj_read <- adj_last & adj_first
  if (isTRUE(verbose)) {
    message("\t", sum(adj_read), " noTSS-noPAS ", mode, " were re-classified as complete;")
  }
  S4Vectors::mcols(all_exons_3)$complete <- S4Vectors::mcols(all_exons_3)$complete | rep(adj_read, times = S4Vectors::elementNROWS(grl))
  if (mode == "reads") {
    S4Vectors::mcols(all_exons_3)$grp_id <- NULL ###
    out <- all_exons_3 %>% S4Vectors::split(S4Vectors::mcols(.)$read_id)
  } else {
    groups <- S4Vectors::mcols(all_exons_3)$read_id
    S4Vectors::mcols(all_exons_3) <- old_mcols
    out <- S4Vectors::split(all_exons_3, groups)
    names(out) <- old_names
  }
  return(out)
}

# -----------------------------------------------------------------------------------------------------------------

add_mcols_to_grl <- function(grl) {
  # GRangesList is assumed to be sorted by coordinates at both levels
  gr <- BiocGenerics::unlist(grl, use.names = FALSE)
  S4Vectors::mcols(gr)$exon_id <- 1:length(gr) # global unique id of the exon
  S4Vectors::mcols(gr)$read_id <- rep(1:length(grl), times = S4Vectors::elementNROWS(grl)) # global unique id of the read
  idx_fw <- BiocGenerics::lapply(S4Vectors::elementNROWS(grl), function(x) { 1:x }) %>% BiocGenerics::unlist()
  idx_rev <- BiocGenerics::lapply(S4Vectors::elementNROWS(grl), function(x) { x:1 }) %>% BiocGenerics::unlist()
  S4Vectors::mcols(gr)$pos_fw <- ifelse(BiocGenerics::strand(gr) == "+", idx_fw, idx_rev) # position of exon relative to the read start
  S4Vectors::mcols(gr)$pos_rev <- ifelse(BiocGenerics::strand(gr) == "+", idx_rev, idx_fw) # position of exon relative to the read end
  S4Vectors::mcols(gr)$terminal <- S4Vectors::mcols(gr)$pos_fw == 1 | S4Vectors::mcols(gr)$pos_rev == 1
  S4Vectors::mcols(gr)$n <- rep(S4Vectors::elementNROWS(grl), times = S4Vectors::elementNROWS(grl)) # total number of exons in the read
  out <- BiocGenerics::relist(gr, grl)
  return(out)
}

# ---------------------------------------------------------------------------------------------------------------

extend_reads <- function(exons, tc, mode, extend_to = "summit", flank_out, flank_in, trunc_limit) {
  stopifnot(mode %in% c("tss", "pas"))
  stopifnot(extend_to %in% c("summit", "outer"))
  stopifnot(class(exons) == "GRanges")
  stopifnot(class(tc) == "GRanges")
  if (mode == "tss") {
    tc_border <- unname(tc) %>% GenomicRanges::resize(1, "start")
    term <- S4Vectors::mcols(exons)$pos_fw == 1
    win <- suppressWarnings(exons[term] %>% GenomicRanges::flank(flank_out) %>% GenomicRanges::trim() %>% GenomicRanges::resize(width(.) + flank_in, "start"))
  } else {
    tc_border <- unname(tc) %>% GenomicRanges::resize(1, "end")
    term <- S4Vectors::mcols(exons)$pos_rev == 1
    win <- suppressWarnings(exons[term] %>% GenomicRanges::flank(flank_out, start = FALSE) %>% GenomicRanges::trim() %>% GenomicRanges::resize(width(.) + flank_in, "end"))
  }
  if (extend_to == "summit") {
    stopifnot(!is.null(S4Vectors::mcols(tc)$thick))
    tc_border <- S4Vectors::mcols(tc)$thick %>% unname()
  }
  out_1 <- exons[!term]
  term_exons <- exons[term] # term_exons is parallel to win
  over <- win %over% tc
  out_2 <- term_exons[!over]
  term_exons <- term_exons[over]
  win <- win[over] # still parallel
  hits <- GenomicRanges::findOverlaps(win, tc)
  tc_par <- tc[S4Vectors::subjectHits(hits)]
  best <- BiocGenerics::score(tc_par) %>% BiocGenerics::tapply(S4Vectors::queryHits(hits), find_uniq_max, simplify = FALSE) %>% BiocGenerics::unlist() %>% unname()
  hits <- hits[best]
  par2 <- tc_border[S4Vectors::subjectHits(hits)]
  if (mode == "tss") {
    par1 <- GenomicRanges::resize(term_exons, 1, "end")
    # Check that par1 is downstream from par2 (to avoid extending terminal exon in the wrong direction):
    valid_pair <- check_gr_up_down(par2, par1, offset = trunc_limit) == "down"
  } else if (mode == "pas") {
    par1 <- GenomicRanges::resize(term_exons, 1, "start")
    # Check that par1 is upstream from par2:
    valid_pair <- check_gr_up_down(par1, par2, offset = trunc_limit) == "down"
  }
  out_3 <- term_exons[!valid_pair]
  par1 <- par1[valid_pair]
  par2 <- par2[valid_pair]
  out_4 <- GenomicRanges::punion(par1, par2, fill.gap = TRUE)
  S4Vectors::mcols(out_4) <- S4Vectors::mcols(par1)
  out <- c(out_1, out_2, out_3, out_4)
  out <- out[BiocGenerics::order(S4Vectors::mcols(out)$exon_id)]
  return(out)
}

#------------------------------------------------------------------------------------------------------------------

extend_terminal_exons <- function(exons, mode, colname_grp = paste0(mode, "_id"), colname_out = "term_adj") {
  stopifnot(mode %in% c("tss", "pas"))
  stopifnot(colname_grp %in% c("tss_id", "pas_id", "grp_id"))
  if (mode == "tss") {
    lgl <- S4Vectors::mcols(exons)$pos_rev == 1 # extend last exons
  } else {
    lgl <- S4Vectors::mcols(exons)$pos_fw == 1 # extend first exons
  }
  target_exons <- exons[lgl]
  out1 <- exons[!lgl]
  if (length(out1) > 0) {
    S4Vectors::mcols(out1)[, colname_out] <- NA
  }
  grp1 <- S4Vectors::mcols(target_exons)[, colname_grp]
  target_exons <- target_exons[BiocGenerics::order(grp1)]
  grp1 <- grp1[BiocGenerics::order(grp1)]
  guide <- BiocGenerics::tapply(S4Vectors::mcols(target_exons)$read_width, grp1, find_uniq_max, simplify = FALSE) %>% BiocGenerics::unlist() %>% unname() # it might be better to choose the shortest candidate...
  guide_exons <- target_exons[guide]
  if (mode == "tss") {
    hits <- GenomicRanges::resize(target_exons, 1, "end") %>% GenomicRanges::findOverlaps(guide_exons)
  } else {
    hits <- GenomicRanges::resize(target_exons, 1, "start") %>% GenomicRanges::findOverlaps(guide_exons)
  }
  grp2 <- S4Vectors::mcols(guide_exons)[, colname_grp]
  valid <- grp1[S4Vectors::queryHits(hits)] == grp2[S4Vectors::subjectHits(hits)] # consider only overlaps within the group
  hits <- hits[valid]
  out2 <- target_exons[-S4Vectors::queryHits(hits)] # target exons which fall outside of the guide exon and therefore cannot be extended
  if (length(out2) > 0) {
    S4Vectors::mcols(out2)[, colname_out] <- FALSE
  }
  par1 <- target_exons[S4Vectors::queryHits(hits)]
  par2 <- guide_exons[S4Vectors::subjectHits(hits)]
  if (mode == "tss") {
    out3 <- GenomicRanges::punion(GenomicRanges::resize(par1, 1, "start"), GenomicRanges::resize(par2, 1, "end"), fill.gap = TRUE)
  } else {
    out3 <- GenomicRanges::punion(GenomicRanges::resize(par1, 1, "end"), GenomicRanges::resize(par2, 1, "start"), fill.gap = TRUE)
  }
  S4Vectors::mcols(out3) <- S4Vectors::mcols(par1)
  if (length(out3) > 0) {
    S4Vectors::mcols(out3)[, colname_out] <- TRUE
  }
  #message(length(out3), " terminal exons were adjusted;")
  out <- c(out1, out2, out3)
  out <- out[BiocGenerics::order(S4Vectors::mcols(out)$exon_id)]
  return(out)
}
