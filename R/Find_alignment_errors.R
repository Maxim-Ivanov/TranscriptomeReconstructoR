#' Detect the most common type of alignment errors in full-length RNA-seq reads aligned by Minimap2
#'
#' @param long_reads \code{GRangesList} object.
#' @param min_read_support Positive integer.
#' @section Details:
#' The input \code{GRangesList} object is returned by the \code{adjust_exons_of_long_reads()} function.\cr
#' The Minimap2 aligner is prone to under-split the full-length RNA-seq reads. In certain reads, two adjacent exons may appear as a single subalignment.
#' Assuming that most of the reads are still aligned correctly, one can detect alignment errors by comparing each subalignment to the set of constitutive exons.
#' The constitutive exons are defined here as exonic intervals supported by more than 50% of reads in given locus (but not less than \code{min_read_support} reads).
#' Subalignment with an alternative 5'- and/or 3-border (relative to a constitutive exon) are considered valid alternative exons, only if the next constitutive exon is also present in the read.
#' Otherwise, the alignment is marked as a possible alignment error.
#' @return \code{GRangesList} object which is identical to the input object, except for the extra logical column \code{aln_error}.
#' @export
detect_alignment_errors <- function(long_reads, min_read_support = 2) {
  stopifnot(BiocGenerics::grepl("GRangesList", class(long_reads)))
  stopifnot(is.numeric(min_read_support) && length(min_read_support) == 1 && min_read_support %% 1 == 0 && min_read_support > 0)
  all_exons <- BiocGenerics::unlist(long_reads, use.names = FALSE)
  ### 1) Call constitutive exons:
  # Exons from complete reads:
  exons_compl <- all_exons[S4Vectors::mcols(all_exons)$complete]
  # Exons from non-complete TSS-only and PAS-only reads (given that TSS/PAS was not used in complete reads):
  exons_with_tss <- all_exons[!S4Vectors::mcols(all_exons)$complete & S4Vectors::mcols(all_exons)$over_tc == "tss_only" & !S4Vectors::mcols(all_exons)$tss_id %in% S4Vectors::mcols(exons_compl)$tss_id]
  exons_with_pas <- all_exons[!S4Vectors::mcols(all_exons)$complete & S4Vectors::mcols(all_exons)$over_tc == "pas_only" & !S4Vectors::mcols(all_exons)$pas_id %in% S4Vectors::mcols(exons_compl)$pas_id]
  # Skip free terminal exons from non-complete reads:
  exons_with_tss <- exons_with_tss %>% `[`(S4Vectors::mcols(.)$pos_rev != 1)
  exons_with_pas <- exons_with_pas %>% `[`(S4Vectors::mcols(.)$pos_fw != 1)
  # Combine complete, TSS-only and PAS-only exons, add grouping info:
  exons_for_const <- c(exons_compl, exons_with_tss, exons_with_pas)
  S4Vectors::mcols(exons_for_const)$grp_id <- BiocGenerics::paste(S4Vectors::mcols(exons_for_const)$tss_id, S4Vectors::mcols(exons_for_const)$pas_id, sep = "|")
  # Add internal exons from reads truncated on both sides:
  exons_trunc <- all_exons[!S4Vectors::mcols(all_exons)$complete & !is.na(S4Vectors::mcols(all_exons)$grp_id)]
  exons_trunc <- exons_trunc %>% `[`(S4Vectors::mcols(.)$pos_fw != 1 & S4Vectors::mcols(.)$pos_rev != 1)
  exons_for_const <- c(exons_for_const, exons_trunc) %>% `[`(BiocGenerics::order(S4Vectors::mcols(.)$exon_id))
  # Find constitutive exons (supported by more than 50% of reads in given TSS-PAS group):
  const_exons <- find_constitutive_exons(exons_for_const, S4Vectors::mcols(exons_for_const)$grp_id) %>% `[`(S4Vectors::mcols(.)$const) # non-unique (one exon per read)
  # Deduplicate constitutive exons, recalculate support:
  const_exons <- BiocGenerics::unique(const_exons) %>% GenomicRanges::granges()
  S4Vectors::mcols(const_exons)$type <- "const"
  S4Vectors::mcols(const_exons)$support <- GenomicRanges::countOverlaps(const_exons, all_exons, type = "equal")
  # Filter constitutive exons by the number of supporting long reads:
  const_exons <- const_exons[S4Vectors::mcols(const_exons)$support >= min_read_support]
  message(length(const_exons), " unique constitutive exons found;")
  ### 2) Use constitutive exons to correct alignment errors in all reads:
  all_exons_out <- find_alignment_errors(const_exons, all_exons)
  # Print stats on the suspected alignment errors:
  message(S4Vectors::mcols(all_exons_out) %>% dplyr::as_tibble() %>% dplyr::filter(aln_error) %>% .$read_id %>% BiocGenerics::unique() %>% length(), " reads contain at least one alignment error;")
  out <- all_exons_out %>% `[`(BiocGenerics::order(S4Vectors::mcols(.)$exon_id)) %>% S4Vectors::split(S4Vectors::mcols(.)$read_id)
  return(out)
}

# -------------------------------------------------------------------------------------------------------------

find_constitutive_exons <- function(exons, grp) {
  S4Vectors::mcols(exons)$grp <- grp
  tbl <- exons %>% dplyr::as_tibble() %>% dplyr::group_by(grp) %>% dplyr::mutate(n_rd = dplyr::n_distinct(read_id)) %>% ungroup()
  tbl <- tbl %>% unite(full_id, start, end, grp, sep = "_") %>% group_by(full_id) %>% mutate(n_ex = n()) %>% dplyr::ungroup()
  tbl <- tbl %>% tidyr::separate(full_id, into = c("start", "end", "grp"), sep = "_") %>% dplyr::mutate(ref = n_ex > n_rd / 2)
  S4Vectors::mcols(exons)$const <- tbl$ref
  return(exons)
}

# ------------------------------------------------------------------------------------------------------------

test_exon_borders <- function(other, const_exons, all_exons, all_id, mode, allow_exitrons) {
  stopifnot(mode %in% c("end", "start"))
  colname <- paste0("good_", mode)
  free_border <- GenomicRanges::resize(other, 1, mode) %outside% GenomicRanges::resize(const_exons, 1, mode)
  other_p1 <- other[!free_border] # reads with correct border
  if (length(other_p1) > 0) {
    S4Vectors::mcols(other_p1)[, colname] <- TRUE
  }
  other <- other[free_border]
  lgl <- GenomicRanges::countOverlaps(other, const_exons) > 1
  other_p2 <- other[lgl] # reads with wrong border which overlap more than 1 constitutive exon
  if (length(other_p2) > 0) {
    S4Vectors::mcols(other_p2)[, colname] <- FALSE
  }
  other <- other[!lgl]
  # Find next subalignment within the same read:
  if (mode == "end") {
    next_subaln_id <- stringr::str_c(S4Vectors::mcols(other)$read_id, S4Vectors::mcols(other)$pos_fw + 1, sep = "_")
  } else {
    next_subaln_id <- stringr::str_c(S4Vectors::mcols(other)$read_id, S4Vectors::mcols(other)$pos_fw - 1, sep = "_")
  }
  idx <- BiocGenerics::match(next_subaln_id, all_id)
  other_p3 <- other[is.na(idx)] # there is no downstream subalignment
  if (length(other_p3) > 0) {
    S4Vectors::mcols(other_p3)[, colname] <- TRUE
  }
  other <- other[!is.na(idx)]
  idx <- idx[!is.na(idx)]
  next_subaln_par <- all_exons[idx]
  # The next subalignment starts within the same constitutive exon:
  # (either exitrons, or homopolymeric errors)
  ref_exon_par <- const_exons[GenomicRanges::findOverlaps(other, const_exons, select = "first")]
  if (mode == "end") {
    lgl_1 <- GenomicRanges::resize(other, 1, "end") %>% poverlaps_gr(ref_exon_par)
    lgl_2 <- GenomicRanges::resize(next_subaln_par, 1, "start") %>% poverlaps_gr(ref_exon_par)
  } else {
    lgl_1 <- GenomicRanges::resize(other, 1, "start") %>% poverlaps_gr(ref_exon_par)
    lgl_2 <- GenomicRanges::resize(next_subaln_par, 1, "end") %>% poverlaps_gr(ref_exon_par)
  }
  lgl <- lgl_1 & lgl_2
  other_p4 <- other[lgl]
  if (length(other_p4) > 0) {
    S4Vectors::mcols(other_p4)[, colname] <- allow_exitrons
  }
  other <- other[!lgl]
  next_subaln_par <- next_subaln_par[!lgl]
  # Check if the next subalignment has a correct start:
  if (mode == "end") {
    lgl <- GenomicRanges::resize(next_subaln_par, 1, "start") %outside% GenomicRanges::resize(const_exons, 1, "start")
  } else {
    lgl <- GenomicRanges::resize(next_subaln_par, 1, "end") %outside% GenomicRanges::resize(const_exons, 1, "end")
  }
  other_p5 <- other[lgl]
  if (length(other_p5) > 0) {
    S4Vectors::mcols(other_p5)[, colname] <- FALSE
  }
  other <- other[!lgl]
  next_subaln_par <- next_subaln_par[!lgl] # its start = start of downstream ref exon
  # Find the upstream ref exon:
  ref_par <- const_exons[GenomicRanges::findOverlaps(other, const_exons, select = "first")]
  # Is there any other ref exon between them?
  gap <- GenomicRanges::pgap(ref_par, next_subaln_par)
  lgl <- BiocGenerics::width(gap) == 0
  other_p6 <- other[lgl]
  if (length(other_p6) > 0) {
    S4Vectors::mcols(other_p6)[, colname] <- TRUE
  }
  other_p7 <- other[!lgl]
  gap <- gap[!lgl]
  hits <- GenomicRanges::findOverlaps(const_exons, gap, type = "within")
  if (length(other_p7) > 0) {
    S4Vectors::mcols(other_p7)[, colname] <- TRUE
    S4Vectors::mcols(other_p7)[BiocGenerics::unique(S4Vectors::subjectHits(hits)), colname] <- FALSE
  }
  # Assemble the exons:
  other <- c(other_p1, other_p2, other_p3, other_p4, other_p5, other_p6, other_p7) %>% BiocGenerics::sort()
  return(other)
}

#------------------------------------------------------------------------------------------------------------------

find_alignment_errors <- function(const_exons, all_exons, allow_exitrons = FALSE) {
  all_id <- stringr::str_c(S4Vectors::mcols(all_exons)$read_id, S4Vectors::mcols(all_exons)$pos_fw, sep = "_")
  # Free terminal exons of non-complete reads:
  first_free <- !S4Vectors::mcols(all_exons)$complete & S4Vectors::mcols(all_exons)$pos_fw == 1 & is.na(S4Vectors::mcols(all_exons)$tss_id)
  last_free <- !S4Vectors::mcols(all_exons)$complete & S4Vectors::mcols(all_exons)$pos_rev == 1 & is.na(S4Vectors::mcols(all_exons)$pas_id)
  lgl_1 <- first_free | last_free
  # Exons which do not overlap with any constitutive exon:
  lgl_2 <- all_exons %outside% const_exons
  # Exons which precisely overlap any constitutive exon:
  lgl_3 <- GenomicRanges::overlapsAny(all_exons, const_exons, type = "equal")
  # Exons which start in one constitutive exon and end in another (IR events):
  lgl_4 <- GenomicRanges::resize(all_exons, 1, "start") %over% GenomicRanges::resize(const_exons, 1, "start") & GenomicRanges::resize(all_exons, 1, "end") %over% GenomicRanges::resize(const_exons, 1, "end")
  lgl <- lgl_1 | lgl_2 | lgl_3 | lgl_4
  out <- all_exons[lgl]
  S4Vectors::mcols(out)$aln_error <- FALSE
  other <- all_exons[!lgl] # read exons which have a partial overlap with at least one constitutive exon
  # Test 3' borders:
  other <- test_exon_borders(other, const_exons, all_exons, all_id, mode = "end", allow_exitrons = allow_exitrons)
  # Test 5' borders:
  other <- test_exon_borders(other, const_exons, all_exons, all_id, mode = "start", allow_exitrons = allow_exitrons)
  # Make decisions on which read exons are not alignment artifacts:
  good <- S4Vectors::mcols(other)$good_start & S4Vectors::mcols(other)$good_end
  S4Vectors::mcols(other)$aln_error <- !good
  S4Vectors::mcols(other)$good_start <- NULL
  S4Vectors::mcols(other)$good_end <- NULL
  message(sum(!good), " exons were suspected to be alignment errors (", length(BiocGenerics::unique(GenomicRanges::granges(other[!good]))), " unique);")
  # Combine all exons:
  out <- c(out, other)
  out <- out[BiocGenerics::order(S4Vectors::mcols(out)$exon_id)]
  return(out)
}
