#' Batch load BAM files
#'
#' @param bamfiles Character vector of BAM filenames.
#' @param mode Character. Must be one of the following: "long_read", "tss", "pas" or "nascent".
#' @param ngs_mode Character. Must be either "SE" or "PE".
#' @param skip_split_aln Logical.
#' @param strand_mode Integer. Must be either 1 or 2.
#' @param filter_tlen Non-negative integer.
#' @section Details:
#' The values of \code{mode} correspond to different NGS library preparation protocols:
#' \itemize{
#'     \item \code{"tss"} for 5' tag sequencing (e.g. CAGE)
#'     \item \code{"pas"} for 3' tag sequencing (e.g. PAT-seq)
#'     \item \code{"long_read"} for full-length RNA sequencing (e.g. Direct RNA-seq by Oxford Nanopore)
#'     \item \code{"nascent"} for nascent RNA sequencing (e.g. NET-seq of GRO-seq)
#' }
#' if \code{mode == "long_reads"}, then each \code{GRanges} element in the returned \code{GRangesList} object corresponds to an individual long read
#' (each genomic range represents an exonic subalignment). The replicate samples are pooled together;\cr
#' If \code{mode %in% c("tss", "pas")}, then each \code{GRanges} element in the returned \code{GRangesList} corresponds to an individual sample (replicates are not pooled).
#' The input short reads are truncated to the first sequenced base and converted to genomic coverage;\cr
#' If \code{mode == "nascent"}, then the whole reads are converted to genomic coverage, and the replicates are pooled.\cr
#' The following arguments are used only \code{mode == "nascent"}:
#' \code{ngs_mode} describes if the sequencing was single-end or paired-end. Most often nascent RNA is sequenced in single-end mode, however some protocols (pNET-seq, plaNET-seq) require PE sequencing.\cr
#' \code{skip_split_aln} set to \code{TRUE} requires to skip split alignments (because nascent RNA molecules usually are not expected to be spliced).
#' However, if nascent RNA is sequenced by a 3rd generation platform and produces single-end long reads (e.g. chrRNA-seq), this argument has to be set to \code{FALSE}.\cr
#' \code{strand_mode} depends on strand orientation of the library preparation protocol. This argument is used only if \code{ngs_mode == "PE"}.\cr
#' \code{filter_tlen} requires to skip paired-end alignments with template length (9th field in SAM/BAM files) below the value provided. Most nascent RNA-seq protocols include an RNA fragmentation step,
#' thus unrealistically long inserts most probably originate either from sample contamination with mature RNA species, or from alignment errors.
#' @return A \code{GRangesList} object, if \code{mode %in% c("long_reads", "tss", "pas"))}. A \code{GRanges} object, if \code{mode == "nascent")}.
#' @export
load_BAM_files <- function(bamfiles, mode, ngs_mode = "SE", skip_split_aln = TRUE, strand_mode = 2, filter_tlen = 500) {
  stopifnot(mode %in% c("long_read", "tss", "pas", "nascent"))
  stopifnot(class(bamfiles) == "character" && length(bamfiles) > 0)
  if (mode == "long_read") {
    out <- bamfiles %>% BiocGenerics::lapply(load_long_read_BAM_file) %>% BiocGenerics::do.call(c, .) %>% sort_grl()
  } else if (mode == "nascent") {
    out <- bamfiles %>% BiocGenerics::lapply(load_nascent_BAM_file, ngs_mode = ngs_mode, skip_split_aln = skip_split_aln, strand_mode = strand_mode, filter_tlen = filter_tlen)
    message("Merging ", length(out), " GRanges...")
    out <- merge_GRanges(out)
  } else {
    if (mode == "tss") {
      out <- BiocGenerics::lapply(bamfiles, load_tag_BAM_file)
    } else {
      out <- BiocGenerics::lapply(bamfiles, load_tag_BAM_file, flip_strands = TRUE)
    }
    names(out) <- bamfiles %>% BiocGenerics::basename() %>% stringr::str_replace(".bam$", "")
    out <- as(out, "CompressedGRangesList")
  }
  return(out)
}

# -----------------------------------------------------------------------------------------------------------------------

extract_first_or_last_element <- function(nested_list, mode = "first") {
  unlisted <- BiocGenerics::unlist(nested_list)
  idx <- rep(1:length(nested_list), times = s4Vectors::elementNROWS(nested_list))
  if (mode == "first") {
    choose <- !BiocGenerics::duplicated(idx)
  } else if (mode == "last") {
    choose <- !BiocGenerics::duplicated(idx, fromLast = TRUE)
  }
  return(unlisted[choose])
}

# ---------------------------------------------------------------------------------------------------------------------

filter_ga_by_terminal_subalignments <- function(ga, abs_threshold, rel_threshold) {
  grl <- GenomicRanges::grglist(ga)
  spl <- S4Vectors::elementNROWS(grl) >= 2
  grl_spl <- grl[spl]
  unal <- GenomicRanges::psetdiff(BiocGenerics::unlist(GenomicRanges::range(grl_spl)), grl_spl)
  w_aln <- BiocGenerics::width(grl_spl)
  w_unal <- BiocGenerics::width(unal)
  first_aln <- extract_first_or_last_element(w_aln)
  last_aln <- extract_first_or_last_element(w_aln, mode = "last")
  first_unal <- extract_first_or_last_element(w_unal)
  last_unal <- extract_first_or_last_element(w_unal, mode = "last")
  if (!is.null(abs_threshold)) {
    good_abs <- first_aln >= abs_threshold & last_aln >= abs_threshold
    good_unspl <- BiocGenerics::width(BiocGenerics::unlist(grl[!spl])) >= abs_threshold
  } else {
    good_abs <- rep(TRUE, sum(spl))
    good_unspl <- rep(TRUE, sum(!spl))
  }
  if (!is.null(rel_threshold)) {
    good_rel <- (first_aln / first_unal >= rel_threshold) & (last_aln / last_unal >= rel_threshold)
  } else (
    good_rel <- rep(TRUE, length(grl_spl))
  )
  good_spl <- good_abs & good_rel
  final_good <- rep(FALSE, length(ga))
  final_good[!spl][good_unspl] <- TRUE
  final_good[spl][good_spl] <- TRUE
  return(ga[final_good])
}

# ----------------------------------------------------------------------------------------------------------------------

skip_duplicated_reads <- function(ga) {
  stopifnot(!is.null(names(ga)))
  nms <- names(ga)
  idx <- BiocGenerics::which(BiocGenerics::duplicated(nms) | BiocGenerics::duplicated(nms, fromLast = TRUE))
  if (length(idx) > 0) {
    message("\tSkipped ", length(idx), " chimeric alignments with duplicated names (", length(BiocGenerics::unique(nms[idx])), " reads);")
    ga <- ga[-idx]
  }
  return(ga)
}

# -----------------------------------------------------------------------------------------------------------------------

load_long_read_BAM_file <- function(bamfile, skip_bad_alignments = TRUE, abs_threshold = 30, rel_threshold = 0.05) {
  message("Loading long read BAM file ", BiocGenerics::basename(bamfile), "...")
  ga <- GenomicAlignments::readGAlignments(bamfile, use.names = TRUE)
  ga <- ga %>% skip_duplicated_reads() %>% GenomeInfoDb::sortSeqlevels() %>% BiocGenerics::sort()
  if (isTRUE(skip_bad_alignments)) {
    ga <- filter_ga_by_terminal_subalignments(ga, abs_threshold = abs_threshold, rel_threshold = rel_threshold)
  }
  return(GenomicRanges::grglist(ga))
}

# ---------------------------------------------------------------------------------------------------------------------

load_tag_BAM_file <- function(bamfile, flip_strands = FALSE) {
  message("Loading tag BAM file ", BiocGenerics::basename(bamfile), "...")
  ga <- GenomicAlignments::readGAlignments(bamfile) %>% GenomeInfoDb::sortSeqlevels()
  gr <- ga %>% GenomicRanges::granges() %>% GenomicRanges::resize(1, "start")
  if (isTRUE(flip_strands)) {
    BiocGenerics::strand(gr) <- ifelse(BiocGenerics::strand(gr) == "+", "-", "+")
  }
  cov_fw <- gr[BiocGenerics::strand(gr) == "+"] %>% GenomicRanges::coverage() %>% GenomicRanges::bindAsGRanges(score = .)
  cov_rev <- gr[BiocGenerics::strand(gr) == "-"] %>% GenomicRanges::coverage() %>% GenomicRanges::bindAsGRanges(score = .)
  BiocGenerics::strand(cov_fw) <- "+"
  BiocGenerics::strand(cov_rev) <- "-"
  out <- c(cov_fw, cov_rev)
  out <- out[BiocGenerics::score(out) > 0]
  return(BiocGenerics::sort(out))
}

# ----------------------------------------------------------------------------------------------------------------------

load_nascent_BAM_file <- function(bamfile, ngs_mode, skip_split_aln, strand_mode, filter_tlen) {
  stopifnot(ngs_mode %in% c("SE", "PE"))
  stopifnot(strand_mode %in% 1:2)
  stopifnot(is.logical(skip_split_aln) && length(skip_split_aln) == 1)
  stopifnot(is.null(filter_tlen) || (is.numeric(filter_tlen) && length(filter_tlen) == 1 && filter_tlen %% 1 == 0 && filter_tlen >= 0))
  message("Loading nascent BAM file ", BiocGenerics::basename(bamfile), "...")
  if (ngs_mode == "SE") {
    data <- GenomicAlignments::readGAlignments(bamfile)
  } else if (ngs_mode == "PE") {
    data <- GenomicAlignments::readGAlignmentPairs(bamfile, strandMode = strand_mode)
    message("\tFiltering by template length...")
    if (!is.null(filter_tlen)) {
      data <- data[BiocGenerics::width(GenomicRanges::granges(data)) <= filter_tlen]
    }
  }
  if (isTRUE(skip_split_aln)) {
    message("\tSkipping split alignments...")
    data <- data[GenomicAlignments::njunc(data) == 0]
  }
  message("\tConverting GAlignments to GRanges...")
  data <- GenomicRanges::granges(data)
  BiocGenerics::strand(data) <- ifelse(BiocGenerics::strand(data) == "+", "-", "+")
  message("\tConverting reads to coverage...")
  data <- convert_GRanges_to_coverage(data)
  message("\tSorting GRanges...")
  data <- data %>% GenomeInfoDb::sortSeqlevels() %>% BiocGenerics::sort()
  return(data)
}

