#' Adjust internal exons of full-length RNA-seq reads by majority vote
#'
#' @param long_reads \code{GRangesList} object.
#' @param max_exon_diff Positive integer.
#' @section Details:
#' The input \code{GRangesList} object is returned by the \code{extend_long_reads_to_TSS_and_PAS()} function.\cr
#' 5'- and 3'-borders of exons in \code{long_reads} are clustered by similarity (within \code{max_exon_diff} bp).
#' Then the exon borders are changed to the median value of genomic coordinates observed within the cluster.
#' This suppresses the "fuzziness" of subalignments which are often observed in full-length RNA-seq reads.
#' @return \code{GRangesList} object with the same order of elements and same mcols as the input object.
#' @export
adjust_exons_of_long_reads <- function(long_reads, max_exon_diff = 10) {
  stopifnot(BiocGenerics::grepl("GRangesList", class(long_reads)))
  stopifnot(is.numeric(max_exon_diff) && length(max_exon_diff) == 1 && max_exon_diff %% 1 == 0 && max_exon_diff > 0)
  all_exons <- BiocGenerics::unlist(long_reads, use.names = FALSE)
  message(length(all_exons), " exons in total (", length(BiocGenerics::unique(GenomicRanges::granges(all_exons))), " unique);")
  # Adjust exon starts:
  message("Adjusting exon starts by majority vote...")
  all_exons_adj <- adjust_splice_sites(all_exons, mode = "start", max_exon_diff = max_exon_diff)
  # Adjust exon ends:
  message("Adjusting exon ends by majority vote...")
  all_exons_adj <- adjust_splice_sites(all_exons_adj, mode = "end", max_exon_diff = max_exon_diff)
  message(sum(BiocGenerics::start(all_exons_adj) != BiocGenerics::start(all_exons) | BiocGenerics::end(all_exons_adj) != BiocGenerics::end(all_exons)), " exons were adjusted at either start or end;")
  message(length(BiocGenerics::unique(GenomicRanges::granges(all_exons))), " -> ", length(BiocGenerics::unique(GenomicRanges::granges(all_exons_adj))), " unique exons after adjustment;")
  out <- S4Vectors::split(all_exons_adj, S4Vectors::mcols(all_exons_adj)$read_id)
  return(out)
}

# ------------------------------------------------------------------------------------------------------------------------

adjust_splice_sites <- function(exons, mode, max_exon_diff) {
  stopifnot(mode %in% c("start", "end"))
  exons <- exons[BiocGenerics::order(S4Vectors::mcols(exons)$exon_id)] # just in case...
  # Exclude borders of terminal exons from the adjustment procedure:
  if (mode == "start") {
    bad <- S4Vectors::mcols(exons)$pos_fw == 1
  } else {
    bad <- S4Vectors::mcols(exons)$pos_rev == 1
  }
  bad_exons <- exons[bad]
  message("\t", length(bad_exons), " terminal exons were excluded (", length(BiocGenerics::unique(GenomicRanges::granges(bad_exons))), " unique);")
  good_exons <- exons[!bad]
  message("\t", length(good_exons), " exons were used (", length(BiocGenerics::unique(GenomicRanges::granges(good_exons))), " unique);")
  borders <- GenomicRanges::resize(good_exons, 0, mode)
  # Cluster exon borders within <max_exon_diff>:
  clusters <- GenomicRanges::reduce(borders, min.gapwidth = max_exon_diff, with.revmap = TRUE)
  revmap <- S4Vectors::mcols(clusters)$revmap
  good_exons <- good_exons[BiocGenerics::unlist(revmap)]
  borders <- borders[BiocGenerics::unlist(revmap)]
  grp <- rep(1:length(revmap), times = S4Vectors::elementNROWS(revmap))
  vote <- BiocGenerics::tapply(BiocGenerics::start(borders), grp, median, simplify = FALSE) %>% as.integer() # median or mode?
  vote_ext <- rep(vote, times = S4Vectors::elementNROWS(revmap))
  vote_gr <- GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(good_exons), IRanges::IRanges(start = vote_ext, width = 0), strand = BiocGenerics::strand(good_exons), seqinfo = GenomeInfoDb::seqinfo(good_exons))
  IRanges::ranges(good_exons) <- GenomicRanges::punion(GenomicRanges::resize(good_exons, 0, ifelse(mode == "start", "end", "start")), vote_gr, fill.gap = TRUE) %>% IRanges::ranges()
  out <- c(good_exons, bad_exons) %>% `[`(BiocGenerics::order(S4Vectors::mcols(.)$exon_id))
  message("\t", sum(BiocGenerics::start(out) != BiocGenerics::start(exons) | BiocGenerics::end(out) != BiocGenerics::end(exons)), " exon ", mode, "s were adjusted;")
  return(out)
}

