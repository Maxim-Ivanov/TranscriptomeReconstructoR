#' Call continuous intervals of nascent transcription from NET-seq or GRO-seq data
#'
#' @param cov \code{GRanges} object.
#' @param min_signal Positive integer.
#' @param max_gapwidth Non-negative integer.
#' @param min_width Positive integer.
#' @section Details:
#' The input file is expected to be returned by the \code{load_BAM_files(mode = "nascent")} function. It contains sequencing coverage
#' of the chromosomes by nascent RNA-seq reads. Continuous intervals of coverage (with values not less than \code{min_signal})
#' are merged if the distance between them does not exceed \code{max_gapwidth} bp.
#' Finally, the merged intervals are filtered to have length not less than \code{min_width}.
#' @return List of length 2:
#' \enumerate{
#'     \item \code{GRanges} object containing merged intervals of nascent transcription.
#'     \item \code{GRanges} object containing low coverage sub-intervals (gaps) within the merged transcribed intervals.
#' }
#' @export
call_transcribed_intervals <- function(cov, min_signal = 3, max_gapwidth = 250, min_width = 500) {
  stopifnot(class(cov) == "GRanges")
  stopifnot(is.numeric(min_signal) && length(min_signal) == 1 && min_signal %% 1 == 0 && min_signal > 0)
  stopifnot(is.numeric(max_gapwidth) && length(max_gapwidth) == 1 && max_gapwidth %% 1 == 0 && max_gapwidth >= 0)
  stopifnot(is.numeric(min_width) && length(min_width) == 1 && min_width %% 1 == 0 && min_width > 0)
  cov_above <- cov[BiocGenerics::score(cov) >= min_signal]
  cov_below <- cov[BiocGenerics::score(cov) < min_signal]
  transcribed <- cov_above %>% GenomicRanges::reduce(min.gapwidth = max_gapwidth) %>% `[`(BiocGenerics::width(.) >= min_width)
  gaps <- cov_below %>% GenomicRanges::reduce() %>% IRanges::subsetByOverlaps(transcribed)
  return(list("transcribed" = transcribed, "gaps" = gaps))
}

# ---------------------------------------------------------------------------------------------------------

#' Classify intervals of nascent transcription and add them to the gene model
#'
#' @param hml_genes,nascent,tss,pas \code{GRanges} objects.
#' @param reads_free \code{GrangesList} object, or NULL.
#' @param gaps \code{GRanges} object, or NULL.
#' @param trim_offset Non-negative integer.
#' @param min_score_2 Non-negative numeric.
#' @param min_lncrna_width Positive integer.
#' @param extend_along_nascent Logical.
#' @param extension_flanks Integer vector of length 2.
#' @section Details:
#' \code{hml_genes} contains the called genes (the first element in the list returned by the \code{call_transcripts_and_genes()} function).
#' \code{nascent} contains continuous intervals of nascent transcription (the first element in the list returned by the \code{call_transcribed_intervals()} function).
#' \code{tss} and \code{pas} are returned by \code{call_TCs()} on 5'- and 3'-tag sequencing data, respectively.
#' \code{gaps} contains intervals of low coverage within the continuous intervals of nascent transcription
#' (the second element in the list returned by the \code{call_transcribed_intervals()} function).
#' /code{reads_free} contains long reads which remain unused during the transcript calling procedure and are located outside of the called genes
#' (the fifth element in the list returned by the \code{call_transcripts_and_genes()} function).\cr
#' The intervals of nascent transcription are first classified by overlap with the called genes (on the same strand). The intervals which start upstream
#' and/or end downstream from a called genes, are considered upstream transcribed intervals and readthrough (RT) tails, respectively.
#' Borders of the called genes are extended to include such gene-associated intervals of nascent transcription. At that, the original gene coordinates
#' (which correspond to the mature RNA molecule) are saved as the "thick" mcols. If \code{extend_along_nascent == TRUE}, then the original ("thick") coordinates
#' of MC and LC genes can be further extended towards strong TSS and/or PAS (with scores not less than \code{min_score_2}) which are found within \code{extension_flanks}
#' windows centered at starts and ends, respectively, of the associated nascent interval.\cr
#' Other nascent intervals which are not associated with any mature transcript, are considered antisense or intergenic lncRNAs. They are filtered to have length
#' not less than \code{min_lncrna_width}.
#' @return List of length 3:
#' \enumerate{
#'     \item \code{GRanges} object which contains intervals covered by mature RNA molecules (i.e. the original coordinates of HC, MC and LC genes,
#'     with the exception that the borders of MC and LC genes could have been extended towards nearby strong TSS and/or PAS).
#'     \item \code{GRanges} object which contains HC, MC and LC genes extended to include the gene-associated intervals of nascent transcription
#'     (the interval of mature transcription moved to the "thick" mcols).
#'     \item \code{GRanges} object which contains lncRNAs called from the nascent-only transcribed intervals.
#' }
#' @export
process_nascent_intervals <- function(hml_genes, nascent, tss, pas, reads_free = NULL, gaps = NULL, trim_offset = 20, min_score_2 = 5, min_lncrna_width = 500,
                                      extend_along_nascent = TRUE, extension_flanks = c("tss" = -50, "pas" = 0)) {
  stopifnot(class(hml_genes) == "GRanges")
  stopifnot(class(nascent) == "GRanges")
  stopifnot(class(tss) == "GRanges")
  stopifnot(class(pas) == "GRanges")
  stopifnot(BiocGenerics::grepl("GRangesList", class(reads_free)) || is.null(reads_free))
  stopifnot(class(gaps) == "GRanges" || is.null(gaps))
  stopifnot(is.numeric(trim_offset) && length(trim_offset) == 1 && trim_offset %% 1 == 0 && trim_offset >= 0)
  stopifnot(is.numeric(min_score_2) && length(min_score_2) == 1 && min_score_2 >= 0)
  stopifnot(is.numeric(min_lncrna_width) && length(min_lncrna_width) == 1 && min_lncrna_width %% 1 == 0 && min_lncrna_width > 0)
  stopifnot(is.logical(extend_along_nascent) && length(extend_along_nascent) == 1)
  stopifnot(is.numeric(extension_flanks) && length(extension_flanks) == 2 && extension_flanks[[1]] %% 1 == 0 && extension_flanks[[2]] %% 1 == 0)
  genic <- GenomicRanges::reduce(hml_genes)
  message("Calling and trimming tails...")
  # Subtract called genes from transcribed intervals:
  nascent_2 <- GenomicRanges::setdiff(nascent, hml_genes)
  # Find nascent tails:
  tails <- nascent_2[suppressWarnings(GenomicRanges::flank(nascent_2, 1)) %over% hml_genes]
  # Trim tails by outer borders of used TSS:
  tss_used <- tss[tss %over% GenomicRanges::resize(genic, 1, "start")]
  tails <- trim_by_down_or_upstream_features(tails, tss_used, mode = "down", offset = trim_offset)
  # Trim tails by strong unused TSS:
  tss_strong <- tss[tss %outside% genic & BiocGenerics::score(tss) >= min_score_2]
  tails <- trim_by_down_or_upstream_features(tails, tss_strong, mode = "down", offset = trim_offset)
  # Subtract tails from the transcribed intervals:
  nascent_3 <- GenomicRanges::setdiff(nascent_2, tails)
  if (!is.null(gaps)) {
    # Also trim tails by the first low coverage interval (without offset):
    tails <- trim_by_down_or_upstream_features(tails, gaps, mode = "down", offset = 0)
  }
  message("Calling and trimming heads...")
  # Find heads (transcribed intervals immediately upstream from called genes):
  heads <- nascent_3[suppressWarnings(GenomicRanges::flank(nascent_3, 1, start = FALSE)) %over% hml_genes]
  # Trim heads by outer borders of used PAS:
  pas_used <- pas[pas %over% GenomicRanges::resize(genic, 1, "end")]
  heads <- trim_by_down_or_upstream_features(heads, pas_used, mode = "up", offset = trim_offset)
  # Trim heads by strong unused PAS:
  pas_strong <- pas[pas %outside% genic & BiocGenerics::score(pas) >= min_score_2]
  heads <- trim_by_down_or_upstream_features(heads, pas_strong, mode = "up", offset = trim_offset)
  # The remaining transcribed intervals can be regarded as lncRNAs:
  lncrna <- GenomicRanges::setdiff(nascent_3, heads) %>% `[`(BiocGenerics::width(.) >= min_lncrna_width)
  if (isTRUE(extend_along_nascent)) {
    message("Extending MC and LC genes along heads and tails towards unused TSS and PAS...")
    # Extend MC and LC genes without TSS along plaNET-seq intervals towards nearby strong TSS (if any):
    wo_tss <- GenomicRanges::resize(hml_genes, 1, "start") %outside% tss
    results <- hml_genes[wo_tss] %>% extend_genes_along_nascent_intervals(tss_strong, heads, mode = "tss", flanks = extension_flanks)
    hml_genes <- c(hml_genes[!wo_tss], results[[1]]) %>% BiocGenerics::sort()
    heads <- results[[2]]
    #tss_strong <- tss_strong[tss_strong %outside% hml_genes] # update the unused strong TSS
    # Extend MC and LC genes without PAS along plaNET-seq intervals towards nearby strong PAS (if any):
    wo_pas <- GenomicRanges::resize(hml_genes, 1, "end") %outside% pas
    results <- hml_genes[wo_pas] %>% extend_genes_along_nascent_intervals(pas_strong, tails, mode = "pas", flanks = extension_flanks)
    hml_genes <- c(hml_genes[!wo_pas], results[[1]]) %>% BiocGenerics::sort()
    tails <- results[[2]]
    #pas_strong <- pas_strong[pas_strong %outside% hml_genes] # update the unused strong PAS
  }
  # Decorate genes with nascent heads and tails:
  message(length(heads), " heads;")
  message(length(tails), " tails;")
  S4Vectors::mcols(hml_genes)$thick <- GenomicRanges::granges(hml_genes) %>% unname()
  hml_genes_with_fl <- hml_genes %>% decorate_genes(heads, mode = "up") %>% decorate_genes(tails, mode = "down")
  if (is.null(reads_free)) {
    # Merge lncRNAs if connected by a free_read:
    range_free <- reads_free %>% range() %>% BiocGenerics::unlist(use.names = FALSE)
    range_over_lncrna <- range_free[GenomicRanges::countOverlaps(range_free, lncrna) >= 2] %>% GenomicRanges::reduce()
    if (length(range_over_lncrna) > 0) {
      hits <- GenomicRanges::findOverlaps(range_over_lncrna, lncrna)
      lncrna_p1 <- lncrna[S4Vectors::subjectHits(hits)] %>% S4Vectors::split(S4Vectors::queryHits(hits)) %>% range() %>% BiocGenerics::unlist() %>% unname()
      lncrna_p2 <- lncrna[-S4Vectors::subjectHits(hits)]
      lncrna <- c(lncrna_p1, lncrna_p2) %>% BiocGenerics::sort()
    }
  }
  # Prepare the output:
  S4Vectors::mcols(lncrna)$score <- 0
  S4Vectors::mcols(lncrna)$type <- "Nascent"
  lncrna_id <- sprintf(paste0("%0", nchar(length(lncrna)), "i"), 1:length(lncrna))
  S4Vectors::mcols(lncrna)$name <- BiocGenerics::paste(S4Vectors::mcols(lncrna)$type, lncrna_id, sep = "_")
  message(length(lncrna), " lncRNAs called;")
  out <- list("hml_genes" = hml_genes, "hml_genes_with_flanks" = hml_genes_with_fl, "lncRNA" = lncrna)
  return(out)
}

# -----------------------------------------------------------------------------------------------------------

extend_genes_along_nascent_intervals <- function(genes, tc, nascent, mode = "tss", flanks) {
  stopifnot(mode %in% c("tss", "pas"))
  # mode == "tss": <nascent> and <TC> are expected to be <heads> and <tss_strong>, respectively;
  # mode == "pas": <nascent> and <TC> are expected to be <tails> and <pas_strong>, respectively;
  # This function searches for unused strong TSS/PAS within the whole head/tail;
  # The <heads> may be additionally extended by <flanks[[1]]>, <tails> = by <flanks[[2]]>;
  h1 <- GenomicRanges::flank(nascent, 1, start = !isTRUE(mode == "tss")) %>% GenomicRanges::findOverlaps(genes)
  if (length(h1) == 0) {
    warning("\tDid not find any ", ifelse(mode == "tss", "heads", "tails"), " among the nascent intervals!")
    return(list("genes" = genes, "nascent" = nascent))
  }
  out_1 <- genes[-S4Vectors::subjectHits(h1)]
  genes <- genes[S4Vectors::subjectHits(h1)]
  win <- nascent[S4Vectors::queryHits(h1)] %>% GenomicRanges::resize(BiocGenerics::width(.) + abs(GenomicRanges::flanks[[mode]]), fix = ifelse(mode == "tss", "end", "start"))
  h2 <- GenomicRanges::findOverlaps(win, tc)
  if (length(h2) == 0) {
    message("\tDid not find any unused strong ", stringr::str_to_upper(mode), " inside of the ", ifelse(mode == "tss", "heads", "tails"), ";")
    return(list("genes" = genes, "nascent" = nascent))
  }
  # In case of multiple TC per window, choose the strongest one:
  scores <- BiocGenerics::score(tc)[S4Vectors::subjectHits(h2)]
  best <- BiocGenerics::tapply(scores, S4Vectors::queryHits(h2), find_uniq_max, simplify = FALSE) %>% BiocGenerics::unlist() %>% unname()
  h2 <- h2[best]
  out_2 <- genes[-S4Vectors::queryHits(h2)]
  genes <- genes[S4Vectors::queryHits(h2)]
  tc_par <- tc[S4Vectors::subjectHits(h2)]
  used <- S4Vectors::queryHits(h1)[S4Vectors::queryHits(h2)]
  # Extend gene border until the TC summit:
  message("\t", length(tc_par), " genes found their missing ", stringr::str_to_upper(mode), ";")
  IRanges::ranges(genes) <- S4Vectors::mcols(tc_par)$thick %>% GenomicRanges::punion(GenomicRanges::resize(genes, 1, fix = ifelse(mode == "tss", "end", "start")), fill.gap = TRUE) %>% IRanges::ranges() %>% unname()
  # Trim the consumed nascent interval:
  nascent_out <- GenomicRanges::setdiff(nascent, nascent[used])
  genes_out <- c(out_1, out_2, genes) %>% BiocGenerics::sort()
  return(list("genes" = genes_out, "nascent" = nascent_out))
}

# --------------------------------------------------------------------------------------------------------------

decorate_genes <- function(genes, planet, mode) {
  stopifnot(mode %in% c("up", "down"))
  if (mode == "up") {
    hits <- suppressWarnings(GenomicRanges::flank(genes, 1) %>% GenomicRanges::trim() %>% GenomicRanges::findOverlaps(planet))
  } else {
    hits <- suppressWarnings(GenomicRanges::flank(genes, 1, start = FALSE) %>% GenomicRanges::trim() %>% GenomicRanges::findOverlaps(planet))
  }
  out_1 <- genes[-S4Vectors::queryHits(hits)]
  par1 <- genes[S4Vectors::queryHits(hits)]
  par2 <- planet[S4Vectors::subjectHits(hits)]
  if (mode == "up") {
    IRanges::ranges(par1) <- GenomicRanges::punion(GenomicRanges::resize(par2, 1, "start"), GenomicRanges::resize(par1, 1, "end"), fill.gap = TRUE) %>% IRanges::ranges() %>% unname()
  } else {
    IRanges::ranges(par1) <- GenomicRanges::punion(GenomicRanges::resize(par1, 1, "start"), GenomicRanges::resize(par2, 1, "end"), fill.gap = TRUE) %>% IRanges::ranges() %>% unname()
  }
  out <- c(out_1, par1) %>% BiocGenerics::sort()
  return(out)
}

