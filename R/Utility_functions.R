#' Save stranded GRanges object as bedGraph file with strand info encoded in the sign of values
#'
#' @param gr \code{GRanges} object with a numeric column in mcols.
#' @param filepath Valid path to the output file.
#' @param colname Character (must match the name of the numeric column).
#' @section Details:
#' This function saves a gzipped bedGraph file which contains both positive and negative values
#' (to visualize the coverage of forward and reverse strands in a genome browser).
#' @export
save_GRanges_as_bedGraph <- function(gr, filepath, colname = "score") {
  stopifnot(class(gr) == "GRanges")
  stopifnot(is.character(filepath) && length(filepath) == 1)
  stopifnot(is.character(colname) && length(colname) == 1 && !is.null(S4Vectors::mcols(gr)[[colname]]))
  S4Vectors::mcols(gr)[[colname]] <- ifelse(BiocGenerics::strand(gr)=="-", -S4Vectors::mcols(gr)[[colname]], S4Vectors::mcols(gr)[[colname]])
  BiocGenerics::strand(gr) <- "*"
  gr <- BiocGenerics::sort(gr)
  if (!stringr::str_detect(filepath, "gz$")) {
    filepath <- stringr::str_c(filepath, ".gz")
  }
  con <- gzfile(filepath, "w")
  writeLines("track type=bedGraph color=0,0,0 altColor=128,128,128", con)
  rtracklayer::export.bedGraph(gr, con)
  close(con)
}

# ------------------------------------------------------------------------------------------------------------------------

#' Save GRangesList object as BED12 file
#'
#' @param grl \code{GRangesList} object.
#' @param filename Valid path to the output file.
#' @param name_in_mcols Logical.
#' @section Details:
#' This function does the same work as \code{rtracklayer::export(format = "BED")},
#' with the added benefit that the \code{S4Vectors::mcols(grl)$score} is properly exported.\cr
#' If \code{name_in_mcols == TRUE}, then the \code{S4Vectors::mcols(grl)$name} is exported as range names (otherwise the \code{names(grl)} is used).
#' @export
write_grl_as_bed12 <- function(grl, filename, name_in_mcols = FALSE) {
  gr <- range(grl) %>% BiocGenerics::unlist(use.names = FALSE)
  grl_unl <- BiocGenerics::unlist(grl, use.names = FALSE)
  idx_fw <- BiocGenerics::lapply(S4Vectors::elementNROWS(grl), function(x) { seq(1, x) }) %>% BiocGenerics::unlist() %>% unname()
  idx_rev <- BiocGenerics::lapply(S4Vectors::elementNROWS(grl), function(x) { seq(x, 1) }) %>% BiocGenerics::unlist() %>% unname()
  idx <- ifelse(BiocGenerics::strand(grl_unl) == "+", idx_fw, idx_rev)
  mc <- grl_unl[idx == 1] %>% S4Vectors::mcols()
  chrom <- GenomeInfoDb::seqnames(gr) %>% as.character()
  chromStart <- BiocGenerics::start(gr)
  chromEnd <- BiocGenerics::end(gr)
  if (isTRUE(name_in_mcols)) {
    if (is.null(mc$name)) {
      name <- "."
    } else {
      name <- mc$name
    }
  } else {
    if (is.null(names(grl))) {
      name <- "."
    } else {
      name <- names(grl)
    }
  }
  if (is.null(mc$score)) {
    score <- "."
  } else {
    score <- mc$score
  }
  strand <- BiocGenerics::strand(gr) %>% as.character()
  if (is.null(mc$thick)) {
    thickStart <- chromStart
    thickEnd <- chromEnd
  } else {
    thickStart <- mc$thick %>% BiocGenerics::start()
    thickEnd <- mc$thick %>% BiocGenerics::end()
  }
  if (is.null(mc$itemRgb)) {
    itemRgb <- "."
  } else {
    itemRgb <- mc$itemRgb
  }
  blockCount <- S4Vectors::elementNROWS(grl)
  blockSizes <- BiocGenerics::width(grl) %>% BiocGenerics::lapply(stringr::str_c, collapse = ",") %>% BiocGenerics::unlist()
  blockStarts <- BiocGenerics::start(grl) %>% `-`(BiocGenerics::start(gr)) %>% BiocGenerics::lapply(stringr::str_c, collapse = ",") %>% BiocGenerics::unlist()
  tbl <- dplyr::tibble(chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts)
  tbl <- dplyr::arrange(tbl, chrom, chromStart, chromEnd)
  readr::write_tsv(tbl, filename, col_names = FALSE)
}

# ----------------------------------------------------------------------------------------------------------------

#' Flip strand orientation of a GenomicRanges object
#'
#' @param gr \code{GRanges} or \code{GRangesList} object
#' @return Object of the same class as input with positive and negative strands flipped
#' @export
flip_strand_info <- function(gr) {
  stopifnot(class(gr) == "GRanges" || BiocGenerics::grepl("GRangesList", class(gr)))
  if (class(gr) == "GRanges") {
    BiocGenerics::strand(gr) <- ifelse(BiocGenerics::strand(gr) == "+", "-", "+")
  } else {
    unl <- BiocGenerics::unlist(gr, use.names = FALSE)
    BiocGenerics::strand(unl) <- ifelse(BiocGenerics::strand(unl) == "+", "-", "+")
    gr <- BiocGenerics::relist(unl, gr)
  }
  return(gr)
}

# -----------------------------------------------------------------------------------------------------------------

merge_GRanges <- function(gr_list) {
  # gr_list is expected to be a list of GRanges with "score" column
  cov_fw <- gr_list %>% BiocGenerics::lapply(function(x) { GenomicRanges::coverage(x[BiocGenerics::strand(x) == "+"], weight = "score") }) %>% BiocGenerics::Reduce(`+`, .) %>% GenomicRanges::bindAsGRanges(score = .)
  cov_rev <- gr_list %>% BiocGenerics::lapply(function(x) { GenomicRanges::coverage(x[BiocGenerics::strand(x) == "-"], weight = "score") }) %>% BiocGenerics::Reduce(`+`, .) %>% GenomicRanges::bindAsGRanges(score = .)
  BiocGenerics::strand(cov_fw) <- "+"
  BiocGenerics::strand(cov_rev) <- "-"
  out <- c(cov_fw, cov_rev)
  out <- out[BiocGenerics::score(out) > 0]
  return(BiocGenerics::sort(out))
}

# -----------------------------------------------------------------------------------------------------------------------

convert_GRanges_to_coverage <- function(gr) {
  # Input: GRanges without score column (e.g. overlapping reads);
  # Output: GRanges with non-overlapping ranges and score;
  cov_fw <- GenomicRanges::coverage(gr[BiocGenerics::strand(gr) == "+"]) %>% GenomicRanges::bindAsGRanges(score = .)
  cov_rev <- GenomicRanges::coverage(gr[BiocGenerics::strand(gr) == "-"]) %>% GenomicRanges::bindAsGRanges(score = .)
  strand(cov_fw) <- "+"
  strand(cov_rev) <- "-"
  out <- c(cov_fw, cov_rev)
  out <- out[BiocGenerics::score(out) > 0] %>% GenomeInfoDb::sortSeqlevels() %>% BiocGenerics::sort()
  return(out)
}

# ----------------------------------------------------------------------------------------------------------------------

find_uniq_max <- function(x) {
  out <- rep(FALSE, times = length(x))
  idx <- BiocGenerics::which.max(x)
  out[[idx]] <- TRUE
  return(out)
}

# ------------------------------------------------------------------------------------------------------------------------

find_uniq_min <- function(x) {
  out <- rep(FALSE, times = length(x))
  idx <- BiocGenerics::which.min(x)
  out[[idx]] <- TRUE
  return(out)
}

# ------------------------------------------------------------------------------------------------------------------------

poverlaps_gr <- function(gr1, gr2, maxgap = 0L, minoverlap = 1L, type = c("any", "start", "end", "within", "equal"), ignore.strand = FALSE) {
  stopifnot(length(gr1) == length(gr2))
  goodSeqnames <- as.character(GenomeInfoDb::seqnames(gr1)) == as.character(GenomeInfoDb::seqnames(gr2))
  if (isTRUE(ignore.strand)) {
    goodStrands <- TRUE
  } else {
    goodStrands <- GenomicRanges:::compatibleStrand(BiocGenerics::strand(gr1), BiocGenerics::strand(gr2))
  }
  goodRanges <- IRanges::poverlaps(IRanges::ranges(gr1), IRanges::ranges(gr2), maxgap, minoverlap, type)
  return(goodSeqnames & goodStrands & goodRanges)
}

# -------------------------------------------------------------------------------------------------------------------------

sort_grl <- function(grl) {
  orig_names <- names(grl)
  names(grl) <- 1:length(grl)
  nms <- grl %>% range() %>% BiocGenerics::unlist() %>% BiocGenerics::sort() %>% names()
  idx <- BiocGenerics::match(nms, names(grl))
  out <- grl[idx]
  names(out) <- orig_names[idx]
  return(out)
}

# -------------------------------------------------------------------------------------------------------------------------

deduplicate_grl <- function(grl, colname = "score") {
  uniq_ranges <- grl %>% BiocGenerics::unlist(use.names = FALSE) %>% BiocGenerics::unique()
  gr <- BiocGenerics::unlist(grl, use.names = FALSE)
  idx <- GenomicRanges::findOverlaps(gr, uniq_ranges, type = "equal", select = "first")
  group <- rep(1:length(grl), times = S4Vectors::elementNROWS(grl))
  id <- S4Vectors::split(idx, group) %>% BiocGenerics::lapply(stringr::str_c, collapse = ",") %>% as.character() %>% unname()
  dupl <- BiocGenerics::duplicated(id)
  out <- grl[!dupl]
  out_unl <- BiocGenerics::unlist(out, use.names = FALSE)
  tbl_1 <- dplyr::tibble(id = unique(id))
  tbl_2 <- suppressMessages(dplyr::tibble(id = id) %>% dplyr::group_by(id) %>% dplyr::summarize(n = n()))
  tbl <- dplyr::left_join(tbl_1, tbl_2, by = "id")
  S4Vectors::mcols(out_unl) <- NULL
  S4Vectors::mcols(out_unl)[[colname]] <- rep(tbl$n, times = S4Vectors::elementNROWS(out))
  out <- BiocGenerics::relist(out_unl, out)
  return(out)
}

# ----------------------------------------------------------------------------------------------------------------

check_gr_up_down <- function(g1, g2, offset = 1) {
  # Returns character vector of "up", "down" or "over"
  # "up": g2 range starts upstream from g1 range
  # "over": g2 range starts in the "grey zone" (g1 start + offset)
  # "down": g2 range starts downstream from g1 range
  stopifnot(class(g1) == "GRanges")
  stopifnot(class(g2) == "GRanges")
  stopifnot(length(g1) == length(g2))
  grey <- GenomicRanges::resize(g1, offset, "start")
  S4Vectors::mcols(g2)$idx <- 1:length(g2)
  g2 <- GenomicRanges::resize(g2, 1, "start")
  over <- g2 %over% grey
  out_1 <- g2[over]
  if (length(out_1) > 0) {
    S4Vectors::mcols(out_1)$decision <- "over"
  }
  grey <- grey[!over]
  g2 <- g2[!over]
  gap <- GenomicRanges::pgap(grey, g2) %>% GenomicRanges::resize(width(.) + 1, "start")
  S4Vectors::mcols(g2)$decision <- ifelse(gap %over% g2, "down", "up")
  out <- c(out_1, g2) %>% `[`(BiocGenerics::order(S4Vectors::mcols(.)$idx)) %>% S4Vectors::mcols() %>% .$decision
  return(out)
}

# --------------------------------------------------------------------------------------------------------------------

cluster_indexes <- function(mat, show_progress = FALSE) {
  stopifnot(all(mat[, 1] < mat[, 2]))
  if (any(BiocGenerics::duplicated(mat))) {
    mat <- BiocGenerics::unique(mat)
  }
  col1 <- as.character(mat[, 1])
  col2 <- as.character(mat[, 2])
  cln <- 1
  d <- collections::dict()
  d$set(col1[[1]], cln)
  d$set(col2[[1]], cln)
  if (isTRUE(show_progress)) {
    pb <- utils::txtProgressBar(min = 2, max = nrow(mat), style = 3)
  }
  if (nrow(mat) >= 2) {
    for (i in 2:nrow(mat)) {
      if (isTRUE(show_progress)) {
        utils::setTxtProgressBar(pb, i)
      }
      key1 <- col1[[i]]
      key2 <- col2[[i]]
      lgl1 <- d$has(key1)
      lgl2 <- d$has(key2)
      if (!lgl1 && !lgl2) {
        cln <- cln + 1
        d$set(key1, cln)
        d$set(key2, cln)
      } else if (lgl1 && !lgl2) {
        d$set(key2, d$get(key1))
      } else if (!lgl1 && lgl2) {
        d$set(key1, d$get(key2))
      } else {
        val1 <- d$get(key1)
        val2 <- d$get(key2)
        if (val1 != val2) {
          d$set(key1, val1)
          d$set(key2, val1)
          bad <- BiocGenerics::which(BiocGenerics::unlist(d$values()) == val2)
          bad_keys <- BiocGenerics::unlist(d$keys())[bad]
          for (j in seq_along(bad_keys)) {
            d$set(bad_keys[[j]], val1)
          }
        }
      }
    }
  }
  if (isTRUE(show_progress)) {
    close(pb)
  }
  lst <- d$as_list()
  indexes <- as.integer(names(lst))
  clusters <- as.integer(lst)
  out <- S4Vectors::split(indexes, clusters)
  out <- BiocGenerics::lapply(out, sort)
  return(unname(out))
}

# -----------------------------------------------------------------------------------------------------------------

# Decisions are interpreted from the perspective of gr2:
# "up" = gr2 overlaps 5' end of gr1;
# "down" = gr2 overlaps 3' end of gr1;
# "inside" = gr2 is within gr1;
# "contains" = gr2 includes gr1;
# "exact" = gr2 and gr1 are equal;
# "no_up" = no overlap, gr2 is upstream of gr1;
# "no_down" = no overlap, gr2 is downstream of gr1;

parallel_overlap_type <- function(gr1, gr2) {
  stopifnot(length(gr1) == length(gr2)) # gr1 and gr2 are expected to be parallel
  stopifnot(all(as.character(BiocGenerics::strand(gr1)) %in% c("+", "-"))) # All intervals in gr1 are expected to have strand info. Strandedness of gr2 is not taken into account.
  out <- vector("character", length(gr1))
  a <- BiocGenerics::start(gr2) <= BiocGenerics::start(gr1) & BiocGenerics::end(gr2) < BiocGenerics::end(gr1) & BiocGenerics::end(gr2) >= BiocGenerics::start(gr1) # gr2 overlaps the beginning of gr1
  b <- BiocGenerics::start(gr2) >= BiocGenerics::start(gr1) & BiocGenerics::end(gr2) <= BiocGenerics::end(gr1) # gr2 is within (inside of) gr1
  c <- BiocGenerics::start(gr2) > BiocGenerics::start(gr1) & BiocGenerics::end(gr2) >= BiocGenerics::end(gr1) & BiocGenerics::start(gr2) <= BiocGenerics::end(gr1) # gr2 overlaps the end of gr1
  d <- BiocGenerics::start(gr2) <= BiocGenerics::start(gr1) & BiocGenerics::end(gr2) >= BiocGenerics::end(gr1) # gr2 includes (contains) gr1
  e <- BiocGenerics::start(gr2) == BiocGenerics::start(gr1) & BiocGenerics::end(gr2) == BiocGenerics::end(gr1)
  b[e] <- FALSE
  d[e] <- FALSE
  f <- BiocGenerics::end(gr2) < BiocGenerics::start(gr1) # no overlap
  g <- BiocGenerics::start(gr2) > BiocGenerics::end(gr1)
  no_up <- ifelse(BiocGenerics::strand(gr1) == "+", f, g)
  no_down <- ifelse(BiocGenerics::strand(gr1) == "+", g, f)
  out[no_up] <- "no_up"
  out[no_down] <- "no_down"
  up <- ifelse(BiocGenerics::strand(gr1) == "+", a, c)
  down <- ifelse(BiocGenerics::strand(gr1) == "+", c, a)
  out[up] <- "up"
  out[down] <- "down"
  out[b] <- "inside"
  out[d] <- "contains"
  out[e] <- "exact"
  return(as.factor(out))
}

# ------------------------------------------------------------------------------------------------------------------

trim_by_down_or_upstream_features <- function(windows, features, mode, offset = 1, ignore.strand = FALSE) {
  stopifnot(mode %in% c("down", "up"))
  S4Vectors::mcols(windows)$orig_order <- 1:length(windows) # enumerate windows
  if (offset > 0) {
    if (mode == "down") {
      features <- suppressWarnings(GenomicRanges::trim(GenomicRanges::resize(features, width = width(features) + offset, fix = "end"))) # extend features upstream by the offset value
    } else if (mode == "up") {
      features <- suppressWarnings(GenomicRanges::trim(GenomicRanges::resize(features, width = width(features) + offset, fix = "start")))
    }
  }
  over_any <- IRanges::overlapsAny(windows, features, ignore.strand = ignore.strand)
  out1 <- windows[!over_any] # exclude and save windows which do not overlap any feature
  win_rem <- windows[over_any]
  if (length(win_rem) == 0) {
    return(windows)
  }
  hits <- GenomicRanges::findOverlaps(win_rem, features, ignore.strand = ignore.strand)
  if (length(hits) == 0) {
    return(windows)
  }
  win_par <- win_rem[S4Vectors::queryHits(hits)]
  feat_par <- features[S4Vectors::subjectHits(hits)]
  over_type <- parallel_overlap_type(win_par, feat_par) # detect the type of overlap
  if (mode == "down") {
    bad <- as.logical(BiocGenerics::tapply(over_type, list(S4Vectors::mcols(win_par)$orig_order), function(x) { any(x == "up") | any(x == "contains") }))
  } else if (mode == "up") {
    bad <- as.logical(BiocGenerics::tapply(over_type, list(S4Vectors::mcols(win_par)$orig_order), function(x) { any(x == "down") | any(x == "contains") }))
  }
  out2 <- win_rem[bad] # exclude and save windows which overlap any feature in undesired orientation
  out2 <- GenomicRanges::resize(out2, width = 0) # trim such windows to zero width
  out3 <- win_rem[!bad] # continue with windows which overlap at least one downstream or upstream feature
  if (length(out3) == 0) {
    out <- c(out1, out2)
    out <- out[BiocGenerics::order(S4Vectors::mcols(out)$orig_order)]
    S4Vectors::mcols(out)$orig_order <- NULL
    return(out)
  }
  if (mode == "down") {
    out3_zw <- GenomicRanges::resize(out3, width = 0, fix = "start") # resize both windows and features to 1 bp width around start
    feat_zw <- GenomicRanges::resize(features, width = 0, fix = "start")
    if (isTRUE(ignore.strand)) {
      BiocGenerics::strand(feat_zw) <- "*"
    }
    feat_near <- features[GenomicRanges::precede(out3_zw, feat_zw)] # extract the nearest downstream feature start
  } else {
    out3_zw <- GenomicRanges::resize(out3, width = 0, fix = "end")
    feat_zw <- GenomicRanges::resize(features, width = 0, fix = "end")
    if (isTRUE(ignore.strand)) {
      BiocGenerics::strand(feat_zw) <- "*"
    }
    feat_near <- features[GenomicRanges::follow(out3_zw, feat_zw)] # extract the nearest downstream feature start
  }
  out3_trimmed <- GenomicRanges::pgap(out3_zw, feat_near, ignore.strand = ignore.strand) # compute truncated windows as gap between the starts
  S4Vectors::mcols(out3_trimmed) <- S4Vectors::mcols(out3) # restore mcols which were removed in the pgap operation
  out <- c(out1, out2, out3_trimmed) # combine all processed windows
  out <- out[BiocGenerics::order(S4Vectors::mcols(out)$orig_order)] # restore the original order
  S4Vectors::mcols(out)$orig_order <- NULL
  return(out)
}


# ------------------------------------------------------------------------------------------------------

pintersect_punion_on_grl <- function(grl_a, grl_b, mode) {
  stopifnot(mode %in% c("pintersect", "punion"))
  stopifnot(BiocGenerics::grepl("GRangesList", class(grl_a)))
  stopifnot(BiocGenerics::grepl("GRangesList", class(grl_b)))
  stopifnot(length(grl_a) == length(grl_b))
  gr_a <- BiocGenerics::unlist(grl_a, use.names = FALSE)
  S4Vectors::mcols(gr_a)$group <- rep(1:length(grl_a), times = S4Vectors::elementNROWS(grl_a))
  gr_b <- BiocGenerics::unlist(grl_b, use.names = FALSE)
  S4Vectors::mcols(gr_b)$group <- rep(1:length(grl_b), times = S4Vectors::elementNROWS(grl_b))
  hits <- GenomicRanges::findOverlaps(gr_a, gr_b)
  par1 <- gr_a[S4Vectors::queryHits(hits)]
  par2 <- gr_b[S4Vectors::subjectHits(hits)]
  same_group <- S4Vectors::mcols(par1)$group == S4Vectors::mcols(par2)$group
  par1 <- par1[same_group]
  par2 <- par2[same_group]
  hits <- hits[same_group]
  if (mode == "pintersect") {
    gr <- GenomicRanges::pintersect(par1, par2)
  } else {
    gr <- GenomicRanges::punion(par1, par2)
  }
  out <- S4Vectors::split(gr, S4Vectors::mcols(par1)$group) %>% GenomicRanges::reduce()
  return(out)
}


