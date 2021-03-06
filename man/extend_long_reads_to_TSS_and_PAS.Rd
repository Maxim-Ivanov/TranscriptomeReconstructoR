% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Extend_long_reads.R
\name{extend_long_reads_to_TSS_and_PAS}
\alias{extend_long_reads_to_TSS_and_PAS}
\title{Adjust outer borders of full-length RNA-seq reads to the nearby TSS and PAS}
\usage{
extend_long_reads_to_TSS_and_PAS(
  grl,
  tss,
  pas,
  read_flanks_up = c(-100, 10),
  read_flanks_down = c(-50, 50),
  trunc_limit = 20,
  mode = "reads",
  verbose = isTRUE(mode == "reads")
)
}
\arguments{
\item{grl}{\code{GRangesList} object.}

\item{tss, pas}{\code{GRanges} objects.}

\item{read_flanks_up, read_flanks_down}{Integer vectors of length 2.}

\item{trunc_limit}{Positive integer.}

\item{mode}{Character. Must be either "reads" (default) or "transcripts".}

\item{verbose}{Logical.}
}
\value{
\code{GRangesList} object containing the adjusted exons + some extra mcols.
}
\description{
Adjust outer borders of full-length RNA-seq reads to the nearby TSS and PAS
}
\section{Details}{

If \code{mode == "reads"}, then \code{grl} is expected to be full-length RNA-seq reads returned by \code{load_BAM_files(mode = "long_read")}.\cr
If \code{mode == "transcripts"}, then \code{grl} is expected to be known transcripts returned by \code{exonsBy(txdb, by = "tx")} (where \code{txdb} is a \code{GenomicFeatures} object).\cr
\code{tss} and \code{pas} are returned by \code{call_TCs()} on 5'- and 3'-tag sequencing data, respectively.\cr
The function searches for called TSS within \code{read_flanks_up} bp around 5' ends of long reads (or known transcripts), and for PAS within \code{read_flanks_down} around their 3' ends.\cr
When neighboring TSS/PAS are found, the long reads are extended towards their summits.\cr
If multiple TSS or PAS are found within the search window, then the one with the highest score is chosen.\cr
The long read may be truncated, if the TSS/PAS summit is located inside of its terminal exon (however truncation to a length smaller than \code{trunc_limit} bp is disallowed).\cr
}

