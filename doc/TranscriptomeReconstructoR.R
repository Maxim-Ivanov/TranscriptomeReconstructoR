## ----knits_options, include = FALSE-------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## ----load_packages_for_knitr, include = FALSE, eval = FALSE-------------------
#  library(GenomicRanges)
#  library(GenomicAlignments)
#  library(rtracklayer)
#  library(tidyverse)
#  library(collections)
#  devtools::load_all(".")

## ----devtools, eval = FALSE---------------------------------------------------
#  if (!"devtools" %in% rownames(installed.packages())) {
#    install.packages("devtools")
#  }
#  devtools::install_github("Maxim-Ivanov/TranscriptomeReconstructoR",
#                           build_vignettes = TRUE)

## ----setup_1, eval = FALSE----------------------------------------------------
#  library(TranscriptomeReconstructoR)

## ----filenames----------------------------------------------------------------
#  pkg <- "TranscriptomeReconstructoR"
#  drs_bamfiles <- system.file("extdata",
#                              paste0("ont_rep", 1:4, "_filt.bam"),
#                              package = pkg)
#  tss_bamfiles <- system.file("extdata",
#                              paste0("cage_rep", 1:3, "_filt.bam"),
#                              package = pkg)
#  pas_bamfiles <- system.file("extdata",
#                              paste0("pat_rep", 1:3, "_filt.bam"),
#                              package = pkg)
#  nascent_bamfiles <- system.file("extdata",
#                                  paste0("planet_rep", 1:2, "_filt.bam"),
#                                  package = pkg)

## ----setup_2, eval = FALSE----------------------------------------------------
#  library(rtracklayer)
#  library(magrittr)

## ----load_bam-----------------------------------------------------------------
#  long_reads <- load_BAM_files(drs_bamfiles, mode = "long_read")
#  tss_data <- load_BAM_files(tss_bamfiles, mode = "tss")
#  pas_data <- load_BAM_files(pas_bamfiles, mode = "pas")
#  nascent_data <- load_BAM_files(nascent_bamfiles,
#                                 mode = "nascent", ngs_mode = "PE")

## ----export_tracks, eval = FALSE----------------------------------------------
#  write_grl_as_bed12(long_reads, "Long_reads.bed")
#  tss_data %>% merge_GRanges() %>%
#    save_GRanges_as_bedGraph("TSS_data_merged.bedgraph.gz")
#  pas_data %>% merge_GRanges() %>%
#    save_GRanges_as_bedGraph("PAS_data_merged.bedgraph.gz")
#  save_GRanges_as_bedGraph(nascent_data,
#                           "Nascent_data_merged.bedgraph.gz")

## ----chrom_sizes, eval = FALSE------------------------------------------------
#  nascent_data %>% seqlengths() %>% as_tibble(rownames = "chr") %>%
#    write_tsv("chrom_sizes", col_names = FALSE)

## ----call_tc------------------------------------------------------------------
#  tss <- call_TCs(tss_data)
#  pas <- call_TCs(pas_data)

## ----export_tc, eval = FALSE--------------------------------------------------
#  rtracklayer::export(tss, "TSS.bed", format = "BED")
#  rtracklayer::export(pas, "PAS.bed", format = "BED")

## ----extend_long_reads--------------------------------------------------------
#  long_reads_2 <- extend_long_reads_to_TSS_and_PAS(long_reads, tss, pas)

## ----adjust_exons-------------------------------------------------------------
#  long_reads_3 <- adjust_exons_of_long_reads(long_reads_2)

## ----detect_errors------------------------------------------------------------
#  long_reads_4 <- detect_alignment_errors(long_reads_3)

## ----call_tx------------------------------------------------------------------
#  out <- call_transcripts_and_genes(long_reads_4)
#  hml_genes <- out[[1]]
#  hml_tx <- out[[2]]
#  fusion_genes <- out[[3]]
#  fusion_tx <- out[[4]]
#  reads_free <- out[[5]]

## ----call_planet--------------------------------------------------------------
#  trans <- call_transcribed_intervals(nascent_data)
#  transcribed <- trans[[1]]
#  gaps <- trans[[2]]

## ----add_planet---------------------------------------------------------------
#  results <- process_nascent_intervals(hml_genes, transcribed,
#                                       tss, pas, reads_free, gaps)
#  hml_genes_v2 <- results[[1]]
#  hml_genes_v2_RT <- results[[2]]
#  lncrna <- results[[3]]

## ----export_results, eval = FALSE---------------------------------------------
#  rtracklayer::export(hml_genes_v2, "Called_genes.bed", format = "BED")
#  rtracklayer::export(hml_genes_v2_RT,
#                      "Called_genes_with_RT_tails.bed", format = "BED")
#  write_grl_as_bed12(hml_tx, "Called_transcripts.bed")
#  rtracklayer::export(fusion_genes, "Fusion_genes.bed", format = "BED")
#  write_grl_as_bed12(fusion_tx, "Fusion_transcripts.bed")
#  rtracklayer::export(lncrna, "Called_lncRNAs.bed", format = "BED")

## ----refine_tx, warning = FALSE-----------------------------------------------
#  library(TxDb.Athaliana.BioMart.plantsmart28)
#  txdb <- TxDb.Athaliana.BioMart.plantsmart28
#  annot_exons <- exonsBy(txdb, by = "tx")
#  
#  ref <- refine_transcripts_by_annotation(hml_tx, annot_exons,
#                                          tss, pas, fusion_tx)
#  hml_genes <- ref[[1]]
#  hml_tx <- ref[[2]]
#  fusion_genes <- ref[[3]]
#  fusion_tx <- ref[[4]]

