# TranscriptomeReconstructoR

*TranscriptomeReconstructoR* is an R package for data-driven annotation of transcriptomes.

## Installation

Install the development version from GitHub using `devtools`:

```{r, eval = FALSE}
if (!"devtools" %in% rownames(installed.packages())) {
    install.packages("devtools")
}
devtools::install_github("Maxim-Ivanov/TranscriptomeReconstructoR", build_vignettes = TRUE, ref = "main")
```
If the Bioconductor dependencies were not resolved automatically by `devtools::install_github()`, please install the required Bioconductor packages manually:

```{r, eval = FALSE}
BiocManager::install(c("GenomicRanges", "GenomicAlignments", "rtracklayer", "BiocStyle"))
```
Linux users may need to resolve some additional dependencies to install `devtools`. For example, on Ubuntu 20.04:

```{bash, eval = FALSE}
sudo apt-get install libssl-dev libcurl4-openssl-dev libxml2-dev
```

Then load the installed package:

```{r, eval = FALSE}
library(TranscriptomeReconstructoR)
```

## User manual

The vignette contains everything you want to know about TranscriptomeReconstructoR:

* A simple usage scenario (the example dataset is provided);
* Detailed description of the underlying algorithm;
* Considerations regarding the choice and pre-processing of the input datasets.

```{r, eval = FALSE}
browseVignettes("TranscriptomeReconstructoR")
```
or

```{r, eval = FALSE}
vignette("TranscriptomeReconstructoR", package = "TranscriptomeReconstructoR")
```
