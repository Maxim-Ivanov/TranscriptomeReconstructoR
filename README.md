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
BiocManager::install(c("GenomicRanges", "GenomicAlignments", "rtracklayer"))
```
Linux users may need to resolve some additional dependencies to install `devtools`. For example, on Ubuntu 20.04:

```{bash, eval = FALSE}
sudo apt-get install libssl-dev libcurl4-openssl-dev libxml2-dev
```
Mac users may need Xcode terminal tools to be installed (otherwise the `brio` dependency may fail to compile). In addition, the `collections` dependency may not compile properly on MacOS X El Capitan 10.11 and should be installed from the binary package on CRAN:

```{r, eval = FALSE}
install.packages("collections", type = "binary")
```

When *TranscriptomeReconstructoR* is finally installed, load it:

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
