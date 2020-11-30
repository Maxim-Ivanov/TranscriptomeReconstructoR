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
vignette("TranscriptomeReconstructoR", package = "TranscriptomeReconstructoR")
```
