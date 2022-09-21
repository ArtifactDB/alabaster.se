# Save `SummarizedExperiment`s to file

The **alabaster.se** package implements methods for saving and loading `SummarizedExperiment` objects under the **alabaster** framework.
It provides a language-agnostic method for serializing experimental data and annotations in these objects, including the genomic coordinates in a `RangedSummarizedExperiment`.
To get started, install the package and its dependencies from GitHub:

```r
devtools::install_github("ArtifactDB/alabaster.schemas")
devtools::install_github("ArtifactDB/alabaster.base")
devtools::install_github("ArtifactDB/alabaster.ranges")
devtools::install_github("ArtifactDB/alabaster.matrix")
devtools::install_github("ArtifactDB/alabaster.se")
```

In the example below, we save a `RangedSummarizedExperiment` object to file:

```r
library(SummarizedExperiment)
example(SummarizedExperiment, echo=FALSE) # can't be bothered to copy it here.
rse
## class: RangedSummarizedExperiment
## dim: 200 6
## metadata(0):
## assays(1): counts
## rownames: NULL
## rowData names(1): feature_id
## colnames(6): A B ... E F
## colData names(1): Treatment

library(alabaster.se)
tmp <- tempfile()
dir.create(tmp)
meta <- stageObject(se, tmp, "se")
meta[["$schema"]]
## [1] "summarized_experiment/v1.json"

roundtrip <- loadObject(meta, tmp)
class(roundtrip)
## [1] "RangedSummarizedExperiment"
## attr(,"package")
## [1] "SummarizedExperiment"
```
