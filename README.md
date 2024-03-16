# Save `SummarizedExperiment`s to file

|Environment|Status|
|---|---|
|[BioC-release](https://bioconductor.org/packages/release/bioc/html/alabaster.se.html)|[![Release OK](https://bioconductor.org/shields/build/release/bioc/alabaster.se.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/alabaster.se/)|
|[BioC-devel](https://bioconductor.org/packages/devel/bioc/html/alabaster.se.html)|[![Devel OK](https://bioconductor.org/shields/build/devel/bioc/alabaster.se.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/alabaster.se/)|

The **alabaster.se** package implements methods for saving and loading `SummarizedExperiment` objects under the **alabaster** framework.
It provides a language-agnostic method for serializing experimental data and annotations in these objects, including the genomic coordinates in a `RangedSummarizedExperiment`.
To get started, install the package and its dependencies from Bioconductor:

```r
# install.packages("BiocManager")
BiocManager::install("alabaster.se")
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
saveObject(rse, tmp)

roundtrip <- readObject(tmp)
class(roundtrip)
## [1] "RangedSummarizedExperiment"
## attr(,"package")
## [1] "SummarizedExperiment"
```
