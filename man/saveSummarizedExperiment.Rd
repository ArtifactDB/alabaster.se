% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saveSummarizedExperiment.R
\name{saveSummarizedExperiment}
\alias{saveSummarizedExperiment}
\alias{stageObject,SummarizedExperiment-method}
\alias{saveObject,SummarizedExperiment-method}
\title{Save a SummarizedExperiment to disk}
\usage{
\S4method{saveObject}{SummarizedExperiment}(x, path, SummarizedExperiment.allow.dataframe.assay = FALSE, ...)
}
\arguments{
\item{x}{A \linkS4class{SummarizedExperiment} object or one of its subclasses.}

\item{path}{String containing the path to a directory in which to save \code{x}.}

\item{SummarizedExperiment.allow.dataframe.assay}{Logical scalar indicating whether to allow data frames as assays of \code{x}.}

\item{...}{Further arguments to pass to internal \code{\link{altSaveObject}} calls.}
}
\value{
\code{x} is saved into \code{path} and \code{NULL} is invisibly returned.
}
\description{
Save a \linkS4class{SummarizedExperiment} to its on-disk representation.
}
\details{
By default, we consider the presence of data frames in the assays to be an error.
Users should coerce these into an appropriate matrix type, e.g., a dense matrix or a sparse dgCMatrix.
If a DataFrame as an assay is truly desired, users may set \code{\link{options}(alabaster.se.reject_data.frames=FALSE)} to skip the error.
Note that this only works for \linkS4class{DataFrame} objects - data.frame objects will not be saved correctly.
}
\examples{
mat <- matrix(rpois(10000, 10), ncol=10)
colnames(mat) <- letters[1:10]
rownames(mat) <- sprintf("GENE_\%i", seq_len(nrow(mat)))

se <- SummarizedExperiment(list(counts=mat))
se$stuff <- LETTERS[1:10]
rowData(se)$blah <- runif(1000)
metadata(se)$whee <- "YAY"

tmp <- tempfile()
saveObject(se, tmp)
list.files(tmp, recursive=TRUE)

}
\seealso{
\code{\link{readSummarizedExperiment}}, to read the SummarizedExperiment back into the R session.
}
\author{
Aaron Lun
}
