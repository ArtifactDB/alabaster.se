#' Is the \code{rowRanges} empty?
#'
#' Check the \code{\link{rowRanges}} of a \linkS4class{RangedSummarizedExperiment} is empty, 
#' i.e., a \linkS4class{GRangesList} with no ranges. 
#' 
#' @param x A \linkS4class{RangedSummarizedExperiment} object or the contents of its \code{\link{rowRanges}}.
#'
#' @return A logical scalar indicating whether \code{x} has empty \code{rowRanges}.
#'
#' @details
#' Metadata in \code{\link{mcols}} is ignored for the purpose of this discussion, 
#' as this can be moved to the \code{\link{rowData}(x)} of the base \linkS4class{SummarizedExperiment} class without loss.
#' In other words, non-empty \code{\link{mcols}} will not be used to determine that the \code{rowRanges} is not empty.
#' However, non-empty fields in the \code{\link{metadata}} or in the inner \code{\link{mcols}} of the \linkS4class{GRanges} will trigger a non-emptiness decision.
#' 
#' @export
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom IRanges PartitioningByEnd
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors mcols<- mcols
#' @importFrom BiocGenerics relist
emptyRowRanges <- function(x) {
    if (is(x, "SummarizedExperiment")) {
        x <- rowRanges(x)
    }
    if (is(x, "GRangesList")) {
        # Creating an empty GRL and comparing it. This re-uses the same
        # code as the SE->RSE coerce method, more or less.
        partitioning <- PartitioningByEnd(integer(length(x)), names = names(x))
        rowRanges <- relist(GRanges(), partitioning)
        mcols(x) <- mcols(rowRanges, use.names = FALSE)
        identical(x, rowRanges)
    } else {
        FALSE
    }
}

