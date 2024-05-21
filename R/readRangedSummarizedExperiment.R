#' Read a RangedSummarizedExperiment from disk
#'
#' Read a \linkS4class{RangedSummarizedExperiment} from its on-disk representation.
#' This is usually not directly called by users, but is instead called by dispatch in \code{\link{readObject}}.
#'
#' @param path String containing a path to a directory, itself created using the \code{\link{saveObject}} method for \linkS4class{RangedSummarizedExperiment} objects.
#' @param metadata Named list of metadata for this object, see \code{\link{readObjectFile}} for details.
#' @param ... Further arguments passed to \code{\link{readSummarizedExperiment}} and internal \code{\link{altReadObject}} calls.
#' 
#' @return A \linkS4class{RangedSummarizedExperiment} object.
#'
#' @author Aaron Lun
#' @seealso
#' \code{"\link{saveObject,RangedSummarizedExperiment-method}"}, to save the RangedSummarizedExperiment to disk.
#'
#' @examples
#' # Mocking up an experiment:
#' mat <- matrix(rpois(10000, 10), ncol=10)
#' colnames(mat) <- letters[1:10]
#' rownames(mat) <- sprintf("GENE_%i", seq_len(nrow(mat)))
#'
#' gr <- GRanges("chrA", IRanges(1:1000, width=10))
#' se <- SummarizedExperiment(list(counts=mat), rowRanges=gr)
#' se$stuff <- LETTERS[1:10]
#' rowData(se)$blah <- runif(1000)
#' metadata(se)$whee <- "YAY"
#' 
#' tmp <- tempfile()
#' saveObject(se, tmp)
#' readObject(tmp)
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment rowData rowData<- rowRanges<-
#' @import alabaster.base
readRangedSummarizedExperiment <- function(path, metadata, ...) {
    metadata$type <- "summarized_experiment"
    se <- altReadObject(path, metadata=metadata, ...)

    rrdir <- file.path(path, "row_ranges")
    if (file.exists(rrdir)) {
        # Avoid overriding the old rowData with the rowRanges's mcols.
        rr <- altReadObject(rrdir, ...)
        old.rd <- rowData(se)
        old.names <- rownames(se)
        rowRanges(se) <- rr
        rowData(se) <- old.rd
        rownames(se) <- old.names
    } else {
        se <- as(se, "RangedSummarizedExperiment")
    }

    se
}
