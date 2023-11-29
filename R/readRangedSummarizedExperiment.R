#' Read a RangedSummarizedExperiment from disk
#'
#' Read a \linkS4class{RangedSummarizedExperiment} from its on-disk representation.
#'
#' @param path String containing a path to a directory, itself created using the \code{\link{stageObject}} method for \linkS4class{RangedSummarizedExperiment} objects.
#' @param ... Further arguments passed to \code{\link{readSummarizedExperiment}} and internal \code{\link{altReadObject}} calls.
#' 
#' @return A \linkS4class{RangedSummarizedExperiment} object.
#'
#' @author Aaron Lun
#' @seealso
#' \code{"\link{saveObject,RangedSummarizedExperiment-method}"}, to save the SummarizedExperiment to disk.
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
#' readRangedSummarizedExperiment(tmp)
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment rowData rowData<-
#' @import alabaster.base
readRangedSummarizedExperiment <- function(path, ...) {
    se <- altReadObject(path, type="summarized_experiment", ...)
    rr <- altReadObject(file.path(path, "row_ranges"), ...)

    # Avoid overriding the old rowData with the rowRanges's mcols.
    old.rd <- rowData(se)
    rowRanges(se) <- rr
    rowData(se) <- old.rd

    se
}
