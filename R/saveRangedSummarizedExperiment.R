#' Save a RangedSummarizedExperiment to disk
#'
#' Save a \linkS4class{RangedSummarizedExperiment} to its on-disk representation.
#' 
#' @param x A \linkS4class{RangedSummarizedExperiment} object or one of its subclasses.
#' @inheritParams alabaster.base::saveObject
#' @param ... Further arguments to pass to \code{"\link{saveObject,SummarizedExperiment-method}"} and internal \code{\link{altSaveObject}} calls.
#'
#' @return \code{x} is saved into \code{path} and \code{NULL} is invisibly returned.
#'
#' @author Aaron Lun
#' @seealso
#' \code{\link{readRangedSummarizedExperiment}}, to read the RangedSummarizedExperiment back into the R session.
#' 
#' @examples
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
#' list.files(tmp, recursive=TRUE)
#' 
#' @aliases 
#' stageObject,RangedSummarizedExperiment-method
#' @name saveRangedSummarizedExperiment
NULL

#' @export
#' @export
#' @rdname saveRangedSummarizedExperiment
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom S4Vectors mcols<-
setMethod("saveObject", "RangedSummarizedExperiment", function(x, path, ...) {
    callNextMethod()

    if (!emptyRowRanges(x)) {
        rr <- rowRanges(x)
        mcols(rr) <- NULL # removing mcols as these are absorbed into the SE's rowData.
        names(rr) <- NULL # names are also absorbed into the SE's rowData.

        tryCatch({
            altSaveObject(rr, file.path(path, "row_ranges"), ...) 
        }, error=function(e) {
            stop("failed to stage 'rowRanges(<", class(x)[1], ">)'\n  - ", e$message)
        })
    }

    write(toJSON(list(dimensions=dim(x), version="1.0"), auto_unbox=TRUE), file=file.path(path, "ranged_summarized_experiment.json"))
    write(file=file.path(path, "OBJECT"), "ranged_summarized_experiment")
    invisible(NULL)
})

##################################
######### OLD STUFF HERE #########
##################################

#' @export
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom alabaster.base .stageObject .writeMetadata
#' @importMethodsFrom alabaster.ranges stageObject
#' @import methods
setMethod("stageObject", "RangedSummarizedExperiment", function(x, dir, path, child=FALSE, ..., skip.ranges=FALSE) {
    dir.create(file.path(dir, path), showWarnings=FALSE)
    meta <- callNextMethod()

    if (!skip.ranges && !emptyRowRanges(x)) {
        rd.processed <- tryCatch({
            rd.info <- .stageObject(rowRanges(x), dir, file.path(path, "rowranges"), mcols.name=NULL, child=TRUE) # skipping the mcols, as this is the row_data.
            .writeMetadata(rd.info, dir=dir)
        }, error=function(e) {
            stop("failed to stage 'rowRanges(<", class(x)[1], ">)'\n  - ", e$message)
        })

        meta$summarized_experiment$row_ranges <- list(resource=rd.processed)
    }

    meta
})
