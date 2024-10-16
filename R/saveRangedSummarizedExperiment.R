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
setMethod("saveObject", "RangedSummarizedExperiment", 
function(x, path, ...) {
    # Note that the use of callNextMethod() means that we cannot respond to
    # application overrides for the SummarizedExperiment base class. Developers
    # should just pretend that saveObject,RSE-method copied all of the code
    # from saveObject,SE-method; callNextMethod() is an implementation detail.
    #
    # This simplifies the dispatch and ensures that an override is only called
    # once. Consider the alternative - namely, casting to the next subclass
    # and then calling altSaveObject to respect the override. This would call
    # the override's SE method repeatedly for every step from the subclass to
    # SE. If the override's behavior is not idempotent, we have a problem.
    # 
    # So, if an application wants to set an override for all SEs, then it
    # should define an altSaveObject,SE-method and then call it. If the
    # override is slightly different for particular SE subclasses, developers
    # should just duplicate the common override logic in the altSaveObject
    # methods for affected subclasses, rather than expecting some injection of
    # the overriding method into the saveObject dispatch hierarchy.
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

    info <- readObjectFile(path)
    info$ranged_summarized_experiment <- list(version="1.0")
    saveObjectFile(path, "ranged_summarized_experiment", info)

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
