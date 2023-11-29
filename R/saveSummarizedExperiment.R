#' Save a SummarizedExperiment to disk
#'
#' Save a \linkS4class{SummarizedExperiment} to its on-disk representation.
#' 
#' @param x A \linkS4class{SummarizedExperiment} object or one of its subclasses.
#' @inheritParams alabaster.base::saveObject
#' @param summarizedexperiment.allow.dataframe.assay Logical scalar indicating whether to allow data frames as assays of \code{x}.
#' @param ... Further arguments to pass to internal \code{\link{altSaveObject}} calls.
#'
#' @return \code{x} is saved into \code{path} and \code{NULL} is invisibly returned.
#'
#' @details
#' If \code{rangedsummarizedexperiment.skip.empty.ranges=TRUE} and \code{x} is a RangedSummarizedExperiment with \dQuote{empty} \code{\link{rowRanges}} (i.e., a \linkS4class{GRangesList} with zero-length entries),
#' \code{stageObject} will save it to file without any genomic range information.
#' This means that any subsequent \code{\link{loadObject}} on the staged assets will return a non-ranged SummarizedExperiment.
#'
#' By default, we consider the presence of data frames in the assays to be an error.
#' Users should coerce these into an appropriate matrix type, e.g., a dense matrix or a sparse dgCMatrix.
#' If a DataFrame as an assay is truly desired, users may set \code{\link{options}(alabaster.se.reject_data.frames=FALSE)} to skip the error.
#' Note that this only works for \linkS4class{DataFrame} objects - data.frame objects will not be saved correctly.
#'
#' @author Aaron Lun
#' @seealso
#' \code{\link{readSummarizedExperiment}}, to read the SummarizedExperiment back into the R session.
#' 
#' @examples
#' mat <- matrix(rpois(10000, 10), ncol=10)
#' colnames(mat) <- letters[1:10]
#' rownames(mat) <- sprintf("GENE_%i", seq_len(nrow(mat)))
#'
#' se <- SummarizedExperiment(list(counts=mat))
#' se$stuff <- LETTERS[1:10]
#' rowData(se)$blah <- runif(1000)
#' metadata(se)$whee <- "YAY"
#' 
#' tmp <- tempfile()
#' saveObject(se, tmp)
#' list.files(tmp, recursive=TRUE)
#' 
#' @export
#' @aliases stageObject,SummarizedExperiment-method
#' @name saveSummarizedExperiment
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom S4Vectors make_zero_col_DFrame
#' @importFrom jsonlite toJSON
#' @import alabaster.base
#' @import methods
setMethod("saveObject", "SummarizedExperiment", function(x, path, summarizedexperiment.allow.dataframe.assay=FALSE, ...) {
    dir.create(path)
    write(file=file.path(path, "OBJECT"), "summarized_experiment")
    write(toJSON(list(dimensions=dim(x), version="1.0"), auto_unbox=TRUE), file=file.path(path, "summarized_experiment.json"))
    args <- list(summarizedexperiment.allow.dataframe.assay=summarizedexperiment.allow.dataframe.assay, ...)

    cd <- colData(x)
    empty.cd <- make_zero_col_DFrame(nrow(cd))
    if (!identical(cd, empty.cd)) { # respect row names, metadata, mcols...
        tryCatch({
            do.call(altSaveObject, c(list(cd, file.path(path, "column_data")), args))
        }, error=function(e) {
            stop("failed to stage 'colData(<", class(x)[1], ">)'\n  - ", e$message)
        })
    }

    rd <- rowData(x)
    empty.rd <- make_zero_col_DFrame(nrow(rd))
    if (!identical(rd, empty.rd)) { # respect row names, metadata, mcols...
        tryCatch({
            do.call(altSaveObject, c(list(rd, file.path(path, "row_data")), args))
        }, error=function(e) {
            stop("failed to stage 'rowData(<", class(x)[1], ">)'\n  - ", e$message)
        })
    }

    adir <- file.path(path, "assays")
    dir.create(adir)
    ass.names <- assayNames(x)
    if (anyDuplicated(ass.names)) {
        stop("assays should be uniquely named")
    }
    write(toJSON(ass.names), file=file.path(adir, "names.json"))

    for (i in seq_along(ass.names)) {
        aname <- as.character(i - 1L)
        curmat <- assay(x, i, withDimnames=FALSE)

        if (is.data.frame(curmat) || (is(curmat, "DataFrame") && !summarizedexperiment.allow.dataframe.assay)) {
            stop("assays should not contain data frames, see ?'saveObject,SummarizedExperiment-method'")
        }

        tryCatch({
            do.call(altSaveObject, c(list(curmat, file.path(adir, aname)), args))
        }, error=function(e) {
            stop("failed to stage 'assay(<", class(x)[1], ">, ", i, ")'\n  - ", e$message)
        })
    }

    saveMetadata(x, metadata.path=file.path(path, "other_data"), mcols.path=NULL)
})

##################################
######### OLD STUFF HERE #########
##################################

#' @export
setMethod("stageObject", "SummarizedExperiment", function(x, dir, path, child=FALSE, meta.name="experiment.json", ...) {
    dir.create(file.path(dir, path), showWarnings=FALSE)

    cd <- colData(x)
    empty.cd <- make_zero_col_DFrame(nrow(cd))
    cd.info <- NULL
    if (!identical(cd, empty.cd)) {
        cd.info <- tryCatch({
            info <- .stageObject(cd, dir, file.path(path, "coldata"), child=TRUE)
            list(resource=.writeMetadata(info, dir=dir))
        }, error=function(e) {
            stop("failed to stage 'colData(<", class(x)[1], ">)'\n  - ", e$message)
        })
    }

    ass.info <- .stage_assays(x, dir, path)
    meta.info <- .processMetadata(x, dir, path, "metadata")

    rd <- rowData(x)
    empty.rd <- make_zero_col_DFrame(nrow(rd))
    rd.info <- NULL
    if (!identical(rd, empty.rd)) {
        rd.info <- tryCatch({
            info <- .stageObject(rd, dir, file.path(path, "rowdata"), child=TRUE)
            list(resource=.writeMetadata(info, dir=dir))
        }, error=function(e) {
            stop("failed to stage 'rowData(<", class(x)[1], ">)'\n  - ", e$message)
        })
    }

    list(
        `$schema`="summarized_experiment/v1.json",
        path=file.path(path, meta.name),
        summarized_experiment=list(
            assays=ass.info,
            column_data=cd.info,
            row_data=rd.info,
            other_data=meta.info,
            dimensions=dim(x)
        ),
        is_child=child
    )
})

#' @importMethodsFrom alabaster.matrix stageObject
#' @importFrom SummarizedExperiment assay assayNames
.stage_assays <- function(x, dir, path) {
    ass.names <- assayNames(x)
    if (length(ass.names) == 0L || anyDuplicated(ass.names)) {
        stop("SummarizedExperiment should contain a non-zero number of uniquely named assays")
    }

    all.meta <- vector("list", length(ass.names))

    for (i in seq_along(all.meta)) {
        curmat <- assay(x, i, withDimnames=FALSE)
        if (is.data.frame(curmat) || (is(curmat, "DataFrame") && getOption("alabaster.se.reject_data.frames", TRUE))) {
            stop("assays should not contain data frames, see ?'stageObject,SummarizedExperiment-method'")
        }

        mat.path <- file.path(path, paste0("assay-", i))
        deets <- tryCatch({
            meta <- .stageObject(curmat, path=mat.path, dir=dir, child=TRUE)
            .writeMetadata(meta, dir=dir)
        }, error=function(e) stop("failed to stage 'assay(<", class(x)[1], ">, ", i, ")'\n  - ", e$message))

        all.meta[[i]] <- list(name=ass.names[i], resource=deets)
    }

    all.meta
}

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
