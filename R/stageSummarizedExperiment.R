#' Stage a SummarizedExperiment
#'
#' Save a \linkS4class{SummarizedExperiment} to file inside the staging directory.
#' 
#' @param x A \linkS4class{SummarizedExperiment} object or one of its subclasses.
#' @inheritParams alabaster.base::stageObject
#' @param meta.name String containing the name of the metadata file.
#' @param ... Further arguments to pass to the \linkS4class{SummarizedExperiment} method.
#' For the SummarizedExperiment itself, all further arguments are just ignored.
#' @param skip.ranges Logical scalar indicating whether to avoid saving the \code{\link{rowRanges}}.
#'
#' @return A named list of metadata that follows the \code{summarized_experiment} schema.
#' The contents of \code{x} are saved into a \code{path} subdirectory inside \code{dir}.
#'
#' @details
#' \code{meta.name} is only needed to set up the output \code{path}, for consistency with the \code{\link{stageObject}} contract.
#' Callers should make sure to write the metadata to the same path by using \code{\link{.writeMetadata}} to create the JSON file.
#'
#' If \code{skip.ranges=TRUE}, the RangedSummarizedExperiment method just calls the SummarizedExperiment method, i.e., \code{\link{rowRanges}} are not saved.
#' This avoids the hassle of switching classes and the associated problems, e.g., \url{https://github.com/Bioconductor/SummarizedExperiment/issues/29}.
#' Note that any subsequent \code{\link{loadObject}} call on the staged assets will return a non-ranged SummarizedExperiment.
#'
#' If \code{x} is a RangedSummarizedExperiment with \dQuote{empty} \code{\link{rowRanges}} (i.e., a \linkS4class{GRangesList} with zero-length entries),
#' \code{stageObject} will save it to file without any genomic range information.
#' This means that any subsequent \code{\link{loadObject}} on the staged assets will return a non-ranged SummarizedExperiment.
#'
#' By default, we consider the presence of data frames in the assays to be an error.
#' Users should coerce these into an appropriate matrix type, e.g., a dense matrix or a sparse dgCMatrix.
#' If a DataFrame as an assay is truly desired, users may set \code{\link{options}(alabaster.se.reject_data.frames=FALSE)} to skip the error.
#' Note that this only works for \linkS4class{DataFrame} objects - data.frame objects will not be saved correctly.
#'
#' @author Aaron Lun
#' 
#' @examples
#' tmp <- tempfile()
#' dir.create(tmp)
#'
#' mat <- matrix(rpois(10000, 10), ncol=10)
#' colnames(mat) <- letters[1:10]
#' rownames(mat) <- sprintf("GENE_%i", seq_len(nrow(mat)))
#'
#' se <- SummarizedExperiment(list(counts=mat))
#' se$stuff <- LETTERS[1:10]
#' rowData(se)$blah <- runif(1000)
#' metadata(se)$whee <- "YAY"
#' 
#' dir.create(tmp)
#' stageObject(se, dir=tmp, "rna-seq") 
#' list.files(file.path(tmp, "rna-seq"))
#' 
#' @export
#' @rdname stageSummarizedExperiment
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom S4Vectors make_zero_col_DFrame
#' @import alabaster.base
#' @import methods
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
#' @importFrom SummarizedExperiment assay assays assayNames
.stage_assays <- function(x, dir, path) {
    ass.names <- assayNames(x)
    if (is.null(ass.names) && length(assays(x)) > 0) {
        stop("assays should be named in a ", class(x)[1], " object")
    }
    if (anyDuplicated(ass.names)) {
        stop("detected duplicate assay names in a ", class(x)[1], " object")
    }
    if (any(ass.names == "")) {
        stop("detected empty assay name in a ", class(x)[1], " object")
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
#' @rdname stageSummarizedExperiment
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
