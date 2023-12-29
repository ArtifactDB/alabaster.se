#' Read a SummarizedExperiment from disk
#'
#' Read a \linkS4class{SummarizedExperiment} from its on-disk representation.
#' This is usually not directly called by users, but is instead called by dispatch in \code{\link{readObject}}.
#'
#' @param path String containing a path to a directory, itself created using the \code{\link{saveObject}} method for \linkS4class{SummarizedExperiment} objects.
#' @param metadata Named list of metadata for this object, see \code{\link{readObjectFile}} for details.
#' @param ... Further arguments passed to internal \code{\link{altReadObject}} calls.
#' 
#' @return A \linkS4class{SummarizedExperiment} object.
#'
#' @author Aaron Lun
#' @seealso
#' \code{"\link{saveObject,SummarizedExperiment-method}"}, to save the SummarizedExperiment to disk.
#'
#' @examples
#' # Mocking up an experiment:
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
#' readObject(tmp)
#'
#' @export
#' @aliases loadSummarizedExperiment
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom jsonlite fromJSON
#' @importFrom S4Vectors make_zero_col_DFrame
#' @import alabaster.base
readSummarizedExperiment <- function(path, metadata, ...) {
    info <- metadata$summarized_experiment

    all.assays <- list()
    names(all.assays) <- character(0)
    if (file.exists(file.path(path, "assays"))) {
        ass.names <- fromJSON(file.path(path, "assays", "names.json"))
        for (y in seq_along(ass.names)) {
            all.assays[[ass.names[y]]] <- altReadObject(file.path(path, "assays", y - 1L), ...)
        }
    }

    cd.path <- file.path(path, "column_data")
    if (file.exists(cd.path)) {
        cd <- altReadObject(cd.path, ...)
    } else {
        cd <- make_zero_col_DFrame(info$dimensions[[2]])
    }

    rd.path <- file.path(path, "row_data")
    if (file.exists(rd.path)) {
        rd <- altReadObject(rd.path, ...)
    } else {
        rd <- make_zero_col_DFrame(info$dimensions[[1]])
    }

    se <- SummarizedExperiment(all.assays, colData=cd, rowData=rd, checkDimnames=FALSE)

    # Need to force the dimnames to match the DFs, because if they're NULL,
    # the dimnames from the assays end up being used instead.
    rownames(se) <- rownames(rd)
    colnames(se) <- rownames(cd)

    readMetadata(se, mcols.path=NULL, metadata.path = file.path(path, "other_data"))
}

##################################
######### OLD STUFF HERE #########
##################################

#' @export
loadSummarizedExperiment <- function(exp.info, project) {
    all.assays <- list()
    for (y in seq_along(exp.info$summarized_experiment$assays)) {
        cur.ass <- exp.info$summarized_experiment$assays[[y]]
        aname <- cur.ass$name
        apath <- cur.ass$resource$path
        ass.info <- acquireMetadata(project, apath)
        all.assays[[aname]] <- .loadObject(ass.info, project=project)
    }

    cd.meta <- exp.info$summarized_experiment$column_data
    if (!is.null(cd.meta)) {
        cd.info <- acquireMetadata(project, cd.meta$resource$path)
        cd <- .loadObject(cd.info, project=project)
    } else {
        cd <- make_zero_col_DFrame(exp.info$summarized_experiment$dimensions[[2]])
    }

    rd.meta <- exp.info$summarized_experiment$row_data
    if (!is.null(rd.meta)) {
        rd.info <- acquireMetadata(project, rd.meta$resource$path)
        rd <- .loadObject(rd.info, project=project)
    } else {
        rd <- make_zero_col_DFrame(exp.info$summarized_experiment$dimensions[[1]])
    }

    range.info <- exp.info$summarized_experiment$row_ranges
    if (is.null(range.info)) {
        se <- SummarizedExperiment(all.assays, colData=cd, rowData=rd, checkDimnames=FALSE)
    } else {
        rr.info <- acquireMetadata(project, range.info$resource$path)
        rr <- .loadObject(rr.info, project=project)
        mcols(rr) <- rd
        se <- SummarizedExperiment(all.assays, colData=cd, rowRanges=rr, checkDimnames=FALSE)
    }

    # Need to force the dimnames to match the DFs, because if they're NULL,
    # the dimnames from the assays end up being used instead.
    rownames(se) <- rownames(rd)
    colnames(se) <- rownames(cd)

    .restoreMetadata(se, NULL, meta.data = exp.info$summarized_experiment$other_data, project = project)
}
