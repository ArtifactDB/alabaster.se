#' Load a SummarizedExperiment
#'
#' Default loading of \linkS4class{SummarizedExperiment}s based on the metadata stored by the corresponding \code{\link{stageObject}} method.
#'
#' @param exp.info Named list containing the metadata for this experiment.
#' @param project Any argument accepted by the acquisition functions, see \code{?\link{acquireFile}}. 
#' By default, this should be a string containing the path to a staging directory.
#' 
#' @return A \linkS4class{SummarizedExperiment} or \linkS4class{RangedSummarizedExperiment} object.
#'
#' @author Aaron Lun
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
#' # Staging it:
#' tmp <- tempfile()
#' dir.create(tmp)
#' info <- stageObject(se, dir=tmp, "rna-seq") 
#'
#' # And loading it back in:
#' loadSummarizedExperiment(info, tmp)
#'
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment assays<-
#' @import alabaster.base
loadSummarizedExperiment <- function(exp.info, project) {
    all.assays <- list()
    for (y in seq_along(exp.info$summarized_experiment$assays)) {
        cur.ass <- exp.info$summarized_experiment$assays[[y]]
        aname <- cur.ass$name
        apath <- cur.ass$resource$path
        ass.info <- acquireMetadata(project, apath)
        all.assays[[aname]] <- .loadObject(ass.info, project=project)
    }

    cd.info <- acquireMetadata(project, exp.info$summarized_experiment$column_data$resource$path)
    cd <- .loadObject(cd.info, project=project)
    rd.info <- acquireMetadata(project, exp.info$summarized_experiment$row_data$resource$path)
    rd <- .loadObject(rd.info, project=project)

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
