.onLoad <- function(libname, pkgname) {
    registerReadObjectFunction("summarized_experiment", readSummarizedExperiment)
    registerReadObjectFunction("ranged_summarized_experiment", readRangedSummarizedExperiment)
}

.onUnload <- function(libname, pkgname) {
    registerReadObjectFunction("summarized_experiment", NULL)
    registerReadObjectFunction("ranged_summarized_experiment", NULL)
}
