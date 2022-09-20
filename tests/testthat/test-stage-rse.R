# This tests the RangedSummarizedExperiment staging.
# library(testthat); library(alabaster.se); source("test-stage-rse.R")

set.seed(100)

# Making an SE and annotating it.
mat <- matrix(rpois(2000, 10), ncol=10)
colnames(mat) <- paste0("SAMPLE_", seq_len(ncol(mat)))

se <- SummarizedExperiment(list(counts=mat, cpm=mat/10), rowRanges=GRanges("chrA", IRanges(1:200, width=1)))
se$stuff <- LETTERS[1:10]
se$blah <- runif(10)
rowData(se)$whee <- runif(nrow(se))

test_that("stageObject works as expected for RSE objects", {
    tmp <- tempfile()
    dir.create(tmp)

    out <- stageObject(se, tmp, "rnaseq")
    expect_identical(vapply(out$summarized_experiment$assays, function(x) x$name, ""), assayNames(se))

    # Make sure that no mcols were saved for the ranges.
    rr.meta <- acquireMetadata(tmp, out$summarized_experiment$row_ranges$resource$path)
    rr.round <- loadObject(rr.meta, tmp)
    expect_identical(ncol(mcols(rr.round)), 0L)

    # Round trip works.
    out2 <- loadSummarizedExperiment(out, tmp)
    expect_equal(colData(se), colData(out2))
    expect_equal(rowRanges(se), rowRanges(out2))
    expect_identical(length(rowData(out2)$whee), nrow(se))
})

test_that("stageObject auto-skips on empty rowRanges", {
    tmp <- tempfile()
    dir.create(tmp)

    copy <- GRangesList(rep(list(GRanges()), nrow(se)))
    names(copy) <- rownames(se)
    rowRanges(se) <- copy
    out <- stageObject(se, tmp, "rnaseq", skip.ranges=NA)

    mpath <- out$summarized_experiment$row_data$resource$path
    expect_match(jsonlite::fromJSON(file.path(tmp, paste0(mpath, ".json")))[["$schema"]], "csv_data_frame")

    out2 <- loadSummarizedExperiment(out, tmp)
    expect_identical(as.character(class(out2)), "SummarizedExperiment")
})

test_that("stageObject allows us to forcibly skip the ranges", {
    tmp <- tempfile()
    dir.create(tmp)

    out <- stageObject(se, tmp, "rnaseq", skip.ranges=TRUE)
    expect_match(jsonlite::fromJSON(file.path(tmp, paste0(out$summarized_experiment$row_data$resource$path, ".json")))[["$schema"]], "csv_data_frame")

    out2 <- loadSummarizedExperiment(out, tmp)
    expect_identical(as.character(class(out2)), "SummarizedExperiment")
})

rowRanges(se) <- splitAsList(rowRanges(se), seq_len(nrow(se)))
rownames(se) <- paste0("FEATURE_", seq_len(nrow(se)))

test_that("stageObject handles GRLs", {
    tmp <- tempfile()
    dir.create(tmp)

    out <- stageObject(se, tmp, "rnaseq")
    rrmeta <- jsonlite::fromJSON(file.path(tmp, paste0(out$summarized_experiment$row_data$resource$path, ".json")), simplifyVector=FALSE)
    expect_true(rrmeta$compressed_list$names)

    grouping <- read.csv(file.path(tmp, rrmeta$path), row.names=1)
    expect_true(all(grouping$number==1L))
    expect_identical(rownames(grouping), rownames(se))

    # Again, make sure that no mcols were saved for the ranges.
    rr.meta <- acquireMetadata(tmp, out$summarized_experiment$row_ranges$resource$path)
    rr.round <- loadObject(rr.meta, tmp)
    expect_identical(ncol(mcols(rr.round)), 0L)

    # Round trip works.
    out2 <- loadSummarizedExperiment(out, tmp)
    ref <- rowRanges(se)
    metadata(ref) <- list()
    metadata(ref@unlistData) <- list()
    expect_equal(rowRanges(out2), ref)
    expect_s4_class(rowRanges(out2), "GenomicRangesList")
})
