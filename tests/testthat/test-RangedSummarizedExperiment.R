# This tests the RangedSummarizedExperiment staging.
# library(testthat); library(alabaster.se); source("test-RangedSummarizedExperiment.R")

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
    expect_error(alabaster.base::.writeMetadata(out, tmp), NA)
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

    # The new world works.
    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_identical(colData(se), colData(out2))
    expect_identical(rowRanges(se), rowRanges(out2))
    expect_true(file.exists(file.path(tmp, "row_ranges", "OBJECT")))
    expect_false(file.exists(file.path(tmp, "row_ranges", "range_annotations")))
})

test_that("saveObject preserves RSE rownames", {
    copy <- se
    rownames(copy) <- sprintf("GENE_%i", seq_len(nrow(copy)))

    tmp <- tempfile()
    saveObject(copy, tmp)
    out2 <- readObject(tmp)

    expect_identical(rownames(out2), rownames(copy))
    expect_identical(rowRanges(out2), rowRanges(copy))
    expect_identical(rowData(out2), rowData(copy))
})

test_that("stageObject auto-skips on empty rowRanges", {
    tmp <- tempfile()
    dir.create(tmp)

    # GRL is empty.
    copy <- GRangesList(rep(list(GRanges()), nrow(se)))
    names(copy) <- rownames(se)
    rowRanges(se) <- copy
    expect_true(emptyRowRanges(se))

    out <- stageObject(se, tmp, "rnaseq", skip.ranges=NA)
    expect_error(alabaster.base::.writeMetadata(out, tmp), NA)

    expect_null(out$summarized_experiment$row_ranges)
    expect_null(out$summarized_experiment$row_data)

    out2 <- loadSummarizedExperiment(out, tmp)
    expect_identical(as.character(class(out2)), "SummarizedExperiment")

    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_s4_class(out2, "RangedSummarizedExperiment")
    expect_true(emptyRowRanges(out2))
    expect_false(file.exists(file.path(tmp, "row_ranges")))

    # Non-empty rowData but GRL is still empty.
    mcols(copy)$FOO <- 2
    rowRanges(se) <- copy

    out <- stageObject(se, tmp, "rnaseq2", skip.ranges=NA)
    expect_error(alabaster.base::.writeMetadata(out, tmp), NA)

    expect_null(out$summarized_experiment$row_ranges)
    mpath <- out$summarized_experiment$row_data$resource$path
    expect_match(jsonlite::fromJSON(file.path(tmp, paste0(mpath, ".json")))[["$schema"]], "csv_data_frame")

    out2 <- loadSummarizedExperiment(out, tmp)
    expect_identical(as.character(class(out2)), "SummarizedExperiment")
    expect_identical(unique(rowData(out2)$FOO), 2)

    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_s4_class(out2, "RangedSummarizedExperiment")
    expect_true(emptyRowRanges(out2))
    expect_identical(unique(rowData(out2)$FOO), 2)
    expect_false(file.exists(file.path(tmp, "row_ranges")))
})

test_that("stageObject allows us to forcibly skip the ranges", {
    tmp <- tempfile()
    dir.create(tmp)

    out <- stageObject(se, tmp, "rnaseq", skip.ranges=TRUE)
    expect_error(alabaster.base::.writeMetadata(out, tmp), NA)
    expect_true("row_data" %in% names(out$summarized_experiment))
    expect_false("row_ranges" %in% names(out$summarized_experiment))

    out2 <- loadSummarizedExperiment(out, tmp)
    expect_identical(as.character(class(out2)), "SummarizedExperiment")
})

rowRanges(se) <- splitAsList(rowRanges(se), seq_len(nrow(se)))
rownames(se) <- paste0("FEATURE_", seq_len(nrow(se)))

test_that("stageObject handles GRLs", {
    tmp <- tempfile()
    dir.create(tmp)

    out <- stageObject(se, tmp, "rnaseq")
    expect_error(alabaster.base::.writeMetadata(out, tmp), NA)
    rrmeta <- jsonlite::fromJSON(file.path(tmp, paste0(out$summarized_experiment$row_ranges$resource$path, ".json")), simplifyVector=FALSE)
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

    # The new world works.
    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_identical(rowRanges(se), ref)
})
