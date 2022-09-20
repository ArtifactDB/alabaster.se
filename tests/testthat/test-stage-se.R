# This tests the stageObject generic for base SE's.
# library(testthat); library(artificer.se); source("test-stage-se.R")

# Making an SE and annotating it.
mat <- matrix(rpois(2000, 10), ncol=10)
colnames(mat) <- paste0("SAMPLE_", seq_len(ncol(mat)))

se <- SummarizedExperiment(list(counts=mat, cpm=mat/10))
se$stuff <- LETTERS[1:10]
se$blah <- runif(10)
rowData(se)$whee <- runif(nrow(se))

rownames(se) <- sprintf("GENE_%i", seq_len(nrow(se)))

test_that("stageObject works as expected for SE objects", {
    tmp <- tempfile()
    dir.create(tmp)

    out <- stageObject(se, tmp, "rnaseq")
    expect_identical(vapply(out$summarized_experiment$assays, function(x) x$name, ""), assayNames(se))

    cdf <- read.csv(file.path(tmp, out$summarized_experiment$column_data$resource$path), row.names=1)
    expect_equal(cdf, as.data.frame(colData(se)))

    rdf <- read.csv(file.path(tmp, out$summarized_experiment$row_data$resource$path), row.names=1)
    expect_equal(rdf, as.data.frame(rowData(se)))

    # Doesn't save empty metadata.
    expect_null(out$summarized_experiment$other_data)

    # Round-tripping it to see that it's okay.
    out2 <- loadSummarizedExperiment(out, tmp)
    expect_equal(colData(se), colData(out2))
    expect_equal(rownames(se), rownames(out2))
    expect_equal(rowData(se)[,1], rowData(out2)$whee)
    expect_equal(sum(assay(se)), sum(assay(out2)))
    expect_equal(sum(assay(se, 2)), sum(assay(out2, 2)))

    # Metadata is succcessfully saved.
    .writeMetadata(out, tmp)
    expect_true(file.exists(file.path(tmp, out$path)))
    deets <- jsonlite::fromJSON(file.path(tmp, out$path))
    expect_null(deets$md5sum)
})

test_that("stageObject works as expected with no row or column names", {
    tmp <- tempfile()
    dir.create(tmp)
    dimnames(se) <- NULL

    out <- stageObject(se, tmp, "rnaseq")
    expect_identical(vapply(out$summarized_experiment$assays, function(x) x$name, ""), assayNames(se))

    cdf <- read.csv(file.path(tmp, out$summarized_experiment$column_data$resource$path))
    expect_equal(cdf, as.data.frame(colData(se)))

    rdf <- read.csv(file.path(tmp, out$summarized_experiment$row_data$resource$path))
    expect_equal(rdf, as.data.frame(rowData(se)))

    # Round-tripping it to see that it's okay.
    # This only works when the row names are NULL,
    # as we can't even form the MAE otherwise.
    out2 <- loadSummarizedExperiment(out, tmp)
    expect_equal(colData(se), colData(out2))
    expect_equal(rowData(se), rowData(out2))
})

test_that("stageObject works as expected with no row or column data", {
    tmp <- tempfile()
    dir.create(tmp)

    colData(se) <- colData(se)[,0]
    rowData(se) <- rowData(se)[,0]

    out <- stageObject(se, tmp, "rnaseq")
    expect_identical(vapply(out$summarized_experiment$assays, function(x) x$name, ""), assayNames(se))

    cdf <- read.csv(file.path(tmp, out$summarized_experiment$column_data$resource$path), row.names=1)
    expect_equal(cdf, as.data.frame(colData(se)))

    rdf <- read.csv(file.path(tmp, out$summarized_experiment$row_data$resource$path), row.names=1)
    expect_equal(rdf, as.data.frame(rowData(se)))

    # Round-tripping to make sure it's okay.
    out2 <- loadSummarizedExperiment(out, tmp)
    expect_equal(colData(se), colData(out2))
    expect_equal(rowData(se), rowData(out2))
})

test_that("stageObject works with the various types of vectors", {
    tmp <- tempfile()
    dir.create(tmp)

    se$stuff <- rep(LETTERS[1:3], length.out=ncol(se))
    se$blah <- factor(se$stuff, LETTERS[10:1])
    se$foo <- seq_len(ncol(se))
    se$whee <- as.numeric(10 + seq_len(ncol(se)))
    se$rabbit <- factor(se$stuff, LETTERS[1:3], ordered=TRUE)
    se$birthday <- rep(Sys.Date(), ncol(se)) - sample(100, ncol(se))

    out <- stageObject(se, tmp, "rnaseq")

    meta <- jsonlite::fromJSON(file.path(tmp, paste0(out$summarized_experiment$column_data$resource$path, ".json")), simplifyVector=FALSE)$data_frame
    expect_identical(unlist(meta$columns[[1]]$values), LETTERS[1:3])
    expect_identical(read.csv(file.path(tmp, meta$columns[[2]]$levels$resource$path))[,1], LETTERS[10:1])
    expect_identical(meta$columns[[3]]$type, "integer")
    expect_identical(meta$columns[[4]]$type, "number")
    expect_identical(meta$columns[[5]]$type, "ordered")
    expect_identical(read.csv(file.path(tmp, meta$columns[[5]]$levels$resource$path))[,1], LETTERS[1:3])
    expect_identical(meta$columns[[6]]$type, "date")

    # Round-tripping to make sure it's okay.
    out2 <- loadSummarizedExperiment(out, tmp)
    expect_identical(out2$stuff, se$stuff)
    expect_identical(out2$blah, se$blah)
    expect_identical(out2$foo, se$foo)
    expect_identical(out2$whee, se$whee)
    expect_identical(out2$rabbit, se$rabbit)
    expect_identical(out2$birthday, se$birthday)
})

test_that("stageObject handles data frames in the assays", {
    tmp <- tempfile()
    dir.create(tmp)

    assay(se) <- as.data.frame(assay(se))
    expect_error(stageObject(se, tmp, "df_fail"), "should not contain data.frame")

    opt <- options(artificer.se.reject_data.frames=FALSE)
    on.exit(options(artificer.se.reject_data.frames=TRUE))

    # still doesn't handle data.frame...
    expect_error(stageObject(se, tmp, "df_fail2"), "should not contain data.frame")

    # but now can handle DataFrame.
    assay(se) <- as(assay(se), "DataFrame")
    expect_error(out <- stageObject(se, tmp, "rnaseq"), NA)

    path <- file.path(tmp, out$summarized_experiment$assays[[1]]$resource$path)
    meta <- jsonlite::fromJSON(paste0(path, ".json"), simplifyVector=FALSE)
    expect_match(meta[["$schema"]], "data_frame/")
})

test_that("stageExperiment saves other metadata when necessary", {
    tmp <- tempfile()
    dir.create(tmp)

    metadata(se)$YAY <- 1
    out <- stageObject(se, tmp, "rnaseq")
    expect_true(file.exists(file.path(tmp, out$summarized_experiment$other_data$resource$path)))

    round <- loadSummarizedExperiment(out, tmp)
    expect_identical(metadata(round), list(YAY=1))
})
