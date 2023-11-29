# This tests the stageObject generic for base SE's.
# library(testthat); library(alabaster.se); source("test-SummarizedExperiment.R")

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

    # Works in the new world.
    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_equal(colData(se), colData(out2))
    expect_equal(rownames(se), rownames(out2))
    expect_equal(rowData(se)[,1], rowData(out2)$whee)
    expect_equal(sum(assay(se)), sum(assay(out2)))
    expect_equal(sum(assay(se, 2)), sum(assay(out2, 2)))
})

test_that("stageObject works as expected with no row or column names", {
    tmp <- tempfile()
    dir.create(tmp)
    dimnames(se) <- NULL

    out <- stageObject(se, tmp, "rnaseq")
    expect_error(alabaster.base::.writeMetadata(out, tmp), NA)
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

    # Works in the new world.
    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_equal(colData(se), colData(out2))
    expect_equal(rowData(se), rowData(out2))
})

test_that("stageObject works as expected with no row or column data, but still names", {
    tmp <- tempfile()
    dir.create(tmp)

    colData(se) <- colData(se)[,0]
    rowData(se) <- rowData(se)[,0]

    out <- stageObject(se, tmp, "rnaseq")
    expect_error(alabaster.base::.writeMetadata(out, tmp), NA)
    expect_identical(vapply(out$summarized_experiment$assays, function(x) x$name, ""), assayNames(se))

    cdf <- read.csv(file.path(tmp, out$summarized_experiment$column_data$resource$path), row.names=1)
    expect_equal(cdf, as.data.frame(colData(se)))

    rdf <- read.csv(file.path(tmp, out$summarized_experiment$row_data$resource$path), row.names=1)
    expect_equal(rdf, as.data.frame(rowData(se)))

    # Round-tripping to make sure it's okay.
    out2 <- loadSummarizedExperiment(out, tmp)
    expect_equal(colData(se), colData(out2))
    expect_equal(rowData(se), rowData(out2))

    # Works in the new world.
    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_equal(colData(se), colData(out2))
    expect_equal(rowData(se), rowData(out2))
})

test_that("stageObject works as expected with no row or column data at all", {
    tmp <- tempfile()
    dir.create(tmp)

    colData(se) <- colData(se)[,0]
    rowData(se) <- rowData(se)[,0]
    dimnames(se) <- list(NULL, NULL)

    out <- stageObject(se, tmp, "rnaseq")
    expect_error(alabaster.base::.writeMetadata(out, tmp), NA)
    expect_identical(vapply(out$summarized_experiment$assays, function(x) x$name, ""), assayNames(se))

    expect_null(out$summarized_experiment$column_data)
    expect_null(out$summarized_experiment$row_data)

    # Round-tripping to make sure it's okay.
    out2 <- loadSummarizedExperiment(out, tmp)
    expect_equal(colData(se), colData(out2))
    expect_equal(rowData(se), rowData(out2))

    # Works in the new world.
    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_equal(colData(se), colData(out2))
    expect_equal(rowData(se), rowData(out2))
    expect_false(file.exists(file.path(tmp, "row_data")))
    expect_false(file.exists(file.path(tmp, "column_data")))
})

test_that("stageObject fails when the assay names are NULL or non-unique", {
    tmp <- tempfile()
    dir.create(tmp)
    ass <- assays(se)
    names(ass) <- NULL
    assays(se) <- ass
    expect_error(stageObject(se, tmp, "rnaseq"), "non-zero number")

    tmp <- tempfile()
    dir.create(tmp)
    assayNames(se) <- rep("FOO", length(assayNames(se)))
    expect_error(stageObject(se, tmp, "rnaseq"), "non-zero number")

    # Fails in the new world.
    tmp <- tempfile()
    expect_error(saveObject(se, tmp), "uniquely named")
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
    expect_error(alabaster.base::.writeMetadata(out, tmp), NA)

    meta <- jsonlite::fromJSON(file.path(tmp, paste0(out$summarized_experiment$column_data$resource$path, ".json")), simplifyVector=FALSE)$data_frame
    expect_identical(meta$columns[[1]]$type, "string")
    expect_identical(read.csv(file.path(tmp, meta$columns[[2]]$levels$resource$path))[,1], LETTERS[10:1])
    expect_identical(meta$columns[[3]]$type, "integer")
    expect_identical(meta$columns[[4]]$type, "number")
    expect_identical(meta$columns[[5]]$type, "factor")
    expect_identical(read.csv(file.path(tmp, meta$columns[[5]]$levels$resource$path))[,1], LETTERS[1:3])
    expect_identical(meta$columns[[6]]$type, "string")

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

    opt <- options(alabaster.se.reject_data.frames=FALSE)
    on.exit(options(alabaster.se.reject_data.frames=TRUE))

    # still doesn't handle data.frame...
    expect_error(stageObject(se, tmp, "df_fail2"), "should not contain data.frame")

    # but now can handle DataFrame.
    se2 <- se
    assay(se2) <- as(assay(se2), "DataFrame")
    expect_error(out <- stageObject(se2, tmp, "rnaseq"), NA)
    expect_error(alabaster.base::.writeMetadata(out, tmp), NA)

    path <- file.path(tmp, out$summarized_experiment$assays[[1]]$resource$path)
    meta <- jsonlite::fromJSON(paste0(path, ".json"), simplifyVector=FALSE)
    expect_match(meta[["$schema"]], "data_frame/")

    # Works in the new world.
    tmp <- tempfile()
    expect_error(saveObject(se, tmp), "should not contain data frames")

    tmp <- tempfile()
    saveObject(se2, tmp, summarizedexperiment.allow.dataframe.assay=TRUE)
    out2 <- readObject(tmp)
    expect_identical(assay(se2), assay(out2))
})

test_that("stageExperiment saves other metadata when necessary", {
    tmp <- tempfile()
    dir.create(tmp)

    metadata(se)$YAY <- 1
    out <- stageObject(se, tmp, "rnaseq")
    expect_error(alabaster.base::.writeMetadata(out, tmp), NA)
    expect_true(file.exists(file.path(tmp, out$summarized_experiment$other_data$resource$path)))

    round <- loadSummarizedExperiment(out, tmp)
    expect_identical(metadata(round), list(YAY=1))

    # Works in the new world.
    tmp <- tempfile()
    saveObject(se, tmp)
    out2 <- readObject(tmp)
    expect_identical(metadata(out2), list(YAY=1))
})
