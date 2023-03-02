testthat::test_that(
  "a named list of dgCMatrix objects is returned after successfully reading in 10x h5 files containing Gene Expression data", {
    
    path <- "../testdata/sample1/raw_feature_bc_matrix.h5"
    rawData <- scReports:::read10xH5(path = path)
    
    testthat::expect_type(object = rawData, type = "list")
    
    testthat::expect_identical(object = names(rawData), expected = "RNA")
    
    testthat::expect_s4_class(object = rawData[["RNA"]], class = "dgCMatrix")
    
  }
)


testthat::test_that(
  "a named list of dgCMatrix objects is returned after successfully reading in 10x Genomics h5 files containing both Gene Expression and Antibody data", {
    
    path <- "../testdata/sample2/raw_feature_bc_matrix.h5"
    rawData <- scReports:::read10xH5(path = path)
    
    testthat::expect_type(object = rawData, type = "list")
    
    testthat::expect_identical(object = names(rawData), 
                               expected = c("RNA", "ADT"))
    
    testthat::expect_s4_class(object = rawData[["RNA"]], class = "dgCMatrix")
    
    testthat::expect_s4_class(object = rawData[["ADT"]], class = "dgCMatrix")
    
  }
)


