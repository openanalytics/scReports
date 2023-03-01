samplePath <- c(
  "../testdata/sample2/raw_feature_bc_matrix.h5ad"
)

rawData = readH5AD(path = samplePath)

testthat::test_that(
  "a dgCMatrix object is returned after successfully reading the h5ad file containing RNA and ADT data", {
    
    testthat::expect_type(object = rawData, type = "list")
    
    testthat::expect_s4_class(
      object = rawData[["RNA"]], 
      class = "dgCMatrix"
    )
    
    testthat::expect_s4_class(
      object = rawData[["ADT"]], 
      class = "dgCMatrix"
    )
    
  }
)
