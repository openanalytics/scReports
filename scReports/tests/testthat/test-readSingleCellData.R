sampleIDs <- sampleNames <- c(
  "sample1_h5",
  "sample1_h5ad",
  "sample2_h5",
  "sample2_h5ad",
  "sample4_h5ad",
  "sample3_h5",
  "sample3_csv",
  "sample3_wrong_path"
)
samplePaths <- c(
  "../testdata/sample1/raw_feature_bc_matrix.h5",
  "../testdata/sample1/raw_feature_bc_matrix.h5ad",
  "../testdata/sample2/raw_feature_bc_matrix.h5",
  "../testdata/sample2/raw_feature_bc_matrix.h5ad",
  "../testdata/sample4/raw_feature_bc_matrix.h5",
  "../testdata/sample3/raw_feature_bc_matrix.h5",
  "../testdata/sample3/metrics_summary.csv",
  "../testdata/missing/raw_feature_bc_matrix.h5"
)


testthat::test_that(
  "an error is returned if the lengths of sampleNames and samplePaths differ", {
    
    testthat::expect_error(
      object = readSingleCellData(sampleNames = sampleNames, 
                                  samplePaths = samplePaths[1:2], 
                                  verbose = FALSE)
    )
    
  }
)


testthat::test_that(
  "an error is returned if samplePaths do not exist", {
    
    testthat::expect_error(
      object = readSingleCellData(sampleNames = sampleNames[8], 
                                  samplePaths = samplePaths[8], 
                                  verbose = FALSE)
    )
    
  }
)


testthat::test_that(
  "an error is returned if samplePaths have an extension other than h5 or h5ad", {
    
    testthat::expect_error(
      object = readSingleCellData(sampleNames = sampleNames[7], 
                                  samplePaths = samplePaths[7], 
                                  verbose = FALSE)
    )
    
  }
)


rawData <- readSingleCellData(
  sampleNames = sampleNames[1:4], 
  samplePaths = samplePaths[1:4], 
  verbose = FALSE
)

testthat::test_that(
  "a named list of dgCMatrix objects is returned after successfully reading h5 and h5ad files containing RNA and/or ADT data", {
    
    testthat::expect_type(object = rawData, type = "list")
    
    testthat::expect_identical(object = names(rawData), 
                               expected = sampleNames[1:4])
    
    testthat::expect_identical(
      object = names(rawData[["sample2_h5ad"]]), 
      expected = c("RNA", "ADT")
    )
    
    testthat::expect_s4_class(
      object = rawData[["sample1_h5"]][["RNA"]], 
      class = "dgCMatrix"
    )
    
    testthat::expect_s4_class(
      object = rawData[["sample2_h5"]][["ADT"]], 
      class = "dgCMatrix"
    )
    
  }
)


testthat::test_that(
  "reading the h5 and h5ad files of the same data produces the same result", {
    
    testthat::expect_equal(
      object = rawData[["sample2_h5ad"]][["RNA"]],
      expected = rawData[["sample2_h5"]][["RNA"]]
    )
    
    testthat::expect_equal(
      object = rawData[["sample1_h5ad"]][["RNA"]], 
      expected = rawData[["sample1_h5"]][["RNA"]]
    )
    
  }
)
