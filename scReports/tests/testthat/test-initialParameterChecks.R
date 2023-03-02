sampleIDs <- c(
  "sample1",
  "sample2",
  "sample4",
  "sample3_wrong_path"
)
metaFile <- c("../testdata/meta/meta.tsv")
samplePaths <- c(
  "../testdata/sample1/raw_feature_bc_matrix.h5", 
  "../testdata/sample2/raw_feature_bc_matrix.h5", 
  "../testdata/sample4/raw_feature_bc_matrix.h5",
  "../testdata/missing/raw_feature_bc_matrix.h5"
)
sampleCRQC <- c(
  "../testdata/sample1/metrics_summary.csv", 
  "../testdata/sample2/metrics_summary.csv", 
  "../testdata/sample4/metrics_summary.csv",
  "../testdata/missing/metrics_summary.csv"
)


testthat::test_that(
  "no output is produced", {
    
  testthat::expect_silent(
    object = initialParameterChecks(
      sampleNames = sampleIDs[1:2], 
      samplePaths = samplePaths[1:2],
      sampleCRQC = sampleCRQC[1:2], 
      metaFile = metaFile))
    
  }
)


testthat::test_that(
  "an error is returned when the lengths of names and paths are not equal", {
    
  testthat::expect_error(
    object = initialParameterChecks(
      sampleNames = sampleIDs[1:2], 
      samplePaths = samplePaths,
      sampleCRQC = sampleCRQC, 
      metaFile = metaFile
    )
  )
    
  }
)


testthat::test_that(
  "an error is returned when paths are incorrectly specified", {
    
  testthat::expect_error(
    object = initialParameterChecks(
      sampleNames = sampleIDs, 
      samplePaths = samplePaths,
      sampleCRQC = sampleCRQC, 
      metaFile = metaFile
    )
  )
    
  }
)
