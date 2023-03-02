sampleNames <- c(
  "sample1_h5",
  "sample2_h5"
)
samplePaths <- c(
  "../testdata/sample1/raw_feature_bc_matrix.h5", 
  "../testdata/sample2/raw_feature_bc_matrix.h5"
)

data <- readSingleCellData(
  sampleNames = sampleNames,
  samplePaths = samplePaths, 
  verbose = FALSE
)
metricsRNA <- computeMetrics(
  X = data[["sample1_h5"]],
  sampleName = "sample1_h5",
  replaceNaN = FALSE
)
metricsRNAandADT <- computeMetrics(
  X = data[["sample2_h5"]],
  sampleName = "sample2_h5",
  replaceNaN = FALSE
)


testthat::test_that(
  "for RNA data, a list containing the required QC metrics is returned", {
    
    testthat::expect_type(object = metricsRNA, type = "list")
    
    expectedNames <- c("sampleName", "nCount_RNA", 
                       "nFeature_RNA", "percentMito")
    testthat::expect_true(
      object = all(colnames(metricsRNA[["cells"]]) %in% expectedNames)
    )
    
  }
)


testthat::test_that(
  "for RNA and ADT data, a list containing the required QC metrics", {
    
    testthat::expect_type(object = metricsRNAandADT, 
                          type = "list")
    
    expectedNames <- c("sampleName", "nCount_RNA", "nCount_ADT", 
                       "nFeature_RNA", "nFeature_ADT",
                       "percentMito")
    testthat::expect_true(
      object = all(colnames(metricsRNAandADT[["cells"]]) %in% expectedNames)
    )
    
  }
)

