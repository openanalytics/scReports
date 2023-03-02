sampleNames <- c(
  "sample1_h5ad",
  "sample2_h5ad"
)
samplePaths <- c(
  "../testdata/sample1/filtered_feature_bc_matrix.h5ad",
  "../testdata/sample2/filtered_feature_bc_matrix.h5ad"
)


metricsRNA <- extractMetrics(
  path = samplePaths[1],
  sampleName = sampleNames[1]
)

testthat::test_that(
  "for RNA data, a list containing the required QC metrics is returned", {
    
    testthat::expect_type(object = metricsRNA, type = "list")
    
    expectedNames <- c("sampleName", "nCount_RNA", 
                       "nFeature_RNA", "percentMito")
    testthat::expect_true(
      object = all(colnames(metricsRNA[["cells"]]) %in% expectedNames))
    
  }
)


# sample2 doesn't have a percent_mito column
testthat::test_that(
  "when percent_mito isn't present in the data, a warning is returned", {
    
    testthat::expect_warning(
      object = extractMetrics(
        path = samplePaths[2], 
        sampleName = sampleNames[2]
      ), 
      regexp = "The following metrics are not present in the AnnData object: percent_mito"
    )
    
  }
)


metricsRNAandADT <- suppressWarnings(
  extractMetrics(path = samplePaths[2],
                 sampleName = sampleNames[2])
)

testthat::test_that(
  "for multi-modal RNA and ADT data, a list containing the required QC metrics is returned", {
    
    testthat::expect_type(object = metricsRNAandADT, type = "list")
    
    expectedNames <- c("nCount_RNA", "nFeature_RNA", "sampleName")
    testthat::expect_true(
      object = all(colnames(metricsRNAandADT[["cells"]]) %in% expectedNames)
    )
    
  }
)

