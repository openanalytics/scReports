sampleNames <- c(
  "sample2_h5",
  "sample2_h5ad",
  "sample4_h5ad"
)
samplePaths <- c(
  "../testdata/sample2/raw_feature_bc_matrix.h5", 
  "../testdata/sample2/raw_feature_bc_matrix.h5ad", 
  "../testdata/sample4/raw_feature_bc_matrix.h5ad"
)

rawData <- readSingleCellData(sampleNames = sampleNames, 
                              samplePaths = samplePaths, 
                              verbose = FALSE)

scoreHeatmapData <- computeADTHeatmapData(rawData = rawData, 
                                          type = "score")
gmHeatmapData <- computeADTHeatmapData(rawData = rawData, 
                                       type = "gm")
countsHeatmapData <- computeADTHeatmapData(rawData = rawData, 
                                           type = "counts")
pctCountsHeatmapData <- computeADTHeatmapData(rawData = rawData, 
                                              type = "pctCounts")


testthat::test_that(
  "a list of three elements is returned", {
    
  testthat::expect_type(object = scoreHeatmapData, type = "list")
    
  testthat::expect_length(object = scoreHeatmapData, n = 3)
  
  testthat::expect_equal(
    object = names(scoreHeatmapData),
    expected = c("heatmapData", "antibodyClust", "sampleClust")
  )
  
  testthat::expect_type(object = gmHeatmapData, type = "list")
  
  testthat::expect_length(object = gmHeatmapData, n = 3)
  
  testthat::expect_equal(
    object = names(gmHeatmapData),
    expected = c("heatmapData", "antibodyClust", "sampleClust")
  )
  
  testthat::expect_type(object = countsHeatmapData, type = "list")
  
  testthat::expect_length(object = countsHeatmapData, n = 3)
  
  testthat::expect_equal(
    object = names(countsHeatmapData),
    expected = c("heatmapData", "antibodyClust", "sampleClust")
  )
  
  testthat::expect_type(object = pctCountsHeatmapData, type = "list")
  
  testthat::expect_length(object = pctCountsHeatmapData, n = 3)
  
  testthat::expect_equal(
    object = names(pctCountsHeatmapData),
    expected = c("heatmapData", "antibodyClust", "sampleClust")
  )
  
  }
)

