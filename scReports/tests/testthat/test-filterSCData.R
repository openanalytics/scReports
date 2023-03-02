sampleNames <- c(
  "sample1",
  "sample2"
)
samplePaths <- c(
  "../testdata/sample1/raw_feature_bc_matrix.h5", 
  "../testdata/sample2/raw_feature_bc_matrix.h5"
)
rawData <- readSingleCellData(
  samplePaths = samplePaths, 
  sampleNames = sampleNames,
  verbose = FALSE
)
metricsRNA <- computeMetrics(
  X = rawData[["sample1"]],
  sampleName = "sample1"
)
metricsRNAandADT <- computeMetrics(
  X = rawData[["sample2"]],
  sampleName = "sample2"
)


out <- filterSCData(
  X = rawData[["sample1"]],
  qcMetrics = metricsRNA,
  sampleName = "sample1",
  verbose = FALSE
)

testthat::test_that(
  "filtering is performed correctly for RNA data", {
    
    testthat::expect_true(
      object = min(out$qcMetricsFiltered[["cells"]][, "nCount_RNA"]) >= 10
    )
    
    testthat::expect_true(
      object = min(out$qcMetricsFiltered[["cells"]][, "nFeature_RNA"]) >= 1
    )
    
    testthat::expect_true(
      object = min(out$qcMetricsFiltered[["features"]][["RNA"]][, "nCells_RNA"]) >= 3
    )
    
  }
)


testthat::test_that(
  "a warning is returned when the RNA and ADT modalities aren't present and applyMultiModalCellFiltering is set to TRUE", {
    
    testthat::expect_warning(
      object = filterSCData(
        X = rawData[["sample1"]],
        qcMetrics = metricsRNA,
        sampleName = "sample1",
        applyMultiModalCellFiltering = TRUE,
        verbose = TRUE
      )
    )
    
  }
)


out <- filterSCData(
  X = rawData[["sample2"]],
  qcMetrics = metricsRNAandADT,
  sampleName = "sample2",
  applyMultiModalCellFiltering = FALSE,
  cellFilteringMinGenes = 2,
  cellFilteringMaxGenes = 10,
  cellFilteringMinCountsRNA = 1,
  cellFilteringMaxCountsRNA = 20,
  cellFilteringMinCountsADT = 1,
  cellFilteringMaxCountsADT = 4,
  cellFilteringMaxPercentMito = 20,
  geneFilteringMinCells = 1,
  antibodyFilteringMinCells = 1
)

testthat::test_that(
  "filtering is performed correctly for multi-modal RNA and ADT data", {
    
    testthat::expect_true(
      object = min(out$qcMetricsFiltered[["cells"]][, "nFeature_RNA"]) >= 2
    )
    testthat::expect_true(
      object = max(out$qcMetricsFiltered[["cells"]][, "nFeature_RNA"]) <= 10
    )
    testthat::expect_true(
      object = min(out$qcMetricsFiltered[["cells"]][, "nCount_RNA"]) >= 1
    )
    testthat::expect_true(
      object = max(out$qcMetricsFiltered[["cells"]][, "nCount_RNA"]) <= 20
    )
    testthat::expect_true(
      object = min(out$qcMetricsFiltered[["cells"]][, "nCount_ADT"]) >= 1
    )
    testthat::expect_true(
      object = max(out$qcMetricsFiltered[["cells"]][, "nCount_ADT"]) <= 4
    )
    testthat::expect_true(
      object = min(out$qcMetricsFiltered[["features"]][["RNA"]][, "nCells_RNA"]) >= 1
    )
    testthat::expect_true(
      object = min(out$qcMetricsFiltered[["features"]][["ADT"]][, "nCells_ADT"]) >= 1
    )
    
  }
)
