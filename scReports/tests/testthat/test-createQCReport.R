# Raw h5 and h5ad files
sampleIDs <- sampleNames <- c(
  "sample1_h5",
  "sample1_h5ad",
  "sample2_h5",
  "sample2_h5ad",
  "sample3_h5",
  "sample3_h5ad",
  "sample4_h5",
  "sample4_h5ad",
  "sample5_h5",
  "sample5_h5ad"
)
samplePaths <- c(
  "../testdata/sample1/raw_feature_bc_matrix.h5",
  "../testdata/sample1/raw_feature_bc_matrix.h5ad",
  "../testdata/sample2/raw_feature_bc_matrix.h5",
  "../testdata/sample2/raw_feature_bc_matrix.h5ad",
  "../testdata/sample3/raw_feature_bc_matrix.h5",
  "../testdata/sample3/raw_feature_bc_matrix.h5ad",
  "../testdata/sample4/raw_feature_bc_matrix.h5",
  "../testdata/sample4/raw_feature_bc_matrix.h5ad",
  "../testdata/sample5/raw_feature_bc_matrix.h5",
  "../testdata/sample5/raw_feature_bc_matrix.h5ad"
)
sampleCRQC <- c(
  "../testdata/sample1/metrics_summary.csv",
  "../testdata/sample1/metrics_summary.csv",
  "../testdata/sample2/metrics_summary.csv",
  "../testdata/sample2/metrics_summary.csv",
  "../testdata/sample3/metrics_summary.csv",
  "../testdata/sample3/metrics_summary.csv",
  "../testdata/sample4/metrics_summary.csv",
  "../testdata/sample4/metrics_summary.csv",
  "../testdata/sample5/metrics_summary.csv",
  "../testdata/sample5/metrics_summary.csv"
)
metaFile <- "../testdata/meta/meta.tsv"

## raw-h5-and-h5ad
directory <- "./wdir-raw-h5-and-h5ad/"
createQCReport(
  TA = "Test",  
  project = "QC",
  experiment = "on raw 10x h5 files (gex/citeseq/crispr) and anndata h5ad files",
  sampleIDs = sampleIDs,
  sampleNames = sampleNames, 
  samplePaths = samplePaths, 
  sampleCRQC = sampleCRQC,
  metaFile = metaFile,
  analyst = "Jane Smith",
  horizontalBarCharts = TRUE,
  horizontalViolinPlots = TRUE,
  saveRObjects = FALSE,
  workingDirectory = directory,
  outputDirectory = directory,
  quietRendering = TRUE
)

testthat::test_that(
  "an html file is returned when input includes raw 10x h5 and anndata h5ad files", {
    
    testthat::expect_true(
      object = file.exists(paste0(directory, "/qualityControl.html"))
    )
    unlink(x = directory, recursive = TRUE)
    
  }
)

## raw-h5-and-h5ad-no-csv
directory <- "./wdir-raw-h5-and-h5ad-no-csv/"
createQCReport(
  TA = "Test",
  project = "QC",
  experiment = "on raw 10x h5 files (gex/citeseq/crispr) and anndata h5ad files with no metrics_summary.csv files",
  sampleIDs = sampleIDs[1:2],
  sampleNames = sampleNames[1:2], 
  samplePaths = samplePaths[1:2], 
  metaFile = metaFile,
  analyst = "Jane Smith",
  horizontalBarCharts = TRUE,
  horizontalViolinPlots = TRUE,
  saveRObjects = FALSE,
  workingDirectory = directory,
  outputDirectory = directory,
  quietRendering = TRUE
)

testthat::test_that(
  "an html file is returned when no metric summary csv files are available", {
    
    testthat::expect_true(
      object = file.exists(paste0(directory, "/qualityControl.html"))
    )
    unlink(x = directory, recursive = TRUE)
    
  }
)

## raw-h5-and-h5ad-no-meta
directory <- "./wdir-raw-h5-and-h5ad-no-meta/"
createQCReport(
  TA = "Test",
  project = "QC",
  experiment = "on raw 10x h5 files (gex/citeseq/crispr) and anndata h5ad files with no meta file",
  sampleIDs = sampleIDs[3:4],
  sampleNames = sampleNames[3:4], 
  samplePaths = samplePaths[3:4], 
  sampleCRQC = sampleCRQC[3:4],
  analyst = "Jane Smith",
  horizontalBarCharts = FALSE,
  horizontalViolinPlots = FALSE,
  saveRObjects = FALSE,
  workingDirectory = directory,
  outputDirectory = directory,
  quietRendering = TRUE
)

testthat::test_that(
  "an html file is returned when a meta file isn't provided", {
    
    testthat::expect_true(
      object = file.exists(paste0(directory, "/qualityControl.html"))
    )
    unlink(x = directory, recursive = TRUE)
    
  }
)

## raw-h5-gex-only
directory <- "./wdir-raw-h5-gex-only/"
createQCReport(
  TA = "Test",
  project = "QC",
  experiment = "on raw 10x h5 files containing gene expression data",
  sampleIDs = sampleIDs[1],
  sampleNames = sampleNames[1], 
  samplePaths = samplePaths[1], 
  sampleCRQC = sampleCRQC[1],
  metaFile = metaFile,
  analyst = "Jane Smith",
  horizontalBarCharts = FALSE,
  horizontalViolinPlots = FALSE,
  saveRObjects = FALSE,
  workingDirectory = directory,
  outputDirectory = directory,
  quietRendering = TRUE
)

testthat::test_that(
  "an html file is returned when there is only gene expression data", {
    
    testthat::expect_true(
      object = file.exists(paste0(directory, "/qualityControl.html"))
    )
    unlink(x = directory, recursive = TRUE)
    
  }
)

# failing-knee-plots
directory <- "./wdir-failing-knee-plots/"
createQCReport(
  TA = "Test",
  project = "QC",
  experiment = "on a raw h5 file with failing knee plot",
  sampleIDs = sampleIDs[5:8],
  sampleNames = sampleNames[5:8], 
  samplePaths = samplePaths[5:8], 
  sampleCRQC = sampleCRQC[5:8],
  metaFile = metaFile,
  analyst = "Jane Smith",
  horizontalBarCharts = FALSE,
  horizontalViolinPlots = FALSE,
  saveRObjects = FALSE,
  workingDirectory = directory,
  outputDirectory = directory,
  quietRendering = TRUE
)

testthat::test_that(
  "a qualityControl html file is returned when the knee plot barcodeRank computation fails", {
     
    testthat::expect_true(
      object = file.exists(paste0(directory, "/qualityControl.html"))
    )
    unlink(x = directory, recursive = TRUE)
    
  }
)


# single-raw-h5ad-with-zero-cells
directory <- "./wdir-single-raw-h5ad-with-zero-cells/"
createQCReport(
  TA = "Test",
  project = "QC",
  experiment = "on raw anndata objects with zero cells after minimal filtering",
  sampleIDs = sampleIDs[6],
  sampleNames = sampleNames[6], 
  samplePaths = samplePaths[6], 
  isFilteredData = FALSE,
  sampleCRQC = NULL,
  metaFile = "",
  analyst = "Jane Smith",
  horizontalBarCharts = FALSE,
  horizontalViolinPlots = FALSE,
  saveRObjects = FALSE,
  workingDirectory = directory,
  outputDirectory = directory,
  quietRendering = TRUE
)


testthat::test_that(
  "a qualityControl html file is returned for a raw h5ad object containing zero cells after minimal filtering", {
    
    testthat::expect_true(
      object = file.exists(paste0(directory, "/qualityControl.html"))
    )
    unlink(x = directory, recursive = TRUE)
    
  }
)


# Filtered h5ad files 
# (1) one without a percent_mito column
# (2) one with zero cells 
sampleIDs <- sampleNames <- c(
  "sample2_h5ad", # gex and adt with no percent_mito column
  "sample3_h5ad" # gex with zero cells 
)
samplePaths <- c(
  "../testdata/sample2/filtered_feature_bc_matrix.h5ad",
  "../testdata/sample3/filtered_feature_bc_matrix.h5ad"
)

# filtered-h5ad
directory <- "./wdir-filtered-h5ad/"
createQCReport(
  TA = "Test",
  project = "QC",
  experiment = "on filtered anndata objects",
  sampleIDs = sampleIDs,
  sampleNames = sampleNames, 
  samplePaths = samplePaths, 
  isFilteredData = TRUE,
  sampleCRQC = NULL,
  metaFile = "",
  analyst = "Jane Smith",
  horizontalBarCharts = FALSE,
  horizontalViolinPlots = FALSE,
  saveRObjects = FALSE,
  workingDirectory = directory,
  outputDirectory = directory,
  quietRendering = TRUE
)


testthat::test_that(
  "a qualityControl html file is returned for filtered h5ad objects", {
    
    testthat::expect_true(
      object = file.exists(paste0(directory, "/qualityControl.html"))
    )
    unlink(x = directory, recursive = TRUE)
    
  }
)


# single-filtered-h5ad-with-zero-cells
directory <- "./wdir-single-filtered-h5ad-with-zero-cells/"
createQCReport(
  TA = "Test",
  project = "QC",
  experiment = "on filtered anndata objects with zero cells",
  sampleIDs = sampleIDs[2],
  sampleNames = sampleNames[2], 
  samplePaths = samplePaths[2], 
  isFilteredData = TRUE,
  sampleCRQC = NULL,
  analyst = "Jane Smith",
  horizontalBarCharts = FALSE,
  horizontalViolinPlots = FALSE,
  saveRObjects = FALSE,
  workingDirectory = directory,
  outputDirectory = directory,
  quietRendering = TRUE
)


testthat::test_that(
  "a qualityControl html file is returned for a filtered h5ad object containing zero cells", {
    
    testthat::expect_true(
      object = file.exists(paste0(directory, "/qualityControl.html"))
    )
    unlink(x = directory, recursive = TRUE)
    
  }
)


# sampleIDs <- sampleNames <- c(
#   "sample1_h5"
# )
# samplePaths <- c(
#   "../testdata/sample1/raw_feature_bc_matrix.h5"
# )
# sampleCRQC <- c(
#   "../testdata/sample1/metrics_summary_multi.csv"
# )
# 
# # h5 and csv from cell ranger multi
# directory <- "./wdir-h5-and-multi-csv/"
# createQCReport(
#   TA = "Test",
#   project = "QC",
#   experiment = "on raw 10x h5 file and cell ranger multi metrics summary file",
#   sampleIDs = sampleIDs,
#   sampleNames = sampleNames,
#   samplePaths = samplePaths,
#   sampleCRQC = sampleCRQC,
#   cellRangerCount = FALSE,
#   analyst = "Jane Smith",
#   workingDirectory = directory,
#   outputDirectory = directory,
#   quietRendering = TRUE
# )
# 
# 
# testthat::test_that(
#   "a qualityControl html file is returned for a Cell Ranger multi csv file", {
#     
#     testthat::expect_true(
#       object = file.exists(paste0(directory, "/qualityControl.html"))
#     )
#     unlink(x = directory, recursive = TRUE)
#     
#   }
# )



