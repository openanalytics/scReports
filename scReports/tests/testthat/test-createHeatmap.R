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

heatmap <- createHeatmap(
  df = scoreHeatmapData$heatmapData,
  obsX = "Sample",
  obsY = "Antibody",
  fillBy = "AntibodyScores",
  palette = "YlOrRd",
  xLabel = "Sample",
  yLabel = "Antibody",
  lLabel = "Percentage of \n Total Scores",
  addText = TRUE,
  makeInteractive = FALSE
)
