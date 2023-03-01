primaryQC <- data.frame(
  "sampleName" = rep(x = c("Sample1", "Sample2", "Sample3"), each = 100),
  "nCount_RNA" = rnbinom(n = 300, size = 0.9, prob = 0.2), 
  "nFeature_RNA" = rnbinom(n = 300, size = 0.5, prob = 0.2),
  "percentMito" = runif(n = 300, min = 0, max = 100),
  "nCount_ADT" = rnbinom(n = 300, size = 2, prob = 0.2),
  "nFeature_ADT" = rnbinom(n = 300, size = 1.5, prob = 0.2), 
  row.names = paste0("Barcode-", 1:300)
) 

scatterPlot <- createCorrelationScatterPlot(
  df = primaryQC,
  xObs = "nCount_RNA", 
  yObs = "nFeature_RNA",
  groupBy = "sampleName",
  xLabel = "UMI Count",
  yLabel = "Number of Genes"
)

