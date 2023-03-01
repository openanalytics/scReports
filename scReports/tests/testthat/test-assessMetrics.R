sampleIDs <- c(
  "sample1",
  "sample2",
  "sample4",
  "sample3"
)
samplePaths <- c(
  "../testdata/sample1/metrics_summary.csv", 
  "../testdata/sample2/metrics_summary.csv", 
  "../testdata/sample4/metrics_summary.csv",
  "../testdata/sample3/metrics_summary.csv"
)

alerts <- c(
  "Error: No Cells Detected<br>Warning: Low Fraction Valid Barcodes<br>Warning: Low Fraction Reads Confidently Mapped To Transcriptome<br>Warning: High Fraction of Reads Mapped Antisense to Genes<br>Warning: Low Barcode Q30 Fraction (Illumina I7 Read for Single Cell 3\\' v1, R1 for Single Cell 3\\' v2/v3 and Single Cell 5\\')<br>Error: Low RNA Read Q30 Fraction (Illumina R1 for Single Cell 3\\' v1 and Single Cell 5\\' paired end, R2 for Single Cell 3\\' v2/v3 and Single Cell 5\\' R2-only)<br>Error: Low Sample Index Q30 Fraction (Illumina I5 Read for Single Cell 3\\' v1, I7 for Single Cell 3\\' v2/v3 and Single Cell 5\\')<br>Warning: Low UMI Q30 Fraction (Illumina R2 Read for Single Cell 3\\' v1, R1 for Single Cell 3\\' v2/v3 and Single Cell 5\\')<br>",
  "", 
  "Warning: Low Number of Cells Detected<br>Error: Low Fraction Valid Barcodes<br>Error: Low Fraction Reads Confidently Mapped To Transcriptome<br>Error: High Fraction of Reads Mapped Antisense to Genes<br>Error: Low Barcode Q30 Fraction (Illumina I7 Read for Single Cell 3\\' v1, R1 for Single Cell 3\\' v2/v3 and Single Cell 5\\')<br>Warning: Low RNA Read Q30 Fraction (Illumina R1 for Single Cell 3\\' v1 and Single Cell 5\\' paired end, R2 for Single Cell 3\\' v2/v3 and Single Cell 5\\' R2-only)<br>Warning: Low Sample Index Q30 Fraction (Illumina I5 Read for Single Cell 3\\' v1, I7 for Single Cell 3\\' v2/v3 and Single Cell 5\\')<br>Error: Low UMI Q30 Fraction (Illumina R2 Read for Single Cell 3\\' v1, R1 for Single Cell 3\\' v2/v3 and Single Cell 5\\')<br>",
  ""
)                                            
alertReference <- data.frame(
  "Sample" = sampleIDs, 
  "Status" = c("Error", "OK", "Error", "OK"),
  "Alerts" = alerts
)

metrics <- readCRQC(sampleIDs, samplePaths)
tableOfAlerts <- assessMetrics(metrics)

testthat::test_that(
  "assessMetrics returns a data.frame equal to the saved alerts file", {
    
    testthat::expect_s3_class(object = tableOfAlerts, 
                              class = "data.frame")
    
    testthat::expect_equal(object = tableOfAlerts, 
                           expected = alertReference, 
                           ignore_attr = TRUE)
    
  }
)
