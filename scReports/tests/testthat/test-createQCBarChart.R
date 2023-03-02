metrics <- data.frame(
  "Sample ID" = c("pbmc_1k_protein_v3", "pbmc_1k_v2", "A375_10k_CRISPR_SI"), 
  "Estimated Number of Cells" = c(713, 996, 11503),
  "Reads Mapped Confidently to Transcriptome" = c(51.3, 60.8, 67.9),
  "CRISPR: Q30 Bases in Barcode" = c(NA, NA, 96.7),
  check.names = FALSE
) 


barChartVertical <- createQCBarChart(
  qcData = metrics,
  x = "Sample ID",
  nameMetric = "Estimated Number of Cells", 
  horizontal = FALSE,
  xLabel = "Sample",
  xAngle = 90,
  makeInteractive = FALSE
)
barChartHorizontal <- createQCBarChart(
  qcData = metrics,
  x = "Sample ID",
  nameMetric = "CRISPR: Q30 Bases in Barcode", 
  horizontal = TRUE,
  xLabel = "Sample",
  makeInteractive = FALSE
)

testthat::test_that(
  "a ggplot object is returned", {
    
    testthat::expect_s3_class(object = barChartVertical, 
                              class = "ggplot")
    
    testthat::expect_s3_class(object = barChartHorizontal, 
                              class = "ggplot")
  
    testthat::expect_s3_class(object = barChartVertical$layers[[1]]$geom, 
                              class = "GeomCol")
  
  }
)


testthat::test_that(
  "the values are correct", {
    
    yValues <- metrics[, "Estimated Number of Cells"]
  
    testthat::expect_equal(object = layer_data(barChartVertical)$y, 
                           expected = yValues)
  
    testthat::expect_identical(object = layer_data(barChartVertical)$text,
                               expected = prettyNumbers(yValues))
  
  }
)


testthat::test_that(
  "a horizontal bar chart is returned when horizontal is TRUE", {
  
    testthat::expect_false(
      object = inherits(barChartVertical$coordinates, "CoordFlip")
    )
    
    testthat::expect_s3_class(object = barChartHorizontal$coordinates, 
                              class = "CoordFlip")
  
  }
)


testthat::test_that(
  "the labels are correctly specified", {
  
    testthat::expect_identical(object = barChartVertical$labels$x, 
                               expected = "Sample")
    
    testthat::expect_identical(object = barChartVertical$labels$y,
                               expected = "Estimated Number of Cells")
  
    testthat::expect_identical(object = barChartVertical$labels$title,
                               expected = "Estimated Number of Cells")
  
    testthat::expect_identical(object = barChartHorizontal$labels$x,
                               expected = "Sample")
  
    testthat::expect_identical(object = barChartHorizontal$labels$y, 
                               expected = "CRISPR: Q30 Bases in Barcode")
  
    testthat::expect_identical(object = barChartHorizontal$labels$title,
                               expected = "CRISPR: Q30 Bases in Barcode")
  
  }
)


testthat::test_that(
  "NA values are displayed as zeros", {
  
    yValues <- metrics[, "CRISPR: Q30 Bases in Barcode"]
    yValues[is.na(yValues)] <- 0
  
    testthat::expect_equal(object = layer_data(barChartHorizontal)$y, 
                           expected = yValues)
  
    testthat::expect_identical(object = layer_data(barChartHorizontal)$text,
                               expected = prettyNumbers(yValues))
  
  }
)


barChartInteractive <- createQCBarChart(
  qcData = metrics,
  x = "Sample ID",
  nameMetric = "Reads Mapped Confidently to Transcriptome", 
  horizontal = TRUE,
  xLabel = "Sample",
  makeInteractive = TRUE
)

testthat::test_that(
  "an interactive plot is returned when makeInteractive is TRUE", {
  
    testthat::expect_s3_class(object = barChartInteractive, 
                              class = "plotly")
  
  }
)




