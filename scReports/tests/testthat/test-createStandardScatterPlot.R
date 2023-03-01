totals <- c(9183, 7413, 6535, 5331, 3898, 2230, 2109, 319, 79, 73, 63, 62,
            61, 60, 59, 54, 53, 52, 50, 49, 49, 49, 48, 47, 46, 46, 44, 44, 
            44, 43, 43, 43, 43, 42, 42, 41, 41, 41, 41, 41, 41, 40, 40, 40, 
            40, 39, 39, 39, 39, 38, 38, 38, 38, 38, 37, 37, 37, 37, 37, 36, 
            36, 36, 36, 36, 36, 36, 36, 36, 36, 35, 35, 35, 35, 35, 35, 
            35, 35, 35, 35, 35, 35, 35, 34, 34, 34, 34, 34, 34, 34, 34, 
            34, 34, 34, 34, 33, 33, 33, 32, 32, 32, 32, 32, 32, 32, 32, 32, 
            31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 30, 30, 30, 29, 29, 
            28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 27, 26, 26, 26, 26, 26, 
            26, 26, 26, 26, 25, 25, 25, 25, 25, 25, 24, 24, 24, 24, 24, 24, 
            23, 23, 23, 23, 23, 22, 22, 22, 22, 22, 22, 21, 21, 21, 21, 21, 
            21, 21, 21, 20, 20, 20, 20, 20, 20, 20, 20, 19, 19, 19, 19, 19,
            18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 17, 17, 17, 
            16, 16, 16, 16, 16, 16, 16, 15, 15, 15, 15, 15, 14, 14, 14, 
            13, 13, 13, 12, 12, 12, 12, 12, 12, 11, 11, 11, 11, 11, 11, 
            10, 10, 10, 10, 10, 10, 10, 10,  5,  4,  2,  2,  2,  1,  0,  
            0, 0,  0,  0,  0)
kneeData <- data.frame(
  "type" = rep(c("RNA", "ADT"), each = 250),
  "rank" = rep(1:250, times = 2),
  "total" = c(totals, totals + 0.2*totals), 
  "fitted" = rep(NA, times = 500),
  row.names = paste0("Barcode-", 1:500)
) 
hlineY = c(1887,1887,2109,2109)

kneePlot <- createStandardScatterPlot(
  df = kneeData,
  obsX = "rank",
  obsY = "total", 
  colourBy = "type",
  colours = c("cornflowerblue", "chocolate"),
  shapeBy = "type",
  shape = c(17,16),
  xLabel = "Barcode Rank", 
  yLabel = "Total UMI Counts",
  xTrans = "log10",
  yTrans = "log10",
  hlineY = hlineY,
  hlineColourBy = c("RNA", "RNA", "ADT", "ADT"),
  hlineTypeBy = c("knee", "inflection", "knee", "inflection"),
  hlineType = c("dotted", "dashed"),
  pointSize = 1
)


testthat::test_that(
  "a ggplot object is returned", {
  
  expect_s3_class(object = kneePlot, class = "ggplot")
  
  expect_s3_class(object = kneePlot$layers[[1]]$geom, class = "GeomPoint")
  
  expect_s3_class(object = kneePlot$layers[[2]]$geom, class = "GeomHline")
  
  expect_s3_class(object = kneePlot$layers[[3]]$geom, class = "GeomText")
  
  }
)


testthat::test_that(
  "labels are specified correctly", {
  
  expect_identical(object = kneePlot$labels$x, expected = "Barcode Rank")
    
  expect_identical(object = kneePlot$labels$y, expected = "Total UMI Counts")
  
  }
)
