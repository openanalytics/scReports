% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readQualityControlMetrics.R
\name{convertCellRangerMultiToCellRangerCount}
\alias{convertCellRangerMultiToCellRangerCount}
\title{Convert Cell Ranger multi metric_summary.csv format to Cell Ranger count format}
\usage{
convertCellRangerMultiToCellRangerCount(path)
}
\arguments{
\item{path}{[\code{character}] \cr 
Path to a metrics_summary.csv file obtained through the Cell Ranger multi 
pipeline.}
}
\value{
A data frame of QC metrics in the layout of the Cell Ranger count
 metrics_summary.csv files.
}
\description{
Convert Cell Ranger multi metric_summary.csv format to Cell Ranger count format
}
