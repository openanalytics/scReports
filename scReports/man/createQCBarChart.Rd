% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createQCBarChart.R
\name{createQCBarChart}
\alias{createQCBarChart}
\title{Create a bar chart of a specific QC metric across samples}
\usage{
createQCBarChart(
  qcData,
  x,
  nameMetric,
  fillBy = x,
  barPosition = "stack",
  legendPosition = "none",
  horizontal = FALSE,
  colours = NULL,
  title = nameMetric,
  xLabel = x,
  yLabel = nameMetric,
  lLabel = NULL,
  xAngle = 0,
  makeInteractive = TRUE,
  saveFigure = FALSE,
  figureDirectory = normalizePath("./figures/")
)
}
\arguments{
\item{qcData}{[\code{data.frame}] \cr
A data frame of quality control metrics with rows corresponding
to samples and columns corresponding to different metrics as returned by 
\code{readCRQC()}.}

\item{x}{[\code{character(1)}] \cr
Name of the ID variable to plot on the x-axis.}

\item{nameMetric}{[\code{character(1)}] \cr
Name of the gene expression metric to be plotted.}

\item{fillBy}{Name of a column in \code{df} that colours should be based on.}

\item{barPosition}{[\code{character(1)}, default: "stack"] \cr
A string indicating the position adjustment of bars that 
would otherwise occupy the same space. For instance, "stack" (default) or 
"dodge".}

\item{legendPosition}{[\code{character(1)}, default: "none"] \cr
Position of the legend.}

\item{horizontal}{[\code{logical(1)}, default: \code{FALSE}] \cr
If \code{TRUE}, samples are plotted on the Y-axis resulting in horizontal 
bars.}

\item{colours}{[\code{character}] \cr
A character vector of colours to be used for plotting.}

\item{title}{[\code{character(1)}] \cr
Plot title.}

\item{xLabel}{[\code{character(1)}] \cr
X-axis label.}

\item{yLabel}{[\code{character(1)}] \cr
Y-axis label.}

\item{lLabel}{\cr
Legend name.}

\item{xAngle}{[\code{numeric(1)}, default: 0] \cr
Angle of text on the x-axis. A value between 0 and 360.}

\item{makeInteractive}{[\code{logical(1)}, default: \code{TRUE}] \cr
Should an interactive plot be returned?}

\item{saveFigure}{[\code{logical(1)}, default: \code{FALSE}] \cr 
Whether to save the figure.}

\item{figureDirectory}{[\code{character(1)}] \cr 
When \code{saveFigure} is \code{TRUE}, this is the directory in which the 
figure is saved.}
}
\value{
[\code{ggplot}] \cr 
 A ggplot object if displayPlot is FALSE
}
\description{
\code{createQCBarChart} creates a bar chart of a quality control metric for 
 multiple samples.
}
\author{
Trishanta Padayachee
}
