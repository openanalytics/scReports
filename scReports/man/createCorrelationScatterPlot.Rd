% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createCorrelationScatterPlot.R
\name{createCorrelationScatterPlot}
\alias{createCorrelationScatterPlot}
\title{Create a correlation scatter plot of two features}
\usage{
createCorrelationScatterPlot(
  df,
  xObs,
  yObs,
  groupBy,
  xLabel = xObs,
  yLabel = yObs,
  xLimits = NULL,
  yLimits = NULL,
  colours = NULL,
  displayAxes = "free_x",
  fixXAxis = TRUE,
  fixYAxis = TRUE,
  nFigPerRow = 4
)
}
\arguments{
\item{df}{[\code{list}] \cr
A \code{data.frame}.}

\item{xObs}{Name of the observation to plot on the X-axis. Should correspond 
to a column in \code{df}.}

\item{yObs}{Name of the observation to plot on the Y-axis. Should correspond 
to a column in \code{df}.}

\item{groupBy}{[\code{character(1)}] \cr
Name of an observation in the meta data slot of the Seurat 
 object(s) that can be used to distinguish one object or sample from another.}

\item{xLabel}{[\code{character(1)}] \cr
X-axis label.}

\item{yLabel}{[\code{character(1)}] \cr
Y-axis label.}

\item{xLimits}{[\code{numeric(2)}] \cr
X-axis limits.}

\item{yLimits}{[\code{numeric(2)}] \cr
Y-axis limits.}

\item{colours}{[\code{character}] \cr
A character vector of colours to be used for plotting.}

\item{displayAxes}{[\code{character(1)}, default: \code{"free_x"}] \cr 
Either "fixed", free", "free_x", or "free_y". 
If "fixed", plots on the extreme left will have a labeled Y axis and
plots on the bottom will have a labeled X axis. "free" results
in labeled X and Y-axes for each plot. "free_x" produces a labeled
X-axis on each plot. "free_y" produces a labeled Y-axis on each plot. 
The scales argument of facet_wrap is set to \code{displayAxes} to 
indicate whether to display separate x and/or y axes per figure. 
To actually obtain free X or Y axes etc., one needs to set \code{displayAxes} to 
"free", "free_x", or "free_y" and set the arguments \code{fixXAxis} and 
\code{fixYAxis} to FALSE.}

\item{fixXAxis}{[\code{logical(1)}, default: \code{TRUE}] \cr
Logical indicating whether to fix the X axis across figures. This argument 
takes priority over scale.}

\item{fixYAxis}{[\code{logical(1)}, default: \code{TRUE}] \cr
Logical indicating whether to fix the y axis across figures. This argument 
takes priority over scale.}

\item{nFigPerRow}{[\code{integer(1)}, default: \code{4}] \cr
Number of figures included per row.}
}
\value{
[\code{ggplot}] \cr
 A ggplot object.
}
\description{
\code{createScatterPlot} Creates a scatter plot of two user specified
features. The Pearson's and Spearman's correlations between the two features 
are indicated above each plot.
}
\author{
Trishanta Padayachee
}
