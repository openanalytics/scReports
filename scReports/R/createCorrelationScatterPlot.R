#' Create a correlation scatter plot of two features 
#'
#' \code{createScatterPlot} Creates a scatter plot of two user specified
#' features. The Pearson's and Spearman's correlations between the two features 
#' are indicated above each plot.
#'
#' @param groupBy [\code{character(1)}] \cr
#' Name of an observation in the meta data slot of the Seurat 
#'  object(s) that can be used to distinguish one object or sample from another.
#' @param displayAxes [\code{character(1)}, default: \code{"free_x"}] \cr 
#'  Either "fixed", free", "free_x", or "free_y". 
#'  If "fixed", plots on the extreme left will have a labeled Y axis and
#'  plots on the bottom will have a labeled X axis. "free" results
#'  in labeled X and Y-axes for each plot. "free_x" produces a labeled
#'  X-axis on each plot. "free_y" produces a labeled Y-axis on each plot. 
#'  The scales argument of facet_wrap is set to \code{displayAxes} to 
#'  indicate whether to display separate x and/or y axes per figure. 
#'  To actually obtain free X or Y axes etc., one needs to set \code{displayAxes} to 
#'  "free", "free_x", or "free_y" and set the arguments \code{fixXAxis} and 
#'  \code{fixYAxis} to FALSE.  
#' @param fixXAxis [\code{logical(1)}, default: \code{TRUE}] \cr
#'  Logical indicating whether to fix the X axis across figures. This argument 
#'  takes priority over scale. 
#' @param fixYAxis [\code{logical(1)}, default: \code{TRUE}] \cr
#'  Logical indicating whether to fix the y axis across figures. This argument 
#'  takes priority over scale. 
#' @param nFigPerRow [\code{integer(1)}, default: \code{4}] \cr
#'  Number of figures included per row.    
#' @inheritParams commonParams
#' 
#' @import ggplot2
#' @importFrom stats cor setNames
#' @importFrom plyr ddply .
#'
#' @return [\code{ggplot}] \cr
#'  A ggplot object.
#' @author Trishanta Padayachee
#' @export
createCorrelationScatterPlot <- function(df, xObs, yObs, groupBy,
                                         xLabel = xObs,
                                         yLabel = yObs,
                                         xLimits = NULL,
                                         yLimits = NULL,
                                         colours = NULL,
                                         displayAxes = "free_x",
                                         fixXAxis = TRUE,
                                         fixYAxis = TRUE,
                                         nFigPerRow = 4){
  
 
  ## Check arguments
  checkmate::assertDataFrame(df,
                             col.names = "unique")
  checkmate::assertSubset(xObs, colnames(df))
  checkmate::assertSubset(yObs, colnames(df))
  checkmate::assertSubset(groupBy, colnames(df))
  checkmate::assertVector(colours,
                          any.missing = FALSE,
                          null.ok = TRUE)
  checkmate::assertCharacter(xLabel,
                             len = 1)
  checkmate::assertCharacter(yLabel,
                             len = 1)
  checkmate::assertSubset(displayAxes, c("fixed", "free", "free_x", "free_y"))
  checkmate::assertLogical(fixXAxis, 
                           any.missing = FALSE,
                           len = 1)
  checkmate::assertLogical(fixYAxis, 
                           any.missing = FALSE,
                           len = 1)
  
  nLevels <- length(unique(df[, groupBy]))
  if (length(colours) != nLevels) {
    
    if (!is.null(colours)) {
      
      message("Length of colours is not equal to the number of Samples. 
              Using the default colour scheme.")
      
    }
    
    colours <- colourUniversalDesign(nLevels)
    
  }
  
    
  ## Calculate the required number of rows and columns of scatter plots 
  nRow <- ceiling(nLevels/nFigPerRow)
  nCol <- nFigPerRow
  
  ## Remove NA and NaN values
  notMissing <- !(is.na(df[, xObs]) | is.na(df[, yObs]) | 
                  is.na(df[, groupBy]) |
                  is.nan(df[,xObs]) | is.nan(df[,yObs]))                               
  dataSubset <- df[notMissing,]
  
  ## Compute correlations and create plot titles
  corr <- plyr::ddply(
    .data = df, 
    .variables = groupBy,
    .fun = function(data) {
      pearson <- stats::cor(data[, xObs], 
                            data[, yObs], 
                            use = "pairwise.complete.obs")
      spearman <- stats::cor(data[, xObs], 
                             data[, yObs], 
                             use = "pairwise.complete.obs", 
                             method = 'spearman')
      return(c(pearson, spearman))
    }
  )
  keys <- corr[, groupBy]
  vals <- as.character(
    paste0(keys,
           "\n Pearson's Corr: ", round(corr[, 2], 2), 
           "\n Spearman's Corr: ", round(corr[, 3], 2))
  )
  lookupCorr <- stats::setNames(vals, keys)
  
  ## Create scatter plots
  gg <- ggplot(data = dataSubset,
               aes(x = .data[[xObs]], 
                   y = .data[[yObs]],
                   color = .data[[groupBy]])) +
    geom_point(shape = 16, size = 1) +
    labs(x = xLabel,
         y = yLabel) +
    facet_wrap(~ .data[[groupBy]], 
               labeller = as_labeller(lookupCorr),
               nrow = nRow,
               ncol = nCol,
               scales = displayAxes) +
    scale_color_manual(values = colours) +
    coord_cartesian(xlim = xLimits, ylim = yLimits) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          axis.text = element_text(size = 6))
  
  if (fixXAxis) {
    
    if (is.null(xLimits)) xLimits <- range(dataSubset[ , xObs])
    gg <- gg + scale_x_continuous(limits = xLimits)
    
  }
  if (fixYAxis) {
    
    if (is.null(yLimits)) yLimits <- range(dataSubset[ , yObs])
    gg <- gg + scale_y_continuous(limits = yLimits)
    
  }
  
  
  return(gg)
}
