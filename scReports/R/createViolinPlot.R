#' Create a violin plot of a feature for a list of Seurat objects or a
#' list of meta data dataframes
#'
#' \code{createViolinPlot} Creates a violin plot for a specific feature for a
#' list of Seurat objects. Violin plots for each object appear in the same
#' figure.
#'
#' @param groupBy [\code{character(1)}] \cr
#'  Name of the variable by which to group the observations (\code{obs}). A 
#'  violin plot will be created for each unique value of the \code{groupBy} 
#'  variable.
#' @param horizontal [\code{logical(1)}, default: \code{NULL}] \cr
#'  Whether to create horizontal violin plots. If \code{TRUE}, the \code{obs} 
#'  variable will be plotted on the X-axis and the \code{groupBy} variable on 
#'  the Y-axis. 
#' @param removeViolin [\code{logical(1)}, default: \code{FALSE}] \cr
#'  If \code{TRUE}, the figure will only include jittered points.
#' @param removeJitter [\code{logical(1)}, default: \code{FALSE}] \cr
#'  If \code{TRUE}, the figure will only contain violins.
#' @param violinScale [\code{character(1)}, default: "area"] \cr
#'  Either "area", "width", or "count". Corresponds to the scale argument in 
#'  \code{ggplot2::geom_violin} which indicates how the violins should be scaled. 
#'  If "area", violins are scale to have the same area (before trimming of the 
#'  tails). If "width", violins will have the same maximum width. If "count", 
#'  areas are scaled in proportion to the number of observations. 
#' @param jitterShape [\code{character or numeric}, default: 16] \cr
#'  Corresponds to the shape aesthetic in \code{ggplot2::geom_point}.
#' @param jitterSize [\code{numeric}, default: 0.01] \cr
#'  Corresponds to the size aesthetic in \code{ggplot2::geom_size}.
#' @param includeDefaultCoordinates [\code{logical(1)}, default: \code{TRUE}] \cr
#'  If \code{TRUE}, a default Cartesian coordinate system is applied using either 
#'  \code{coord_cartesian} or \code{coord_flip}. If \code{FALSE}, the default 
#'  coordinate system is not applied. 
#' @param includeDefaultPositionalScales [\code{logical(1)}, default: \code{TRUE}] \cr
#'  If \code{TRUE}, the default positional scale is applied (i.e., 
#'  \code{scale_x_discrete}. If \code{FALSE}, the default positional scale is 
#'  not applied.
#' @inheritParams commonParams
#'  
#' @import ggplot2  
#' @importFrom rlang .data
#'
#' @return [\code{ggplot}] \cr
#'  A ggplot object.
#' @author Trishanta Padayachee
#' @export
createViolinPlot <- function(df, obs, groupBy,
                             horizontal = FALSE,
                             colours = NULL,
                             title = obs,
                             xLabel = groupBy,
                             yLabel = obs,
                             xAxisTextAngle = 0,
                             yLimits = NULL,
                             xLimits = NULL,
                             removeViolin = FALSE,
                             removeJitter = FALSE,
                             violinScale = "area",
                             jitterShape = 16,
                             jitterSize = 0.01,
                             includeDefaultCoordinates = FALSE,
                             includeDefaultPositionalScales = TRUE,
                             performChecks = TRUE,
                             verbose = FALSE){
  
  if (verbose) cat("Checking the arguments...")
  
  if (performChecks) {
    
    checkmate::assertDataFrame(df,
                               col.names = "unique")
    checkmate::assertSubset(obs, colnames(df))
    checkmate::assertSubset(groupBy, colnames(df))
    checkmate::assertVector(colours,
                            any.missing = FALSE,
                            null.ok = TRUE)
    checkmate::assertCharacter(title,
                               len = 1)
    checkmate::assertCharacter(xLabel,
                               len = 1)
    checkmate::assertCharacter(yLabel,
                               len = 1)
    checkmate::assertNumeric(xAxisTextAngle, 
                             lower = 0,
                             upper = 360,
                             len = 1)
    checkmate::assertNumeric(yLimits,
                             any.missing = FALSE,
                             len = 2,
                             null.ok = TRUE)
    checkmate::assertNumeric(xLimits,
                             any.missing = FALSE,
                             len = 2,
                             null.ok = TRUE)
    checkmate::assertLogical(horizontal, 
                             any.missing = FALSE,
                             len = 1)
    checkmate::assertLogical(removeViolin, 
                             any.missing = FALSE,
                             len = 1)
    checkmate::assertLogical(removeJitter, 
                             any.missing = FALSE,
                             len = 1)
    checkmate::assertChoice(violinScale,
                            choices = c("area", "width", "count"))
    checkmate::assertNumeric(jitterSize,
                             any.missing = FALSE,
                             null.ok = FALSE)
    checkmate::assertLogical(includeDefaultCoordinates, 
                             any.missing = FALSE,
                             len = 1)
    checkmate::assertLogical(includeDefaultPositionalScales, 
                             any.missing = FALSE,
                             len = 1)
    checkmate::assertLogical(performChecks, 
                             any.missing = FALSE,
                             len = 1)
    
  }

  nLevels <- length(unique(df[, groupBy]))
  if (length(colours) != nLevels) {
    
    if (!is.null(colours)) {
      
      message("Length of colours is not equal to the number of Samples. 
              Using the default colour scheme.")
      
    }
    
    colours <- colourUniversalDesign(nLevels)
    
  }
  
  if (verbose) cat("Generating the plot...")
  
  dataSubset <- df[!(is.nan(df[, obs]) | is.na(df[, obs]) | 
                       is.na(df[, groupBy])), ]
  
  gg <- ggplot(
    data = dataSubset,
    aes(x = .data[[groupBy]],
        y = .data[[obs]],
        fill = .data[[groupBy]])
  )
  
  if (!removeViolin) {
    
    if (verbose) cat(" Adding a violin geom...")
    gg <- gg + geom_violin(scale = violinScale)
    
  }
  
  if (!removeJitter) {
    
    if (verbose) cat(" Adding a jitter geom...")
    gg <- gg + geom_jitter(shape = jitterShape,
                           size = jitterSize)
    
  }
  
  if (includeDefaultCoordinates) {
    
    if (!horizontal) {
      
      if (verbose) cat(" Including a cartesian coordinate system...")
      gg <- gg + 
       coord_cartesian(ylim = yLimits, 
                       xlim = xLimits, 
                       expand = FALSE,
                       default = TRUE)
      
    } else if (horizontal) {
      
      if (verbose) cat(" Including a flipped cartesian coordinate system...")
      gg <- gg + 
        coord_flip(ylim = yLimits, 
                   xlim = xLimits, 
                   expand = FALSE)
      
    }
     
  }
  
  gg <- gg + 
    labs(x = xLabel, y = yLabel, title = title) + 
    scale_fill_manual(values = colours)
  
  if (includeDefaultPositionalScales) {
    
    if (verbose) cat(" Including a positional scale...")
    gg <- gg + 
      scale_x_discrete(drop = FALSE)
    
  }
    
  gg <- gg + 
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(family = "monospace"),
          axis.text.x = element_text(angle = xAxisTextAngle, 
                                     hjust = 1),
          legend.position = "none",
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "in")) 
  
  
  return(gg)
}



