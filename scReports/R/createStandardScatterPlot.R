#' Create a scatter plot of two features
#'
#' \code{createScatterPlot} creates a scatter plot of two user specified
#' features.
#'
#' @param colourBy [\code{character(1)}] \cr
#'  Name of the observation by which to colour the points.
#' @param shapeBy [\code{character(1)}] \cr
#'  Name of the observation by which to shape the points.
#' @param shapes [\code{numeric}] \cr
#'  A vector of integers representing the shapes that should be used when 
#'  plotting points.
#' @param xTrans [\code{character(1)}] \cr
#'  X-axis transformation.
#' @param yTrans [\code{character(1)}] \cr
#'  Y-axis transformation.
#' @param hlineY [\code{numeric}] \cr
#'  Y intercept for horizontal line.
#' @param hlineColourBy [\code{character}] \cr
#'  A character vector of length equal to that of hlineY providing the 
#'  categories for the hline colours.
#' @param hlineTypeBy [\code{character(1)}] \cr
#'  A character vector of length equal to that of hlineY providing the 
#'  categories for the hline types.
#' @param hlineType [\code{character(1)}] \cr
#'  A character vector of line types (e.g., a combination of 'solid','dotted', 
#'  'dashed', 'dotdash'). 
#' @param pointSize [\code{numeric(1)}] \cr
#'  Size of the points.
#' @param flipCoordinates [\code{logical(1)}, default: \code{FALSE}] \cr
#'  Whether to flip coordinates such that \code{obsX} appears on the y-axis and 
#'  \code{obsY} appears on the x-axis.
#' @inheritParams commonParams
#' 
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#'
#' @return [\code{ggplot}] \cr
#'  A ggplot object.
#' @author Trishanta Padayachee
#' @export
createStandardScatterPlot <- function(
  df, 
  obsX, obsY, 
  colourBy = NULL,
  colours = NULL,
  shapeBy = NULL,
  shapes = NULL,
  xLabel = NA, yLabel = NA, lLabel = NULL,
  xLimits = NULL, 
  yLimits = NULL,
  xTrans = "identity",
  yTrans = "identity",
  hlineY = NULL,
  hlineColourBy = NULL, 
  hlineTypeBy = NULL,
  hlineType = NULL,
  xAngle = 0,
  pointSize = 3,
  flipCoordinates = FALSE){
  
  # Palette
  
  if (is.null(colours)) {
    n = length(unique(df[, colourBy]))
    colours <- colourUniversalDesign(n)
  }
  if (is.null(shapes)) {
    n = length(unique(df[, shapeBy]))
    if (n <= 6) {
      shapes <- c(16,17,15,8,3,18)[1:n]
    } else {
      shapes <- 19
    }
  }
  
  
  # Create ggplot object
  gg <- ggplot(data = df,
               aes(x = .data[[obsX]], 
                   y = .data[[obsY]],
                   colour = switch(is.null(colourBy) + 1, .data[[colourBy]], NULL),
                   shape = switch(is.null(shapeBy) + 1, .data[[shapeBy]], NULL))) +
    geom_point(size = pointSize) +
    labs(x = xLabel,
         y = yLabel) +
    theme_minimal() + 
    theme(axis.text = element_text(family = "monospace"),
          axis.text.x = element_text(angle = xAngle)) +
    scale_color_manual(name = lLabel, values = colours) + 
    scale_shape_manual(name = lLabel, values = shapes) 
  
  
  if (is.numeric(df[,obsX])) {
    gg <- gg +
      scale_x_continuous(trans = xTrans, limits = xLimits)
  }
  
  if (is.numeric(df[,obsY])) {
    gg <- gg +
      scale_y_continuous(trans = yTrans, limits = yLimits) 
  }
  
  if (!is.null(hlineY)) {
    gg <- gg +
      geom_hline(data = as.data.frame(hlineY),
                 aes(yintercept = hlineY, color = hlineColourBy, 
                     linetype = hlineTypeBy, shape = NULL)) +
      geom_text(data = as.data.frame(hlineY),
                aes(x = max(df[,obsX]), y = hlineY, 
                    color = hlineColourBy, label = hlineY,
                    shape = NULL), 
                vjust = -1, size = 3) + 
      scale_linetype_manual(name = NULL, values = hlineType) 
  }
  
  if (flipCoordinates == TRUE) {
    
    gg <- gg + coord_flip(ylim = yLimits)
    
  }
  
  # Display plot
  return(gg)
  
}

