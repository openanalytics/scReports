#' Create a bar chart of a specific QC metric across samples 
#' 
#' \code{createQCBarChart} creates a bar chart of a quality control metric for 
#'  multiple samples. 
#'
#' @param qcData [\code{data.frame}] \cr
#'  A data frame of quality control metrics with rows corresponding
#'  to samples and columns corresponding to different metrics as returned by 
#'  \code{readCRQC()}.
#' @param x [\code{character(1)}] \cr
#'  Name of the ID variable to plot on the x-axis.
#' @param nameMetric [\code{character(1)}] \cr
#'  Name of the gene expression metric to be plotted.
#' @param barPosition [\code{character(1)}, default: "stack"] \cr
#'  A string indicating the position adjustment of bars that 
#'  would otherwise occupy the same space. For instance, "stack" (default) or 
#'  "dodge".
#' @param legendPosition [\code{character(1)}, default: "none"] \cr
#'  Position of the legend. 
#' @param horizontal [\code{logical(1)}, default: \code{FALSE}] \cr
#'  If \code{TRUE}, samples are plotted on the Y-axis resulting in horizontal 
#'  bars. 
#' @param saveFigure [\code{logical(1)}, default: \code{FALSE}] \cr 
#'  Whether to save the figure.
#' @param figureDirectory [\code{character(1)}] \cr 
#'  When \code{saveFigure} is \code{TRUE}, this is the directory in which the 
#'  figure is saved.
#' @inheritParams commonParams
#'
#' @import ggplot2
#' @importFrom plotly ggplotly style
#' @importFrom htmlwidgets saveWidget
#' @importFrom magrittr %>%
#'
#' @return [\code{ggplot}] \cr 
#'  A ggplot object if displayPlot is FALSE
#' @author Trishanta Padayachee
#' @export 
createQCBarChart <- function(qcData, 
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
                             figureDirectory = normalizePath("./figures/")) {
  
  nLevels <- length(unique(qcData[, fillBy]))
  
  # Colours
  if (length(colours) != nLevels) {
    
    if (!is.null(colours)) {
      
      message("Length of colours is not equal to the number of Samples. 
              Using the default colour scheme.")
      
    }
    
    colours <- colourUniversalDesign(nLevels)
    
  }
  
  # Replace NA values with 0
  qcData[is.na(qcData)] <- 0
  
  # Format numbers
  prettyNumbers <- function(x) {
    
    int <- findInterval(x, c(0, 1e4,1e6, 1e9, 1e12))
    num <- paste0(sprintf("%.4g",x/10^(3*(int-1))), 
                  c("","K","M", "B", "T")[int])
    
  }
  
  # GGPlot
  
  
  gg <- ggplot(data = qcData, 
               aes(x = get(x),
                   y = get(nameMetric), 
                   fill = get(fillBy),
                   color = get(fillBy),
                   text = prettyNumbers(get(nameMetric)))) +
    geom_col(position = barPosition) +
    labs(title = title, 
         x = xLabel, 
         y = yLabel) +
    scale_colour_manual(name = NULL, values = colours) + 
    scale_fill_manual(name = lLabel, values = alpha(colours, 0.8)) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = xAngle, 
                                     hjust = 1),
          legend.position = legendPosition)
  
  if (horizontal) gg <- gg + coord_flip()
  
  if (makeInteractive) {
    
    gg <- plotly::ggplotly(gg, tooltip = c("text") ) %>%
      plotly::style(textposition = "auto")
    
    if (saveFigure == TRUE) {
      
      if (dir.exists(figureDirectory)) {
        
        htmlwidgets::saveWidget(widget = gg,
                                file = paste0(figureDirectory, "/",
                                              make.names(nameMetric), ".html"),
                                selfcontained = TRUE)
        
      } else {
        
        message("HTML of figure not saved. Output directory does not exist.")
        
      }
    }
    
  } 
  
  
  
  return(gg)
  
}