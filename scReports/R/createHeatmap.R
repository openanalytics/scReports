#' Create a heatmap
#'
#' @param limits [\code{numeric}, default: \code{NULL}] \cr
#'  A vector of length two representing the range of values associated with 
#'  the colour palette.
#' @inheritParams commonParams
#'  
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal brewer.pal.info  
#'
#' @return [\code{ggplot}] \cr 
#'  A ggplot object.
#' @author Trishanta Padayachee
#' @export
createHeatmap <- function(df,
                          obsX, obsY, fillBy,
                          palette = "YlOrRd",
                          xLabel = NULL, 
                          yLabel = NULL,
                          lLabel = NULL,
                          addText = TRUE,
                          textSize = 2.5,
                          xAngle = 90,
                          limits = NULL,
                          makeInteractive = TRUE) {
  
  
  if (is.null(limits)) {
    
    limits <- c(min(df[, fillBy]), max(df[, fillBy]))
    
  }
  
  gg <- ggplot(data = df,
               aes(x = !!sym(obsX),
                   y = !!sym(obsY),
                   fill = !!sym(fillBy),
                   text = round(get(fillBy),1))) + 
    geom_tile(stat = "identity",
              color = "white",
              width = 1) +
    labs(x = xLabel,
         y = yLabel) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = xAngle,
                                     hjust = 1)) + 
    scale_fill_gradientn(name = lLabel,
                         colours = 
                           RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[palette,]$maxcolors,
                                                    palette),
                         na.value = 'white',
                         guide = guide_colorbar(frame.colour = "white", ticks.colour = "white"),
                         
                         limits = limits)
  
  if (addText == TRUE){
    
    gg <- gg + geom_text(aes(label = ifelse(is.na(.data[[fillBy]]), "", paste0(round(.data[[fillBy]],1)))),
                         size = textSize)
    
  }
  
  if (makeInteractive) {
    
    gg <- plotly::ggplotly(gg, 
                           tooltip = c("text")) 
    
  }
  
  
  return(gg)
}
