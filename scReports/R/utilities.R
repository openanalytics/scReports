#' Get path to R markdown child document
#'
#' @param child [\code{character(1)}] \cr
#'  Name of the R markdown child document. 
#' @param report [\code{character(1)}] \cr
#'  Name of the report folder in the package that contains the R markdown 
#'  child documents. For example, 'qc' or 'integration'. 
#'
#' @return [\code{character(1)}] \cr
#'  Path to the R markdown child document.
getPathChild <- function(child, report = NA){
  
  if (is.na(report)) {
    
    path <- system.file("template", "rmarkdown", child,
                        package = "scReports")
    
  } else {
    
    path <- system.file("template", "rmarkdown", report, child,
                        package = "scReports")
    
  }
  
  
  return(path)
}



#' Extract colours from the colour universal design colour-blind friendly palette
#' 
#' @param n [\code{numeric(1)}] \cr
#'  Number of colours to extract
#' 
#' @importFrom grDevices colorRampPalette 
#'
#' @return [\code{character(1)}] \cr
#'  A vector of length \code{n} of colours.
#' @author Trishanta Padayachee
colourUniversalDesign <- function(n) {
  
  cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  palette <- grDevices::colorRampPalette(cbp1)
  colours <- palette(n)
  
  
  return(colours)
} 



#' Convert to pretty numbers
#'
#' @param x [\code{numeric}] \cr
#'  A single numeric value or a numeric vector.
#'
#' @return [\code{character}] \cr
#'  A character vector of prettified numbers.
prettyNumbers <- function(x) {
  
  int <- findInterval(x, c(0, 1e4,1e6, 1e9, 1e12))
  num <- paste0(sprintf("%.4g",x/10^(3*(int-1))), 
                c("","K","M", "B", "T")[int])
  
  
  return(num)
}



#' Determine the number of samples included in each figure and the 
#' corresponding figure heights 
#'
#' @param geom [\code{character(1)}] \cr
#'  Type of plot; either "violin", "scatter", "faceted_scatter", or "heatmap".
#' @param nBreaks [\code{integer(1)}] \cr 
#'  Number of samples that need to be plotted. 
#' @param nBreaksInStandardHeight [\code{integer(1)}] \cr
#'  Number of samples that fit in the standard figure width and height.
#' @param standardFigHeight [\code{numeric(1)}, default: \code{8.5}] \cr
#'  Standard figure height in inches.
#' @param standardFigWidth [\code{numeric(1)}, default: \code{8.5}] \cr 
#'  Standard figure width in inches.
#' @param splitFigure [\code{logical(1)}, default: \code{TRUE}] \cr
#'  Whether to split samples across multiple figures if necessary. 
#' @inheritParams commonParams
#'
#' @importFrom grid convertUnit
#' @return A list containing the following values:
#'  geom: [\code{character(1)}] the type of figure (e.g., "violin", "scatter", etc.)
#'  multipleFigures: [\code{logical(1)}] logical indicating whether figures will be 
#'   split
#'  nFigures: [\code{numeric(1)}] number of figures  
#'  width: [\code{numeric(1)}] width of plots
#'  height: [\code{numeric(1)}] height of full-size plots
#'  heightRemainder: [\code{numeric(1)}] height of the last plot if 
#'   \code{remainingSamples} if \code{TRUE}
#'  remainder: [\code{logical(1)}] whether there is a remainder when 
#'   nSamples is divided by the sample breakLimit
#'  seqStart: [\code{numeric}] sequence of start positions 
#'  seqEnd: [\code{numeric}] sequence of end positions 
#'  nFullFigs: [\code{integer(1)}] number of full-sized plots
#'
#' @export
#' @author Trishanta Padayachee
computeFigParams <- function(geom,
                             nBreaks, 
                             nBreaksInStandardHeight, 
                             standardFigWidth = 8.5,
                             standardFigHeight = 6.5,
                             splitFigure = TRUE,
                             performChecks = TRUE) {
  
  
  if (performChecks) {
    
    # Check arguments
    checkmate::assertChoice(x = geom, 
                            choices = c("violin", "scatter", 
                                        "faceted_scatter", "heatmap"))
    checkmate::assertNumber(x = nBreaks,
                             lower = 1, 
                             finite = TRUE)
    checkmate::assertNumber(x = nBreaksInStandardHeight,
                            lower = 1, 
                            finite = TRUE)
    checkmate::assertNumber(x = standardFigWidth,
                            lower = 1, 
                            finite = TRUE)
    checkmate::assertNumber(x = standardFigHeight,
                            lower = 1, 
                            finite = TRUE)
    checkmate::assertLogical(x = splitFigure,
                             any.missing = FALSE, 
                             len = 1)
    
  }
  
  
  height <- standardFigHeight
  width <- standardFigWidth
  heightRemainder <- standardFigHeight
  remainder <- FALSE
  seqStart <- 1
  seqEnd <- nBreaks
  nFullFigs <- 1
  nFigures <- 1
  multipleFigures <- FALSE
  breakLimit <- Inf
  
  # approximate calculation of break limit
  # note: graphics device errors occur for values around 175
  if (geom %in% c("scatter", "violin")) {
    
    # note: violin plot uses theme_minimal() with plot.margins of zero
    heightMarginPoints <- 
      as.numeric(theme_minimal()$plot.title$margin)[3] + 
      as.numeric(theme_minimal()$axis.title.x$margin)[1] + 
      as.numeric(theme_minimal()$axis.ticks.length) + 
      as.numeric(theme_minimal()$axis.text.x$margin)[1] + 
      theme_minimal()$text$size*(
        as.numeric(theme_minimal()$plot.title$size) +
          as.numeric(theme_minimal()$axis.text$size) +
          as.numeric(theme_minimal()$panel.grid.minor$linewidth)
      )
    heightMarginInches <- as.numeric(
      grid::convertUnit(unit(heightMarginPoints, "pt"), 
                        unitTo = "inches")
    )
    
    breakLimit <- floor((170 - heightMarginInches)/0.33)
    
  } else if (geom == "faceted_scatter") {
    
    breakLimit <- floor(170*nBreaksInStandardHeight/(4*standardFigHeight))*4
    
  } else if (geom == "heatmap") {
   
    breakLimit <- floor(170/(standardFigHeight/nBreaksInStandardHeight)) 

  }

  if (splitFigure & nBreaks >= breakLimit) {
    
    nFullFigs <- nBreaks %/% breakLimit
    nRemainingBreaks <- nBreaks %% breakLimit
    remainder <- nRemainingBreaks > 0
    nFigures <- nFullFigs + (nRemainingBreaks > 0)
    multipleFigures <- TRUE
    
    if (geom %in% c("scatter", "violin")) {
      
      # height = 0.5631659 + 0.33*breakLimit
      height <- heightMarginInches + 0.3*breakLimit + 0.1*0.3*(breakLimit)
      
      if (remainder) {
        
        seqStart <- c(seq(from = nBreaks-breakLimit+1, to = 1, by = -breakLimit), 1)
        seqEnd <- seq(from = nBreaks, to = 1, by = -breakLimit)
        
        heightRemainder <- heightMarginInches + 0.3*nRemainingBreaks + 
          0.1*0.3*(nRemainingBreaks)
        
      } else {
        
        # note: assumes that breaks correspond to the y-axis (after coord-flip) 
        # and first figure includes samples that lie furthest from the origin
        seqStart <- seq(from = nBreaks-breakLimit+1, to = 1, by = -breakLimit)
        seqEnd <- seq(from = nBreaks, to = 1, by = -breakLimit)
        
      }
      
    } else if (geom == "faceted_scatter") {
      
      # note: faceted scatter plots use four columns 
      height <- (standardFigHeight/(ceiling(nBreaksInStandardHeight/4)))*
        (ceiling(breakLimit/4))
      
      if (remainder) {
        
        seqStart <- seq(from = 1, to = nBreaks, by = breakLimit)
        seqEnd <- c(seq(from = breakLimit, to = nBreaks, by = breakLimit), nBreaks)
        
        # note: faceted scatter plots use four columns 
        heightRemainder <- (standardFigHeight/(ceiling(nBreaksInStandardHeight/4)))*
          (ceiling(nRemainingBreaks/4))
        
      } else {
        
        seqStart <- seq(from = 1, to = nBreaks, by = breakLimit)
        seqEnd <- seq(from = breakLimit, to = nBreaks, by = breakLimit)
        
      }
      
    } else if (geom %in% c("heatmap")) {
      
      height <- breakLimit*(standardFigHeight/nBreaksInStandardHeight)
      
      if (remainder) {
        
        heightRemainder <- nRemainingBreaks*(standardFigHeight/nBreaksInStandardHeight)
        
      }
      
    }
    
  } else {
    
    height <- standardFigHeight + 
      (nBreaks-nBreaksInStandardHeight)*
      (standardFigHeight/nBreaksInStandardHeight)*
      (nBreaks > nBreaksInStandardHeight)
    height <- min(height, 170)
    
  }
  
  return(list("geom" = geom,
              "multipleFigures" = multipleFigures,
              "nFigures" = nFigures,
              "width" = width,
              "height" = height, 
              "heightRemainder" = heightRemainder,
              "remainder" = remainder,
              "seqStart" = seqStart,
              "seqEnd" = seqEnd,
              "nFullFigs" = nFullFigs))
}
  


#' Display plot
#' 
#' \code{displayPlot} is used to extract a part of specific ggplot. It is 
#' assumed that the y-axis corresponds to a continuous variable and the x-axis 
#' to a discrete variable. The positional arguments \code{startPos} and 
#' \code{endPos} corresponds to the discrete x-axis variable. Zooming applies 
#' only to the y-axis. An option is included to expand the y-axis only 
#' (i.e., \code{xpansionFactor}).
#'
#' @param gg [\code{ggplot}] \cr
#'  A ggplot object as returned from the function \code{createViolinPlot}.
#' @param startPos [\code{numeric}] \cr
#'  A numeric vector representing the point along the x-axis at which each 
#'  figure begins as returned from the function \code{computeFigParams}. 
#' @param endPos [\code{numeric}] \cr
#'  A numeric vector representing the point along the x-axis at which each 
#'  figure ends as returned from the function \code{computeFigParams}.
#' @param flip [\code{logical{1}}, default: \code{FALSE}] \cr
#'  Whether to flip the Cartesian coordinates such that the discrete \code{X} 
#'  observation appears on the y-axis and the continuous \code{Y} observation 
#'  appears on the x-axis.
#' @param zoom [\code{logical{1}}, default: \code{FALSE}] \cr
#'  Whether to display a zoomed in version of the plot. Zooming is applied to 
#'  the continuous \code{Y} observation only. 
#' @param expansionFactor [\code{numeric{1}}, default: \code{0}] \cr
#'  A numeric value between \code{0} and \code{1} representing the fraction of 
#'  the y-axis range (i.e., before flipping the coordinates if \code{flip} is 
#'  \code{TRUE}) to be subtracted from the lower y-axis limit and added to the 
#'  upper y-axis limit.   
#' @param splitFigure [\code{logical{1}}, default: \code{FALSE}] \cr
#'  Whether \code{gg} should be displayed as a single figure or split into 
#'  multiple figures.
#' @inheritParams commonParams
#' 
#' @return A ggplot object.
#' @export
#' @author Trishanta Padayachee
displayPlot <- function(gg, 
                        startPos, 
                        endPos,
                        flip = FALSE,
                        zoom = FALSE,
                        expansionFactor = 0,
                        splitFigure = FALSE,
                        performChecks = TRUE) {
  
  if (performChecks) {
    
    # Check arguments
    checkmate::assertClass(x = gg,
                           classes = "ggplot")
    checkmate::assertNumeric(x = startPos,
                             lower = 1, 
                             finite = TRUE,
                             any.missing = FALSE, 
                             unique = TRUE,
                             null.ok = FALSE)
    checkmate::assertNumeric(x = endPos,
                             lower = 1, 
                             finite = TRUE,
                             any.missing = FALSE, 
                             unique = TRUE,
                             null.ok = FALSE)
    checkmate::assertLogical(x = flip,
                             any.missing = FALSE, 
                             len = 1)
    checkmate::assertLogical(x = zoom,
                             any.missing = FALSE, 
                             len = 1)
    checkmate::assertNumber(x = expansionFactor,
                            lower = 0,
                            upper = 1)
    checkmate::assertLogical(x = splitFigure,
                             any.missing = FALSE, 
                             len = 1)
    
  }
  
  # Compute Y-axis limits
  yRange <- layer_scales(gg)$y$get_limits()
  buffer <- expansionFactor*(yRange[2] - 0)
  if (zoom) {
    
    zoomRatio <- 0.2
    upperLimit <- max(yRange[1], 0) + zoomRatio*diff(yRange)
    yLimits <- c(0-buffer, upperLimit+buffer)
    
  } else {
    
    yLimits <- c(0-buffer, yRange[2]+buffer)
    
  }
  
  # Generate and display plot
  if (!splitFigure) {
    
    if (flip) {
      
      gg + 
        coord_flip(ylim = yLimits)
       
    } else {
      
      gg + 
        coord_cartesian(ylim = yLimits)
      
    }

    
  } else {
    
    for (i in 1:length(startPos)) {
      
      allBreaks <- layer_scales(gg)$x$get_labels()
      breaks <- allBreaks[startPos[i]:endPos[i]]
      labelWidth <- max(sapply(allBreaks, nchar))
      
      ggSnippet <- gg +
        coord_flip(xlim = c(startPos[i] - 0.5, endPos[i] + 0.5),
                   ylim = yLimits, 
                   expand = FALSE) +
        scale_x_discrete(breaks = breaks,
                         labels = function(x) sprintf(paste0("% ", labelWidth, "s"), x),
                         drop = FALSE) + 
        theme(axis.text.x = element_text(angle = 0))

      print(ggSnippet)
      
    }
  }
} 


#' Extract file extensions
#'
#' @param filePaths [\code{character}] \cr 
#'  One or more file paths.
#' @importFrom checkmate assertCharacter
#'   
#' @return A vector of file extensions. 
#' @export
#' @author Trishanta Padayachee
getFileExtensions <- function(filePaths) {
  
  checkmate::assertCharacter(
    x = filePaths,
    any.missing = FALSE
  )
  
  fileFormats <- lapply(
    filePaths, 
    FUN = function(x) {
      splitPath <- strsplit(basename(filePaths), split="[.]")[[1]]
      ext <- splitPath[length(splitPath)]
    }
  )
  
  fileFormats <- unlist(fileFormats)
  
  return(fileFormats)
}


